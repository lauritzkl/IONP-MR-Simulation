/**
 * @file MRsequence.cpp
 *
 * Simulation kernel. Reading order: random streams and pulse timing; spatial
 * queries and material selection; segment/boundary walking; initial sampling;
 * Sequence orchestration, proton evolution and complex-signal aggregation.
 * See docs/source-guide.md for the data flow and numerical invariants.
 */

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <omp.h>

#include "mptmacros.h"

#include "RandomList.hpp"
#include "Point3D.hpp"
#include "Vector3D.hpp"
#include "Particle.hpp"
#include "Proton.hpp"
#include "Ionp.hpp"
#include "Aggregate.hpp"
#include "SimMatrix.hpp"
#include "Data.hpp"
#include "CorrectedHybridField.hpp"
#include "MRsequence.hpp"
#include "OutputUtils.hpp"
#include "SimSpace.hpp"

// Helper functions that are non-member beacuse they don't need access to MRsequence internals
namespace {

// Existing numerical controls: MIN_TIME_STEP is the base recording cadence,
// not a lower bound on diffusion integration steps. Variable steps subdivide
// a base step by VARSTEP_FACTOR; keep these constants and rounding rules stable.
struct SimConsts {
	static constexpr double MIN_TIME_STEP = 1e-5;   // time threshold for data saving
	static constexpr double VARSTEP_FACTOR = 64.0;  // factor to reduce time step for variable stepping
	static constexpr double FIELD_RANGE_FACTOR = 8.0; // distance check for variable stepping
	static constexpr unsigned int INIT_RETRY_LIMIT = 1000000;
};

// Stable stream tags separate placement, diffusion and relaxation draws.
// Values feed seed mixing; changing/reordering tags changes seeded trajectories.
enum class ProtonRandomStream : std::uint64_t {
	PositionX = 0,
	PositionY = 1,
	PositionZ = 2,
	DiffusionAndBoundaries = 3,
	RelaxationOutside = 4,
	RelaxationInside = 5,
	SubstepRelaxationOutside = 6,
	SubstepRelaxationInside = 7,
	ImportanceInitialization = 8,
	ImportanceFallbackPositionX = 9,
	ImportanceFallbackPositionY = 10,
	ImportanceFallbackPositionZ = 11,
	ImportanceFallbackInitialization = 12
};

// Derive streams from the run seed and proton index, not the OpenMP worker.
// Each proton owns its generators; scheduling does not choose its random stream.
std::uint32_t proton_random_seed(std::uint32_t run_seed,
								 std::size_t proton_index,
								 ProtonRandomStream stream) {
	std::uint64_t value =
		static_cast<std::uint64_t>(run_seed) +
		0x9e3779b97f4a7c15ULL *
			(static_cast<std::uint64_t>(proton_index) + 1ULL);
	value ^= 0xbf58476d1ce4e5b9ULL *
		(static_cast<std::uint64_t>(stream) + 1ULL);
	value += 0x9e3779b97f4a7c15ULL;
	value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
	value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
	value ^= value >> 31U;
	return static_cast<std::uint32_t>(value ^ (value >> 32U));
}

void place_proton_random_sphere(Proton& proton,
								double center_x,
								double center_y,
								double center_z,
								double minimum_radius,
								double maximum_radius,
								RandomList& random_unit) {
	const double theta = TWO_PI * random_unit.getNumber();
	const double phi = std::acos(1.0 - 2.0 * random_unit.getNumber());
	double radius = 0.0;
	do {
		radius = maximum_radius * std::cbrt(random_unit.getNumber());
	} while (radius < minimum_radius);

	proton.placeRandom(center_x + std::sin(phi) * std::cos(theta) * radius,
					   center_y + std::sin(phi) * std::sin(theta) * radius,
					   center_z + std::cos(phi) * radius,
					   true);
}

// Map SE/MSE half-echo pulse times onto the base-step indices using round().
// XOR preserves the parity of pulses that map to the same discrete step.
std::vector<unsigned char> build_refocusing_schedule(int sequence_num,
													 double echo_spacing,
													 double time_step,
													 double echo_time) {
	if (time_step <= 0.0 || echo_time <= 0.0) {
		return {};
	}
	const std::size_t step_count =
		static_cast<std::size_t>(std::ceil(echo_time / time_step)) + 2;
	std::vector<unsigned char> schedule(step_count, 0);
	auto toggle_step = [&](unsigned int step) {
		if (step < schedule.size()) {
			schedule[step] ^= 1U;
		}
	};

	if (sequence_num == 1) {
		toggle_step(SC_UI(round(echo_spacing / (2 * time_step))));
	}
	else if (sequence_num == 2) {
		for (unsigned int k = 0;
			 k < SC_UI(round(echo_time / echo_spacing));
			 ++k) {
			toggle_step(
				SC_UI(round((2 * k + 1) * echo_spacing / (2 * time_step))));
		}
	}
	return schedule;
}

inline void apply_spin_echo_refocusing(
	double& phase,
	unsigned int step,
	const std::vector<unsigned char>& refocusing_schedule) {
	if (step < refocusing_schedule.size() && refocusing_schedule[step] != 0U) {
		phase *= -1.;
	}
}

// Update phase including dephasing and refocusing for SE/MSE sequences
inline void update_phase(double& phase,
						 double bfield_tot,
						 double time_step,
						 double agar_relax,
						 unsigned int step,
						 bool smallerTimeStep,
						 const std::vector<unsigned char>& refocusing_schedule) {
	double dephase = D_GI* bfield_tot * time_step + agar_relax;
	phase += dephase; //phase at t
	if (!smallerTimeStep){
		apply_spin_echo_refocusing(phase, step, refocusing_schedule);
	}
}

} // namespace

std::size_t IonpSpatialIndex::CellKeyHash::operator()(const CellKey& key) const {
	std::size_t h = static_cast<std::size_t>(key.x);
	h = h * 73856093u ^ (static_cast<std::size_t>(key.y) * 19349663u);
	h ^= (static_cast<std::size_t>(key.z) * 83492791u);
	return h;
}

IonpSpatialIndex::IonpSpatialIndex(
	const std::vector<Ionp>& ionpInstances,
	double minimum_cell_size) {
	rebuild(ionpInstances, minimum_cell_size);
}

// Snapshot source coordinates into bins. Candidate lists accelerate exact
// geometry tests; they are not a cutoff for the full magnetic-source sum.
void IonpSpatialIndex::rebuild(
	const std::vector<Ionp>& ionpInstances,
	double minimum_cell_size) {
	m_cells.clear();
	m_variable_step_cells.clear();
	m_all_indices.clear();
	m_cell_keys.clear();
	m_all_indices.reserve(ionpInstances.size());
	m_cell_keys.reserve(ionpInstances.size());
	m_max_outer_radius = 0.0;
	m_max_core_radius = 0.0;
	m_has_coatings = false;
	m_translation_x = 0.0;
	m_translation_y = 0.0;
	m_translation_z = 0.0;

	for (const Ionp& ionp : ionpInstances) {
		m_max_outer_radius = std::max(m_max_outer_radius, ionp.outer_radius());
		m_max_core_radius = std::max(m_max_core_radius, ionp.radius());
		m_has_coatings = m_has_coatings || ionp.coating() > 0.0;
	}
	m_cell_size = std::max(m_max_outer_radius, minimum_cell_size);
	if (!(m_cell_size > 0.0) || !std::isfinite(m_cell_size)) {
		m_cell_size = 1.0;
	}
	m_variable_step_cell_size =
		(m_max_outer_radius > 0.0)
			? SimConsts::FIELD_RANGE_FACTOR * m_max_outer_radius
			: 1.0;

	for (unsigned int idx = 0; idx < SC_UI(ionpInstances.size()); ++idx) {
		m_all_indices.push_back(idx);
		const CellKey key = key_for_point(ionpInstances[idx].xt());
		m_cell_keys.push_back(key);
		m_cells[key].push_back(idx);
		const Point3D position = ionpInstances[idx].xt();
		const CellKey variable_key{
			static_cast<int>(std::floor(position.x() / m_variable_step_cell_size)),
			static_cast<int>(std::floor(position.y() / m_variable_step_cell_size)),
			static_cast<int>(std::floor(position.z() / m_variable_step_cell_size))
		};
		m_variable_step_cells[variable_key].push_back(idx);
	}
}

void IonpSpatialIndex::translate(double delta_x, double delta_y, double delta_z) {
	m_translation_x += delta_x;
	m_translation_y += delta_y;
	m_translation_z += delta_z;
}

bool IonpSpatialIndex::empty() const {
	return m_all_indices.empty();
}

double IonpSpatialIndex::max_outer_radius() const {
	return m_max_outer_radius;
}

double IonpSpatialIndex::max_core_radius() const {
	return m_max_core_radius;
}

bool IonpSpatialIndex::has_coatings() const {
	return m_has_coatings;
}

const std::vector<unsigned int>& IonpSpatialIndex::all_indices() const {
	return m_all_indices;
}

IonpSpatialIndex::CellKey IonpSpatialIndex::key_for_point(const Point3D& point) const {
	return {
		static_cast<int>(std::floor((point.x() - m_translation_x) / m_cell_size)),
		static_cast<int>(std::floor((point.y() - m_translation_y) / m_cell_size)),
		static_cast<int>(std::floor((point.z() - m_translation_z) / m_cell_size))
	};
}

void IonpSpatialIndex::append_candidates_in_box(double min_x,
												double max_x,
												double min_y,
												double max_y,
												double min_z,
												double max_z,
												std::vector<unsigned int>& candidates,
												double translation_x,
												double translation_y,
												double translation_z,
												bool sort_candidates) const {
	if (empty()) {
		return;
	}

	const CellKey min_key{
		static_cast<int>(std::floor((min_x - m_translation_x - translation_x) / m_cell_size)),
		static_cast<int>(std::floor((min_y - m_translation_y - translation_y) / m_cell_size)),
		static_cast<int>(std::floor((min_z - m_translation_z - translation_z) / m_cell_size))
	};
	const CellKey max_key{
		static_cast<int>(std::floor((max_x - m_translation_x - translation_x) / m_cell_size)),
		static_cast<int>(std::floor((max_y - m_translation_y - translation_y) / m_cell_size)),
		static_cast<int>(std::floor((max_z - m_translation_z - translation_z) / m_cell_size))
	};
	const std::size_t key_count =
		static_cast<std::size_t>(max_key.x - min_key.x + 1) *
		static_cast<std::size_t>(max_key.y - min_key.y + 1) *
		static_cast<std::size_t>(max_key.z - min_key.z + 1);
	if (key_count >= m_all_indices.size()) {
		for (unsigned int index : m_all_indices) {
			const CellKey& key = m_cell_keys[index];
			if (key.x >= min_key.x && key.x <= max_key.x &&
				key.y >= min_key.y && key.y <= max_key.y &&
				key.z >= min_key.z && key.z <= max_key.z) {
				candidates.push_back(index);
			}
		}
		return;
	}

	for (int z = min_key.z; z <= max_key.z; ++z) {
		for (int x = min_key.x; x <= max_key.x; ++x) {
			for (int y = min_key.y; y <= max_key.y; ++y) {
				const auto it = m_cells.find(CellKey{x, y, z});
				if (it == m_cells.end()) {
					continue;
				}
				candidates.insert(candidates.end(), it->second.begin(), it->second.end());
			}
		}
	}

	if (sort_candidates) {
		std::sort(candidates.begin(), candidates.end());
	}
}

std::vector<unsigned int> IonpSpatialIndex::point_candidates(const Point3D& point, double search_radius) const {
	std::vector<unsigned int> candidates;
	point_candidates(point, search_radius, candidates);
	return candidates;
}

void IonpSpatialIndex::point_candidates(const Point3D& point,
										double search_radius,
										std::vector<unsigned int>& candidates,
										double translation_x,
										double translation_y,
										double translation_z) const {
	candidates.clear();
	const double radius = std::max(0.0, search_radius);
	append_candidates_in_box(point.x() - radius,
							 point.x() + radius,
							 point.y() - radius,
							 point.y() + radius,
							 point.z() - radius,
							 point.z() + radius,
							 candidates,
							 translation_x,
							 translation_y,
							 translation_z,
							 true);
}

std::vector<unsigned int> IonpSpatialIndex::segment_candidates(const Point3D& start, const Point3D& end, double search_radius) const {
	std::vector<unsigned int> candidates;
	segment_candidates(start, end, search_radius, candidates);
	return candidates;
}

void IonpSpatialIndex::segment_candidates(const Point3D& start,
										  const Point3D& end,
										  double search_radius,
										  std::vector<unsigned int>& candidates,
										  double translation_x,
										  double translation_y,
										  double translation_z,
										  bool sort_candidates) const {
	candidates.clear();
	const double radius = std::max(0.0, search_radius);
	append_candidates_in_box(std::min(start.x(), end.x()) - radius,
							 std::max(start.x(), end.x()) + radius,
							 std::min(start.y(), end.y()) - radius,
							 std::max(start.y(), end.y()) + radius,
							 std::min(start.z(), end.z()) - radius,
							 std::max(start.z(), end.z()) + radius,
							 candidates,
							 translation_x,
							 translation_y,
							 translation_z,
							 sort_candidates);
}

std::vector<unsigned int> IonpSpatialIndex::variable_step_candidates(const Point3D& point) const {
	std::vector<unsigned int> candidates;
	variable_step_candidates(point, candidates);
	return candidates;
}

void IonpSpatialIndex::variable_step_candidates(
	const Point3D& point,
	std::vector<unsigned int>& candidates,
	double translation_x,
	double translation_y,
	double translation_z,
	bool sort_candidates) const {
	candidates.clear();
	if (empty() || m_max_outer_radius <= 0.0) {
		return;
	}

	const double search_radius =
		SimConsts::FIELD_RANGE_FACTOR * m_max_outer_radius;
	const CellKey min_key{
		static_cast<int>(std::floor((point.x() - search_radius - m_translation_x - translation_x) / m_variable_step_cell_size)),
		static_cast<int>(std::floor((point.y() - search_radius - m_translation_y - translation_y) / m_variable_step_cell_size)),
		static_cast<int>(std::floor((point.z() - search_radius - m_translation_z - translation_z) / m_variable_step_cell_size))
	};
	const CellKey max_key{
		static_cast<int>(std::floor((point.x() + search_radius - m_translation_x - translation_x) / m_variable_step_cell_size)),
		static_cast<int>(std::floor((point.y() + search_radius - m_translation_y - translation_y) / m_variable_step_cell_size)),
		static_cast<int>(std::floor((point.z() + search_radius - m_translation_z - translation_z) / m_variable_step_cell_size))
	};

	const std::size_t key_count =
		static_cast<std::size_t>(max_key.x - min_key.x + 1) *
		static_cast<std::size_t>(max_key.y - min_key.y + 1) *
		static_cast<std::size_t>(max_key.z - min_key.z + 1);
	if (key_count >= m_all_indices.size()) {
		candidates.assign(m_all_indices.begin(), m_all_indices.end());
		return;
	}

	for (int z = min_key.z; z <= max_key.z; ++z) {
		for (int x = min_key.x; x <= max_key.x; ++x) {
			for (int y = min_key.y; y <= max_key.y; ++y) {
				const auto it = m_variable_step_cells.find(CellKey{x, y, z});
				if (it != m_variable_step_cells.end()) {
					candidates.insert(candidates.end(), it->second.begin(), it->second.end());
				}
			}
		}
	}
	if (sort_candidates) {
		std::sort(candidates.begin(), candidates.end());
	}
}

const std::vector<unsigned int>& point_ionp_candidates(const IonpSpatialIndex& ionp_index,
													   IonpCheckMode ionp_check_mode,
													   const Point3D& point,
													   double search_radius,
													   std::vector<unsigned int>& spatial_candidates,
													   const Point3D* ionp_translation = nullptr) {
	if (ionp_check_mode == IonpCheckMode::All) {
		return ionp_index.all_indices();
	}
	if (ionp_translation) {
		ionp_index.point_candidates(point,
									search_radius,
									spatial_candidates,
									ionp_translation->x(),
									ionp_translation->y(),
									ionp_translation->z());
	}
	else {
		ionp_index.point_candidates(point, search_radius, spatial_candidates);
	}
	return spatial_candidates;
}

const std::vector<unsigned int>& segment_ionp_candidates(const IonpSpatialIndex& ionp_index,
														 IonpCheckMode ionp_check_mode,
														 bool force_all_ionp_candidates,
														 bool allow_unordered_spatial_candidates,
														 const Point3D& start,
														 const Point3D& end,
														 double search_radius,
														 std::vector<unsigned int>& spatial_candidates,
														 const Point3D* ionp_translation = nullptr) {
	if (force_all_ionp_candidates || ionp_check_mode == IonpCheckMode::All) {
		return ionp_index.all_indices();
	}
	if (ionp_translation) {
		ionp_index.segment_candidates(start,
									  end,
									  search_radius,
									  spatial_candidates,
									  ionp_translation->x(),
									  ionp_translation->y(),
									  ionp_translation->z(),
									  !allow_unordered_spatial_candidates);
	}
	else {
		ionp_index.segment_candidates(start,
									  end,
									  search_radius,
									  spatial_candidates,
									  0.0,
									  0.0,
									  0.0,
									  !allow_unordered_spatial_candidates);
	}
	return spatial_candidates;
}

Point3D translated_ionp_center(const Ionp& ionp, const Point3D* ionp_translation) {
	const Point3D center = ionp.xt();
	if (!ionp_translation) {
		return center;
	}
	return Point3D(center.x() + ionp_translation->x(),
				   center.y() + ionp_translation->y(),
				   center.z() + ionp_translation->z());
}

static double squared_distance_to_center(const Point3D& point, const Point3D& center) {
	const double dx = point.x() - center.x();
	const double dy = point.y() - center.y();
	const double dz = point.z() - center.z();
	return dx * dx + dy * dy + dz * dz;
}

bool should_use_variable_time_steps(bool var_steps,
									bool importance_sampling_enabled,
									ProtonSampleRegion sample_region) {
	if (!var_steps) {
		return false;
	}
	(void) importance_sampling_enabled;
	return sample_region != ProtonSampleRegion::BoundIntracellular &&
		   sample_region != ProtonSampleRegion::BoundExtracellular;
}

bool should_use_field_grid_for_step(bool normal_grid_selection,
									bool variable_step_active) {
	return normal_grid_selection && !variable_step_active;
}

// Output cadence is distinct from integration dt. Keep the historical
// short/long echo-time branches aligned with record scheduling below.
double phase_record_interval(double time_step, double echo_time) {
	if (time_step > SimConsts::MIN_TIME_STEP && echo_time <= 0.1) {
		return time_step;
	}
	return SimConsts::MIN_TIME_STEP;
}

bool reached_record_time(double time, double next_record_time) {
	const double tolerance = std::max(1e-12, 1e-9 * std::max(1.0, std::max(std::abs(time), std::abs(next_record_time))));
	return time + tolerance >= next_record_time;
}

double diffusion_coefficient_for_position(const Point3D& position,
										  const std::vector<Ionp>& ionpInstances,
										  const IonpSpatialIndex& ionp_index,
										  IonpCheckMode ionp_check_mode,
										  double coating_diffusion_coefficient,
										  bool in_cell,
										  double cell_diffusion_coefficient) {
	if (!ionp_index.has_coatings()) {
		return in_cell ? cell_diffusion_coefficient : D_D_CONST;
	}
	thread_local std::vector<unsigned int> spatial_candidates;
	const std::vector<unsigned int>& candidate_indices =
		point_ionp_candidates(ionp_index, ionp_check_mode, position, ionp_index.max_outer_radius(), spatial_candidates);
	for (unsigned int idx : candidate_indices) {
		const Ionp& ionp = ionpInstances[idx];
		if (ionp.coating() <= 0.) {
			continue;
		}
		const double radius_squared =
			squared_distance_to_center(position, ionp.xt());
		const double core_radius = ionp.radius();
		const double outer_radius = ionp.outer_radius();
		if (radius_squared >= core_radius * core_radius &&
			radius_squared < outer_radius * outer_radius) {
			return coating_diffusion_coefficient;
		}
	}
	return in_cell ? cell_diffusion_coefficient : D_D_CONST;
}

// Locate the coating containing this point; the cached index lets subsequent
// steps choose material diffusion without rescanning all particles.
int coating_index_for_position_translated(const Point3D& position,
										  const std::vector<Ionp>& ionpInstances,
										  const IonpSpatialIndex& ionp_index,
										  IonpCheckMode ionp_check_mode,
										  const Point3D* ionp_translation) {
	if (!ionp_index.has_coatings()) {
		return -1;
	}
	thread_local std::vector<unsigned int> spatial_candidates;
	const std::vector<unsigned int>& candidate_indices =
		point_ionp_candidates(ionp_index,
							  ionp_check_mode,
							  position,
							  ionp_index.max_outer_radius(),
							  spatial_candidates,
							  ionp_translation);
	for (unsigned int idx : candidate_indices) {
		const Ionp& ionp = ionpInstances[idx];
		if (ionp.coating() <= 0.) {
			continue;
		}
		const double radius_squared = squared_distance_to_center(
			position,
			translated_ionp_center(ionp, ionp_translation));
		const double core_radius = ionp.radius();
		const double outer_radius = ionp.outer_radius();
		if (radius_squared >= core_radius * core_radius &&
			radius_squared < outer_radius * outer_radius) {
			return static_cast<int>(idx);
		}
	}
	return -1;
}

int coating_index_for_position(const Point3D& position,
							   const std::vector<Ionp>& ionpInstances,
							   const IonpSpatialIndex& ionp_index,
							   IonpCheckMode ionp_check_mode) {
	return coating_index_for_position_translated(position,
												 ionpInstances,
												 ionp_index,
												 ionp_check_mode,
												 nullptr);
}

int coating_index_for_position_candidates(const Point3D& position,
										  const std::vector<Ionp>& ionpInstances,
										  const std::vector<unsigned int>& candidate_indices,
										  const Point3D* ionp_translation = nullptr) {
	for (unsigned int idx : candidate_indices) {
		const Ionp& ionp = ionpInstances[idx];
		if (ionp.coating() <= 0.) {
			continue;
		}
		const double radius_squared = squared_distance_to_center(
			position,
			translated_ionp_center(ionp, ionp_translation));
		const double core_radius = ionp.radius();
		const double outer_radius = ionp.outer_radius();
		if (radius_squared >= core_radius * core_radius &&
			radius_squared < outer_radius * outer_radius) {
			return static_cast<int>(idx);
		}
	}
	return -1;
}

double diffusion_coefficient_from_coating_state(const Point3D& position,
												const std::vector<Ionp>& ionpInstances,
												const IonpSpatialIndex& ionp_index,
												IonpCheckMode ionp_check_mode,
												int& coating_ionp_index,
												bool& coating_state_known,
												double coating_diffusion_coefficient,
												bool in_cell,
												double cell_diffusion_coefficient,
												const Point3D* ionp_translation = nullptr) {
	if (!ionp_index.has_coatings()) {
		coating_ionp_index = -1;
		coating_state_known = true;
		return in_cell ? cell_diffusion_coefficient : D_D_CONST;
	}
	if (!coating_state_known) {
		coating_ionp_index = coating_index_for_position_translated(position,
																 ionpInstances,
																 ionp_index,
																 ionp_check_mode,
																 ionp_translation);
		coating_state_known = true;
	}
	if (coating_ionp_index >= 0) {
		return coating_diffusion_coefficient;
	}
	return in_cell ? cell_diffusion_coefficient : D_D_CONST;
}

bool has_near_ionp_for_variable_step_translated(const Point3D& position,
												const std::vector<Ionp>& ionpInstances,
												const IonpSpatialIndex& ionp_index,
												IonpCheckMode ionp_check_mode,
												const Point3D* ionp_translation) {
	if (ionp_index.empty()) {
		return false;
	}
	thread_local std::vector<unsigned int> spatial_candidates;
	const std::vector<unsigned int>* candidate_indices = &ionp_index.all_indices();
	if (ionp_check_mode == IonpCheckMode::Spatial) {
		if (ionp_translation) {
			ionp_index.variable_step_candidates(position,
												spatial_candidates,
												ionp_translation->x(),
												ionp_translation->y(),
												ionp_translation->z(),
												false);
		}
		else {
			ionp_index.variable_step_candidates(
				position,
				spatial_candidates,
				0.0,
				0.0,
				0.0,
				false);
		}
		candidate_indices = &spatial_candidates;
	}
	for (unsigned int idx : *candidate_indices) {
		const Ionp& ionp = ionpInstances[idx];
		const double range =
			SimConsts::FIELD_RANGE_FACTOR * ionp.outer_radius();
		const Point3D center = translated_ionp_center(ionp, ionp_translation);
		const double dx = position.x() - center.x();
		const double dy = position.y() - center.y();
		const double dz = position.z() - center.z();
		if (dx * dx + dy * dy + dz * dz < range * range) {
			return true;
		}
	}
	return false;
}

bool has_near_ionp_for_variable_step(const Point3D& position,
									 const std::vector<Ionp>& ionpInstances,
									 const IonpSpatialIndex& ionp_index,
									 IonpCheckMode ionp_check_mode) {
	return has_near_ionp_for_variable_step_translated(position,
													  ionpInstances,
													  ionp_index,
													  ionp_check_mode,
													  nullptr);
}

bool has_ionp_within_grid_influence(
	const Point3D& position,
	const std::vector<Ionp>& ionpInstances,
	const IonpSpatialIndex& ionp_index,
	IonpCheckMode ionp_check_mode,
	const Point3D* ionp_translation) {
	if (ionp_index.empty()) {
		return false;
	}
	const double search_radius =
		SimConsts::FIELD_RANGE_FACTOR * ionp_index.max_outer_radius();
	thread_local std::vector<unsigned int> spatial_candidates;
	const std::vector<unsigned int>& candidate_indices =
		point_ionp_candidates(ionp_index,
							  ionp_check_mode,
							  position,
							  search_radius,
							  spatial_candidates,
							  ionp_translation);
	for (unsigned int index : candidate_indices) {
		const Ionp& ionp = ionpInstances[index];
		const double influence_radius =
			SimConsts::FIELD_RANGE_FACTOR *
			ionp.outer_radius();
		const Point3D center =
			translated_ionp_center(ionp, ionp_translation);
		const double dx = position.x() - center.x();
		const double dy = position.y() - center.y();
		const double dz = position.z() - center.z();
		if (dx * dx + dy * dy + dz * dz <
			influence_radius * influence_radius) {
			return true;
		}
	}
	return false;
}

namespace {

enum class IonpBoundaryKind {
	Core,
	Coating
};

struct IonpBoundaryHit {
	bool found = false;
	double path_distance = std::numeric_limits<double>::infinity();
	unsigned int ionp_index = 0;
	IonpBoundaryKind kind = IonpBoundaryKind::Core;
	Point3D center;
	double radius = 0.0;
	Point3D point;
};

// Select the earliest segment event. Near-equal hits use a deterministic
// particle-index/type tie break, preserving boundary decisions and random draws.
bool ionp_boundary_precedes(const IonpBoundaryHit& candidate,
							const IonpBoundaryHit& current,
							double segment_length) {
	if (!current.found) {
		return true;
	}
	const double tolerance =
		std::max(1e-15, 1e-12 * std::max(segment_length, std::max(candidate.radius, current.radius)));
	if (candidate.path_distance < current.path_distance - tolerance) {
		return true;
	}
	if (std::abs(candidate.path_distance - current.path_distance) > tolerance) {
		return false;
	}
	if (candidate.ionp_index != current.ionp_index) {
		return candidate.ionp_index < current.ionp_index;
	}
	return candidate.kind == IonpBoundaryKind::Core && current.kind == IonpBoundaryKind::Coating;
}

// Solve segment/sphere intersections for entry and exit events. Roots must
// lie on the travelled segment; tangent contacts are not membrane crossings.
void consider_sphere_boundaries(const Point3D& start,
								const Point3D& end,
								double direction_x,
								double direction_y,
								double direction_z,
								double segment_length,
								const Point3D& center,
								double radius,
								unsigned int ionp_index,
								IonpBoundaryKind kind,
								double minimum_path_distance,
								IonpBoundaryHit& best) {
	if (radius <= 0.0 || segment_length <= 0.0) {
		return;
	}

	const double offset_x = start.x() - center.x();
	const double offset_y = start.y() - center.y();
	const double offset_z = start.z() - center.z();
	const double projection =
		offset_x * direction_x + offset_y * direction_y + offset_z * direction_z;
	const double radius_squared = radius * radius;
	const double start_delta =
		offset_x * offset_x + offset_y * offset_y + offset_z * offset_z - radius_squared;

	const double end_offset_x = end.x() - center.x();
	const double end_offset_y = end.y() - center.y();
	const double end_offset_z = end.z() - center.z();
	const double end_delta =
		end_offset_x * end_offset_x + end_offset_y * end_offset_y +
		end_offset_z * end_offset_z - radius_squared;
	const bool start_inside = start_delta < 0.0;
	const bool end_inside = end_delta < 0.0;

	double discriminant = projection * projection - start_delta;
	const double tangent_tolerance = std::max(1e-15, 1e-12 * radius);
	if (discriminant <= tangent_tolerance * tangent_tolerance) {
		return;
	}
	const double root = std::sqrt(discriminant);
	const double roots[2] = {-projection - root, -projection + root};
	const double path_tolerance =
		std::max(1e-15, 1e-12 * std::max(radius, segment_length));

	for (double candidate_distance : roots) {
		if (candidate_distance < -path_tolerance ||
			candidate_distance > segment_length + path_tolerance) {
			continue;
		}
		candidate_distance = std::max(0.0, std::min(candidate_distance, segment_length));
		if (candidate_distance < minimum_path_distance) {
			continue;
		}

		const double normal_projection = projection + candidate_distance;
		if (std::abs(normal_projection) <= tangent_tolerance) {
			continue;
		}
		const bool entering = normal_projection < 0.0;

		bool is_crossing = true;
		if (candidate_distance <= path_tolerance) {
			is_crossing = (start_inside && !entering) || (!start_inside && entering);
		}
		else if (candidate_distance >= segment_length - path_tolerance) {
			is_crossing = (end_inside && entering) || (!end_inside && !entering);
		}
		if (!is_crossing) {
			continue;
		}

		IonpBoundaryHit candidate;
		candidate.found = true;
		candidate.path_distance = candidate_distance;
		candidate.ionp_index = ionp_index;
		candidate.kind = kind;
		candidate.center = center;
		candidate.radius = radius;
		candidate.point.set(start.x() + direction_x * candidate_distance,
							start.y() + direction_y * candidate_distance,
							start.z() + direction_z * candidate_distance);
		if (ionp_boundary_precedes(candidate, best, segment_length)) {
			best = candidate;
		}
	}
}

IonpBoundaryHit find_next_ionp_boundary(
	const Point3D& start,
	const Point3D& end,
	double direction_x,
	double direction_y,
	double direction_z,
	double segment_length,
	const std::vector<Ionp>& ionp_instances,
	const IonpSpatialIndex& ionp_index,
	IonpCheckMode ionp_check_mode,
	bool check_coating_permeability,
	double coating_permeability,
	double minimum_path_distance,
	std::vector<unsigned int>& spatial_candidates,
	const Point3D* ionp_translation = nullptr) {
	IonpBoundaryHit next_hit;
	const std::vector<unsigned int>& candidate_indices =
			segment_ionp_candidates(ionp_index,
									ionp_check_mode,
									false,
									true,
									start,
									end,
									ionp_index.max_outer_radius(),
									spatial_candidates,
									ionp_translation);
	for (unsigned int ionp_index_value : candidate_indices) {
		const Ionp& ionp = ionp_instances[ionp_index_value];
		const Point3D center = translated_ionp_center(ionp, ionp_translation);
		consider_sphere_boundaries(start,
								 end,
								 direction_x,
								 direction_y,
								 direction_z,
								 segment_length,
								 center,
								 ionp.radius(),
								 ionp_index_value,
								 IonpBoundaryKind::Core,
								 minimum_path_distance,
								 next_hit);
		if (check_coating_permeability &&
			coating_permeability < 1.0 &&
			ionp.coating() > 0.0) {
			consider_sphere_boundaries(start,
									 end,
									 direction_x,
									 direction_y,
									 direction_z,
									 segment_length,
									 center,
									 ionp.outer_radius(),
									 ionp_index_value,
									 IonpBoundaryKind::Coating,
									 minimum_path_distance,
									 next_hit);
		}
	}
	return next_hit;
}

// Reflect the remaining displacement about the sphere normal at the hit;
// core and coating events keep their separate physical rules.
bool reflect_at_ionp_boundary(const IonpBoundaryHit& hit,
							  double segment_length,
							  double direction_x,
							  double direction_y,
							  double direction_z,
							  Proton& proton,
							  Point3D& last_position,
							  double current_time,
							  std::vector<TrajectoryEvent>* event_log) {
	double normal_x = hit.point.x() - hit.center.x();
	double normal_y = hit.point.y() - hit.center.y();
	double normal_z = hit.point.z() - hit.center.z();
	const double normal_length =
		std::sqrt(normal_x * normal_x + normal_y * normal_y + normal_z * normal_z);
	if (normal_length == 0.0) {
		return false;
	}
	normal_x /= normal_length;
	normal_y /= normal_length;
	normal_z /= normal_length;

	const double projection =
		direction_x * normal_x + direction_y * normal_y + direction_z * normal_z;
	const double reflected_x = direction_x - 2.0 * projection * normal_x;
	const double reflected_y = direction_y - 2.0 * projection * normal_y;
	const double reflected_z = direction_z - 2.0 * projection * normal_z;
	const double remaining = std::max(0.0, segment_length - hit.path_distance);

	const double epsilon = std::max(1e-12, 1e-6 * hit.radius);
	const Point3D restart(hit.point.x() + reflected_x * epsilon,
						  hit.point.y() + reflected_y * epsilon,
						  hit.point.z() + reflected_z * epsilon);
	const double reflected_travel = remaining > epsilon ? remaining - epsilon : 0.0;
	proton.setXt(restart.x() + reflected_x * reflected_travel,
				 restart.y() + reflected_y * reflected_travel,
				 restart.z() + reflected_z * reflected_travel);
	last_position = restart;
	if (event_log) {
		event_log->push_back({current_time, hit.point, 1});
	}
	return true;
}

} // namespace

bool enforce_ionp_boundaries(Proton& proton,
							 Point3D& lastPosition,
							 const std::vector<Ionp>& ionpInstances,
							 const IonpSpatialIndex& ionp_index,
							 IonpCheckMode ionp_check_mode,
							 bool force_all_ionp_candidates,
							 int& coating_ionp_index,
							 bool& coating_state_known,
							 double current_time,
							 bool check_coating_permeability,
							 double coating_permeability,
							 RandomList& random_unit,
							 std::vector<TrajectoryEvent>* event_log) {
	std::size_t processed_events = 0;
	while (processed_events < 4096) {
		const Point3D start = lastPosition;
		const Point3D end = proton.xt();
		const double delta_x = end.x() - start.x();
		const double delta_y = end.y() - start.y();
		const double delta_z = end.z() - start.z();
		const double segment_length =
			std::sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);
		const double direction_x = segment_length > 0.0 ? delta_x / segment_length : 0.0;
		const double direction_y = segment_length > 0.0 ? delta_y / segment_length : 0.0;
		const double direction_z = segment_length > 0.0 ? delta_z / segment_length : 0.0;
		thread_local std::vector<unsigned int> spatial_candidates;
		const std::vector<unsigned int>& candidate_indices =
				segment_ionp_candidates(ionp_index,
										ionp_check_mode,
										force_all_ionp_candidates,
										false,
										start,
									end,
									ionp_index.max_outer_radius(),
									spatial_candidates);

		double minimum_path_distance = -std::numeric_limits<double>::infinity();
		while (processed_events < 4096) {
			IonpBoundaryHit next_hit;
			for (unsigned int k : candidate_indices) {
				const Ionp& ionp = ionpInstances[k];
				const Point3D center = ionp.xt();
				consider_sphere_boundaries(start,
										 end,
										 direction_x,
										 direction_y,
										 direction_z,
										 segment_length,
										 center,
										 ionp.radius(),
										 k,
										 IonpBoundaryKind::Core,
										 minimum_path_distance,
										 next_hit);

				if (check_coating_permeability &&
					coating_permeability < 1.0 &&
					ionp.coating() > 0.0) {
					consider_sphere_boundaries(start,
											 end,
											 direction_x,
											 direction_y,
											 direction_z,
											 segment_length,
											 center,
											 ionp.outer_radius(),
											 k,
											 IonpBoundaryKind::Coating,
											 minimum_path_distance,
											 next_hit);
				}
			}

			if (!next_hit.found) {
				for (unsigned int k : candidate_indices) {
					const Ionp& ionp = ionpInstances[k];
					if (squared_distance_to_center(end, ionp.xt()) <
						ionp.radius() * ionp.radius()) {
						coating_state_known = false;
						return true;
					}
				}
				coating_ionp_index =
					coating_index_for_position_candidates(end, ionpInstances, candidate_indices);
				coating_state_known = true;
				return false;
			}

			++processed_events;
			if (next_hit.kind == IonpBoundaryKind::Coating) {
				const bool crossing_allowed =
					random_unit.getNumber() < coating_permeability;
				if (crossing_allowed) {
					minimum_path_distance =
						std::nextafter(next_hit.path_distance,
									   std::numeric_limits<double>::infinity());
					continue;
				}
			}

			if (!reflect_at_ionp_boundary(next_hit,
										 segment_length,
										 direction_x,
										 direction_y,
										 direction_z,
										 proton,
										 lastPosition,
										 current_time,
										 event_log)) {
				coating_state_known = false;
				return true;
			}
			break;
		}
	}
	throw std::runtime_error("Exceeded IONP boundary-event limit in one diffusion step");
}

// Compute a specular reflection for a step that crosses a spherical boundary.
// last_pos is updated to the intersection point so multiple bounces in one step work.
bool reflect_from_sphere(
    const Point3D& sphere_center,
    double sphere_rad,
    Proton& proton,
    Point3D& last_pos,
    double current_time,
    std::vector<TrajectoryEvent>* event_log) {
	const Point3D start = last_pos;
	const Point3D end = proton.xt();
	const double delta_x = end.x() - start.x();
	const double delta_y = end.y() - start.y();
	const double delta_z = end.z() - start.z();
	const double segment_length =
		std::sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);
	if (segment_length == 0.0) {
		return false;
	}
	const double direction_x = delta_x / segment_length;
	const double direction_y = delta_y / segment_length;
	const double direction_z = delta_z / segment_length;
	IonpBoundaryHit hit;
	consider_sphere_boundaries(start,
							 end,
							 direction_x,
							 direction_y,
							 direction_z,
							 segment_length,
							 sphere_center,
							 sphere_rad,
							 0,
							 IonpBoundaryKind::Core,
							 -std::numeric_limits<double>::infinity(),
							 hit);
	if (!hit.found) {
		return false;
	}
	return reflect_at_ionp_boundary(hit,
								   segment_length,
								   direction_x,
								   direction_y,
								   direction_z,
								   proton,
								   last_pos,
								   current_time,
								   event_log);
}

// Gaussian phase-noise scale for the requested material T2 and integration dt.
double calc_std_T2( double time_step, double back_T2 ){
	return (0.14218444427715876 * pow((time_step/0.00001), 0.5014107086422922)) * pow(back_T2*1e3, -0.5014107088655116);
}

std::size_t available_cell_count(const std::vector<Point3D>& cell_origins, double cell_r) {
	if (cell_r <= 0.) {
		return 0;
	}
	return cell_origins.empty() ? 1u : cell_origins.size();
}

Point3D cell_center_at(const std::vector<Point3D>& cell_origins, const Point3D& cell_origin, std::size_t index) {
	if (cell_origins.empty()) {
		return cell_origin;
	}
	return cell_origins[index];
}

namespace {

struct UnitDirection {
	double x;
	double y;
	double z;
};

struct CellBoundaryHit {
	bool found = false;
	double path_distance = 0.;
	std::size_t cell_index = 0;
	Point3D image_center;
	bool entering = false;
};

struct VoxelFaceHit {
	bool found = false;
	double path_distance = 0.;
	bool wrap_x = false;
	bool wrap_y = false;
	bool wrap_z = false;
};

struct CellWalkResult {
	bool reflected = false;
	bool wrapped = false;
	bool retry = false;
};

struct IonpBoundaryContext {
	const std::vector<Ionp>& ionp_instances;
	const IonpSpatialIndex& ionp_index;
	IonpCheckMode ionp_check_mode;
	int& coating_ionp_index;
	bool& coating_state_known;
	bool check_coating_permeability;
	double coating_permeability;
	const Point3D* ionp_translation;
};

struct CellImageRange {
	int x_min = 0;
	int x_max = 0;
	int y_min = 0;
	int y_max = 0;
	int z_min = 0;
	int z_max = 0;
};

constexpr std::size_t kCellBoundaryEventLimit = 100000;

Point3D advance_point(const Point3D& point, const UnitDirection& direction, double length) {
	return Point3D(point.x() + direction.x * length,
				   point.y() + direction.y * length,
				   point.z() + direction.z * length);
}

double cell_surface_epsilon(double radius) {
	return std::max(1e-12, 1e-6 * radius);
}

double path_event_tolerance(double radius, double remaining_length) {
	return std::max(1e-12, 1e-12 * std::max({1.0, radius, remaining_length}));
}

CellImageRange cell_image_range(const Point3D& center,
								double cell_radius,
								double voxel_x_length,
								double voxel_y_length,
								double voxel_z_length,
								bool periodic) {
	CellImageRange range;
	if (!periodic) {
		return range;
	}

	if (voxel_x_length > 0.) {
		range.x_min = center.x() + cell_radius >= voxel_x_length / 2. ? -1 : 0;
		range.x_max = center.x() - cell_radius <= -voxel_x_length / 2. ? 1 : 0;
	}
	if (voxel_y_length > 0.) {
		range.y_min = center.y() + cell_radius >= voxel_y_length / 2. ? -1 : 0;
		range.y_max = center.y() - cell_radius <= -voxel_y_length / 2. ? 1 : 0;
	}
	if (voxel_z_length > 0.) {
		range.z_min = center.z() + cell_radius >= voxel_z_length / 2. ? -1 : 0;
		range.z_max = center.z() - cell_radius <= -voxel_z_length / 2. ? 1 : 0;
	}
	return range;
}

} // namespace

std::size_t CellSpatialIndex::CellKeyHash::operator()(const CellKey& key) const {
	std::size_t hash = static_cast<std::size_t>(key.x);
	hash = hash * 73856093u ^ (static_cast<std::size_t>(key.y) * 19349663u);
	hash ^= static_cast<std::size_t>(key.z) * 83492791u;
	return hash;
}

// Create periodic image proxies for cell spheres near voxel faces. A proxy
// has its own centre but retains the physical cell's index.
CellSpatialIndex::CellSpatialIndex(const std::vector<Point3D>& cell_origins,
								   const Point3D& fallback_origin,
								   double cell_radius,
								   double voxel_x_length,
								   double voxel_y_length,
								   double voxel_z_length,
								   bool periodic)
	: m_cell_radius(cell_radius),
	  m_bin_size(cell_radius > 0.0 ? 2.0 * cell_radius : 1.0) {
	const std::size_t cell_count = available_cell_count(cell_origins, cell_radius);
	m_proxy_indices_by_cell.resize(cell_count);
	for (std::size_t cell_index = 0; cell_index < cell_count; ++cell_index) {
		const Point3D center = cell_center_at(cell_origins, fallback_origin, cell_index);
		const CellImageRange image_range =
			cell_image_range(center,
							 cell_radius,
							 voxel_x_length,
							 voxel_y_length,
							 voxel_z_length,
							 periodic);
		for (int image_x = image_range.x_min; image_x <= image_range.x_max; ++image_x) {
			for (int image_y = image_range.y_min; image_y <= image_range.y_max; ++image_y) {
				for (int image_z = image_range.z_min; image_z <= image_range.z_max; ++image_z) {
					const Point3D image_center(center.x() + image_x * voxel_x_length,
											 center.y() + image_y * voxel_y_length,
											 center.z() + image_z * voxel_z_length);
					const std::size_t proxy_index = m_proxies.size();
						m_proxies.push_back({
							image_center,
							cell_index,
							image_x == 0 && image_y == 0 && image_z == 0
						});
					m_all_proxy_indices.push_back(proxy_index);
					const CellKey key{
						static_cast<int>(std::floor(image_center.x() / m_bin_size)),
						static_cast<int>(std::floor(image_center.y() / m_bin_size)),
						static_cast<int>(std::floor(image_center.z() / m_bin_size))
					};
					if (m_disjoint_interiors) {
						for (int dz = -1; dz <= 1 && m_disjoint_interiors; ++dz) {
							for (int dx = -1; dx <= 1 && m_disjoint_interiors; ++dx) {
								for (int dy = -1; dy <= 1 && m_disjoint_interiors; ++dy) {
									const auto nearby = m_bins.find(
										CellKey{key.x + dx, key.y + dy, key.z + dz});
									if (nearby == m_bins.end()) {
										continue;
									}
										for (std::size_t other_proxy_index : nearby->second) {
											const Proxy& other = m_proxies[other_proxy_index];
											const double diameter = 2.0 * cell_radius;
											if (squared_distance_to_center(image_center, other.center) <
												diameter * diameter) {
											m_disjoint_interiors = false;
											break;
										}
									}
								}
							}
						}
					}
					m_bins[key].push_back(proxy_index);
					m_proxy_indices_by_cell[cell_index].push_back(proxy_index);
				}
			}
		}
	}
}

bool CellSpatialIndex::empty() const {
	return m_proxies.empty();
}

bool CellSpatialIndex::has_disjoint_interiors() const {
	return m_disjoint_interiors;
}

double CellSpatialIndex::cell_radius() const {
	return m_cell_radius;
}

const CellSpatialIndex::Proxy& CellSpatialIndex::proxy(std::size_t proxy_index) const {
	return m_proxies[proxy_index];
}

const std::vector<std::size_t>& CellSpatialIndex::all_proxy_indices() const {
	return m_all_proxy_indices;
}

const std::vector<std::size_t>& CellSpatialIndex::proxies_for_cell(
	std::size_t cell_index) const {
	static const std::vector<std::size_t> empty_indices;
	return cell_index < m_proxy_indices_by_cell.size()
		? m_proxy_indices_by_cell[cell_index]
		: empty_indices;
}

void CellSpatialIndex::append_candidates_in_box(
	double min_x,
	double max_x,
	double min_y,
	double max_y,
	double min_z,
	double max_z,
	std::vector<std::size_t>& candidates) const {
	candidates.clear();
	if (empty()) {
		return;
	}
	const CellKey min_key{
		static_cast<int>(std::floor(min_x / m_bin_size)),
		static_cast<int>(std::floor(min_y / m_bin_size)),
		static_cast<int>(std::floor(min_z / m_bin_size))
	};
	const CellKey max_key{
		static_cast<int>(std::floor(max_x / m_bin_size)),
		static_cast<int>(std::floor(max_y / m_bin_size)),
		static_cast<int>(std::floor(max_z / m_bin_size))
	};
	const std::size_t key_count =
		static_cast<std::size_t>(max_key.x - min_key.x + 1) *
		static_cast<std::size_t>(max_key.y - min_key.y + 1) *
		static_cast<std::size_t>(max_key.z - min_key.z + 1);
	if (key_count >= m_all_proxy_indices.size()) {
		candidates.assign(m_all_proxy_indices.begin(), m_all_proxy_indices.end());
		return;
	}
	for (int z = min_key.z; z <= max_key.z; ++z) {
		for (int x = min_key.x; x <= max_key.x; ++x) {
			for (int y = min_key.y; y <= max_key.y; ++y) {
				const auto it = m_bins.find(CellKey{x, y, z});
				if (it != m_bins.end()) {
					candidates.insert(candidates.end(), it->second.begin(), it->second.end());
				}
			}
		}
	}
	std::sort(candidates.begin(), candidates.end());
}

void CellSpatialIndex::segment_candidates(
	const Point3D& start,
	const Point3D& end,
	double expansion,
	std::vector<std::size_t>& candidates) const {
	const double radius = std::max(0.0, expansion);
	append_candidates_in_box(std::min(start.x(), end.x()) - radius,
							 std::max(start.x(), end.x()) + radius,
							 std::min(start.y(), end.y()) - radius,
							 std::max(start.y(), end.y()) + radius,
							 std::min(start.z(), end.z()) - radius,
							 std::max(start.z(), end.z()) + radius,
							 candidates);
}

void CellSpatialIndex::point_candidates(
	const Point3D& point,
	double expansion,
	std::vector<std::size_t>& candidates) const {
	const double radius = std::max(0.0, expansion);
	append_candidates_in_box(point.x() - radius,
							 point.x() + radius,
							 point.y() - radius,
							 point.y() + radius,
							 point.z() - radius,
							 point.z() + radius,
							 candidates);
}

namespace {

// Refresh inCell/Cellnum together from exact containment after placement
// or wrapping; a proxy index is never stored as a physical cell number.
void synchronize_cell_membership(Proton& proton,
								 const CellSpatialIndex& cell_index,
								 const Point3D& position) {
	if (cell_index.empty()) {
		proton.inCell(false);
		proton.Cellnum(0);
		return;
	}

	const double cell_radius = cell_index.cell_radius();
	const double tolerance = cell_surface_epsilon(cell_radius);
	thread_local std::vector<std::size_t> point_candidates;
	cell_index.point_candidates(position, cell_radius + tolerance, point_candidates);

	if (proton.inCell() && proton.Cellnum() > 0) {
		const std::size_t active_cell_index =
			static_cast<std::size_t>(proton.Cellnum() - 1);
		const double active_limit = cell_radius + tolerance;
		for (std::size_t proxy_index : point_candidates) {
			const CellSpatialIndex::Proxy& proxy = cell_index.proxy(proxy_index);
			if (proxy.cell_index == active_cell_index &&
				squared_distance_to_center(position, proxy.center) <=
					active_limit * active_limit) {
				return;
			}
		}
	}

	const double inside_limit = std::max(0.0, cell_radius - tolerance);
	for (std::size_t proxy_index : point_candidates) {
		const CellSpatialIndex::Proxy& proxy = cell_index.proxy(proxy_index);
		if (squared_distance_to_center(position, proxy.center) <
			inside_limit * inside_limit) {
			proton.inCell(true);
			proton.Cellnum(SC_I(proxy.cell_index + 1));
			return;
		}
	}

	proton.inCell(false);
	proton.Cellnum(0);
}

CellBoundaryHit find_next_cell_boundary(const Point3D& start,
										const UnitDirection& direction,
										double remaining_length,
										const CellSpatialIndex& spatial_index,
										int active_cell_number) {
	CellBoundaryHit best;
	if (spatial_index.empty() || remaining_length < 0.) {
		return best;
	}

	const double cell_radius = spatial_index.cell_radius();
	const double tolerance = path_event_tolerance(cell_radius, remaining_length);
	const double tangent_tolerance = std::max(1e-15, 1e-12 * cell_radius);
	const Point3D end = advance_point(start, direction, remaining_length);
	thread_local std::vector<std::size_t> proxy_candidates;
	bool use_active_cell_only = false;
	if (active_cell_number > 0 && spatial_index.has_disjoint_interiors()) {
		const std::size_t active_cell_index =
			static_cast<std::size_t>(active_cell_number - 1);
		const double active_limit = cell_radius + tolerance;
		const double endpoint_limit = std::max(0.0, cell_radius - tolerance);
		for (std::size_t proxy_index : spatial_index.proxies_for_cell(active_cell_index)) {
			const CellSpatialIndex::Proxy& active_proxy = spatial_index.proxy(proxy_index);
			if (squared_distance_to_center(start, active_proxy.center) <=
				active_limit * active_limit) {
				use_active_cell_only = true;
				if (squared_distance_to_center(end, active_proxy.center) <
					endpoint_limit * endpoint_limit) {
					return best;
				}
				break;
			}
		}
		if (use_active_cell_only) {
			const std::vector<std::size_t>& active_proxies =
				spatial_index.proxies_for_cell(active_cell_index);
			proxy_candidates.assign(active_proxies.begin(), active_proxies.end());
		}
	}
	if (!use_active_cell_only) {
		spatial_index.segment_candidates(start,
										 end,
										 cell_radius + tolerance,
										 proxy_candidates);
	}

	for (std::size_t proxy_index : proxy_candidates) {
		const CellSpatialIndex::Proxy& proxy = spatial_index.proxy(proxy_index);
		const Point3D& image_center = proxy.center;
		const std::size_t cell_index = proxy.cell_index;
		const double ox = start.x() - image_center.x();
		const double oy = start.y() - image_center.y();
		const double oz = start.z() - image_center.z();
		const double projection = ox * direction.x + oy * direction.y + oz * direction.z;
		const double radial_delta = ox * ox + oy * oy + oz * oz - cell_radius * cell_radius;
		const double discriminant = projection * projection - radial_delta;

		// A zero discriminant is a tangential contact, not a membrane crossing.
		if (discriminant <= tangent_tolerance * tangent_tolerance) {
			continue;
		}

		const double root = std::sqrt(discriminant);
		const double roots[2] = {-projection - root, -projection + root};
		for (double candidate : roots) {
			if (candidate < -tolerance || candidate > remaining_length + tolerance) {
				continue;
			}
			candidate = std::max(0., std::min(candidate, remaining_length));
			const double normal_projection = projection + candidate;
			if (std::abs(normal_projection) <= tangent_tolerance) {
				continue;
			}
			const bool entering = normal_projection < 0.;
			const bool currently_inside_this_cell = active_cell_number == SC_I(cell_index + 1);
			if (candidate <= tolerance && entering == currently_inside_this_cell) {
				// Starting on a membrane and moving farther into the current
				// compartment is not another crossing of that membrane.
				continue;
			}

			const bool candidate_is_active_exit = !entering && currently_inside_this_cell;
			const bool best_is_active_exit = best.found && !best.entering &&
				active_cell_number == SC_I(best.cell_index + 1);
			if (!best.found || candidate < best.path_distance - tolerance ||
				(std::abs(candidate - best.path_distance) <= tolerance &&
				 candidate_is_active_exit && !best_is_active_exit)) {
				best.found = true;
				best.path_distance = candidate;
				best.cell_index = cell_index;
				best.image_center = image_center;
				best.entering = entering;
			}
		}
	}

	return best;
}

VoxelFaceHit find_next_voxel_face(const Point3D& start,
								  const UnitDirection& direction,
								  double remaining_length,
								  double voxel_x_length,
								  double voxel_y_length,
								  double voxel_z_length) {
	VoxelFaceHit hit;
	const double tolerance = path_event_tolerance(0., remaining_length);

	auto consider_axis = [&](double coordinate, double component, double voxel_length, bool& selected_axis) {
		if (voxel_length <= 0. || component == 0.) {
			return;
		}
		const double half_length = voxel_length / 2.;
		const double face = component > 0. ? half_length : -half_length;
		double candidate = (face - coordinate) / component;
		if (candidate < -tolerance || candidate > remaining_length + tolerance) {
			return;
		}
		candidate = std::max(0., std::min(candidate, remaining_length));
		if (!hit.found || candidate < hit.path_distance - tolerance) {
			hit.found = true;
			hit.path_distance = candidate;
			hit.wrap_x = false;
			hit.wrap_y = false;
			hit.wrap_z = false;
			selected_axis = true;
		}
		else if (std::abs(candidate - hit.path_distance) <= tolerance) {
			selected_axis = true;
		}
	};

	consider_axis(start.x(), direction.x, voxel_x_length, hit.wrap_x);
	consider_axis(start.y(), direction.y, voxel_y_length, hit.wrap_y);
	consider_axis(start.z(), direction.z, voxel_z_length, hit.wrap_z);
	return hit;
}

void set_cell_state_after_allowed_crossing(Proton& proton, const CellBoundaryHit& hit) {
	if (hit.entering) {
		proton.inCell(true);
		proton.Cellnum(SC_I(hit.cell_index + 1));
		return;
	}

	if (proton.inCell() && proton.Cellnum() == SC_I(hit.cell_index + 1)) {
		proton.inCell(false);
		proton.Cellnum(0);
	}
}

// Walk one proposed displacement in travelled-path order. Compete cell
// membranes, IONP core/coating spheres and periodic voxel faces for the next
// event. Reflection/wrapping preserves remaining path length; small geometric
// offsets avoid immediately rediscovering the same surface. If the event
// budget is exhausted, report retry so the caller can restore pre-step state.
CellWalkResult walk_cell_boundaries(Proton& proton,
									 Point3D& last_position,
									 const CellSpatialIndex& cell_index,
								 IonpBoundaryContext* ionp_context,
								 double cell_radius,
								 double permeability,
								 double voxel_x_length,
								 double voxel_y_length,
								 double voxel_z_length,
								 bool periodic,
								 RandomList& random_unit,
								 double current_time,
								 std::vector<TrajectoryEvent>* event_log) {
	CellWalkResult result;
	auto finalize_ionp_state = [&](const Point3D& endpoint) {
		if (!ionp_context) {
			return false;
		}
		thread_local std::vector<unsigned int> endpoint_candidates;
			const std::vector<unsigned int>& candidate_indices =
				point_ionp_candidates(ionp_context->ionp_index,
								  ionp_context->ionp_check_mode,
									  endpoint,
									  ionp_context->ionp_index.has_coatings()
										  ? ionp_context->ionp_index.max_outer_radius()
										  : ionp_context->ionp_index.max_core_radius(),
									  endpoint_candidates,
									  ionp_context->ionp_translation);
			for (unsigned int index : candidate_indices) {
				const Ionp& ionp = ionp_context->ionp_instances[index];
				if (squared_distance_to_center(
						endpoint,
						translated_ionp_center(ionp, ionp_context->ionp_translation)) <
					ionp.radius() * ionp.radius()) {
				ionp_context->coating_state_known = false;
					return true;
				}
			}
			if (!ionp_context->ionp_index.has_coatings()) {
				ionp_context->coating_ionp_index = -1;
				ionp_context->coating_state_known = true;
				return false;
			}
			ionp_context->coating_ionp_index =
				coating_index_for_position_candidates(endpoint,
													 ionp_context->ionp_instances,
													 candidate_indices,
													 ionp_context->ionp_translation);
		ionp_context->coating_state_known = true;
		return false;
	};
	Point3D requested_end = proton.xt();
	Point3D current = last_position;
	if (periodic) {
		double current_x = current.x();
		double current_y = current.y();
		double current_z = current.z();
		double endpoint_x = requested_end.x();
		double endpoint_y = requested_end.y();
		double endpoint_z = requested_end.z();
		auto canonicalize_axis = [&](double& coordinate, double& endpoint, double voxel_length) {
			if (voxel_length <= 0. ||
				(coordinate >= -voxel_length / 2. && coordinate <= voxel_length / 2.)) {
				return;
			}
			const double period_count = std::floor((coordinate + voxel_length / 2.) / voxel_length);
			coordinate -= period_count * voxel_length;
			endpoint -= period_count * voxel_length;
			result.wrapped = true;
		};
		canonicalize_axis(current_x, endpoint_x, voxel_x_length);
		canonicalize_axis(current_y, endpoint_y, voxel_y_length);
		canonicalize_axis(current_z, endpoint_z, voxel_z_length);
		current.set(current_x, current_y, current_z);
		requested_end.set(endpoint_x, endpoint_y, endpoint_z);
		if (result.wrapped) {
			last_position = current;
		}
	}
	const double full_length = distance(current, requested_end);

	synchronize_cell_membership(proton,
								cell_index,
								current);
	if (full_length == 0.) {
		proton.setXt(current);
		result.retry = finalize_ionp_state(current);
		return result;
	}

	UnitDirection direction{
		(requested_end.x() - current.x()) / full_length,
		(requested_end.y() - current.y()) / full_length,
		(requested_end.z() - current.z()) / full_length
	};
	double remaining_length = full_length;
	const bool skip_cell_boundary_events =
		permeability >= 1.0 && cell_index.has_disjoint_interiors();
	auto finalize_walk_endpoint = [&](const Point3D& endpoint) {
		proton.setXt(endpoint);
		if (skip_cell_boundary_events) {
			synchronize_cell_membership(proton, cell_index, endpoint);
		}
		result.retry = finalize_ionp_state(endpoint);
	};

	if (permeability >= 1.0 &&
		cell_index.has_disjoint_interiors() &&
		(!ionp_context || ionp_context->ionp_index.empty())) {
		Point3D fast_current = current;
		Point3D fast_last_position = last_position;
		double fast_remaining_length = remaining_length;
		bool fast_wrapped = result.wrapped;
		Point3D fast_endpoint;
		bool endpoint_found = false;

		for (std::size_t event_count = 0;
			 event_count < kCellBoundaryEventLimit;
			 ++event_count) {
			const VoxelFaceHit face_hit = periodic
				? find_next_voxel_face(fast_current,
									  direction,
									  fast_remaining_length,
									  voxel_x_length,
									  voxel_y_length,
									  voxel_z_length)
				: VoxelFaceHit{};
			if (!face_hit.found) {
				fast_endpoint =
					advance_point(fast_current, direction, fast_remaining_length);
				endpoint_found = true;
				break;
			}

			const Point3D face_point =
				advance_point(fast_current, direction, face_hit.path_distance);
			fast_remaining_length =
				std::max(0.0, fast_remaining_length - face_hit.path_distance);
			double wrapped_x = face_point.x();
			double wrapped_y = face_point.y();
			double wrapped_z = face_point.z();
			if (face_hit.wrap_x) {
				wrapped_x =
					direction.x > 0.0 ? -voxel_x_length / 2.0 : voxel_x_length / 2.0;
			}
			if (face_hit.wrap_y) {
				wrapped_y =
					direction.y > 0.0 ? -voxel_y_length / 2.0 : voxel_y_length / 2.0;
			}
			if (face_hit.wrap_z) {
				wrapped_z =
					direction.z > 0.0 ? -voxel_z_length / 2.0 : voxel_z_length / 2.0;
			}
			fast_current.set(wrapped_x, wrapped_y, wrapped_z);
			fast_last_position = fast_current;
			fast_wrapped = true;

			const double event_tolerance =
				path_event_tolerance(cell_radius, fast_remaining_length);
			if (fast_remaining_length <= event_tolerance) {
				fast_endpoint =
					advance_point(fast_current, direction, fast_remaining_length);
				endpoint_found = true;
				break;
			}
		}
		if (!endpoint_found) {
			throw std::runtime_error(
				"Exceeded the periodic-face event limit during one proton diffusion step");
		}

		const double surface_tolerance = cell_surface_epsilon(cell_radius);
		thread_local std::vector<std::size_t> endpoint_candidates;
		cell_index.point_candidates(fast_endpoint,
									cell_radius + surface_tolerance,
									endpoint_candidates);
		bool endpoint_on_membrane = false;
		int endpoint_cell_number = 0;
		for (std::size_t proxy_index : endpoint_candidates) {
			const CellSpatialIndex::Proxy& proxy = cell_index.proxy(proxy_index);
			const double center_distance = distance(fast_endpoint, proxy.center);
			if (std::abs(center_distance - cell_radius) <= surface_tolerance) {
				endpoint_on_membrane = true;
				break;
			}
			if (endpoint_cell_number == 0 &&
				center_distance < cell_radius - surface_tolerance) {
				endpoint_cell_number = SC_I(proxy.cell_index + 1);
			}
		}

		if (!endpoint_on_membrane) {
			proton.setXt(fast_endpoint);
			proton.inCell(endpoint_cell_number > 0);
			proton.Cellnum(endpoint_cell_number);
			last_position = fast_last_position;
			result.wrapped = fast_wrapped;
			result.retry = finalize_ionp_state(fast_endpoint);
			return result;
		}
	}

	double minimum_ionp_path_distance =
		-std::numeric_limits<double>::infinity();
	thread_local std::vector<unsigned int> ionp_spatial_candidates;
	for (std::size_t event_count = 0; event_count < kCellBoundaryEventLimit; ++event_count) {
		const CellBoundaryHit cell_hit =
			skip_cell_boundary_events
				? CellBoundaryHit{}
				: find_next_cell_boundary(current,
										 direction,
										 remaining_length,
										 cell_index,
										 proton.inCell() ? proton.Cellnum() : 0);
		const VoxelFaceHit face_hit = periodic
			? find_next_voxel_face(current, direction, remaining_length, voxel_x_length, voxel_y_length, voxel_z_length)
			: VoxelFaceHit{};
		const double event_tolerance = path_event_tolerance(cell_radius, remaining_length);
		IonpBoundaryHit ionp_hit;
		if (ionp_context && !ionp_context->ionp_index.empty()) {
			const Point3D segment_end =
				advance_point(current, direction, remaining_length);
			ionp_hit = find_next_ionp_boundary(current,
											 segment_end,
											 direction.x,
											 direction.y,
											 direction.z,
											 remaining_length,
											 ionp_context->ionp_instances,
											 ionp_context->ionp_index,
											 ionp_context->ionp_check_mode,
											 ionp_context->check_coating_permeability,
											 ionp_context->coating_permeability,
											 minimum_ionp_path_distance,
											 ionp_spatial_candidates,
											 ionp_context->ionp_translation);
		}

		if (!cell_hit.found && !face_hit.found && !ionp_hit.found) {
			const Point3D endpoint =
				advance_point(current, direction, remaining_length);
			finalize_walk_endpoint(endpoint);
			return result;
		}

		const bool ionp_event_is_first =
			ionp_hit.found &&
			(!cell_hit.found ||
			 ionp_hit.path_distance <= cell_hit.path_distance + event_tolerance) &&
			(!face_hit.found ||
			 ionp_hit.path_distance <= face_hit.path_distance + event_tolerance);
		if (ionp_event_is_first) {
			if (ionp_hit.kind == IonpBoundaryKind::Coating) {
				const bool crossing_allowed =
					random_unit.getNumber() < ionp_context->coating_permeability;
				if (crossing_allowed) {
					minimum_ionp_path_distance =
						std::nextafter(ionp_hit.path_distance,
									   std::numeric_limits<double>::infinity());
					continue;
				}
			}

			remaining_length =
				std::max(0.0, remaining_length - ionp_hit.path_distance);
			double normal_x = ionp_hit.point.x() - ionp_hit.center.x();
			double normal_y = ionp_hit.point.y() - ionp_hit.center.y();
			double normal_z = ionp_hit.point.z() - ionp_hit.center.z();
			const double normal_length =
				std::sqrt(normal_x * normal_x +
						  normal_y * normal_y +
						  normal_z * normal_z);
			if (normal_length == 0.0) {
				ionp_context->coating_state_known = false;
				result.retry = true;
				return result;
			}
			normal_x /= normal_length;
			normal_y /= normal_length;
			normal_z /= normal_length;
			const double projection =
				direction.x * normal_x +
				direction.y * normal_y +
				direction.z * normal_z;
			direction.x -= 2.0 * projection * normal_x;
			direction.y -= 2.0 * projection * normal_y;
			direction.z -= 2.0 * projection * normal_z;

			const double epsilon = std::max(1e-12, 1e-6 * ionp_hit.radius);
			current.set(ionp_hit.point.x() + direction.x * epsilon,
						ionp_hit.point.y() + direction.y * epsilon,
						ionp_hit.point.z() + direction.z * epsilon);
			remaining_length =
				remaining_length > epsilon ? remaining_length - epsilon : 0.0;
			last_position = current;
			result.reflected = true;
			minimum_ionp_path_distance =
				-std::numeric_limits<double>::infinity();
			if (event_log) {
				event_log->push_back({current_time, ionp_hit.point, 1});
			}
			if (remaining_length <= event_tolerance) {
				const Point3D endpoint =
					advance_point(current, direction, remaining_length);
				finalize_walk_endpoint(endpoint);
				return result;
			}
			continue;
		}

		if (cell_hit.found &&
			(!face_hit.found || cell_hit.path_distance <= face_hit.path_distance + event_tolerance)) {
			const Point3D hit_point = advance_point(current, direction, cell_hit.path_distance);
			remaining_length = std::max(0., remaining_length - cell_hit.path_distance);

			const bool crossing_allowed = permeability >= 1. ||
				(permeability > 0. && random_unit.getNumber() < permeability);
			if (crossing_allowed) {
				set_cell_state_after_allowed_crossing(proton, cell_hit);
			}
			else {
				result.reflected = true;
				double nx = hit_point.x() - cell_hit.image_center.x();
				double ny = hit_point.y() - cell_hit.image_center.y();
				double nz = hit_point.z() - cell_hit.image_center.z();
				const double normal_length = std::sqrt(nx * nx + ny * ny + nz * nz);
				if (normal_length == 0.) {
					throw std::runtime_error("Cannot reflect a proton at a degenerate cell boundary");
				}
				nx /= normal_length;
				ny /= normal_length;
				nz /= normal_length;
				const double projection = direction.x * nx + direction.y * ny + direction.z * nz;
				direction.x -= 2. * projection * nx;
				direction.y -= 2. * projection * ny;
				direction.z -= 2. * projection * nz;
				if (event_log) {
					event_log->push_back({current_time, hit_point, 1});
				}
			}

			// Stay on the event point.  The state-aware zero-root filter prevents
			// reprocessing this membrane without skipping a nearby cell or voxel face.
			current = hit_point;
			minimum_ionp_path_distance =
				-std::numeric_limits<double>::infinity();
			if (!crossing_allowed) {
				last_position = current;
			}
			if (remaining_length <= event_tolerance) {
				const Point3D endpoint =
					advance_point(current, direction, remaining_length);
				finalize_walk_endpoint(endpoint);
				return result;
			}
			continue;
		}

		const Point3D face_point = advance_point(current, direction, face_hit.path_distance);
		remaining_length = std::max(0., remaining_length - face_hit.path_distance);
		double wrapped_x = face_point.x();
		double wrapped_y = face_point.y();
		double wrapped_z = face_point.z();
		if (face_hit.wrap_x) {
			wrapped_x = direction.x > 0. ? -voxel_x_length / 2. : voxel_x_length / 2.;
		}
		if (face_hit.wrap_y) {
			wrapped_y = direction.y > 0. ? -voxel_y_length / 2. : voxel_y_length / 2.;
		}
		if (face_hit.wrap_z) {
			wrapped_z = direction.z > 0. ? -voxel_z_length / 2. : voxel_z_length / 2.;
		}
		current.set(wrapped_x, wrapped_y, wrapped_z);
		last_position = current;
		result.wrapped = true;
		minimum_ionp_path_distance =
			-std::numeric_limits<double>::infinity();
		if (remaining_length <= event_tolerance) {
			const Point3D endpoint =
				advance_point(current, direction, remaining_length);
			finalize_walk_endpoint(endpoint);
			return result;
		}
	}

	throw std::runtime_error("Exceeded the cell-boundary event limit during one proton diffusion step");
}

} // namespace

bool enforce_cell_boundaries(Proton& proton,
							 Point3D& last_position,
							 const std::vector<Point3D>& cell_origins,
							 const Point3D& fallback_origin,
							 double cell_radius,
							 double permeability,
							 RandomList& random_unit,
								 double current_time,
								 std::vector<TrajectoryEvent>* event_log) {
	const CellSpatialIndex cell_index(cell_origins,
									 fallback_origin,
									 cell_radius,
									 0.,
									 0.,
									 0.,
									 false);
	return walk_cell_boundaries(proton,
								last_position,
								cell_index,
								nullptr,
								cell_radius,
								permeability,
								0.,
								0.,
								0.,
								false,
								random_unit,
								current_time,
								event_log).reflected;
}

bool enforce_cell_boundaries_periodic(Proton& proton,
								  Point3D& last_position,
								  const std::vector<Point3D>& cell_origins,
								  const Point3D& fallback_origin,
								  double cell_radius,
								  double permeability,
								  double voxel_x_length,
								  double voxel_y_length,
								  double voxel_z_length,
								  RandomList& random_unit,
								  double current_time,
								  std::vector<TrajectoryEvent>* event_log) {
	const CellSpatialIndex cell_index(cell_origins,
									 fallback_origin,
									 cell_radius,
									 voxel_x_length,
									 voxel_y_length,
									 voxel_z_length,
									 true);
	return walk_cell_boundaries(proton,
								last_position,
								cell_index,
								nullptr,
								cell_radius,
								permeability,
								voxel_x_length,
								voxel_y_length,
								voxel_z_length,
								true,
								random_unit,
								current_time,
								event_log).reflected;
}

double importance_shell_volume(const Ionp& ionp,
							   IonpDistanceRegion region,
							   double near_cutoff,
							   double far_cutoff) {
	const double inner = ionp.outer_radius() +
		(region == IonpDistanceRegion::Intermediate ? near_cutoff : 0.);
	const double outer = ionp.outer_radius() +
		(region == IonpDistanceRegion::Near ? near_cutoff : far_cutoff);
	return (4. / 3.) * D_PI *
		(outer * outer * outer - inner * inner * inner);
}

void importance_shell_bounds(const Ionp& ionp,
							 IonpDistanceRegion region,
							 double near_cutoff,
							 double far_cutoff,
							 double& inner,
							 double& outer) {
	inner = ionp.outer_radius() +
		(region == IonpDistanceRegion::Intermediate ? near_cutoff : 0.);
	outer = ionp.outer_radius() +
		(region == IonpDistanceRegion::Near ? near_cutoff : far_cutoff);
}

std::size_t importance_shell_coverage(
	const Point3D& point,
	const std::vector<Ionp>& ionps,
	const IonpPopulationMap& populations,
	const std::vector<unsigned int>& candidates,
	IonpPopulation population,
	IonpDistanceRegion region,
	double near_cutoff,
	double far_cutoff) {
	std::size_t coverage = 0;
	for (unsigned int candidate : candidates) {
		const std::size_t index = static_cast<std::size_t>(candidate);
		if (populations.populations[index] != population) {
			continue;
		}
		double inner = 0.;
		double outer = 0.;
		importance_shell_bounds(
			ionps[index], region, near_cutoff, far_cutoff, inner, outer);
		const double radius = distance(point, ionps[index].xt());
		if (radius >= inner && radius <= outer) {
			++coverage;
		}
	}
	return coverage;
}

// Volume-weighted choice of source shells. Subsequent inverse-coverage
// acceptance makes the union uniform even where individual shells overlap.
struct ImportanceShellProposal {
	const std::vector<std::size_t>* indices = nullptr;
	IonpPopulation population = IonpPopulation::Extracellular;
	IonpDistanceRegion region = IonpDistanceRegion::Far;
	double near_cutoff = 0.;
	double far_cutoff = 0.;
	double total_volume = 0.;
	std::vector<double> cumulative_volumes;
};

ImportanceShellProposal make_importance_shell_proposal(
	const std::vector<Ionp>& ionps,
	const std::vector<std::size_t>& indices,
	IonpPopulation population,
	IonpDistanceRegion region,
	double near_cutoff,
	double far_cutoff) {
	ImportanceShellProposal proposal;
	if (region == IonpDistanceRegion::Far) {
		return proposal;
	}
	proposal.indices = &indices;
	proposal.population = population;
	proposal.region = region;
	proposal.near_cutoff = near_cutoff;
	proposal.far_cutoff = far_cutoff;
	proposal.cumulative_volumes.reserve(indices.size());
	for (std::size_t index : indices) {
		proposal.total_volume += importance_shell_volume(
			ionps[index], region, near_cutoff, far_cutoff);
		proposal.cumulative_volumes.push_back(proposal.total_volume);
	}
	return proposal;
}

void synchronize_initial_cell_membership(
	Proton& proton,
	const std::vector<Point3D>& cell_origins,
	const Point3D& cell_origin,
	double cell_r) {
	proton.inCell(false);
	proton.Cellnum(0);
	const std::size_t cell_count = available_cell_count(cell_origins, cell_r);
	for (std::size_t cell_index = 0; cell_index < cell_count; ++cell_index) {
		if (distance(
				proton.xt(),
				cell_center_at(cell_origins, cell_origin, cell_index)) <= cell_r) {
			proton.inCell(true);
			proton.Cellnum(SC_I(cell_index + 1));
			return;
		}
	}
}

// Draw a start in one exclusive legacy joint region, rejecting inaccessible
// points. The assigned region remains the estimator label during diffusion.
void initialize_proton_in_region(Proton& proton,
								 ProtonSampleRegion sample_region,
								 const ImportanceSamplingConfig& importance_sampling,
								 const std::vector<Point3D>& cell_origins,
								 const Point3D& cell_origin,
								 double cell_r,
								 const std::vector<Ionp>& ionps,
								 const IonpPopulationMap& populations,
								 const IonpSpatialIndex& ionp_index,
								 double sim_x_length,
								 double sim_y_length,
								 double sim_z_length,
								 RandomList& random_x,
								 RandomList& random_y,
								 RandomList& random_z,
								 RandomList& random_unit) {
	if (sample_region == ProtonSampleRegion::BoundIntracellular ||
		sample_region == ProtonSampleRegion::BoundExtracellular) {
		throw std::runtime_error("Bound proton sampling is not valid for free-proton initialization");
	}

	const JointSampleRegion target =
		decode_joint_sample_region(sample_region);
	const ImportanceShellProposal inside_proposal =
		make_importance_shell_proposal(
			ionps,
			populations.intracellular_indices,
			IonpPopulation::Intracellular,
			target.intracellular_ionp,
			importance_sampling.intracellular_ionp_near_cutoff,
			importance_sampling.intracellular_ionp_far_cutoff);
	const ImportanceShellProposal outside_proposal =
		make_importance_shell_proposal(
			ionps,
			populations.extracellular_indices,
			IonpPopulation::Extracellular,
			target.extracellular_ionp,
			importance_sampling.extracellular_ionp_near_cutoff,
			importance_sampling.extracellular_ionp_far_cutoff);

	const ImportanceShellProposal* shell_proposal = nullptr;
	if (inside_proposal.total_volume > 0.) {
		shell_proposal = &inside_proposal;
	}
	if (outside_proposal.total_volume > 0. &&
		(!shell_proposal ||
		 outside_proposal.total_volume < shell_proposal->total_volume)) {
		shell_proposal = &outside_proposal;
	}
	if ((target.intracellular_ionp != IonpDistanceRegion::Far ||
		 target.extracellular_ionp != IonpDistanceRegion::Far) &&
		!shell_proposal) {
		throw std::runtime_error(
			"No IONPs are available for requested region " +
			sample_region_name(sample_region));
	}

	const double half_x = sim_x_length / 2.;
	const double half_y = sim_y_length / 2.;
	const double half_z = sim_z_length / 2.;
	const double maximum_distance_cutoff = std::max(
		importance_sampling.intracellular_ionp_sampling_enabled()
			? importance_sampling.intracellular_ionp_far_cutoff
			: 0.,
		importance_sampling.extracellular_ionp_sampling_enabled()
			? importance_sampling.extracellular_ionp_far_cutoff
			: 0.);
	const double ionp_search_radius =
		ionp_index.max_outer_radius() + maximum_distance_cutoff;
	for (unsigned int attempt = 0; attempt < SimConsts::INIT_RETRY_LIMIT; ++attempt) {
		if (shell_proposal) {
			const double draw =
				random_unit.getNumber() * shell_proposal->total_volume;
			const auto chosen_position = std::lower_bound(
				shell_proposal->cumulative_volumes.begin(),
				shell_proposal->cumulative_volumes.end(),
				draw);
			const std::size_t chosen_offset =
				chosen_position == shell_proposal->cumulative_volumes.end()
					? shell_proposal->cumulative_volumes.size() - 1
					: static_cast<std::size_t>(
						chosen_position -
						shell_proposal->cumulative_volumes.begin());
			const std::size_t chosen =
				(*shell_proposal->indices)[chosen_offset];
			double inner = 0.;
			double outer = 0.;
			importance_shell_bounds(
				ionps[chosen],
				shell_proposal->region,
				shell_proposal->near_cutoff,
				shell_proposal->far_cutoff,
				inner,
				outer);
			place_proton_random_sphere(
				proton,
				ionps[chosen].xt().x(),
				ionps[chosen].xt().y(),
				ionps[chosen].xt().z(),
				inner,
				outer,
				random_unit);
		}
		else if (target.cell == CellSampleRegion::Inside) {
			const std::size_t cell_count =
				available_cell_count(cell_origins, cell_r);
			if (cell_count == 0) {
				throw std::runtime_error(
					"Inside-cell importance sampling requires at least one cell");
			}
			const std::size_t selected = std::min(
				static_cast<std::size_t>(
					random_unit.getNumber() *
					static_cast<double>(cell_count)),
				cell_count - 1);
			const Point3D center =
				cell_center_at(cell_origins, cell_origin, selected);
			place_proton_random_sphere(
				proton,
				center.x(),
				center.y(),
				center.z(),
				0.,
				cell_r,
				random_unit);
		}
		else {
			proton.placeRandom(
				random_x.getNumber(),
				random_y.getNumber(),
				random_z.getNumber(),
				true);
		}

		if (proton.xt().x() < -half_x || proton.xt().x() > half_x ||
			proton.xt().y() < -half_y || proton.xt().y() > half_y ||
			proton.xt().z() < -half_z || proton.xt().z() > half_z) {
			continue;
		}

		const std::vector<unsigned int> ionp_candidates =
			ionp_index.point_candidates(
				proton.xt(), ionp_search_radius);
		if (shell_proposal) {
			const std::size_t coverage = importance_shell_coverage(
				proton.xt(),
				ionps,
				populations,
				ionp_candidates,
				shell_proposal->population,
				shell_proposal->region,
				shell_proposal->near_cutoff,
				shell_proposal->far_cutoff);
			if (coverage == 0 ||
				random_unit.getNumber() >
					1.0 / static_cast<double>(coverage)) {
				continue;
			}
		}

		ProtonSampleRegion classified;
		if (classify_joint_sample_region_from_ionp_candidates(
				proton.xt(),
				cell_origins,
				cell_origin,
				cell_r,
				ionps,
				populations,
				ionp_candidates,
				importance_sampling,
				classified) &&
			classified == sample_region) {
			synchronize_initial_cell_membership(
				proton, cell_origins, cell_origin, cell_r);
			return;
		}
	}

	throw std::runtime_error(
		"Failed to place a proton in requested importance-sampling region " +
		sample_region_name(sample_region) + " after " +
		std::to_string(SimConsts::INIT_RETRY_LIMIT) + " attempts");
}

bool proposal_matches_membership(
	FreeSamplingProposal proposal,
	const JointSampleRegion& membership) {
	switch (proposal) {
		case FreeSamplingProposal::Uniform:
			return true;
		case FreeSamplingProposal::CellInside:
			return membership.cell == CellSampleRegion::Inside;
		case FreeSamplingProposal::CellNear:
			return membership.cell == CellSampleRegion::Near;
		case FreeSamplingProposal::CellFar:
			return membership.cell == CellSampleRegion::Far;
		case FreeSamplingProposal::IntracellularIonpNear:
			return membership.intracellular_ionp == IonpDistanceRegion::Near;
		case FreeSamplingProposal::IntracellularIonpIntermediate:
			return membership.intracellular_ionp == IonpDistanceRegion::Intermediate;
		case FreeSamplingProposal::IntracellularIonpFar:
			return membership.intracellular_ionp == IonpDistanceRegion::Far;
		case FreeSamplingProposal::ExtracellularIonpNear:
			return membership.extracellular_ionp == IonpDistanceRegion::Near;
		case FreeSamplingProposal::ExtracellularIonpIntermediate:
			return membership.extracellular_ionp == IonpDistanceRegion::Intermediate;
		case FreeSamplingProposal::ExtracellularIonpFar:
			return membership.extracellular_ionp == IonpDistanceRegion::Far;
	}
	return false;
}

IonpDistanceRegion classify_importance_surface_distance(
	double surface_distance,
	bool enabled,
	double near_cutoff,
	double far_cutoff) {
	if (!enabled || !std::isfinite(surface_distance)) {
		return IonpDistanceRegion::Far;
	}
	if (surface_distance <= near_cutoff) {
		return IonpDistanceRegion::Near;
	}
	if (surface_distance <= far_cutoff) {
		return IonpDistanceRegion::Intermediate;
	}
	return IonpDistanceRegion::Far;
}

CellSampleRegion classify_importance_cell_region(
	const Point3D& point,
	const CellSpatialIndex& cell_index,
	double cell_radius,
	double near_cutoff,
	int& inside_cell_number) {
	inside_cell_number = 0;
	if (cell_index.empty()) {
		return CellSampleRegion::Far;
	}

	thread_local std::vector<std::size_t> candidates;
	cell_index.point_candidates(
		point, std::max(cell_radius, near_cutoff), candidates);
	double nearest_squared = std::numeric_limits<double>::infinity();
	std::size_t nearest_cell = 0;
	for (std::size_t proxy_index : candidates) {
		const CellSpatialIndex::Proxy& proxy = cell_index.proxy(proxy_index);
		const double squared = squared_distance_to_center(point, proxy.center);
		if (squared < nearest_squared) {
			nearest_squared = squared;
			nearest_cell = proxy.cell_index;
		}
	}

	if (nearest_squared <= cell_radius * cell_radius) {
		inside_cell_number = static_cast<int>(nearest_cell + 1u);
		return CellSampleRegion::Inside;
	}
	if (near_cutoff > 0. && nearest_squared <= near_cutoff * near_cutoff) {
		return CellSampleRegion::Near;
	}
	return CellSampleRegion::Far;
}

bool classify_overlapping_start(
	const Point3D& point,
	const ImportanceSamplingConfig& config,
	const std::vector<Ionp>& ionps,
	const IonpPopulationMap& populations,
	const std::vector<unsigned int>& ionp_candidates,
	const CellSpatialIndex& cell_index,
	double cell_radius,
	JointSampleRegion& membership,
	int& inside_cell_number) {
	double intracellular_distance = std::numeric_limits<double>::infinity();
	double extracellular_distance = std::numeric_limits<double>::infinity();
	for (unsigned int candidate : ionp_candidates) {
		const std::size_t index = static_cast<std::size_t>(candidate);
		const double surface_distance =
			distance(point, ionps[index].xt()) - ionps[index].outer_radius();
		if (surface_distance < 0.) {
			return false;
		}
		if (populations.populations[index] == IonpPopulation::Intracellular) {
			intracellular_distance =
				std::min(intracellular_distance, surface_distance);
		}
		else {
			extracellular_distance =
				std::min(extracellular_distance, surface_distance);
		}
	}

	membership.cell = classify_importance_cell_region(
		point,
		cell_index,
		cell_radius,
		config.cell_near_cutoff,
		inside_cell_number);
	membership.intracellular_ionp = classify_importance_surface_distance(
		intracellular_distance,
		config.intracellular_ionp_sampling_enabled(),
		config.intracellular_ionp_near_cutoff,
		config.intracellular_ionp_far_cutoff);
	membership.extracellular_ionp = classify_importance_surface_distance(
		extracellular_distance,
		config.extracellular_ionp_sampling_enabled(),
		config.extracellular_ionp_near_cutoff,
		config.extracellular_ionp_far_cutoff);
	return true;
}

std::array<ImportanceShellProposal, kFreeSamplingProposalCount>
make_overlapping_shell_proposals(
	const ImportanceSamplingConfig& config,
	const std::vector<Ionp>& ionps,
	const IonpPopulationMap& populations) {
	std::array<ImportanceShellProposal, kFreeSamplingProposalCount> proposals{};
	proposals[free_sampling_proposal_index(
		FreeSamplingProposal::IntracellularIonpNear)] =
		make_importance_shell_proposal(
			ionps,
			populations.intracellular_indices,
			IonpPopulation::Intracellular,
			IonpDistanceRegion::Near,
			config.intracellular_ionp_near_cutoff,
			config.intracellular_ionp_far_cutoff);
	proposals[free_sampling_proposal_index(
		FreeSamplingProposal::IntracellularIonpIntermediate)] =
		make_importance_shell_proposal(
			ionps,
			populations.intracellular_indices,
			IonpPopulation::Intracellular,
			IonpDistanceRegion::Intermediate,
			config.intracellular_ionp_near_cutoff,
			config.intracellular_ionp_far_cutoff);
	proposals[free_sampling_proposal_index(
		FreeSamplingProposal::ExtracellularIonpNear)] =
		make_importance_shell_proposal(
			ionps,
			populations.extracellular_indices,
			IonpPopulation::Extracellular,
			IonpDistanceRegion::Near,
			config.extracellular_ionp_near_cutoff,
			config.extracellular_ionp_far_cutoff);
	proposals[free_sampling_proposal_index(
		FreeSamplingProposal::ExtracellularIonpIntermediate)] =
		make_importance_shell_proposal(
			ionps,
			populations.extracellular_indices,
			IonpPopulation::Extracellular,
			IonpDistanceRegion::Intermediate,
			config.extracellular_ionp_near_cutoff,
			config.extracellular_ionp_far_cutoff);
	return proposals;
}

std::vector<Point3D> importance_cell_centers(
	const std::vector<Point3D>& cell_origins,
	const Point3D& fallback_origin,
	double cell_radius) {
	std::vector<Point3D> centers;
	const std::size_t count = available_cell_count(cell_origins, cell_radius);
	centers.reserve(count);
	for (std::size_t index = 0; index < count; ++index) {
		centers.push_back(cell_center_at(cell_origins, fallback_origin, index));
	}
	return centers;
}

std::size_t importance_cell_shell_coverage(
	const Point3D& point,
	const CellSpatialIndex& cell_index,
	FreeSamplingProposal proposal,
	double cell_radius,
	double near_cutoff) {
	const double outer = proposal == FreeSamplingProposal::CellInside
		? cell_radius
		: near_cutoff;
	const double inner = proposal == FreeSamplingProposal::CellInside
		? 0.
		: cell_radius;
	thread_local std::vector<std::size_t> candidates;
	cell_index.point_candidates(point, outer, candidates);
	std::size_t coverage = 0;
	for (std::size_t proxy_index : candidates) {
		const double radius = distance(point, cell_index.proxy(proxy_index).center);
		if (radius >= inner && radius <= outer) {
			++coverage;
		}
	}
	return coverage;
}

// Draw from one overlapping marginal proposal, then test accessible-water
// membership on every axis. Return false if its existing retry budget expires.
bool initialize_proton_from_proposal(
	Proton& proton,
	FreeSamplingProposal proposal,
	JointSampleRegion& membership,
	const ImportanceSamplingConfig& config,
	const std::vector<Ionp>& ionps,
	const IonpPopulationMap& populations,
	const IonpSpatialIndex& ionp_index,
	const CellSpatialIndex& cell_index,
	const std::array<ImportanceShellProposal, kFreeSamplingProposalCount>& shell_proposals,
	const std::vector<Point3D>& cell_centers,
	double cell_radius,
	double sim_x_length,
	double sim_y_length,
	double sim_z_length,
	RandomList& random_x,
	RandomList& random_y,
	RandomList& random_z,
	RandomList& random_unit) {
	const std::size_t proposal_index = free_sampling_proposal_index(proposal);
	const ImportanceShellProposal& shell_proposal = shell_proposals[proposal_index];
	const bool use_ionp_shell = shell_proposal.total_volume > 0.;
	const bool use_cell_shell =
		proposal == FreeSamplingProposal::CellInside ||
		proposal == FreeSamplingProposal::CellNear;
	if (use_cell_shell && cell_centers.empty()) {
		return false;
	}
	if ((proposal == FreeSamplingProposal::IntracellularIonpNear ||
		 proposal == FreeSamplingProposal::IntracellularIonpIntermediate ||
		 proposal == FreeSamplingProposal::ExtracellularIonpNear ||
		 proposal == FreeSamplingProposal::ExtracellularIonpIntermediate) &&
		!use_ionp_shell) {
		return false;
	}

	const double half_x = sim_x_length / 2.;
	const double half_y = sim_y_length / 2.;
	const double half_z = sim_z_length / 2.;
	const double maximum_distance_cutoff = std::max(
		config.intracellular_ionp_sampling_enabled()
			? config.intracellular_ionp_far_cutoff
			: 0.,
		config.extracellular_ionp_sampling_enabled()
			? config.extracellular_ionp_far_cutoff
			: 0.);
	const double ionp_search_radius =
		ionp_index.max_outer_radius() + maximum_distance_cutoff;

	for (unsigned int attempt = 0; attempt < SimConsts::INIT_RETRY_LIMIT; ++attempt) {
		if (use_ionp_shell) {
			const double draw = random_unit.getNumber() * shell_proposal.total_volume;
			const auto chosen_position = std::lower_bound(
				shell_proposal.cumulative_volumes.begin(),
				shell_proposal.cumulative_volumes.end(),
				draw);
			const std::size_t chosen_offset =
				chosen_position == shell_proposal.cumulative_volumes.end()
					? shell_proposal.cumulative_volumes.size() - 1u
					: static_cast<std::size_t>(
						chosen_position - shell_proposal.cumulative_volumes.begin());
			const std::size_t chosen = (*shell_proposal.indices)[chosen_offset];
			double inner = 0.;
			double outer = 0.;
			importance_shell_bounds(
				ionps[chosen],
				shell_proposal.region,
				shell_proposal.near_cutoff,
				shell_proposal.far_cutoff,
				inner,
				outer);
			place_proton_random_sphere(
				proton,
				ionps[chosen].xt().x(),
				ionps[chosen].xt().y(),
				ionps[chosen].xt().z(),
				inner,
				outer,
				random_unit);
		}
		else if (use_cell_shell) {
			const std::size_t selected = std::min(
				static_cast<std::size_t>(
					random_unit.getNumber() *
					static_cast<double>(cell_centers.size())),
				cell_centers.size() - 1u);
			const Point3D& center = cell_centers[selected];
			place_proton_random_sphere(
				proton,
				center.x(),
				center.y(),
				center.z(),
				proposal == FreeSamplingProposal::CellInside ? 0. : cell_radius,
				proposal == FreeSamplingProposal::CellInside
					? cell_radius
					: config.cell_near_cutoff,
				random_unit);
		}
		else {
			proton.placeRandom(
				random_x.getNumber(),
				random_y.getNumber(),
				random_z.getNumber(),
				true);
		}

		if (proton.xt().x() < -half_x || proton.xt().x() > half_x ||
			proton.xt().y() < -half_y || proton.xt().y() > half_y ||
			proton.xt().z() < -half_z || proton.xt().z() > half_z) {
			continue;
		}

		const std::vector<unsigned int> ionp_candidates =
			ionp_index.point_candidates(proton.xt(), ionp_search_radius);
		if (use_ionp_shell) {
			const std::size_t coverage = importance_shell_coverage(
				proton.xt(),
				ionps,
				populations,
				ionp_candidates,
				shell_proposal.population,
				shell_proposal.region,
				shell_proposal.near_cutoff,
				shell_proposal.far_cutoff);
			if (coverage == 0 ||
				random_unit.getNumber() > 1. / static_cast<double>(coverage)) {
				continue;
			}
		}
		if (use_cell_shell) {
			const std::size_t coverage = importance_cell_shell_coverage(
				proton.xt(),
				cell_index,
				proposal,
				cell_radius,
				config.cell_near_cutoff);
			if (coverage == 0 ||
				random_unit.getNumber() > 1. / static_cast<double>(coverage)) {
				continue;
			}
		}

		int inside_cell_number = 0;
		if (!classify_overlapping_start(
				proton.xt(),
				config,
				ionps,
				populations,
				ionp_candidates,
				cell_index,
				cell_radius,
				membership,
				inside_cell_number) ||
			!proposal_matches_membership(proposal, membership)) {
			continue;
		}

		proton.inCell(membership.cell == CellSampleRegion::Inside);
		proton.Cellnum(inside_cell_number);
		return true;
	}
	return false;
}

// Extract the legacy MSE samples from the full recorded signal. The k=0
// sample is retained and the final nominal echo is excluded. These historical
// index formulas can disagree with recording cadence; see the simulation
// reference and use the written time column when interpreting sampled output.
void build_mse_signal(const std::vector<std::complex<double>>& full_signal,
					  const std::vector<double>& timeV,
					  int sequence_num,
					  double echo_time,
					  double echo_spacing,
					  double time_step,
					  std::vector<std::complex<double>>& signalVectorMSE,
					  std::vector<double>& timeVMSE) {
	signalVectorMSE.clear();
	timeVMSE.clear();
	if (sequence_num != 2) {
		return;
	}

	for (unsigned int k = 0; k < SC_UI(round(echo_time / echo_spacing)); ++k) {
		std::size_t signal_index = 0;
		if (time_step <= 0.00001 && echo_time >= 0.1) {
			signal_index = SC_UI(round(k * echo_spacing / 0.00001 / echo_time * 0.1));
		} else {
			signal_index = SC_UI(round(k * echo_spacing / time_step));
		}
		if (signal_index >= full_signal.size() || signal_index >= timeV.size()) {
			continue;
		}
		signalVectorMSE.push_back(full_signal[signal_index]);
		timeVMSE.push_back(timeV[signal_index]);
	}
}

// -- MRsequence member helper functions --

// Initialize a proton at j=0, or propose a diffusion step and resolve its boundaries.
// Failed walks restore the cached pre-step material state before another attempt.
void MRsequence::replaceProton(
    std::vector<Proton>& protonInstances,
    const std::vector<Ionp>& ionpInstances,
    const IonpSpatialIndex& ionp_index,
    const CellSpatialIndex& cell_index,
    int& coating_ionp_index,
    bool& coating_state_known,
    const std::vector<Point3D>& cell_origins,
    double dt,
    unsigned int i,
    unsigned int j,
    double current_time,
    ProtonSampleRegion sample_region,
    RandomList& RandomNormalX,
    RandomList& RandomNormalY,
    RandomList& RandomNormalZ,
    RandomList& RandomDiff,
    std::vector<TrajectoryEvent>* event_log,
    const Point3D* ionp_translation){

		/* Cell membranes, periodic faces, and IONP core/coating boundaries are
		   processed together in traveled-path order. */

	Proton& proton = protonInstances[i];
	const Proton proton_before_step = proton;
	const int coating_ionp_index_before_step = coating_ionp_index;
	const bool coating_state_known_before_step = coating_state_known;

	for (unsigned int attempt = 0; attempt < SimConsts::INIT_RETRY_LIMIT; ++attempt) {
		if (attempt > 0 && j != 0) {
			proton = proton_before_step;
		}
		Point3D lastPosition;

		if (!j) {
			if (importance_sampling.enabled &&
				!importance_sampling.overlapping_proposals_enabled()) {
				initialize_proton_in_region(
					proton,
					sample_region,
					importance_sampling,
					cell_origins,
					cell_o,
					cell_r,
					ionpInstances,
					ionp_populations,
					ionp_index,
					simMatrix.xLength(),
					simMatrix.yLength(),
					simMatrix.zLength(),
					RandomNormalX,
					RandomNormalY,
					RandomNormalZ,
					RandomDiff);
			}
			else if (!importance_sampling.enabled) {
				proton.placeRandom(RandomNormalX.getNumber(), RandomNormalY.getNumber(), RandomNormalZ.getNumber(), true);
			}

			lastPosition = proton.xt();
			coating_state_known = false;
		}
		else {
			lastPosition = proton.xt();
			const double diffusion_coefficient = diffusion_coefficient_from_coating_state(proton.xt(),
																						  ionpInstances,
																						  ionp_index,
																						  ionp_check_mode,
																						  coating_ionp_index,
																						  coating_state_known,
																  coating_diffusion_coefficient,
																  proton.inCell(),
																  cell_diffusion_coefficient,
																  ionp_translation);
			double radial_step = 0.0;
			const bool base_step = dt == time_step;
			const bool variable_substep =
				dt == time_step / SimConsts::VARSTEP_FACTOR;
			if (base_step || variable_substep) {
				if (coating_ionp_index >= 0) {
					radial_step = base_step
						? coating_diffusion_step
						: coating_diffusion_substep;
				}
				else if (proton.inCell()) {
					radial_step = base_step
						? cell_diffusion_step
						: cell_diffusion_substep;
				}
				else {
					radial_step = base_step
						? outside_diffusion_step
						: outside_diffusion_substep;
				}
			}
			else {
				radial_step = std::sqrt(6 * diffusion_coefficient * dt);
			}
			proton.diffuseWithStep(
				radial_step,
				RandomDiff.getNumber(),
				RandomDiff.getNumber());
		}

		// Cell membranes and periodic voxel faces share one path-ordered event
		// loop, so a seam never turns into a chord across the whole voxel.
		IonpBoundaryContext ionp_context{
			ionpInstances,
			ionp_index,
			ionp_check_mode,
			coating_ionp_index,
			coating_state_known,
			check_coating_permeability,
			coating_permeability,
			ionp_translation
		};
		const CellWalkResult cell_walk = walk_cell_boundaries(proton,
															 lastPosition,
															 cell_index,
															 &ionp_context,
															 cell_r,
															 cell_permeability,
															 simMatrix.xLength(),
															 simMatrix.yLength(),
															 simMatrix.zLength(),
															 true,
															 RandomDiff,
															 current_time,
															 event_log);

		if (!cell_walk.retry) {
			proton.setDiffTime(current_time);
			return;
		}

		coating_ionp_index = coating_ionp_index_before_step;
		coating_state_known = coating_state_known_before_step;
	}

	throw std::runtime_error("Failed to find a valid proton position after " +
							 std::to_string(SimConsts::INIT_RETRY_LIMIT) +
							 " attempts; check geometry, time_step, and excluded IONP/cell volumes");
}

// Simulated movement of IONP particles
// Translate every source by the same drift step and update the optional
// spatial index by that translation; relative source geometry is unchanged.
void MRsequence::replaceIONP(std::vector<Ionp>& ionpInstances, IonpSpatialIndex* ionp_index) {
	Point3D previous_first_position;
	if (ionp_index && !ionpInstances.empty()) {
		previous_first_position = ionpInstances.front().xt();
	}
	for (unsigned int k = 0; k < SC_UI(ionpInstances.size()); ++k) {
		ionpInstances[k].move(directionIONP, vIONP, time_step);
	}
	if (ionp_index && !ionpInstances.empty()) {
		const Point3D current_first_position = ionpInstances.front().xt();
		ionp_index->translate(current_first_position.x() - previous_first_position.x(),
							  current_first_position.y() - previous_first_position.y(),
							  current_first_position.z() - previous_first_position.z());
	}
}

// Calculates the B-field at the current position of an proton (with interpolation if precomputed grid is used)
// Dispatch the legacy grid/direct path. Direct evaluation sums all supplied
// sources, optionally with precomputed amplitudes and a common translation.
double MRsequence::computeBFieldAt(
    const Point3D& proton_pos,
    const std::vector<Ionp>& ionpInstances,
    const std::vector<double>& BfieldGrid,
    bool useGrid,
    const Point3D* ionp_translation,
    const std::vector<double>* field_amplitudes) const {
	double bfield_tot = 0.;
	// If grid is used and not empty, interpolate B-field from grid
	if ( useGrid && !BfieldGrid.empty() ) {
		return simMatrix.interpolate( proton_pos, BfieldGrid );
	}
	if (field_amplitudes && field_amplitudes->size() == ionpInstances.size()) {
		for (unsigned int l = 0; l < SC_UI(ionpInstances.size()); ++l) {
			bfield_tot += simMatrix.calculateBfieldWithAmplitude(
				proton_pos,
				translated_ionp_center(ionpInstances[l], ionp_translation),
				(*field_amplitudes)[l]);
		}
		return bfield_tot;
	}
	// Otherwise, compute B-field contributions from all IONP particles directly - using the IONP magnetic moment or the magnetization
	if ( magnetic_moment != 0. ) {
		for ( unsigned l = 0; l < SC_UI( ionpInstances.size() ); ++l ) {
			bfield_tot += simMatrix.calculateBfield( proton_pos, translated_ionp_center(ionpInstances[l], ionp_translation), ionpInstances[l].radius(), magnetic_moment );
		}
	} else {
		for ( unsigned l = 0; l < SC_UI( ionpInstances.size() ); ++l ) {
			bfield_tot += simMatrix.calculateBfield( proton_pos, translated_ionp_center(ionpInstances[l], ionp_translation), ionpInstances[l].radius() );
		}
	}
	return bfield_tot;
}

void MRsequence::set_output_context(std::filesystem::path output_dir) {
	output_tmp_dir = std::move(output_dir);
}

void MRsequence::set_bfield_solver(BFieldSolver solver) {
	bfield_solver = solver;
}

void MRsequence::set_surrounding_voxels(bool enabled) {
	surrounding_voxels = enabled;
}

void MRsequence::set_surrounding_voxel_superparticles(std::size_t count) {
	surrounding_voxel_superparticles = count;
}

void MRsequence::writeData ( const std::vector<std::complex<double>>& signalShellProtons, const std::vector<double>& timeV, const std::vector<Ionp>& ionpInstances_orig ) const
{
	if (output_tmp_dir.empty()) {
		throw std::runtime_error("MRsequence output directory is not configured");
	}

	std::filesystem::create_directories(output_tmp_dir);

	const std::size_t file_index = output::next_file_index(output_tmp_dir);
	const std::string prefix = std::to_string(file_index);

	std::string protonData_out_name;

	if ( sequence_num == 0 ){
		protonData_out_name = prefix + "_SignalCube_FID_" + std::to_string(int (r_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionpInstances_orig.size()) + "IONPs_" + std::to_string(protons_n) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te.dat";
	}
	else if ( sequence_num == 1 ){
		protonData_out_name = prefix + "_SignalCube_SE_" + std::to_string(int (r_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionpInstances_orig.size()) + "IONPs_" + std::to_string(protons_n) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te.dat";
	}
	else if ( sequence_num == 2 ){
		protonData_out_name = prefix + "_SignalCube_MSE_" + std::to_string(int (r_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionpInstances_orig.size()) + "IONPs_" + std::to_string(protons_n) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te_" + std::to_string(int (echo_spacing*pow(10., 3.))) + "ms_spacing.dat";
	}

	const auto out_path = output_tmp_dir / protonData_out_name;

	std::cout << out_path.string() << std::endl;

	std::ofstream outfile(out_path, std::ios::out | std::ios::binary);
	outfile.precision(15);
	if (!outfile.is_open()) {
		throw std::runtime_error("Unable to open output file: " + out_path.string());
	}
	for (unsigned int i = 0; i < signalShellProtons.size(); ++i) {
		outfile << real(signalShellProtons[i]) << "\t" << timeV[i] << std::endl;
	}
	outfile.close();
}

// Flatten each proton's event history to the existing trajectory columns;
// hit and post-reflection events supplement the regular sampling points.
void MRsequence::writeProtonPositions(const std::vector<std::vector<TrajectoryEvent>>& events, std::size_t ionp_count) const {
	if (!record_positions) {
		return;
	}
	if (output_tmp_dir.empty()) {
		throw std::runtime_error("MRsequence output directory is not configured");
	}
	if (events.empty()) {
		return;
	}

	std::filesystem::create_directories(output_tmp_dir);
	const std::size_t file_index = output::next_file_index(output_tmp_dir);
	const std::string prefix = std::to_string(file_index);

	std::string seq_tag = "FID";
	if (sequence_num == 1) {
		seq_tag = "SE";
	} else if (sequence_num == 2) {
		seq_tag = "MSE";
	}

	const std::size_t proton_count = events.size();

	std::string filename = prefix + "_ProtonTrajectories_" + seq_tag + "_" + std::to_string(int(r_ionp * 1e9)) + "nm_" + std::to_string(ionp_count) + "IONPs_" + std::to_string(proton_count) + "Protons_" + std::to_string(int(time_step * 1e9)) + "ns_dt_" + std::to_string(int(echo_time * 1e3)) + "ms_te";
	if (sequence_num == 2) {
		filename += "_" + std::to_string(int(echo_spacing * 1e3)) + "ms_spacing";
	}
	filename += ".dat";

	const auto out_path = output_tmp_dir / filename;
	std::cout << out_path.string() << std::endl;

	std::ofstream outfile(out_path, std::ios::out | std::ios::binary);
	outfile.precision(15);
	if (!outfile.is_open()) {
		throw std::runtime_error("Unable to open output file: " + out_path.string());
	}

	outfile << "# time[s]\tproton_index\tkind\tx\ty\tz\n";
	for (std::size_t p = 0; p < events.size(); ++p) {
		for (const auto& ev : events[p]) {
			outfile << ev.time << "\t" << p << "\t" << ev.kind << "\t" << ev.position.x() << "\t" << ev.position.y() << "\t" << ev.position.z() << std::endl;
		}
	}
}

// Orchestrate one run: prepare sources/grid, resolve sample allocations,
// initialize and propagate free/bound spins, then aggregate recorded phases.
// Optional diagnostics share the same time axis as their corresponding signal.
void MRsequence::Sequence(std::vector<Proton>& protonInstances,
							  std::vector<Ionp>& ionpInstances_orig,
							  std::vector<std::complex<double>>& signalVector,
							  std::vector<double>& timeV,
							  std::vector<std::complex<double>>& signalVectorMSE,
							  std::vector<double>& timeVMSE,
							  std::vector<Point3D>& cell_origins,
							  ImportanceRegionSignals* regionSignalVectors,
							  ImportanceRegionSignals* regionSignalVectorsMSE,
							  ImportanceProposalSignalDiagnostics* proposalSignalDiagnostics,
							  ImportanceProposalSignalDiagnostics* proposalSignalDiagnosticsMSE) {
	if (regionSignalVectors) {
		for (auto& signal : *regionSignalVectors) {
			signal.clear();
		}
	}
	if (regionSignalVectorsMSE) {
		for (auto& signal : *regionSignalVectorsMSE) {
			signal.clear();
		}
	}
	auto clear_proposal_diagnostics =
		[](ImportanceProposalSignalDiagnostics* diagnostics) {
			if (!diagnostics) {
				return;
			}
			for (auto& signal : diagnostics->means) {
				signal.clear();
			}
			for (auto& signal : diagnostics->contributions) {
				signal.clear();
			}
		};
	clear_proposal_diagnostics(proposalSignalDiagnostics);
	clear_proposal_diagnostics(proposalSignalDiagnosticsMSE);

	if (surrounding_voxel_superparticles > 0 && !surrounding_voxels) {
		throw std::invalid_argument(
			"Surrounding_voxel_superparticles requires Surrounding_voxels = Yes");
	}
	if (surrounding_voxel_superparticles > 0 && magnetic_moment != 0.0) {
		throw std::invalid_argument(
			"Surrounding_voxel_superparticles requires Magnetic Moment = 0");
	}
	if (echo_time > 0. && time_step <= 0.) {
		throw std::runtime_error("MRsequence requires a positive time_step for non-zero echo_time");
	}

	std::vector<Ionp> surrounding_field_sources;
	const std::vector<Ionp>* magnetic_field_sources_ptr = &ionpInstances_orig;
	if (surrounding_voxels) {
		surrounding_field_sources =
			SimSpace::surroundingSpace(
				ionpInstances_orig,
				simMatrix,
				surrounding_voxel_superparticles,
				cell_origins,
				cell_o,
				cell_r);
		magnetic_field_sources_ptr = &surrounding_field_sources;
	}
	const std::vector<Ionp>& magnetic_field_sources =
		*magnetic_field_sources_ptr;
	if (magnetic_field_sources.size() < ionpInstances_orig.size()) {
		throw std::logic_error(
			"Surrounding magnetic source collection lost physical IONPs");
	}
	const std::size_t replica_source_count =
		magnetic_field_sources.size() - ionpInstances_orig.size();
	constexpr std::size_t surrounding_voxel_replica_count = 26;
	if (surrounding_voxels &&
		replica_source_count % surrounding_voxel_replica_count != 0) {
		throw std::logic_error(
			"Surrounding magnetic source collection does not contain 26 identical templates");
	}
	const std::size_t template_source_count = surrounding_voxels
		? replica_source_count / surrounding_voxel_replica_count
		: 0;
	auto core_radius_cubed_sum = [](const std::vector<Ionp>& ionps,
									  std::size_t begin,
									  std::size_t count) {
		double sum = 0.0;
		for (std::size_t index = begin; index < begin + count; ++index) {
			const double radius = ionps[index].radius();
			sum += radius * radius * radius;
		}
		return sum;
	};
	const double physical_core_volume = core_radius_cubed_sum(
		ionpInstances_orig, 0, ionpInstances_orig.size());
	const double template_core_volume = surrounding_voxels
		? core_radius_cubed_sum(
			  magnetic_field_sources,
			  ionpInstances_orig.size(),
			  template_source_count)
		: 0.0;
	const char* surrounding_source_mode = !surrounding_voxels
		? "disabled"
		: (surrounding_voxel_superparticles == 0 ? "exact" : "superparticles");
	std::cout << "Physical IONP count: " << ionpInstances_orig.size()
		<< std::endl;
	std::cout << "Surrounding voxel source mode: "
		<< surrounding_source_mode << std::endl;
	std::cout << "Surrounding voxel superparticle cap: "
		<< surrounding_voxel_superparticles << std::endl;
	std::cout << "Surrounding voxel template source count: "
		<< template_source_count << std::endl;
	std::cout << "Physical IONP sum(r^3): " << physical_core_volume
		<< std::endl;
	if (surrounding_voxels) {
		std::cout << "Surrounding template sum(r^3): "
			<< template_core_volume << std::endl;
		std::cout << "Surrounding template iron difference sum(r^3): "
			<< (template_core_volume - physical_core_volume) << std::endl;
	}
	std::cout << "Surrounding voxel replica count: "
		<< (surrounding_voxels ? surrounding_voxel_replica_count : 0)
		<< std::endl;
	std::cout << "Surrounding magnetic replica count: "
		<< replica_source_count << std::endl;
	std::cout << "Total magnetic field-source count: "
		<< magnetic_field_sources.size() << std::endl;

	// Collision geometry uses physical central IONPs. The field-source index
	// may additionally contain surrounding replicas and magnetic proxies.
	const IonpSpatialIndex static_ionp_index(ionpInstances_orig);
	std::optional<IonpSpatialIndex> surrounding_field_source_index;
	const IonpSpatialIndex* magnetic_field_source_index = &static_ionp_index;
	if (surrounding_voxels) {
		surrounding_field_source_index.emplace(magnetic_field_sources);
		magnetic_field_source_index = &*surrounding_field_source_index;
	}

	//Precompute grid if requested
	std::vector<double> BfieldGrid;
	if (Grid == true){
		const unsigned int nx = SC_UI(simMatrix.nxNodes());
		const unsigned int ny = SC_UI(simMatrix.nyNodes());
		const unsigned int nz = SC_UI(simMatrix.nzNodes());
		BfieldGrid.resize(static_cast<std::size_t>(nx) * ny * nz);
		simMatrix.calculateBfield(
			BfieldGrid, magnetic_field_sources, magnetic_moment);
	}
	std::vector<double> direct_field_amplitudes;
	direct_field_amplitudes.reserve(magnetic_field_sources.size());
	if (magnetic_moment != 0.0) {
		const double amplitude =
			D_MO / (4.0 * D_PI) * magnetic_moment;
		direct_field_amplitudes.assign(
			magnetic_field_sources.size(), amplitude);
	}
	else {
		for (const Ionp& ionp : magnetic_field_sources) {
			const double radius_cubed =
				ionp.radius() * ionp.radius() * ionp.radius();
			direct_field_amplitudes.push_back(B_fixed * radius_cubed);
		}
	}
	std::optional<CorrectedHybridField> corrected_hybrid_field;
	if (bfield_solver == BFieldSolver::CorrectedHybrid) {
		if (!Grid || BfieldGrid.empty()) {
			throw std::runtime_error(
				"BField_solver=corrected_hybrid requires a non-empty B-field grid");
		}
		if (vIONP != 0.0) {
			throw std::runtime_error(
				"BField_solver=corrected_hybrid requires stationary IONPs");
		}
		corrected_hybrid_field.emplace(
			simMatrix,
			magnetic_field_sources,
			direct_field_amplitudes,
			BfieldGrid);
		std::cout << "Corrected-hybrid B-field evaluation enabled" << std::endl;
		std::cout << "Hybrid grid spacing: "
			<< simMatrix.maxGridSpacing() << std::endl;
		std::cout << "Hybrid maximum correction radius: "
			<< corrected_hybrid_field->maximumCorrectionEnd()
			<< std::endl;
	}

	std::cout << "Real Time Steps: " << time_step << std::endl;
	if (!ionpInstances_orig.empty()) {
		std::cout << "Realized outer radius (first IONP if radius distribution): " << ionpInstances_orig[0].outer_radius() << std::endl;
	} else {
		std::cout << "No IONPs present for this sequence" << std::endl;
	}

	//==================== Proton diffusion ==========================//

	//Use maximum number available threads if no number is given - for parallelization
	if (numThreads == 0){
    	numThreads = omp_get_max_threads( );
		std::cout << "Number of max Threads (used): " << omp_get_max_threads( ) << std::endl;
  	}

	omp_set_num_threads( numThreads );

	const std::complex<double> i_compl(0, 1); // i
	ionp_populations = classify_ionp_populations(
		ionpInstances_orig, cell_origins, cell_o, cell_r);
	if (importance_sampling.enabled &&
		importance_sampling.overlapping_proposals_enabled() &&
		!importance_sampling.proposal_volumes_estimated) {
		const OverlappingProposalVolumeEstimate estimate =
			estimate_overlapping_proposal_volumes(
				importance_sampling,
				cell_origins,
				cell_o,
				cell_r,
				ionpInstances_orig,
				ionp_populations,
				simMatrix.xLength(),
				simMatrix.yLength(),
				simMatrix.zLength());
		if (!(estimate.accessible_volume > 0.) ||
			!(estimate.box_volume > 0.)) {
			throw std::runtime_error(
				"Overlapping importance sampling found no accessible free-water volume");
		}
		const EffectiveProposalCounts effective =
			resolve_effective_proposal_counts(
				importance_sampling.proposal_protons,
				estimate.region_fractions);
		for (std::size_t proposal_index = 0;
			 proposal_index < kFreeSamplingProposalCount;
			 ++proposal_index) {
			if (effective.counts[proposal_index] >
				static_cast<std::size_t>(std::numeric_limits<int>::max())) {
				throw std::runtime_error(
					"Overlapping proposal count exceeds supported integer range");
			}
			importance_sampling.proposal_protons[proposal_index] =
				static_cast<int>(effective.counts[proposal_index]);
			if (effective.redistributed[proposal_index]) {
				std::cout << "Warning: proposal "
					<< free_sampling_proposal_name(
						static_cast<FreeSamplingProposal>(proposal_index))
					<< " has zero estimated support; its samples were moved to uniform"
					<< std::endl;
			}
		}
		importance_sampling.proposal_volume_fractions =
			estimate.region_fractions;
		importance_sampling.free_water_volume_fraction =
			std::clamp(
				estimate.accessible_volume / estimate.box_volume,
				0.0,
				1.0);
		importance_sampling.proposal_volumes_estimated = true;
	}

	// Bound counts are per coated particle. Record jobs in physical source
	// order so each bound spin receives a stable index and stream seed.
	struct BoundSampleJob {
		std::size_t ionp_index;
		ProtonSampleRegion region;
	};
	std::vector<BoundSampleJob> bound_jobs;
	std::size_t bound_intracellular_sample_count = 0;
	std::size_t bound_extracellular_sample_count = 0;
	for (std::size_t ionp_index = 0;
		 ionp_index < ionpInstances_orig.size();
		 ++ionp_index) {
		if (ionpInstances_orig[ionp_index].coating() <= 0.) {
			continue;
		}
		const bool intracellular =
			ionp_populations.populations[ionp_index] ==
			IonpPopulation::Intracellular;
		const int samples_per_particle =
			importance_sampling.enabled &&
				importance_sampling.split_bound_configured
					? (intracellular
					? importance_sampling.bound_protons_intracellular
					: importance_sampling.bound_protons_extracellular)
				: bound_protons_n;
		const ProtonSampleRegion bound_region =
			intracellular
				? ProtonSampleRegion::BoundIntracellular
				: ProtonSampleRegion::BoundExtracellular;
		for (int sample = 0; sample < samples_per_particle; ++sample) {
			bound_jobs.push_back({ionp_index, bound_region});
			if (intracellular) {
				++bound_intracellular_sample_count;
			}
			else {
				++bound_extracellular_sample_count;
			}
		}
	}
	const std::size_t bound_sample_count = bound_jobs.size();
	unsigned int all_proton_num = SC_UI(protons_n) + SC_UI(bound_sample_count);

	//Function to determine ideal std for specific T2 of Water/Agar
	double standarddeviation_back = calc_std_T2( time_step, background_T2 );
	double standarddeviation_cell = calc_std_T2( time_step, incell_T2 );

	for (unsigned int i = SC_UI(protonInstances.size()); i < all_proton_num; ++i) protonInstances.emplace_back(Proton(0, 0));

	const bool ionpMoves = (vIONP != 0. && Grid == false);
	Point3D ionp_step_displacement;
	if (ionpMoves) {
		const double direction_norm =
			std::sqrt(square(directionIONP.x()) +
					  square(directionIONP.y()) +
					  square(directionIONP.z()));
		if (std::fabs(direction_norm - 1.0) > 1.e-10) {
			const double inverse_norm = 1.0 / direction_norm;
			directionIONP.set(inverse_norm * directionIONP.x(),
							  inverse_norm * directionIONP.y(),
							  inverse_norm * directionIONP.z());
		}
		ionp_step_displacement.set(directionIONP.x() * vIONP * time_step,
								   directionIONP.y() * vIONP * time_step,
								   directionIONP.z() * vIONP * time_step);
	}

    std::vector<std::vector<double>> phaseVector(all_proton_num); // per-proton phase history
    std::vector<std::vector<int>> protonPosition(all_proton_num);  // per-proton position mask
	std::vector<std::vector<TrajectoryEvent>> protonTrajectory;
	if (record_positions) {
		protonTrajectory.resize(all_proton_num); // per-proton position history
	}

	// Sampling labels describe initialization; material membership evolves
	// separately in Proton. Keep these labels fixed for the estimator weights.
	std::vector<ProtonSampleRegion> protonSampleRegions(
		all_proton_num, ProtonSampleRegion::Far);
	const bool overlapping_proposals =
		importance_sampling.enabled &&
		importance_sampling.overlapping_proposals_enabled();
	std::vector<FreeSamplingProposal> protonProposalSources;
	std::vector<JointSampleRegion> protonInitialMemberships;
	std::vector<double> protonImportanceWeights;
	if (overlapping_proposals) {
		protonProposalSources.reserve(SC_UI(protons_n));
		for (std::size_t proposal_index = 0;
			 proposal_index < kFreeSamplingProposalCount;
			 ++proposal_index) {
			const int count = importance_sampling.proposal_protons[proposal_index];
			if (count < 0) {
				throw std::runtime_error(
					"Overlapping importance-sampling counts must be non-negative");
			}
			for (int sample = 0; sample < count; ++sample) {
				protonProposalSources.push_back(
					static_cast<FreeSamplingProposal>(proposal_index));
			}
		}
		if (protonProposalSources.size() != SC_UI(protons_n)) {
			throw std::runtime_error(
				"Overlapping proposal counts do not match configured free proton total");
		}
		protonInitialMemberships.resize(SC_UI(protons_n));
	}
	else if (importance_sampling.enabled) {
		std::size_t proton_index = 0;
		for (std::size_t region_index = 0;
			 region_index < kFreeSampleRegionCount;
			 ++region_index) {
			for (int sample = 0;
				 sample < importance_sampling.joint_protons[region_index];
				 ++sample) {
				if (proton_index >= SC_UI(protons_n)) {
					throw std::runtime_error(
						"Joint importance-sampling counts exceed configured free proton total");
				}
				protonSampleRegions[proton_index++] =
					static_cast<ProtonSampleRegion>(region_index);
			}
		}
		if (proton_index != SC_UI(protons_n)) {
			throw std::runtime_error(
				"Joint importance-sampling counts do not match configured free proton total");
		}
	}
	for (std::size_t bound_index = 0;
		 bound_index < bound_jobs.size();
		 ++bound_index) {
		protonSampleRegions[SC_UI(protons_n) + bound_index] =
			bound_jobs[bound_index].region;
	}

	std::array<std::size_t, kSampleRegionCount> regionSampleCounts{};
	std::array<double, kSampleRegionCount> normalizedWeights{};
	if (importance_sampling.enabled && !overlapping_proposals) {
		regionSampleCounts = sample_region_counts(
			importance_sampling,
			bound_intracellular_sample_count,
			bound_extracellular_sample_count);
		normalizedWeights = normalized_sample_weights(
			importance_sampling,
			bound_intracellular_sample_count,
			bound_extracellular_sample_count);
		std::cout << "Normalized active importance-sampling weights:" << std::endl;
		for (std::size_t region_index = 0;
			 region_index < kSampleRegionCount;
			 ++region_index) {
			if (regionSampleCounts[region_index] == 0) {
				continue;
			}
			std::cout << "  "
				<< sample_region_name(
					static_cast<ProtonSampleRegion>(region_index))
				<< ": " << normalizedWeights[region_index]
				<< std::endl;
		}
	}

	if ( cube_s.x() != 0. ){
		if ( std::abs(cube_o.x() - cube_s.x()/2.) > simMatrix.xLength() / 2. || std::abs(cube_o.x() + cube_s.x()/2.) > simMatrix.xLength() / 2. || std::abs(cube_o.y() - cube_s.y()/2.) > simMatrix.yLength() / 2. || std::abs(cube_o.y() + cube_s.y()/2.) > simMatrix.yLength() / 2. || std::abs(cube_o.z() - cube_s.z()/2.) > simMatrix.zLength() / 2. || std::abs(cube_o.z() + cube_s.z()/2.) > simMatrix.zLength() / 2.) {
			std::cout << "Cube is outside of Simulation Volume" << std::endl;
		}
		if (importance_sampling.enabled) {
			std::cout << "Cube signal output is disabled in importance-sampling mode" << std::endl;
		}
	}

	//Proton loop (every proton is simulated separately)

	std::cout << "Now protons are propagating..." << std::endl;

	// helper lambda to check and optionally record position in cube
    auto recordPositionInCube = [this](std::vector<int> &posVec, Proton const &p) {
        if (cube_s.x() == 0.) return;
        posVec.push_back(0);
        if (p.xt().x() > (cube_o.x() - cube_s.x() / 2.) && p.xt().x() < (cube_o.x() + cube_s.x() / 2.) &&
            p.xt().y() > (cube_o.y() - cube_s.y() / 2.) && p.xt().y() < (cube_o.y() + cube_s.y() / 2.) &&
            p.xt().z() > (cube_o.z() - cube_s.z() / 2.) && p.xt().z() < (cube_o.z() + cube_s.z() / 2.)) {
            posVec.back() = 1;
        }
    };
	// helper lambda to optionally record full proton trajectory
	auto recordTrajectory = [this](std::vector<TrajectoryEvent>& events, Proton const &p, double current_time, int kind) {
		if (!record_positions) return;
		events.push_back({current_time, p.xt(), kind});
		};
		const double output_time_step = phase_record_interval(time_step, echo_time);
		const std::size_t expected_phase_samples =
			output_time_step > 0.0
				? static_cast<std::size_t>(
					  std::floor(echo_time / output_time_step + 1e-9)) +
					  1
				: 1;
		timeV.reserve(timeV.size() + expected_phase_samples);
		const std::vector<unsigned char> refocusing_schedule =
			build_refocusing_schedule(
				sequence_num,
				echo_spacing,
				time_step,
				echo_time);
	auto recordPhaseSample = [&](std::vector<double>& phaseOneProt,
								 std::vector<int>& positionOneProt,
								 std::vector<TrajectoryEvent>& trajectoryOneProt,
								 Proton const& proton,
								 unsigned int proton_index,
								 double phase,
								 double current_time,
								 double& next_record_time) {
		if (next_record_time > echo_time + output_time_step * 1e-9 || !reached_record_time(current_time, next_record_time)) {
			return;
		}

		phaseOneProt.push_back(phase);
		if (!proton_index) {
			timeV.push_back(current_time);
		}
		recordPositionInCube(positionOneProt, proton);
		recordTrajectory(trajectoryOneProt, proton, current_time, 0);

		do {
			next_record_time += output_time_step;
		} while (reached_record_time(current_time, next_record_time));
		};
			IonpSpatialIndex importance_ionp_index;
			if (overlapping_proposals) {
				const double importance_maximum_cutoff = std::max(
					importance_sampling.intracellular_ionp_sampling_enabled()
						? importance_sampling.intracellular_ionp_far_cutoff
						: 0.,
					importance_sampling.extracellular_ionp_sampling_enabled()
						? importance_sampling.extracellular_ionp_far_cutoff
						: 0.);
				importance_ionp_index.rebuild(
					ionpInstances_orig,
					static_ionp_index.max_outer_radius() +
						importance_maximum_cutoff);
			}
				const CellSpatialIndex static_cell_index(cell_origins,
										cell_o,
										cell_r,
											simMatrix.xLength(),
											simMatrix.yLength(),
											simMatrix.zLength(),
											true);
				std::optional<CellSpatialIndex> importance_cell_index;
				if (overlapping_proposals) {
					importance_cell_index.emplace(
						cell_origins,
						cell_o,
						cell_r,
						simMatrix.xLength(),
						simMatrix.yLength(),
						simMatrix.zLength(),
						false);
				}
			auto should_use_field_grid =
				[&](const Proton& proton,
					const Point3D* ionp_translation) {
					if (!Grid) {
						return false;
					}
					if (!proton.inCell()) {
						return true;
					}
					return !has_ionp_within_grid_influence(
						proton.xt(),
						magnetic_field_sources,
						*magnetic_field_source_index,
						ionp_check_mode,
						ionp_translation);
				};
			auto compute_field_for_step =
				[&](const Proton& proton,
					const Point3D* ionp_translation,
					bool variable_step_active) {
					if (corrected_hybrid_field) {
						return corrected_hybrid_field->evaluate(proton.xt());
					}
					const bool use_field_grid = should_use_field_grid_for_step(
						should_use_field_grid(proton, ionp_translation),
						variable_step_active);
					return computeBFieldAt(
						proton.xt(),
						magnetic_field_sources,
						BfieldGrid,
						use_field_grid,
						ionp_translation,
						&direct_field_amplitudes);
				};
				const std::uint32_t simulation_random_seed = rng_config::next_seed();
				if (overlapping_proposals) {
					const auto proposal_initialization_start =
						std::chrono::steady_clock::now();
					const auto shell_proposals =
						make_overlapping_shell_proposals(
							importance_sampling,
							ionpInstances_orig,
							ionp_populations);
					const std::vector<Point3D> cell_centers =
						importance_cell_centers(cell_origins, cell_o, cell_r);
					std::array<std::size_t, kFreeSamplingProposalCount> effective_counts{};
					std::array<std::size_t, kFreeSamplingProposalCount> proposal_starts{};
					std::size_t next_start = 0;
					for (std::size_t proposal_index = 0;
						 proposal_index < kFreeSamplingProposalCount;
						 ++proposal_index) {
						proposal_starts[proposal_index] = next_start;
						effective_counts[proposal_index] = static_cast<std::size_t>(
							importance_sampling.proposal_protons[proposal_index]);
						next_start += effective_counts[proposal_index];
					}

					auto initialize_start =
						[&](std::size_t proton_index,
							FreeSamplingProposal proposal,
							bool fallback_stream) {
							const ProtonRandomStream position_x_stream =
								fallback_stream
									? ProtonRandomStream::ImportanceFallbackPositionX
									: ProtonRandomStream::PositionX;
							const ProtonRandomStream position_y_stream =
								fallback_stream
									? ProtonRandomStream::ImportanceFallbackPositionY
									: ProtonRandomStream::PositionY;
							const ProtonRandomStream position_z_stream =
								fallback_stream
									? ProtonRandomStream::ImportanceFallbackPositionZ
									: ProtonRandomStream::PositionZ;
							const ProtonRandomStream initialization_stream =
								fallback_stream
									? ProtonRandomStream::ImportanceFallbackInitialization
									: ProtonRandomStream::ImportanceInitialization;
							RandomList random_x(
								-simMatrix.xLength() / 2.,
								simMatrix.xLength() / 2.,
								false,
								proton_random_seed(
									simulation_random_seed,
									proton_index,
									position_x_stream));
							RandomList random_y(
								-simMatrix.yLength() / 2.,
								simMatrix.yLength() / 2.,
								false,
								proton_random_seed(
									simulation_random_seed,
									proton_index,
									position_y_stream));
							RandomList random_z(
								-simMatrix.zLength() / 2.,
								simMatrix.zLength() / 2.,
								false,
								proton_random_seed(
									simulation_random_seed,
									proton_index,
									position_z_stream));
							RandomList random_unit(
								0.,
								1.,
								false,
								proton_random_seed(
									simulation_random_seed,
									proton_index,
									initialization_stream));
							return initialize_proton_from_proposal(
								protonInstances[proton_index],
								proposal,
								protonInitialMemberships[proton_index],
								importance_sampling,
								ionpInstances_orig,
								ionp_populations,
								importance_ionp_index,
								*importance_cell_index,
								shell_proposals,
								cell_centers,
								cell_r,
								simMatrix.xLength(),
								simMatrix.yLength(),
								simMatrix.zLength(),
								random_x,
								random_y,
								random_z,
								random_unit);
						};

					for (std::size_t proposal_index = 0;
						 proposal_index < kFreeSamplingProposalCount;
						 ++proposal_index) {
						const std::size_t count = effective_counts[proposal_index];
						if (count == 0) {
							continue;
						}
						const FreeSamplingProposal proposal =
							static_cast<FreeSamplingProposal>(proposal_index);
						int placement_success = 1;
						#pragma omp parallel for schedule(dynamic, 1) reduction(&:placement_success)
						for (std::size_t offset = 0; offset < count; ++offset) {
							const std::size_t proton_index =
								proposal_starts[proposal_index] + offset;
							placement_success &=
								initialize_start(proton_index, proposal, false) ? 1 : 0;
						}
						if (placement_success != 0) {
							continue;
						}
						if (proposal == FreeSamplingProposal::Uniform) {
							throw std::runtime_error(
								"Failed to initialize the mandatory uniform free-water proposal");
						}

						std::cout << "Warning: proposal "
							<< free_sampling_proposal_name(proposal)
							<< " could not fill its quota; all " << count
							<< " samples were moved to uniform" << std::endl;
						int uniform_success = 1;
						#pragma omp parallel for schedule(dynamic, 1) reduction(&:uniform_success)
						for (std::size_t offset = 0; offset < count; ++offset) {
							const std::size_t proton_index =
								proposal_starts[proposal_index] + offset;
							protonProposalSources[proton_index] =
								FreeSamplingProposal::Uniform;
							uniform_success &= initialize_start(
								proton_index,
								FreeSamplingProposal::Uniform,
								true) ? 1 : 0;
						}
						if (uniform_success == 0) {
							throw std::runtime_error(
								"Failed to replace an unusable proposal with uniform free-water samples");
						}
						effective_counts[
							free_sampling_proposal_index(FreeSamplingProposal::Uniform)] += count;
						effective_counts[proposal_index] = 0;
						importance_sampling.proposal_protons[
							free_sampling_proposal_index(FreeSamplingProposal::Uniform)] +=
							static_cast<int>(count);
						importance_sampling.proposal_protons[proposal_index] = 0;
					}

					protonImportanceWeights =
						overlapping_mis_self_normalized_weights(
							protonInitialMemberships,
							effective_counts,
							importance_sampling.proposal_volume_fractions);
					protonInitialMemberships.clear();
					protonInitialMemberships.shrink_to_fit();
					std::array<double, kFreeSamplingProposalCount> proposal_weight_mass{};
					double inverse_ess = 0.;
					for (std::size_t proton_index = 0;
						 proton_index < protonImportanceWeights.size();
						 ++proton_index) {
						const double weight = protonImportanceWeights[proton_index];
						proposal_weight_mass[free_sampling_proposal_index(
							protonProposalSources[proton_index])] += weight;
						inverse_ess += weight * weight;
					}
					std::cout << "Overlapping proposal MIS weight mass:" << std::endl;
					for (std::size_t proposal_index = 0;
						 proposal_index < kFreeSamplingProposalCount;
						 ++proposal_index) {
						if (effective_counts[proposal_index] == 0) {
							continue;
						}
						std::cout << "  "
							<< free_sampling_proposal_name(
								static_cast<FreeSamplingProposal>(proposal_index))
							<< ": samples=" << effective_counts[proposal_index]
							<< ", normalized_weight_mass="
							<< proposal_weight_mass[proposal_index]
							<< std::endl;
					}
					std::cout << "Overlapping proposal effective sample size: "
						<< (inverse_ess > 0. ? 1. / inverse_ess : 0.)
						<< " / " << protonImportanceWeights.size() << std::endl;
					const std::chrono::duration<double>
						proposal_initialization_elapsed =
							std::chrono::steady_clock::now() -
							proposal_initialization_start;
					std::cout << "Overlapping proposal initialization seconds: "
						<< proposal_initialization_elapsed.count() << std::endl;
				}

				auto simulate_single_proton = [&](unsigned int index, const std::vector<Ionp>& ionpInstances_orig_local, bool ionpMovesLocal) {
				const ProtonSampleRegion sample_region = protonSampleRegions[index];
				const bool use_variable_steps = should_use_variable_time_steps(var_Steps, importance_sampling.enabled, sample_region);
				RandomList RandomNormalX(
					-simMatrix.xLength() / 2.,
					simMatrix.xLength() / 2.,
					false,
					proton_random_seed(simulation_random_seed, index, ProtonRandomStream::PositionX));
				RandomList RandomNormalY(
					-simMatrix.yLength() / 2.,
					simMatrix.yLength() / 2.,
					false,
					proton_random_seed(simulation_random_seed, index, ProtonRandomStream::PositionY));
				RandomList RandomNormalZ(
					-simMatrix.zLength() / 2.,
					simMatrix.zLength() / 2.,
					false,
					proton_random_seed(simulation_random_seed, index, ProtonRandomStream::PositionZ));
				RandomList RandomDiff(
					0.,
					1.,
					false,
					proton_random_seed(
						simulation_random_seed,
						index,
						ProtonRandomStream::DiffusionAndBoundaries));

			Point3D ionp_translation;
			Point3D ionp_reference_position =
				ionpInstances_orig_local.empty()
					? Point3D()
					: ionpInstances_orig_local.front().xt();
			const Point3D* ionp_translation_for_calc =
				ionpMovesLocal ? &ionp_translation : nullptr;
			auto advance_ionp_translation = [&]() {
				const Point3D previous_reference_position = ionp_reference_position;
				ionp_reference_position.set(
					ionp_reference_position.x() + ionp_step_displacement.x(),
					ionp_reference_position.y() + ionp_step_displacement.y(),
					ionp_reference_position.z() + ionp_step_displacement.z());
				ionp_translation.set(
					ionp_translation.x() +
						(ionp_reference_position.x() - previous_reference_position.x()),
					ionp_translation.y() +
						(ionp_reference_position.y() - previous_reference_position.y()),
					ionp_translation.z() +
						(ionp_reference_position.z() - previous_reference_position.z()));
			};

		//================================= Start Initialization =================================//

		//Initiate Proton specific variables
		double phase = 0., time = 0.;
		unsigned int j = 0, leapfrog = 0;
		int coating_ionp_index = -1;
		bool coating_state_known = false;
		double next_record_time = output_time_step;

				std::vector<double> phaseOneProt;
				std::vector<int> positionOneProt;
				std::vector<TrajectoryEvent> trajectoryOneProt;
				phaseOneProt.reserve(expected_phase_samples);
				if (cube_s.x() != 0.0) {
					positionOneProt.reserve(expected_phase_samples);
				}
				if (record_positions) {
					trajectoryOneProt.reserve(expected_phase_samples);
				}

			//Initiate Gaussian distribution for background T2
				RandNumberGaussian rand_phase_back(
					0,
					standarddeviation_back,
					proton_random_seed(
						simulation_random_seed,
						index,
						ProtonRandomStream::RelaxationOutside));
				RandNumberGaussian rand_phase_cell(
					0,
					standarddeviation_cell,
					proton_random_seed(
						simulation_random_seed,
						index,
						ProtonRandomStream::RelaxationInside));
				const double substep_dt = time_step / SimConsts::VARSTEP_FACTOR;
				RandNumberGaussian rand_phase_back_sub(
					0,
					background_T2 != 0. ? calc_std_T2(substep_dt, background_T2) : 0.,
					proton_random_seed(
						simulation_random_seed,
						index,
						ProtonRandomStream::SubstepRelaxationOutside));
				RandNumberGaussian rand_phase_cell_sub(
					0,
					incell_T2 != 0. ? calc_std_T2(substep_dt, incell_T2) : 0.,
					proton_random_seed(
						simulation_random_seed,
						index,
						ProtonRandomStream::SubstepRelaxationInside));

			//Initial random placement
				replaceProton(protonInstances, ionpInstances_orig_local, static_ionp_index, static_cell_index, coating_ionp_index, coating_state_known, cell_origins, time_step, index, j, time, sample_region, RandomNormalX, RandomNormalY, RandomNormalZ, RandomDiff, record_positions ? &trajectoryOneProt : nullptr, ionp_translation_for_calc);

				recordPositionInCube(positionOneProt, protonInstances[index]);
				recordTrajectory(trajectoryOneProt, protonInstances[index], time, 0);

			//Start Signal at t=0
			phaseOneProt.push_back(phase);
			if (!index) {
				timeV.push_back(time);
			}

		//================================= End Initialization =================================//

		//================================= Start Proton (IONP) propagation =================================//
		while(time < echo_time){

			//Reset Time step size for every iteration - for if a variable Time step is used
			double dt = time_step;
			bool smallerTimeStep = false;

			if (leapfrog % 2 == 0){
				if (use_variable_steps) {
					if (has_near_ionp_for_variable_step_translated(
							protonInstances[index].xt(),
							magnetic_field_sources,
							*magnetic_field_source_index,
							ionp_check_mode,
							ionp_translation_for_calc)) {
							dt /= SimConsts::VARSTEP_FACTOR;
							smallerTimeStep = true;
								j++;
								time += time_step;
								for (unsigned int s = 1; s <= SC_UI(SimConsts::VARSTEP_FACTOR); ++s){
									const double sub_time = time - time_step + s * dt;
									replaceProton(protonInstances, ionpInstances_orig_local, static_ionp_index, static_cell_index, coating_ionp_index, coating_state_known, cell_origins, dt, index, j, sub_time, sample_region, RandomNormalX, RandomNormalY, RandomNormalZ, RandomDiff, record_positions ? &trajectoryOneProt : nullptr, ionp_translation_for_calc);
									double bfield_tot = compute_field_for_step(
										protonInstances[index],
										ionp_translation_for_calc,
										smallerTimeStep);
									double agar_relax = protonInstances[index].inCell() ? (incell_T2 != 0. ? rand_phase_cell_sub() : 0.) : (background_T2 != 0. ? rand_phase_back_sub() : 0.);

								update_phase(phase, bfield_tot, dt, agar_relax, j, smallerTimeStep, refocusing_schedule);
							}

							// after subdivision, apply refocusing if needed
							apply_spin_echo_refocusing(phase, j, refocusing_schedule);

								recordPhaseSample(phaseOneProt, positionOneProt, trajectoryOneProt, protonInstances[index], index, phase, time, next_record_time);
					}
				}

					if (smallerTimeStep == false){
						time += time_step;
						j++;
						replaceProton(protonInstances, ionpInstances_orig_local, static_ionp_index, static_cell_index, coating_ionp_index, coating_state_known, cell_origins, time_step, index, j, time, sample_region, RandomNormalX, RandomNormalY, RandomNormalZ, RandomDiff, record_positions ? &trajectoryOneProt : nullptr, ionp_translation_for_calc);
					//============= Calculation of the Bfield at proton position =================================//
						double bfield_tot = compute_field_for_step(
							protonInstances[index],
							ionp_translation_for_calc,
						false);
					double agar_relax = protonInstances[index].inCell() ? (incell_T2 != 0. ? rand_phase_cell() : 0.) : (background_T2 != 0. ? rand_phase_back() : 0.);

					update_phase(phase, bfield_tot, time_step, agar_relax, j, smallerTimeStep, refocusing_schedule);

						recordPhaseSample(phaseOneProt, positionOneProt, trajectoryOneProt, protonInstances[index], index, phase, time, next_record_time);
					}
			}
			else{
				//IONP movement can only be simulated if the B-Field is calculated directly
				if (ionpMovesLocal){
					advance_ionp_translation();
					coating_state_known = false;
				}
			}

			if (ionpMovesLocal){
				leapfrog ++;
			}

		} // end time loop

			phaseVector[index] = std::move(phaseOneProt);
				protonPosition[index] = std::move(positionOneProt);
				if (record_positions) {
					protonTrajectory[index] = std::move(trajectoryOneProt);
				}

			};

			auto simulate_single_bound_proton = [&](unsigned int index, unsigned int ionp_index, const std::vector<Ionp>& ionpInstances_orig_local, bool ionpMovesLocal) {
				const bool use_variable_steps = should_use_variable_time_steps(var_Steps, importance_sampling.enabled, ProtonSampleRegion::Bound);
				RandomList random_unit(
					0.,
					1.,
					false,
					proton_random_seed(
						simulation_random_seed,
						index,
						ProtonRandomStream::DiffusionAndBoundaries));
				place_proton_random_sphere(
					protonInstances[index],
					ionpInstances_orig_local[ionp_index].x0().x(),
					ionpInstances_orig_local[ionp_index].x0().y(),
					ionpInstances_orig_local[ionp_index].x0().z(),
					ionpInstances_orig_local[ionp_index].radius(),
					ionpInstances_orig_local[ionp_index].outer_radius(),
					random_unit);

			Point3D ionp_translation;
			Point3D ionp_reference_position =
				ionpInstances_orig_local.empty()
					? Point3D()
					: ionpInstances_orig_local.front().xt();
			const Point3D* ionp_translation_for_calc =
				ionpMovesLocal ? &ionp_translation : nullptr;
			auto advance_ionp_translation = [&]() {
				const Point3D previous_reference_position = ionp_reference_position;
				ionp_reference_position.set(
					ionp_reference_position.x() + ionp_step_displacement.x(),
					ionp_reference_position.y() + ionp_step_displacement.y(),
					ionp_reference_position.z() + ionp_step_displacement.z());
				ionp_translation.set(
					ionp_translation.x() +
						(ionp_reference_position.x() - previous_reference_position.x()),
					ionp_translation.y() +
						(ionp_reference_position.y() - previous_reference_position.y()),
					ionp_translation.z() +
						(ionp_reference_position.z() - previous_reference_position.z()));
			};

		//================================= Start Initialization =================================//

		//Initiate Proton specific variables
		double phase = 0., time = 0.;
		unsigned int j = 0, leapfrog = 0;
		double next_record_time = output_time_step;

			std::vector<double> phaseOneProt;
			std::vector<int> positionOneProt;
			std::vector<TrajectoryEvent> trajectoryOneProt;
			phaseOneProt.reserve(expected_phase_samples);
			if (cube_s.x() != 0.0) {
				positionOneProt.reserve(expected_phase_samples);
			}
			if (record_positions) {
				trajectoryOneProt.reserve(expected_phase_samples);
			}

		//Initiate Gaussian distribution for background T2
			RandNumberGaussian rand_phase_back(
				0,
				standarddeviation_back,
				proton_random_seed(
					simulation_random_seed,
					index,
					ProtonRandomStream::RelaxationOutside));
			RandNumberGaussian rand_phase_cell(
				0,
				standarddeviation_cell,
				proton_random_seed(
					simulation_random_seed,
					index,
					ProtonRandomStream::RelaxationInside));
			const double substep_dt = time_step / SimConsts::VARSTEP_FACTOR;
			RandNumberGaussian rand_phase_back_sub(
				0,
				background_T2 != 0. ? calc_std_T2(substep_dt, background_T2) : 0.,
				proton_random_seed(
					simulation_random_seed,
					index,
					ProtonRandomStream::SubstepRelaxationOutside));
			RandNumberGaussian rand_phase_cell_sub(
				0,
				incell_T2 != 0. ? calc_std_T2(substep_dt, incell_T2) : 0.,
				proton_random_seed(
					simulation_random_seed,
					index,
					ProtonRandomStream::SubstepRelaxationInside));
			const bool static_bound_near_ionp =
				!ionpMovesLocal &&
				use_variable_steps &&
					has_near_ionp_for_variable_step_translated(
						protonInstances[index].xt(),
						magnetic_field_sources,
						*magnetic_field_source_index,
						ionp_check_mode,
						nullptr);
			const double static_bound_bfield =
				!ionpMovesLocal
						? compute_field_for_step(
							  protonInstances[index],
							  nullptr,
						  false)
					: 0.0;

		recordPositionInCube(positionOneProt, protonInstances[index]);
		recordTrajectory(trajectoryOneProt, protonInstances[index], time, 0);

		//Start Signal at t=0
		phaseOneProt.push_back(phase);
		if (!index) {
			timeV.push_back(time);
		}

		//================================= End Initialization =================================//

		//================================= Start Proton (IONP) propagation =================================//
		while(time < echo_time){

			//Reset Time step size for every iteration - for if a variable Time step is used
			double dt = time_step;
			bool smallerTimeStep = false;

			if (leapfrog % 2 == 0){
				if (use_variable_steps) {
					const bool near_ionp =
						ionpMovesLocal
							? has_near_ionp_for_variable_step_translated(
								  protonInstances[index].xt(),
								  magnetic_field_sources,
								  *magnetic_field_source_index,
								  ionp_check_mode,
								  ionp_translation_for_calc)
							: static_bound_near_ionp;
					if (near_ionp) {
							dt /= SimConsts::VARSTEP_FACTOR;
							smallerTimeStep = true;
							j++;
							time += time_step;
							const double bfield_tot =
								ionpMovesLocal
									? compute_field_for_step(
										  protonInstances[index],
										  ionp_translation_for_calc,
										  smallerTimeStep)
									: static_bound_bfield;
							for (unsigned int s = 1; s <= SC_UI(SimConsts::VARSTEP_FACTOR); ++s){
								double agar_relax = protonInstances[index].inCell() ? (incell_T2 != 0. ? rand_phase_cell_sub() : 0.) : (background_T2 != 0. ? rand_phase_back_sub() : 0.);

								update_phase(phase, bfield_tot, dt, agar_relax, j, smallerTimeStep, refocusing_schedule);
							}

							// after subdivision, apply refocusing if needed
							apply_spin_echo_refocusing(phase, j, refocusing_schedule);

								recordPhaseSample(phaseOneProt, positionOneProt, trajectoryOneProt, protonInstances[index], index, phase, time, next_record_time);
					}
				}

				if (smallerTimeStep == false){
					time += time_step;
					j++;
					//============= Calculation of the Bfield at proton position =================================//
					const double bfield_tot =
						ionpMovesLocal
							? compute_field_for_step(
								  protonInstances[index],
								  ionp_translation_for_calc,
								  false)
							: static_bound_bfield;
					double agar_relax = protonInstances[index].inCell() ? (incell_T2 != 0. ? rand_phase_cell() : 0.) : (background_T2 != 0. ? rand_phase_back() : 0.);

					update_phase(phase, bfield_tot, time_step, agar_relax, j, smallerTimeStep, refocusing_schedule);

						recordPhaseSample(phaseOneProt, positionOneProt, trajectoryOneProt, protonInstances[index], index, phase, time, next_record_time);
					}
			}
			else{
				//IONP movement can only be simulated if the B-Field is calculated directly
				if (ionpMovesLocal){
					advance_ionp_translation();
					protonInstances[index].diffuse(ionp_step_displacement.x(),
												   ionp_step_displacement.y(),
												   ionp_step_displacement.z());
				}
			}

			if (ionpMovesLocal){
				leapfrog ++;
			}

		} // end time loop

			phaseVector[index] = std::move(phaseOneProt);
			protonPosition[index] = std::move(positionOneProt);
			if (record_positions) {
				protonTrajectory[index] = std::move(trajectoryOneProt);
			}

		};

	//================================= Dephasing of free protons =================================//

		const auto proposal_propagation_start =
			std::chrono::steady_clock::now();
		#pragma omp parallel for schedule(dynamic, 1)
		for (unsigned int i = 0; i < SC_UI(protons_n); ++i) {
			simulate_single_proton(i, ionpInstances_orig, ionpMoves);
		}

	//================================= Now dephasing of bound protons at the shell of the IONPs =================================//

	if (!bound_jobs.empty()) {
		std::cout << "There are bound protons" << std::endl;
	}

			#pragma omp parallel for schedule(dynamic, 1)
			for (std::size_t bound_index = 0;
				 bound_index < bound_sample_count;
				 ++bound_index) {
				const unsigned int index =
					SC_UI(protons_n) + SC_UI(bound_index);
				simulate_single_bound_proton(
					index,
					SC_UI(bound_jobs[bound_index].ionp_index),
					ionpInstances_orig,
					ionpMoves);
			}
		if (overlapping_proposals) {
			const std::chrono::duration<double> proposal_propagation_elapsed =
				std::chrono::steady_clock::now() - proposal_propagation_start;
			std::cout << "Overlapping proposal propagation seconds: "
				<< proposal_propagation_elapsed.count() << std::endl;
		}

		std::cout << "Magnetization calculation start!" << std::endl;

	if ( phaseVector.empty() ){
		return;
	}
	// Propagation writes independent per-proton histories. Validate their
	// common time axis before summing complex signals in the existing order.
	const std::size_t phase_time_points = phaseVector.front().size();
	for (unsigned int i = 0; i < all_proton_num; ++i) {
		if (phaseVector[i].size() != phase_time_points) {
			throw std::runtime_error("Phase history length mismatch before magnetization for proton " + std::to_string(i) +
									 ": got " + std::to_string(phaseVector[i].size()) +
									 ", expected " + std::to_string(phase_time_points));
		}
	}
	if (timeV.size() != phase_time_points) {
		throw std::runtime_error("Time vector length mismatch before magnetization: got " + std::to_string(timeV.size()) +
								 ", expected " + std::to_string(phase_time_points));
	}

	ImportanceRegionSignals regionalComplexSignals;
	ImportanceProposalSignalDiagnostics proposalComplexSignals;
	signalVector.clear();
	const auto proposal_aggregation_start =
		std::chrono::steady_clock::now();
	if (overlapping_proposals) {
		OverlappingWeightedSignalAggregation aggregation =
			aggregate_overlapping_weighted_signals(
				phaseVector,
				protonProposalSources,
				protonImportanceWeights);
		proposalComplexSignals.means = std::move(aggregation.proposal_means);
		proposalComplexSignals.contributions =
			std::move(aggregation.proposal_contributions);

		std::vector<std::complex<double>> bound_intracellular_signal(
			phase_time_points, std::complex<double>(0., 0.));
		std::vector<std::complex<double>> bound_extracellular_signal(
			phase_time_points, std::complex<double>(0., 0.));
		for (std::size_t proton_index = SC_UI(protons_n);
			 proton_index < phaseVector.size();
			 ++proton_index) {
			std::vector<std::complex<double>>& bound_signal =
				protonSampleRegions[proton_index] ==
					ProtonSampleRegion::BoundIntracellular
					? bound_intracellular_signal
					: bound_extracellular_signal;
			const std::size_t count =
				protonSampleRegions[proton_index] ==
					ProtonSampleRegion::BoundIntracellular
					? bound_intracellular_sample_count
					: bound_extracellular_sample_count;
			if (count == 0) {
				continue;
			}
			for (std::size_t time_index = 0;
				 time_index < phase_time_points;
				 ++time_index) {
				bound_signal[time_index] +=
					std::exp(i_compl * phaseVector[proton_index][time_index]) /
					static_cast<double>(count);
			}
		}

		double bound_intracellular_weight = 0.;
		double bound_extracellular_weight = 0.;
		if (importance_sampling.split_bound_configured) {
			bound_intracellular_weight =
				bound_intracellular_sample_count > 0
					? importance_sampling.weight_bound_intracellular
					: 0.;
			bound_extracellular_weight =
				bound_extracellular_sample_count > 0
					? importance_sampling.weight_bound_extracellular
					: 0.;
		}
		else {
			const std::size_t total_bound =
				bound_intracellular_sample_count +
				bound_extracellular_sample_count;
			if (total_bound > 0) {
				bound_intracellular_weight = importance_sampling.weight_bound *
					static_cast<double>(bound_intracellular_sample_count) /
					static_cast<double>(total_bound);
				bound_extracellular_weight = importance_sampling.weight_bound *
					static_cast<double>(bound_extracellular_sample_count) /
					static_cast<double>(total_bound);
			}
		}

		const double free_weight =
			importance_sampling.free_water_volume_fraction;
		const double total_compartment_weight =
			free_weight + bound_intracellular_weight +
			bound_extracellular_weight;
		if (!(free_weight > 0.) || !(total_compartment_weight > 0.)) {
			throw std::runtime_error(
				"Overlapping importance-sampling compartment weights must be positive");
		}
		std::cout << "Overlapping normalized compartment weights:" << std::endl;
		std::cout << "  free: "
			<< free_weight / total_compartment_weight << std::endl;
		if (bound_intracellular_weight > 0.) {
			std::cout << "  bound_intracellular: "
				<< bound_intracellular_weight / total_compartment_weight
				<< std::endl;
		}
		if (bound_extracellular_weight > 0.) {
			std::cout << "  bound_extracellular: "
				<< bound_extracellular_weight / total_compartment_weight
				<< std::endl;
		}
		signalVector.assign(
			phase_time_points, std::complex<double>(0., 0.));
		for (std::size_t time_index = 0;
			 time_index < phase_time_points;
			 ++time_index) {
			signalVector[time_index] =
				(free_weight * aggregation.combined_signal[time_index] +
				 bound_intracellular_weight *
					 bound_intracellular_signal[time_index] +
				 bound_extracellular_weight *
					 bound_extracellular_signal[time_index]) /
				total_compartment_weight;
		}
	}
	else if (importance_sampling.enabled) {
		WeightedSignalAggregation aggregation = aggregate_weighted_signals(phaseVector, protonSampleRegions, regionSampleCounts, normalizedWeights);
		if (regionSignalVectors || regionSignalVectorsMSE) {
			regionalComplexSignals = std::move(aggregation.region_signals);
		}
		signalVector = std::move(aggregation.combined_signal);
	}
	else {
		signalVector.assign(phaseVector[0].size(), std::complex<double>(0., 0.));
		for (unsigned int i = 0; i < all_proton_num; ++i) {
			for (unsigned int j = 0; j < phaseVector[0].size(); ++j){
				signalVector[j] += std::exp(i_compl * phaseVector[i][j]) / SC_D(all_proton_num);
			}
		}
	}
	if (overlapping_proposals) {
		const std::chrono::duration<double> proposal_aggregation_elapsed =
			std::chrono::steady_clock::now() - proposal_aggregation_start;
		std::cout << "Overlapping proposal aggregation seconds: "
			<< proposal_aggregation_elapsed.count() << std::endl;
	}
	build_mse_signal(signalVector, timeV, sequence_num, echo_time, echo_spacing, time_step, signalVectorMSE, timeVMSE);
	if (importance_sampling.enabled && !overlapping_proposals &&
		(regionSignalVectors || regionSignalVectorsMSE)) {
		for (std::size_t region_index = 0;
			 region_index < kSampleRegionCount;
			 ++region_index) {
			if (regionSampleCounts[region_index] == 0) {
				continue;
			}
			if (regionSignalVectorsMSE) {
				std::vector<double> regionalTimeVMSE;
				build_mse_signal(
					regionalComplexSignals[region_index],
					timeV,
					sequence_num,
					echo_time,
					echo_spacing,
					time_step,
					(*regionSignalVectorsMSE)[region_index],
					regionalTimeVMSE);
				if (regionalTimeVMSE != timeVMSE) {
					throw std::runtime_error(
						"Regional MSE time points do not match the combined signal");
				}
			}
			if (regionSignalVectors) {
				(*regionSignalVectors)[region_index] =
					std::move(regionalComplexSignals[region_index]);
			}
		}
	}
	if (overlapping_proposals &&
		(proposalSignalDiagnostics || proposalSignalDiagnosticsMSE)) {
		for (std::size_t proposal_index = 0;
			 proposal_index < kFreeSamplingProposalCount;
			 ++proposal_index) {
			if (proposalComplexSignals.means[proposal_index].empty()) {
				continue;
			}
			if (proposalSignalDiagnosticsMSE) {
				std::vector<double> proposal_time_mse;
				build_mse_signal(
					proposalComplexSignals.means[proposal_index],
					timeV,
					sequence_num,
					echo_time,
					echo_spacing,
					time_step,
					proposalSignalDiagnosticsMSE->means[proposal_index],
					proposal_time_mse);
				if (proposal_time_mse != timeVMSE) {
					throw std::runtime_error(
						"Proposal-mean MSE time points do not match the combined signal");
				}
				proposal_time_mse.clear();
				build_mse_signal(
					proposalComplexSignals.contributions[proposal_index],
					timeV,
					sequence_num,
					echo_time,
					echo_spacing,
					time_step,
					proposalSignalDiagnosticsMSE->contributions[proposal_index],
					proposal_time_mse);
				if (proposal_time_mse != timeVMSE) {
					throw std::runtime_error(
						"Proposal-contribution MSE time points do not match the combined signal");
				}
			}
			if (proposalSignalDiagnostics) {
				proposalSignalDiagnostics->means[proposal_index] =
					std::move(proposalComplexSignals.means[proposal_index]);
				proposalSignalDiagnostics->contributions[proposal_index] =
					std::move(proposalComplexSignals.contributions[proposal_index]);
			}
		}
	}

	//get the Magnitude, so calculate abs of complex signal
	std::for_each(signalVector.begin(), signalVector.end(), [](std::complex<double>& c) { c = std::abs(c); });
	std::for_each(signalVectorMSE.begin(), signalVectorMSE.end(), [](std::complex<double>& c) { c = std::abs(c); });
	std::cout << "Magnetization calculated." << std::endl;

	unsigned int inCell = 0;
	for (unsigned int k = 0; k < SC_UI(protonInstances.size()); ++k) {
		if (protonInstances[k].inCell() == true){
			inCell++;
		}
	}
	std::cout << "Protons inside Cell: " << inCell << std::endl;

	if (!importance_sampling.enabled && cube_s.x() != 0. ){
		std::vector<std::complex<double>> signalShellProtons(phaseVector[0].size(), 0.);
		long int numberShellProtons = 0;
		for (const auto& row : protonPosition) {
			numberShellProtons += std::count(row.begin(), row.end(), 1);
		}

		std::cout << "Number of Protons*TimeStep in defined Voxel: " << numberShellProtons << std::endl;

		for (unsigned int i = 0; i < all_proton_num; ++i) {
			for (unsigned int j = 0; j < phaseVector[0].size(); ++j){
				if ( protonPosition[i][j] == 1 ){
					signalShellProtons[j] += std::exp(i_compl * phaseVector[i][j])/SC_D(all_proton_num);
				}
			}
		}

		std::for_each(signalShellProtons.begin(), signalShellProtons.end(), [](std::complex<double>& c) { c = std::abs(c); });

		writeData( signalShellProtons, timeV, ionpInstances_orig );
	}

	if (record_positions) {
		writeProtonPositions(protonTrajectory, ionpInstances_orig.size());
	}

}
