/**
 * @file ImportanceSampling.cpp
 *
 * Importance-sampling support: classify physical IONP populations, estimate
 * accessible proposal volumes, validate allocations and combine complex
 * signals. Sampling changes where starts are drawn and how they are weighted;
 * MRsequence still applies the same motion and phase evolution to each spin.
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include "RandNumberBetween.hpp"

#include "ImportanceSampling.hpp"

namespace {

constexpr double kGeometryTolerance = 1e-12;
constexpr double kFourThirdsPi = 4.1887902047863909846;

std::size_t available_cell_count(const std::vector<Point3D>& cell_origins,
								 double cell_r) {
	if (cell_r <= 0.) {
		return 0;
	}
	return cell_origins.empty() ? 1u : cell_origins.size();
}

Point3D cell_center_at(const std::vector<Point3D>& cell_origins,
					   const Point3D& cell_origin,
					   std::size_t index) {
	return cell_origins.empty() ? cell_origin : cell_origins[index];
}

double nearest_cell_center_distance(const Point3D& point,
									const std::vector<Point3D>& cell_origins,
									const Point3D& cell_origin,
									double cell_r) {
	const std::size_t cell_count = available_cell_count(cell_origins, cell_r);
	if (cell_count == 0) {
		return std::numeric_limits<double>::infinity();
	}

	double nearest_distance = std::numeric_limits<double>::infinity();
	for (std::size_t cell_index = 0; cell_index < cell_count; ++cell_index) {
		nearest_distance = std::min(
			nearest_distance,
			distance(point, cell_center_at(cell_origins, cell_origin, cell_index)));
	}
	return nearest_distance;
}

CellSampleRegion classify_cell_region(const Point3D& point,
									  const std::vector<Point3D>& cell_origins,
									  const Point3D& cell_origin,
									  double cell_r,
									  double near_cutoff) {
	const double nearest_distance =
		nearest_cell_center_distance(point, cell_origins, cell_origin, cell_r);
	if (nearest_distance <= cell_r) {
		return CellSampleRegion::Inside;
	}
	if (near_cutoff > 0. && nearest_distance <= near_cutoff) {
		return CellSampleRegion::Near;
	}
	return CellSampleRegion::Far;
}

double nearest_surface_distance(const Point3D& point,
								const std::vector<Ionp>& ionps,
								const std::vector<std::size_t>& indices) {
	double nearest = std::numeric_limits<double>::infinity();
	for (std::size_t index : indices) {
		nearest = std::min(
			nearest,
			distance(point, ionps[index].xt()) - ionps[index].outer_radius());
	}
	return nearest;
}

IonpDistanceRegion classify_ionp_distance(double surface_distance,
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

const char* cell_region_name(CellSampleRegion region) {
	switch (region) {
		case CellSampleRegion::Inside: return "inside";
		case CellSampleRegion::Near: return "near";
		case CellSampleRegion::Far: return "far";
	}
	return "unknown";
}

const char* ionp_region_name(IonpDistanceRegion region) {
	switch (region) {
		case IonpDistanceRegion::Near: return "near";
		case IonpDistanceRegion::Intermediate: return "intermediate";
		case IonpDistanceRegion::Far: return "far";
	}
	return "unknown";
}

double shell_inner_radius(const Ionp& ionp,
						  IonpDistanceRegion region,
						  double near_cutoff) {
	if (region == IonpDistanceRegion::Near) {
		return ionp.outer_radius();
	}
	return ionp.outer_radius() + near_cutoff;
}

double shell_outer_radius(const Ionp& ionp,
						  IonpDistanceRegion region,
						  double near_cutoff,
						  double far_cutoff) {
	if (region == IonpDistanceRegion::Near) {
		return ionp.outer_radius() + near_cutoff;
	}
	return ionp.outer_radius() + far_cutoff;
}

double sphere_shell_volume(double inner_radius, double outer_radius) {
	return kFourThirdsPi *
		(outer_radius * outer_radius * outer_radius -
		 inner_radius * inner_radius * inner_radius);
}

Point3D random_point_in_spherical_shell(
	const Point3D& center,
	double inner,
	double outer,
	RandNumberBetween& random_unit) {
	const double inner_cubed = inner * inner * inner;
	const double outer_cubed = outer * outer * outer;
	const double radius =
		std::cbrt(inner_cubed + random_unit() * (outer_cubed - inner_cubed));
	const double theta = 2.0 * 3.14159265358979323846 * random_unit();
	const double phi = std::acos(1.0 - 2.0 * random_unit());
	return Point3D(
		center.x() + radius * std::sin(phi) * std::cos(theta),
		center.y() + radius * std::sin(phi) * std::sin(theta),
		center.z() + radius * std::cos(phi));
}

Point3D random_point_in_shell(const Ionp& ionp,
							  IonpDistanceRegion region,
							  double near_cutoff,
							  double far_cutoff,
							  RandNumberBetween& random_unit) {
	return random_point_in_spherical_shell(
		ionp.xt(),
		shell_inner_radius(ionp, region, near_cutoff),
		shell_outer_radius(ionp, region, near_cutoff, far_cutoff),
		random_unit);
}

struct ShellFamily {
	IonpPopulation population;
	IonpDistanceRegion region;
	const std::vector<std::size_t>* particle_indices;
	double near_cutoff;
	double far_cutoff;
	double total_shell_volume;
	double proposal_probability;
	std::vector<double> cumulative_shell_volumes;
};

struct CellShellFamily {
	CellSampleRegion region;
	double inner_radius;
	double outer_radius;
	double total_shell_volume;
	double proposal_probability;
};

bool any_active_region_for_axis(const ImportanceSamplingConfig& config,
								IonpPopulation population,
								IonpDistanceRegion distance_region) {
	for (std::size_t index = 0; index < kFreeSampleRegionCount; ++index) {
		if (config.joint_protons[index] <= 0) {
			continue;
		}
		const JointSampleRegion joint =
			decode_joint_sample_region(static_cast<ProtonSampleRegion>(index));
		const IonpDistanceRegion selected =
			population == IonpPopulation::Intracellular
				? joint.intracellular_ionp
				: joint.extracellular_ionp;
		if (selected == distance_region) {
			return true;
		}
	}
	return false;
}

struct ParticleGridKey {
	int x;
	int y;
	int z;

	bool operator==(const ParticleGridKey& other) const {
		return x == other.x && y == other.y && z == other.z;
	}
};

struct ParticleGridKeyHash {
	std::size_t operator()(const ParticleGridKey& key) const {
		std::size_t hash = static_cast<std::size_t>(key.x);
		hash ^= static_cast<std::size_t>(key.y) +
			0x9e3779b9u + (hash << 6u) + (hash >> 2u);
		hash ^= static_cast<std::size_t>(key.z) +
			0x9e3779b9u + (hash << 6u) + (hash >> 2u);
		return hash;
	}
};

class ImportanceParticleGrid {
public:
	ImportanceParticleGrid(const std::vector<Ionp>& ionps, double cell_size)
		: m_cell_size(cell_size) {
		if (m_cell_size <= 0.) {
			return;
		}
		for (std::size_t index = 0; index < ionps.size(); ++index) {
			m_cells[key(ionps[index].xt())].push_back(index);
		}
	}

	void candidates(const Point3D& point,
					double search_radius,
					std::vector<std::size_t>& result) const {
		result.clear();
		if (m_cell_size <= 0. || search_radius < 0.) {
			return;
		}
		const ParticleGridKey minimum = key(Point3D(
			point.x() - search_radius,
			point.y() - search_radius,
			point.z() - search_radius));
		const ParticleGridKey maximum = key(Point3D(
			point.x() + search_radius,
			point.y() + search_radius,
			point.z() + search_radius));
		for (int z = minimum.z; z <= maximum.z; ++z) {
			for (int y = minimum.y; y <= maximum.y; ++y) {
				for (int x = minimum.x; x <= maximum.x; ++x) {
					const auto found =
						m_cells.find(ParticleGridKey{x, y, z});
					if (found != m_cells.end()) {
						result.insert(
							result.end(),
							found->second.begin(),
							found->second.end());
					}
				}
			}
		}
	}

private:
	ParticleGridKey key(const Point3D& point) const {
		return {
			static_cast<int>(std::floor(point.x() / m_cell_size)),
			static_cast<int>(std::floor(point.y() / m_cell_size)),
			static_cast<int>(std::floor(point.z() / m_cell_size))
		};
	}

	double m_cell_size;
	std::unordered_map<
		ParticleGridKey,
		std::vector<std::size_t>,
		ParticleGridKeyHash> m_cells;
};

class ImportanceCellGrid {
public:
	ImportanceCellGrid(const std::vector<Point3D>& cell_origins,
					   const Point3D& fallback_origin,
					   double cell_radius,
					   double near_cutoff)
		: m_cell_radius(cell_radius),
		  m_near_cutoff(near_cutoff),
		  m_cell_size(std::max(cell_radius, near_cutoff)) {
		const std::size_t count =
			available_cell_count(cell_origins, cell_radius);
		m_centers.reserve(count);
		for (std::size_t index = 0; index < count; ++index) {
			m_centers.push_back(
				cell_center_at(cell_origins, fallback_origin, index));
		}
		if (m_cell_size <= 0. || m_centers.size() <= 1) {
			return;
		}
		for (std::size_t index = 0; index < m_centers.size(); ++index) {
			m_cells[key(m_centers[index])].push_back(index);
		}
	}

	CellSampleRegion classify(const Point3D& point) const {
		if (m_centers.empty()) {
			return CellSampleRegion::Far;
		}
		double nearest = std::numeric_limits<double>::infinity();
		if (m_centers.size() == 1) {
			nearest = distance(point, m_centers.front());
		}
		else {
			const ParticleGridKey center_key = key(point);
			for (int z = center_key.z - 1; z <= center_key.z + 1; ++z) {
				for (int y = center_key.y - 1; y <= center_key.y + 1; ++y) {
					for (int x = center_key.x - 1; x <= center_key.x + 1; ++x) {
						const auto found =
							m_cells.find(ParticleGridKey{x, y, z});
						if (found == m_cells.end()) {
							continue;
						}
						for (std::size_t index : found->second) {
							nearest = std::min(
								nearest,
								distance(point, m_centers[index]));
						}
					}
				}
			}
		}
		if (nearest <= m_cell_radius) {
			return CellSampleRegion::Inside;
		}
		if (m_near_cutoff > 0. && nearest <= m_near_cutoff) {
			return CellSampleRegion::Near;
		}
		return CellSampleRegion::Far;
	}

	std::size_t coverage(const Point3D& point,
					 double inner_radius,
					 double outer_radius) const {
		if (m_centers.empty() || outer_radius <= 0.) {
			return 0;
		}

		std::size_t result = 0;
		auto consider = [&](std::size_t index) {
			const double radius = distance(point, m_centers[index]);
			if (radius >= inner_radius && radius <= outer_radius) {
				++result;
			}
		};
		if (m_centers.size() == 1 || m_cell_size <= 0.) {
			for (std::size_t index = 0; index < m_centers.size(); ++index) {
				consider(index);
			}
			return result;
		}

		const ParticleGridKey center_key = key(point);
		const int cell_range = std::max(
			1,
			static_cast<int>(std::ceil(outer_radius / m_cell_size)));
		for (int z = center_key.z - cell_range;
			 z <= center_key.z + cell_range;
			 ++z) {
			for (int y = center_key.y - cell_range;
				 y <= center_key.y + cell_range;
				 ++y) {
				for (int x = center_key.x - cell_range;
					 x <= center_key.x + cell_range;
					 ++x) {
					const auto found =
						m_cells.find(ParticleGridKey{x, y, z});
					if (found == m_cells.end()) {
						continue;
					}
					for (std::size_t index : found->second) {
						consider(index);
					}
				}
			}
		}
		return result;
	}

private:
	ParticleGridKey key(const Point3D& point) const {
		return {
			static_cast<int>(std::floor(point.x() / m_cell_size)),
			static_cast<int>(std::floor(point.y() / m_cell_size)),
			static_cast<int>(std::floor(point.z() / m_cell_size))
		};
	}

	double m_cell_radius;
	double m_near_cutoff;
	double m_cell_size;
	std::vector<Point3D> m_centers;
	std::unordered_map<
		ParticleGridKey,
		std::vector<std::size_t>,
		ParticleGridKeyHash> m_cells;
};

bool classify_joint_sample_region_candidates(
	const Point3D& point,
	CellSampleRegion cell_region,
	const std::vector<Ionp>& ionps,
	const IonpPopulationMap& populations,
	const std::vector<std::size_t>& candidates,
	const ImportanceSamplingConfig& config,
	ProtonSampleRegion& region) {
	double intracellular_distance = std::numeric_limits<double>::infinity();
	double extracellular_distance = std::numeric_limits<double>::infinity();
	for (std::size_t index : candidates) {
		const double surface_distance =
			distance(point, ionps[index].xt()) - ionps[index].outer_radius();
		if (populations.populations[index] == IonpPopulation::Intracellular) {
			intracellular_distance =
				std::min(intracellular_distance, surface_distance);
		}
		else {
			extracellular_distance =
				std::min(extracellular_distance, surface_distance);
		}
	}
	if (intracellular_distance < 0. || extracellular_distance < 0.) {
		return false;
	}
	region = make_joint_sample_region(
		cell_region,
		classify_ionp_distance(
			intracellular_distance,
			config.intracellular_ionp_sampling_enabled(),
			config.intracellular_ionp_near_cutoff,
			config.intracellular_ionp_far_cutoff),
		classify_ionp_distance(
			extracellular_distance,
			config.extracellular_ionp_sampling_enabled(),
			config.extracellular_ionp_near_cutoff,
			config.extracellular_ionp_far_cutoff));
	return true;
}

std::size_t shell_coverage(const Point3D& point,
						   const std::vector<Ionp>& ionps,
						   const IonpPopulationMap& populations,
						   const std::vector<std::size_t>& candidates,
						   const ShellFamily& family) {
	std::size_t coverage = 0;
	for (std::size_t index : candidates) {
		if (populations.populations[index] != family.population) {
			continue;
		}
		const double radius = distance(point, ionps[index].xt());
		const double inner = shell_inner_radius(
			ionps[index], family.region, family.near_cutoff);
		const double outer = shell_outer_radius(
			ionps[index],
			family.region,
			family.near_cutoff,
			family.far_cutoff);
		if (radius >= inner && radius <= outer) {
			++coverage;
		}
	}
	return coverage;
}

void validate_overlapping_fraction_bounds(
	const std::array<double, kFreeSamplingProposalCount>& fractions) {
	constexpr double tolerance = 1e-9;
	const double uniform = fractions[free_sampling_proposal_index(
		FreeSamplingProposal::Uniform)];
	if (!std::isfinite(uniform) || std::abs(uniform - 1.0) > tolerance) {
		throw std::runtime_error(
			"The overlapping uniform proposal volume fraction must equal one");
	}
	for (double fraction : fractions) {
		if (!std::isfinite(fraction) || fraction < 0. ||
			fraction > 1.0 + tolerance) {
			throw std::runtime_error(
				"Overlapping proposal volume fractions must be finite and lie in [0, 1]");
		}
	}
}

} // namespace

int ImportanceSamplingConfig::total_free_protons() const {
	validate_mode_exclusivity();
	if (overlapping_proposals_enabled()) {
		return total_proposal_protons();
	}
	const int joint_total =
		std::accumulate(joint_protons.begin(), joint_protons.end(), 0);
	if (joint_regions_configured || joint_total != 0) {
		return joint_total;
	}
	return protons_inside + protons_near + protons_far;
}

int ImportanceSamplingConfig::total_proposal_protons() const {
	validate_mode_exclusivity();
	long long total = 0;
	for (int count : proposal_protons) {
		if (count < 0) {
			throw std::runtime_error(
				"Overlapping importance-sampling proposal counts must be non-negative");
		}
		total += static_cast<long long>(count);
	}
	if (total > static_cast<long long>(std::numeric_limits<int>::max())) {
		throw std::runtime_error(
			"Total overlapping importance-sampling proton count exceeds the supported integer range");
	}
	return static_cast<int>(total);
}

bool ImportanceSamplingConfig::overlapping_proposals_enabled() const {
	return mode == ImportanceSamplingMode::OverlappingProposal;

}

void ImportanceSamplingConfig::validate_mode_exclusivity() const {
	const bool has_legacy_free_values =
		protons_inside != 0 || protons_near != 0 || protons_far != 0 ||
		weight_inside != 0. || weight_near != 0. || weight_far != 0. ||
		joint_regions_configured ||
		std::any_of(
			joint_protons.begin(),
			joint_protons.end(),
			[](int count) { return count != 0; }) ||
		std::any_of(
			joint_weights.begin(),
			joint_weights.end(),
			[](double weight) { return weight != 0.; });
	const bool has_proposal_values = std::any_of(
		proposal_protons.begin(),
		proposal_protons.end(),
		[](int count) { return count != 0; });
	if (mode == ImportanceSamplingMode::OverlappingProposal &&
		has_legacy_free_values) {
		throw std::runtime_error(
			"Overlapping and legacy/joint free importance-sampling controls cannot be mixed");
	}
	if (mode == ImportanceSamplingMode::LegacyJoint &&
		has_proposal_values) {
		throw std::runtime_error(
			"Overlapping proposal counts require OverlappingProposal importance-sampling mode");
	}
}

bool ImportanceSamplingConfig::intracellular_ionp_sampling_enabled() const {
	return intracellular_ionp_near_cutoff != 0. ||
		intracellular_ionp_far_cutoff != 0.;
}

bool ImportanceSamplingConfig::extracellular_ionp_sampling_enabled() const {
	return extracellular_ionp_near_cutoff != 0. ||
		extracellular_ionp_far_cutoff != 0.;
}

void ImportanceSamplingConfig::materialize_legacy_regions() {
	validate_mode_exclusivity();
	if (overlapping_proposals_enabled()) {
		return;
	}
	if (joint_regions_configured) {
		if (protons_inside != 0 || protons_near != 0 || protons_far != 0 ||
			weight_inside != 0. || weight_near != 0. || weight_far != 0.) {
			throw std::runtime_error(
				"Legacy IS_protons_inside/near/far or IS_weight_inside/near/far "
				"cannot be mixed with joint cell/IONP importance-sampling keys");
		}
		return;
	}
	if (std::any_of(
			joint_protons.begin(),
			joint_protons.end(),
			[](int count) { return count != 0; })) {
		// Legacy regions have already been materialized. In particular, keep
		// automatic weights that may have since been stored in joint_weights.
		return;
	}

	joint_protons.fill(0);
	joint_weights.fill(0.);
	const std::size_t inside_index = joint_sample_region_index(
		CellSampleRegion::Inside,
		IonpDistanceRegion::Far,
		IonpDistanceRegion::Far);
	const std::size_t near_index = joint_sample_region_index(
		CellSampleRegion::Near,
		IonpDistanceRegion::Far,
		IonpDistanceRegion::Far);
	const std::size_t far_index = joint_sample_region_index(
		CellSampleRegion::Far,
		IonpDistanceRegion::Far,
		IonpDistanceRegion::Far);
	joint_protons[inside_index] = protons_inside;
	joint_protons[near_index] = protons_near;
	joint_protons[far_index] = protons_far;
	joint_weights[inside_index] = weight_inside;
	joint_weights[near_index] = weight_near;
	joint_weights[far_index] = weight_far;
}

std::string free_sampling_proposal_name(FreeSamplingProposal proposal) {
	switch (proposal) {
		case FreeSamplingProposal::Uniform:
			return "uniform";
		case FreeSamplingProposal::CellInside:
			return "cell_inside";
		case FreeSamplingProposal::CellNear:
			return "cell_near";
		case FreeSamplingProposal::CellFar:
			return "cell_far";
		case FreeSamplingProposal::IntracellularIonpNear:
			return "intracellular_ionp_near";
		case FreeSamplingProposal::IntracellularIonpIntermediate:
			return "intracellular_ionp_intermediate";
		case FreeSamplingProposal::IntracellularIonpFar:
			return "intracellular_ionp_far";
		case FreeSamplingProposal::ExtracellularIonpNear:
			return "extracellular_ionp_near";
		case FreeSamplingProposal::ExtracellularIonpIntermediate:
			return "extracellular_ionp_intermediate";
		case FreeSamplingProposal::ExtracellularIonpFar:
			return "extracellular_ionp_far";
	}
	throw std::runtime_error("Invalid free importance-sampling proposal");
}

FreeSamplingProposal cell_sampling_proposal(CellSampleRegion region) {
	switch (region) {
		case CellSampleRegion::Inside:
			return FreeSamplingProposal::CellInside;
		case CellSampleRegion::Near:
			return FreeSamplingProposal::CellNear;
		case CellSampleRegion::Far:
			return FreeSamplingProposal::CellFar;
	}
	throw std::runtime_error("Invalid cell importance-sampling region");
}

FreeSamplingProposal ionp_sampling_proposal(
	IonpPopulation population,
	IonpDistanceRegion region) {
	if (population == IonpPopulation::Intracellular) {
		switch (region) {
			case IonpDistanceRegion::Near:
				return FreeSamplingProposal::IntracellularIonpNear;
			case IonpDistanceRegion::Intermediate:
				return FreeSamplingProposal::IntracellularIonpIntermediate;
			case IonpDistanceRegion::Far:
				return FreeSamplingProposal::IntracellularIonpFar;
		}
	}
	else if (population == IonpPopulation::Extracellular) {
		switch (region) {
			case IonpDistanceRegion::Near:
				return FreeSamplingProposal::ExtracellularIonpNear;
			case IonpDistanceRegion::Intermediate:
				return FreeSamplingProposal::ExtracellularIonpIntermediate;
			case IonpDistanceRegion::Far:
				return FreeSamplingProposal::ExtracellularIonpFar;
		}
	}
	throw std::runtime_error("Invalid IONP importance-sampling population or region");
}

std::array<FreeSamplingProposal, 3> marginal_sampling_proposals(
	const JointSampleRegion& membership) {
	return {
		cell_sampling_proposal(membership.cell),
		ionp_sampling_proposal(
			IonpPopulation::Intracellular,
			membership.intracellular_ionp),
		ionp_sampling_proposal(
			IonpPopulation::Extracellular,
			membership.extracellular_ionp)
	};
}

std::size_t sample_region_index(ProtonSampleRegion region) {
	return static_cast<std::size_t>(region);
}

ProtonSampleRegion make_joint_sample_region(
	CellSampleRegion cell_region,
	IonpDistanceRegion intracellular_ionp_region,
	IonpDistanceRegion extracellular_ionp_region) {
	return static_cast<ProtonSampleRegion>(
		joint_sample_region_index(
			cell_region,
			intracellular_ionp_region,
			extracellular_ionp_region));
}

JointSampleRegion decode_joint_sample_region(ProtonSampleRegion region) {
	const std::size_t index = sample_region_index(region);
	if (index >= kFreeSampleRegionCount) {
		throw std::runtime_error("Bound sample region has no joint free-region coordinates");
	}
	return {
		static_cast<CellSampleRegion>(index / 9u),
		static_cast<IonpDistanceRegion>((index % 9u) / 3u),
		static_cast<IonpDistanceRegion>(index % 3u)
	};
}

std::string sample_region_name(ProtonSampleRegion region) {
	if (region == ProtonSampleRegion::BoundIntracellular) {
		return "bound_intracellular_ionp";
	}
	if (region == ProtonSampleRegion::BoundExtracellular) {
		return "bound_extracellular_ionp";
	}
	const JointSampleRegion joint = decode_joint_sample_region(region);
	return std::string("cell_") + cell_region_name(joint.cell) +
		"__intracellular_ionp_" + ionp_region_name(joint.intracellular_ionp) +
		"__extracellular_ionp_" + ionp_region_name(joint.extracellular_ionp);
}

bool has_explicit_free_region_weights(const ImportanceSamplingConfig& config) {
	const double joint_total = std::accumulate(
		config.joint_weights.begin(), config.joint_weights.end(), 0.0);
	if (config.joint_regions_configured || joint_total != 0.) {
		return joint_total > 0.;
	}
	return config.weight_inside + config.weight_near + config.weight_far > 0.;
}

std::array<std::size_t, kSampleRegionCount> sample_region_counts(
	const ImportanceSamplingConfig& config,
	std::size_t bound_intracellular_sample_count,
	std::size_t bound_extracellular_sample_count) {
	std::array<std::size_t, kSampleRegionCount> counts{};
	for (std::size_t index = 0; index < kFreeSampleRegionCount; ++index) {
		counts[index] = static_cast<std::size_t>(config.joint_protons[index]);
	}
	counts[kBoundIntracellularSampleRegionIndex] =
		bound_intracellular_sample_count;
	counts[kBoundExtracellularSampleRegionIndex] =
		bound_extracellular_sample_count;
	return counts;
}

std::array<std::size_t, kSampleRegionCount> sample_region_counts(
	const ImportanceSamplingConfig& config,
	std::size_t bound_sample_count) {
	ImportanceSamplingConfig materialized = config;
	materialized.materialize_legacy_regions();
	return sample_region_counts(materialized, 0, bound_sample_count);
}

// A particle belongs to an intracellular population only when its outer
// sphere fits in a cell. Straddling a membrane is rejected instead of silently
// assigning the source from its centre alone.
IonpPopulationMap classify_ionp_populations(
	const std::vector<Ionp>& ionps,
	const std::vector<Point3D>& cell_origins,
	const Point3D& cell_origin,
	double cell_radius) {
	IonpPopulationMap result;
	result.populations.resize(ionps.size(), IonpPopulation::Extracellular);

	const std::size_t cell_count =
		available_cell_count(cell_origins, cell_radius);
	for (std::size_t ionp_index = 0; ionp_index < ionps.size(); ++ionp_index) {
		const Ionp& ionp = ionps[ionp_index];
		bool intracellular = false;
		bool overlaps_membrane = false;

		for (std::size_t cell_index = 0; cell_index < cell_count; ++cell_index) {
			const double center_distance = distance(
				ionp.xt(),
				cell_center_at(cell_origins, cell_origin, cell_index));
			if (center_distance + ionp.outer_radius() <=
				cell_radius + kGeometryTolerance) {
				intracellular = true;
				continue;
			}
			if (center_distance <
				cell_radius + ionp.outer_radius() - kGeometryTolerance) {
				overlaps_membrane = true;
			}
		}

		if (overlaps_membrane) {
			throw std::runtime_error(
				"IONP " + std::to_string(ionp_index) +
				" straddles a cell membrane when classified by its full outer radius");
		}
		if (intracellular) {
			result.populations[ionp_index] = IonpPopulation::Intracellular;
			result.intracellular_indices.push_back(ionp_index);
			if (ionp.coating() > 0.) {
				++result.coated_intracellular_count;
			}
		}
		else {
			result.extracellular_indices.push_back(ionp_index);
			if (ionp.coating() > 0.) {
				++result.coated_extracellular_count;
			}
		}
	}
	return result;
}

ProtonSampleRegion classify_sample_region(
	const Point3D& point,
	const std::vector<Point3D>& cell_origins,
	const Point3D& cell_origin,
	double cell_r,
	const ImportanceSamplingConfig& importance_sampling) {
	return make_joint_sample_region(
		classify_cell_region(
			point,
			cell_origins,
			cell_origin,
			cell_r,
			importance_sampling.cell_near_cutoff),
		IonpDistanceRegion::Far,
		IonpDistanceRegion::Far);
}

bool classify_joint_sample_region(
	const Point3D& point,
	const std::vector<Point3D>& cell_origins,
	const Point3D& cell_origin,
	double cell_r,
	const std::vector<Ionp>& ionps,
	const IonpPopulationMap& populations,
	const ImportanceSamplingConfig& importance_sampling,
	ProtonSampleRegion& region) {
	const double intracellular_distance =
		nearest_surface_distance(point, ionps, populations.intracellular_indices);
	const double extracellular_distance =
		nearest_surface_distance(point, ionps, populations.extracellular_indices);

	// All realized particle volumes are excluded from free-proton initialization,
	// including populations whose distance axis is not actively stratified.
	if (intracellular_distance < 0. || extracellular_distance < 0.) {
		return false;
	}

	region = make_joint_sample_region(
		classify_cell_region(
			point,
			cell_origins,
			cell_origin,
			cell_r,
			importance_sampling.cell_near_cutoff),
		classify_ionp_distance(
			intracellular_distance,
			importance_sampling.intracellular_ionp_sampling_enabled(),
			importance_sampling.intracellular_ionp_near_cutoff,
			importance_sampling.intracellular_ionp_far_cutoff),
		classify_ionp_distance(
			extracellular_distance,
			importance_sampling.extracellular_ionp_sampling_enabled(),
			importance_sampling.extracellular_ionp_near_cutoff,
			importance_sampling.extracellular_ionp_far_cutoff));
	return true;
}

bool classify_joint_sample_region_from_ionp_candidates(
	const Point3D& point,
	const std::vector<Point3D>& cell_origins,
	const Point3D& cell_origin,
	double cell_r,
	const std::vector<Ionp>& ionps,
	const IonpPopulationMap& populations,
	const std::vector<unsigned int>& ionp_candidates,
	const ImportanceSamplingConfig& importance_sampling,
	ProtonSampleRegion& region) {
	double intracellular_distance = std::numeric_limits<double>::infinity();
	double extracellular_distance = std::numeric_limits<double>::infinity();
	for (unsigned int candidate : ionp_candidates) {
		const std::size_t index = static_cast<std::size_t>(candidate);
		const double surface_distance =
			distance(point, ionps[index].xt()) - ionps[index].outer_radius();
		if (populations.populations[index] == IonpPopulation::Intracellular) {
			intracellular_distance =
				std::min(intracellular_distance, surface_distance);
		}
		else {
			extracellular_distance =
				std::min(extracellular_distance, surface_distance);
		}
	}
	if (intracellular_distance < 0. || extracellular_distance < 0.) {
		return false;
	}
	region = make_joint_sample_region(
		classify_cell_region(
			point,
			cell_origins,
			cell_origin,
			cell_r,
			importance_sampling.cell_near_cutoff),
		classify_ionp_distance(
			intracellular_distance,
			importance_sampling.intracellular_ionp_sampling_enabled(),
			importance_sampling.intracellular_ionp_near_cutoff,
			importance_sampling.intracellular_ionp_far_cutoff),
		classify_ionp_distance(
			extracellular_distance,
			importance_sampling.extracellular_ionp_sampling_enabled(),
			importance_sampling.extracellular_ionp_near_cutoff,
			importance_sampling.extracellular_ionp_far_cutoff));
	return true;
}

// Deterministic pilot integration estimates accessible joint-region volumes.
// Shell proposals help resolve small near-source regions; coverage correction
// prevents overlapping shells from counting their intersection multiple times.
std::array<double, kFreeSampleRegionCount>
estimate_active_joint_region_volume_weights(
	const ImportanceSamplingConfig& config,
	const std::vector<Point3D>& cell_origins,
	const Point3D& cell_origin,
	double cell_r,
	const std::vector<Ionp>& ionps,
	const IonpPopulationMap& populations,
	double sim_x_length,
	double sim_y_length,
	double sim_z_length,
	std::size_t sample_count) {
	if (sample_count == 0) {
		throw std::runtime_error(
			"Importance-sampling volume-weight estimation requires a positive sample count");
	}

	std::vector<ShellFamily> families;
	auto add_family = [&](IonpPopulation population,
						  IonpDistanceRegion region,
						  const std::vector<std::size_t>& indices,
						  double near_cutoff,
						  double far_cutoff) {
		if (indices.empty() ||
			!any_active_region_for_axis(config, population, region)) {
			return;
		}
		ShellFamily family{
			population,
			region,
			&indices,
			near_cutoff,
			far_cutoff,
			0.,
			0.,
			{}
		};
		family.cumulative_shell_volumes.reserve(indices.size());
		for (std::size_t index : indices) {
			family.total_shell_volume += sphere_shell_volume(
				shell_inner_radius(ionps[index], region, near_cutoff),
				shell_outer_radius(
					ionps[index], region, near_cutoff, far_cutoff));
			family.cumulative_shell_volumes.push_back(
				family.total_shell_volume);
		}
		if (family.total_shell_volume > 0.) {
			families.push_back(std::move(family));
		}
	};

	add_family(
		IonpPopulation::Intracellular,
		IonpDistanceRegion::Near,
		populations.intracellular_indices,
		config.intracellular_ionp_near_cutoff,
		config.intracellular_ionp_far_cutoff);
	add_family(
		IonpPopulation::Intracellular,
		IonpDistanceRegion::Intermediate,
		populations.intracellular_indices,
		config.intracellular_ionp_near_cutoff,
		config.intracellular_ionp_far_cutoff);
	add_family(
		IonpPopulation::Extracellular,
		IonpDistanceRegion::Near,
		populations.extracellular_indices,
		config.extracellular_ionp_near_cutoff,
		config.extracellular_ionp_far_cutoff);
	add_family(
		IonpPopulation::Extracellular,
		IonpDistanceRegion::Intermediate,
		populations.extracellular_indices,
		config.extracellular_ionp_near_cutoff,
		config.extracellular_ionp_far_cutoff);

	const double uniform_probability = families.empty() ? 1.0 : 0.5;
	const double family_probability =
		families.empty() ? 0.0 : 0.5 / static_cast<double>(families.size());
	for (ShellFamily& family : families) {
		family.proposal_probability = family_probability;
	}

	const double half_x = sim_x_length / 2.;
	const double half_y = sim_y_length / 2.;
	const double half_z = sim_z_length / 2.;
	const double box_volume =
		sim_x_length * sim_y_length * sim_z_length;
	if (box_volume <= 0.) {
		throw std::runtime_error(
			"Importance-sampling volume-weight estimation requires a positive simulation volume");
	}

	double maximum_outer_radius = 0.;
	for (const Ionp& ionp : ionps) {
		maximum_outer_radius =
			std::max(maximum_outer_radius, ionp.outer_radius());
	}
	double maximum_distance_cutoff = 0.;
	if (config.intracellular_ionp_sampling_enabled()) {
		maximum_distance_cutoff =
			std::max(maximum_distance_cutoff, config.intracellular_ionp_far_cutoff);
	}
	if (config.extracellular_ionp_sampling_enabled()) {
		maximum_distance_cutoff =
			std::max(maximum_distance_cutoff, config.extracellular_ionp_far_cutoff);
	}
	const double candidate_search_radius =
		maximum_outer_radius + maximum_distance_cutoff;
	const ImportanceParticleGrid particle_grid(
		ionps, candidate_search_radius);
	const ImportanceCellGrid cell_grid(
		cell_origins, cell_origin, cell_r, config.cell_near_cutoff);
	std::vector<std::size_t> particle_candidates;

	// Weight estimation is deterministic and does not consume the trajectory
	// RNG sequence, so adding a geometry-support check cannot change a run.
	RandNumberBetween x_random(-half_x, half_x, 0x49a6f3d1u);
	RandNumberBetween y_random(-half_y, half_y, 0x8b17c2e5u);
	RandNumberBetween z_random(-half_z, half_z, 0xd4e9037bu);
	RandNumberBetween random_unit(0., 1., 0x2f6c81a9u);

	std::array<double, kFreeSampleRegionCount> accumulated_volume{};
	for (std::size_t sample_index = 0; sample_index < sample_count; ++sample_index) {
		Point3D point;
		const double component_draw = random_unit();
		if (component_draw < uniform_probability) {
			point.set(x_random(), y_random(), z_random());
		}
		else {
			const std::size_t family_index = std::min(
				static_cast<std::size_t>(
					(component_draw - uniform_probability) /
					family_probability),
				families.size() - 1);
			const ShellFamily& family = families[family_index];
			const double particle_draw =
				random_unit() * family.total_shell_volume;
			const auto chosen_position = std::lower_bound(
				family.cumulative_shell_volumes.begin(),
				family.cumulative_shell_volumes.end(),
				particle_draw);
			const std::size_t chosen_offset =
				chosen_position == family.cumulative_shell_volumes.end()
					? family.cumulative_shell_volumes.size() - 1
					: static_cast<std::size_t>(
						chosen_position -
						family.cumulative_shell_volumes.begin());
			const std::size_t chosen =
				(*family.particle_indices)[chosen_offset];
			point = random_point_in_shell(
				ionps[chosen],
				family.region,
				family.near_cutoff,
				family.far_cutoff,
				random_unit);
		}

		const bool in_box =
			point.x() >= -half_x && point.x() <= half_x &&
			point.y() >= -half_y && point.y() <= half_y &&
			point.z() >= -half_z && point.z() <= half_z;
		if (!in_box) {
			continue;
		}

		particle_grid.candidates(
			point, candidate_search_radius, particle_candidates);
		double proposal_density = uniform_probability / box_volume;
		for (const ShellFamily& family : families) {
			proposal_density += family.proposal_probability *
				static_cast<double>(
					shell_coverage(
						point,
						ionps,
						populations,
						particle_candidates,
						family)) /
				family.total_shell_volume;
		}
		if (proposal_density <= 0.) {
			continue;
		}

		ProtonSampleRegion region;
		if (!classify_joint_sample_region_candidates(
				point,
				cell_grid.classify(point),
				ionps,
				populations,
				particle_candidates,
				config,
				region)) {
			continue;
		}
		const std::size_t region_index = sample_region_index(region);
		if (config.joint_protons[region_index] > 0) {
			accumulated_volume[region_index] += 1.0 / proposal_density;
		}
	}

	std::array<double, kFreeSampleRegionCount> weights{};
	for (std::size_t index = 0; index < kFreeSampleRegionCount; ++index) {
		if (config.joint_protons[index] <= 0) {
			continue;
		}
		weights[index] =
			accumulated_volume[index] /
			static_cast<double>(sample_count) /
			box_volume;
		if (weights[index] <= 0.) {
			throw std::runtime_error(
				"Automatic importance-sampling weight is zero for requested region " +
				sample_region_name(static_cast<ProtonSampleRegion>(index)));
		}
	}
	return weights;
}

std::array<double, 3> estimate_active_free_region_volume_weights(
	const ImportanceSamplingConfig& config,
	const std::vector<Point3D>& cell_origins,
	const Point3D& cell_origin,
	double cell_r,
	double sim_x_length,
	double sim_y_length,
	double sim_z_length,
	std::size_t sample_count) {
	ImportanceSamplingConfig materialized = config;
	materialized.materialize_legacy_regions();
	const IonpPopulationMap populations;
	const auto joint = estimate_active_joint_region_volume_weights(
		materialized,
		cell_origins,
		cell_origin,
		cell_r,
		{},
		populations,
		sim_x_length,
		sim_y_length,
		sim_z_length,
		sample_count);
	return {
		joint[joint_sample_region_index(
			CellSampleRegion::Inside,
			IonpDistanceRegion::Far,
			IonpDistanceRegion::Far)],
		joint[joint_sample_region_index(
			CellSampleRegion::Near,
			IonpDistanceRegion::Far,
			IonpDistanceRegion::Far)],
		joint[joint_sample_region_index(
			CellSampleRegion::Far,
			IonpDistanceRegion::Far,
			IonpDistanceRegion::Far)]
	};
}

std::array<double, kSampleRegionCount> normalized_sample_weights(
	const ImportanceSamplingConfig& config,
	std::size_t bound_intracellular_sample_count,
	std::size_t bound_extracellular_sample_count) {
	std::array<double, kSampleRegionCount> weights{};
	for (std::size_t index = 0; index < kFreeSampleRegionCount; ++index) {
		weights[index] = config.joint_weights[index];
	}

	if (config.split_bound_configured) {
		weights[kBoundIntracellularSampleRegionIndex] =
			bound_intracellular_sample_count > 0
				? config.weight_bound_intracellular
				: 0.;
		weights[kBoundExtracellularSampleRegionIndex] =
			bound_extracellular_sample_count > 0
				? config.weight_bound_extracellular
				: 0.;
	}
	else {
		const std::size_t total_bound =
				bound_intracellular_sample_count +
				bound_extracellular_sample_count;
		if (total_bound > 0) {
			weights[kBoundIntracellularSampleRegionIndex] =
				config.weight_bound *
				static_cast<double>(bound_intracellular_sample_count) /
				static_cast<double>(total_bound);
			weights[kBoundExtracellularSampleRegionIndex] =
				config.weight_bound *
				static_cast<double>(bound_extracellular_sample_count) /
				static_cast<double>(total_bound);
		}
	}

	double total_weight = 0.;
	for (std::size_t index = 0; index < kSampleRegionCount; ++index) {
		if (weights[index] < 0.) {
			throw std::runtime_error("Importance-sampling weights must be non-negative");
		}
		const std::size_t sample_count =
			index == kBoundIntracellularSampleRegionIndex
				? bound_intracellular_sample_count
				: index == kBoundExtracellularSampleRegionIndex
					? bound_extracellular_sample_count
					: static_cast<std::size_t>(config.joint_protons[index]);
		if (weights[index] > 0. && sample_count == 0) {
			throw std::runtime_error(
				"Positive importance-sampling weight has no samples for region " +
				sample_region_name(static_cast<ProtonSampleRegion>(index)));
		}
		total_weight += weights[index];
	}
	if (total_weight <= 0.) {
		throw std::runtime_error(
			"Importance-sampling weights must sum to a positive value");
	}
	for (double& weight : weights) {
		weight /= total_weight;
	}
	return weights;
}

std::array<double, kSampleRegionCount> normalized_sample_weights(
	const ImportanceSamplingConfig& config,
	bool include_bound) {
	ImportanceSamplingConfig materialized = config;
	materialized.materialize_legacy_regions();
	return normalized_sample_weights(
		materialized,
		0,
		include_bound ? 1u : 0u);
}

// Legacy estimator: first average exp(i*phase) within each exclusive region,
// then combine those regional means with normalized region weights.
WeightedSignalAggregation aggregate_weighted_signals(
	const std::vector<std::vector<double>>& phase_vector,
	const std::vector<ProtonSampleRegion>& sample_regions,
	const std::array<std::size_t, kSampleRegionCount>& sample_counts,
	const std::array<double, kSampleRegionCount>& normalized_weights) {
	if (phase_vector.size() != sample_regions.size()) {
		throw std::runtime_error(
			"Phase histories and sample-region labels must have matching sizes");
	}

	WeightedSignalAggregation aggregation;
	if (phase_vector.empty()) {
		return aggregation;
	}

	const std::size_t time_points = phase_vector.front().size();
	for (const auto& phases : phase_vector) {
		if (phases.size() != time_points) {
			throw std::runtime_error(
				"All proton phase histories must have the same length");
		}
	}

	for (auto& signal : aggregation.region_signals) {
		signal.assign(time_points, std::complex<double>(0., 0.));
	}
	aggregation.combined_signal.assign(
		time_points, std::complex<double>(0., 0.));

	const std::complex<double> i_compl(0., 1.);
	std::array<std::size_t, kSampleRegionCount> seen_counts{};
	for (std::size_t proton_index = 0;
		 proton_index < phase_vector.size();
		 ++proton_index) {
		const std::size_t region_index =
			sample_region_index(sample_regions[proton_index]);
		if (region_index >= kSampleRegionCount) {
			throw std::runtime_error("Invalid importance-sampling region label");
		}
		++seen_counts[region_index];
		if (sample_counts[region_index] == 0) {
			continue;
		}
		for (std::size_t time_index = 0;
			 time_index < time_points;
			 ++time_index) {
			aggregation.region_signals[region_index][time_index] +=
				std::exp(i_compl * phase_vector[proton_index][time_index]) /
				static_cast<double>(sample_counts[region_index]);
		}
	}

	for (std::size_t region_index = 0;
		 region_index < kSampleRegionCount;
		 ++region_index) {
		if (seen_counts[region_index] != sample_counts[region_index]) {
			throw std::runtime_error(
				"Sample-region labels do not match the expected importance-sampling counts for " +
				sample_region_name(
					static_cast<ProtonSampleRegion>(region_index)));
		}
		for (std::size_t time_index = 0;
			 time_index < time_points;
			 ++time_index) {
			aggregation.combined_signal[time_index] +=
				normalized_weights[region_index] *
				aggregation.region_signals[region_index][time_index];
		}
	}
	return aggregation;
}

// Estimate free-water support and marginal proposal fractions with a
// separately seeded pilot, not the proton propagation streams. Fractions for
// overlapping proposals need not sum to one.
OverlappingProposalVolumeEstimate estimate_overlapping_proposal_volumes(
	const ImportanceSamplingConfig& config,
	const std::vector<Point3D>& cell_origins,
	const Point3D& cell_origin,
	double cell_r,
	const std::vector<Ionp>& ionps,
	const IonpPopulationMap& populations,
	double sim_x_length,
	double sim_y_length,
	double sim_z_length,
	std::size_t sample_count) {
	config.validate_mode_exclusivity();
	if (!config.overlapping_proposals_enabled()) {
		throw std::runtime_error(
			"Overlapping proposal volume estimation requires OverlappingProposal mode");
	}
	if (sample_count == 0) {
		throw std::runtime_error(
			"Overlapping importance-sampling volume estimation requires a positive sample count");
	}

	const double box_volume =
		sim_x_length * sim_y_length * sim_z_length;
	if (!(box_volume > 0.) || !std::isfinite(box_volume)) {
		throw std::runtime_error(
			"Overlapping importance sampling requires a finite positive simulation volume");
	}

	std::vector<ShellFamily> shell_families;
	auto add_shell_family = [&](IonpPopulation population,
								IonpDistanceRegion region,
								FreeSamplingProposal proposal,
								const std::vector<std::size_t>& indices,
								double near_cutoff,
								double far_cutoff) {
		if (config.proposal_protons[free_sampling_proposal_index(proposal)] <= 0 ||
			indices.empty() || region == IonpDistanceRegion::Far) {
			return;
		}
		ShellFamily family{
			population,
			region,
			&indices,
			near_cutoff,
			far_cutoff,
			0.,
			0.,
			{}
		};
		family.cumulative_shell_volumes.reserve(indices.size());
		for (std::size_t index : indices) {
			const double volume = sphere_shell_volume(
				shell_inner_radius(ionps[index], region, near_cutoff),
				shell_outer_radius(
					ionps[index], region, near_cutoff, far_cutoff));
			if (volume > 0. && std::isfinite(volume)) {
				family.total_shell_volume += volume;
			}
			family.cumulative_shell_volumes.push_back(
				family.total_shell_volume);
		}
		if (family.total_shell_volume > 0.) {
			shell_families.push_back(std::move(family));
		}
	};

	add_shell_family(
		IonpPopulation::Intracellular,
		IonpDistanceRegion::Near,
		FreeSamplingProposal::IntracellularIonpNear,
		populations.intracellular_indices,
		config.intracellular_ionp_near_cutoff,
		config.intracellular_ionp_far_cutoff);
	add_shell_family(
		IonpPopulation::Intracellular,
		IonpDistanceRegion::Intermediate,
		FreeSamplingProposal::IntracellularIonpIntermediate,
		populations.intracellular_indices,
		config.intracellular_ionp_near_cutoff,
		config.intracellular_ionp_far_cutoff);
	add_shell_family(
		IonpPopulation::Extracellular,
		IonpDistanceRegion::Near,
		FreeSamplingProposal::ExtracellularIonpNear,
		populations.extracellular_indices,
		config.extracellular_ionp_near_cutoff,
		config.extracellular_ionp_far_cutoff);
	add_shell_family(
		IonpPopulation::Extracellular,
		IonpDistanceRegion::Intermediate,
		FreeSamplingProposal::ExtracellularIonpIntermediate,
		populations.extracellular_indices,
		config.extracellular_ionp_near_cutoff,
		config.extracellular_ionp_far_cutoff);

	std::vector<CellShellFamily> cell_families;
	const std::size_t cell_count =
		available_cell_count(cell_origins, cell_r);
	auto add_cell_family = [&](CellSampleRegion region,
							   FreeSamplingProposal proposal,
							   double inner_radius,
							   double outer_radius) {
		if (config.proposal_protons[free_sampling_proposal_index(proposal)] <= 0 ||
			cell_count == 0 || !(outer_radius > inner_radius)) {
			return;
		}
		const double primitive_volume =
			sphere_shell_volume(inner_radius, outer_radius);
		const double total_volume =
			static_cast<double>(cell_count) * primitive_volume;
		if (total_volume > 0. && std::isfinite(total_volume)) {
			cell_families.push_back(
				{region, inner_radius, outer_radius, total_volume, 0.});
		}
	};
	add_cell_family(
		CellSampleRegion::Inside,
		FreeSamplingProposal::CellInside,
		0.,
		cell_r);
	add_cell_family(
		CellSampleRegion::Near,
		FreeSamplingProposal::CellNear,
		cell_r,
		config.cell_near_cutoff);

	const std::size_t rare_family_count =
		shell_families.size() + cell_families.size();
	const double uniform_probability =
		rare_family_count == 0 ? 1.0 : 0.5;
	const double rare_family_probability =
		rare_family_count == 0
			? 0.0
			: 0.5 / static_cast<double>(rare_family_count);
	for (ShellFamily& family : shell_families) {
		family.proposal_probability = rare_family_probability;
	}
	for (CellShellFamily& family : cell_families) {
		family.proposal_probability = rare_family_probability;
	}

	double maximum_outer_radius = 0.;
	for (const Ionp& ionp : ionps) {
		maximum_outer_radius =
			std::max(maximum_outer_radius, ionp.outer_radius());
	}
	double maximum_distance_cutoff = 0.;
	if (config.intracellular_ionp_sampling_enabled()) {
		maximum_distance_cutoff = std::max(
			maximum_distance_cutoff,
			config.intracellular_ionp_far_cutoff);
	}
	if (config.extracellular_ionp_sampling_enabled()) {
		maximum_distance_cutoff = std::max(
			maximum_distance_cutoff,
			config.extracellular_ionp_far_cutoff);
	}
	const double candidate_search_radius =
		maximum_outer_radius + maximum_distance_cutoff;
	const ImportanceParticleGrid particle_grid(
		ionps, candidate_search_radius);
	const ImportanceCellGrid cell_grid(
		cell_origins, cell_origin, cell_r, config.cell_near_cutoff);
	std::vector<std::size_t> particle_candidates;

	const double half_x = sim_x_length / 2.;
	const double half_y = sim_y_length / 2.;
	const double half_z = sim_z_length / 2.;
	RandNumberBetween x_random(-half_x, half_x, 0x7c18b42du);
	RandNumberBetween y_random(-half_y, half_y, 0x35e91af7u);
	RandNumberBetween z_random(-half_z, half_z, 0xb84276c3u);
	RandNumberBetween random_unit(0., 1., 0x61d94e0bu);

	std::array<double, kFreeSamplingProposalCount> accumulated_volumes{};
	for (std::size_t sample_index = 0;
		 sample_index < sample_count;
		 ++sample_index) {
		Point3D point;
		const double component_draw = random_unit();
		if (component_draw < uniform_probability) {
			point.set(x_random(), y_random(), z_random());
		}
		else {
			const std::size_t family_index = std::min(
				static_cast<std::size_t>(
					(component_draw - uniform_probability) /
					rare_family_probability),
				rare_family_count - 1);
			if (family_index < shell_families.size()) {
				const ShellFamily& family = shell_families[family_index];
				const double particle_draw =
					random_unit() * family.total_shell_volume;
				const auto chosen_position = std::lower_bound(
					family.cumulative_shell_volumes.begin(),
					family.cumulative_shell_volumes.end(),
					particle_draw);
				const std::size_t chosen_offset =
					chosen_position == family.cumulative_shell_volumes.end()
						? family.cumulative_shell_volumes.size() - 1
						: static_cast<std::size_t>(
							chosen_position -
							family.cumulative_shell_volumes.begin());
				const std::size_t chosen =
					(*family.particle_indices)[chosen_offset];
				point = random_point_in_shell(
					ionps[chosen],
					family.region,
					family.near_cutoff,
					family.far_cutoff,
					random_unit);
			}
			else {
				const CellShellFamily& family =
					cell_families[family_index - shell_families.size()];
				const std::size_t selected = std::min(
					static_cast<std::size_t>(
						random_unit() * static_cast<double>(cell_count)),
					cell_count - 1);
				point = random_point_in_spherical_shell(
					cell_center_at(
						cell_origins, cell_origin, selected),
					family.inner_radius,
					family.outer_radius,
					random_unit);
			}
		}

		const bool in_box =
			point.x() >= -half_x && point.x() <= half_x &&
			point.y() >= -half_y && point.y() <= half_y &&
			point.z() >= -half_z && point.z() <= half_z;
		if (!in_box) {
			continue;
		}

		particle_grid.candidates(
			point, candidate_search_radius, particle_candidates);
		double pilot_density = uniform_probability / box_volume;
		for (const ShellFamily& family : shell_families) {
			pilot_density += family.proposal_probability *
				static_cast<double>(shell_coverage(
					point,
					ionps,
					populations,
					particle_candidates,
					family)) /
				family.total_shell_volume;
		}
		for (const CellShellFamily& family : cell_families) {
			pilot_density += family.proposal_probability *
				static_cast<double>(cell_grid.coverage(
					point,
					family.inner_radius,
					family.outer_radius)) /
				family.total_shell_volume;
		}
		if (!(pilot_density > 0.) || !std::isfinite(pilot_density)) {
			continue;
		}

		ProtonSampleRegion joint_region;
		if (!classify_joint_sample_region_candidates(
				point,
				cell_grid.classify(point),
				ionps,
				populations,
				particle_candidates,
				config,
				joint_region)) {
			continue;
		}

		const double volume_contribution = 1.0 / pilot_density;
		accumulated_volumes[free_sampling_proposal_index(
			FreeSamplingProposal::Uniform)] += volume_contribution;
		const JointSampleRegion membership =
			decode_joint_sample_region(joint_region);
		for (FreeSamplingProposal proposal :
			 marginal_sampling_proposals(membership)) {
			accumulated_volumes[free_sampling_proposal_index(proposal)] +=
				volume_contribution;
		}
	}

	OverlappingProposalVolumeEstimate result;
	result.box_volume = box_volume;
	result.pilot_samples = sample_count;
	for (std::size_t index = 0;
		 index < kFreeSamplingProposalCount;
		 ++index) {
		result.region_volumes[index] =
			accumulated_volumes[index] /
			static_cast<double>(sample_count);
	}
	result.accessible_volume =
		result.region_volumes[free_sampling_proposal_index(
			FreeSamplingProposal::Uniform)];
	if (!(result.accessible_volume > 0.) ||
		!std::isfinite(result.accessible_volume)) {
		throw std::runtime_error(
			"Overlapping importance-sampling pilot found no accessible free-proton volume");
	}
	for (std::size_t index = 0;
		 index < kFreeSamplingProposalCount;
		 ++index) {
		result.region_fractions[index] =
			result.region_volumes[index] / result.accessible_volume;
	}
	result.region_fractions[free_sampling_proposal_index(
		FreeSamplingProposal::Uniform)] = 1.0;
	return result;
}

// Move allocations with zero estimated support into the uniform proposal,
// preserving the total requested number of free samples.
EffectiveProposalCounts resolve_effective_proposal_counts(
	const std::array<int, kFreeSamplingProposalCount>& requested_counts,
	const std::array<double, kFreeSamplingProposalCount>& region_fractions) {
	validate_overlapping_fraction_bounds(region_fractions);
	constexpr double partition_tolerance = 1e-8;
	for (std::size_t family_start : {std::size_t{1}, std::size_t{4}, std::size_t{7}}) {
		const double family_total =
			region_fractions[family_start] +
			region_fractions[family_start + 1] +
			region_fractions[family_start + 2];
		if (std::abs(family_total - 1.0) > partition_tolerance) {
			throw std::runtime_error(
				"Each overlapping proposal family must partition the accessible free-water volume");
		}
	}
	EffectiveProposalCounts result;
	const std::size_t uniform_index = free_sampling_proposal_index(
		FreeSamplingProposal::Uniform);
	for (std::size_t index = 0;
		 index < kFreeSamplingProposalCount;
		 ++index) {
		if (requested_counts[index] < 0) {
			throw std::runtime_error(
				"Overlapping importance-sampling proposal counts must be non-negative");
		}
		result.counts[index] =
			static_cast<std::size_t>(requested_counts[index]);
	}
	if (requested_counts[uniform_index] <= 0) {
		throw std::runtime_error(
			"Overlapping importance sampling requires a positive uniform proposal count");
	}
	if (!(region_fractions[uniform_index] > 0.) ||
		!std::isfinite(region_fractions[uniform_index])) {
		throw std::runtime_error(
			"The uniform importance-sampling proposal has no accessible volume");
	}

	for (std::size_t index = 1;
		 index < kFreeSamplingProposalCount;
		 ++index) {
		if (result.counts[index] == 0) {
			continue;
		}
		if (region_fractions[index] > 0. &&
			std::isfinite(region_fractions[index])) {
			continue;
		}
		result.redistributed[index] = true;
		result.redistributed_to_uniform += result.counts[index];
		result.counts[index] = 0;
	}
	result.counts[uniform_index] += result.redistributed_to_uniform;
	return result;
}

// For a start x, each covering proposal k contributes n_k/f_k to the
// deterministic-mixture density (f_k is its accessible volume fraction).
// The reciprocal is normalized across sampled starts before aggregation.
double overlapping_mis_raw_weight(
	const JointSampleRegion& membership,
	const std::array<std::size_t, kFreeSamplingProposalCount>& effective_counts,
	const std::array<double, kFreeSamplingProposalCount>& region_fractions) {
	validate_overlapping_fraction_bounds(region_fractions);
	const std::size_t uniform_index = free_sampling_proposal_index(
		FreeSamplingProposal::Uniform);
	if (effective_counts[uniform_index] == 0) {
		throw std::runtime_error(
			"Overlapping MIS weights require a positive uniform proposal count");
	}
	if (!(region_fractions[uniform_index] > 0.) ||
		!std::isfinite(region_fractions[uniform_index])) {
		throw std::runtime_error(
			"Overlapping MIS weights require positive uniform support");
	}

	double denominator =
		static_cast<double>(effective_counts[uniform_index]) /
		region_fractions[uniform_index];
	for (FreeSamplingProposal proposal :
		 marginal_sampling_proposals(membership)) {
		const std::size_t index = free_sampling_proposal_index(proposal);
		if (effective_counts[index] == 0) {
			continue;
		}
		const double fraction = region_fractions[index];
		if (!(fraction > 0.) || !std::isfinite(fraction)) {
			throw std::runtime_error(
				"Active importance-sampling proposal " +
				free_sampling_proposal_name(proposal) +
				" has no positive finite volume fraction");
		}
		denominator +=
			static_cast<double>(effective_counts[index]) / fraction;
	}
	if (!(denominator > 0.) || !std::isfinite(denominator)) {
		throw std::runtime_error(
			"Overlapping MIS proposal density is not positive and finite");
	}
	return 1.0 / denominator;
}

std::vector<double> overlapping_mis_raw_weights(
	const std::vector<JointSampleRegion>& memberships,
	const std::array<std::size_t, kFreeSamplingProposalCount>& effective_counts,
	const std::array<double, kFreeSamplingProposalCount>& region_fractions) {
	std::vector<double> result;
	result.reserve(memberships.size());
	for (const JointSampleRegion& membership : memberships) {
		result.push_back(overlapping_mis_raw_weight(
			membership, effective_counts, region_fractions));
	}
	return result;
}

// Divide each positive raw weight by the total so the free-sample weights
// sum to one. Bound-water mixture weights are handled separately by MRsequence.
std::vector<double> self_normalize_importance_weights(
	const std::vector<double>& raw_weights) {
	if (raw_weights.empty()) {
		return {};
	}
	double total = 0.;
	for (double weight : raw_weights) {
		if (!(weight > 0.) || !std::isfinite(weight)) {
			throw std::runtime_error(
				"Raw importance-sampling weights must be positive and finite");
		}
		total += weight;
	}
	if (!(total > 0.) || !std::isfinite(total)) {
		throw std::runtime_error(
			"Raw importance-sampling weights cannot be normalized");
	}

	std::vector<double> result = raw_weights;
	for (double& weight : result) {
		weight /= total;
	}
	return result;
}

std::vector<double> overlapping_mis_self_normalized_weights(
	const std::vector<JointSampleRegion>& memberships,
	const std::array<std::size_t, kFreeSamplingProposalCount>& effective_counts,
	const std::array<double, kFreeSamplingProposalCount>& region_fractions) {
	return self_normalize_importance_weights(
		overlapping_mis_raw_weights(
			memberships, effective_counts, region_fractions));
}

// Histories are indexed [proton][recorded time]. Bin ordinary means and
// weighted contributions by the proposal that GENERATED each start. Only the
// complex weighted contributions add up to the combined free-water signal.
OverlappingWeightedSignalAggregation aggregate_overlapping_weighted_signals(
	const std::vector<std::vector<double>>& phase_vector,
	const std::vector<FreeSamplingProposal>& sample_proposals,
	const std::vector<double>& normalized_sample_weights) {
	if (sample_proposals.size() != normalized_sample_weights.size() ||
		phase_vector.size() < sample_proposals.size()) {
		throw std::runtime_error(
			"Overlapping importance-sampling proposal and weight counts must match and cannot exceed phase histories");
	}

	OverlappingWeightedSignalAggregation result;
	if (sample_proposals.empty()) {
		return result;
	}
	const std::size_t time_points = phase_vector.front().size();
	for (std::size_t sample_index = 0;
		 sample_index < sample_proposals.size();
		 ++sample_index) {
		if (phase_vector[sample_index].size() != time_points) {
			throw std::runtime_error(
				"All free overlapping importance-sampling phase histories must have the same length");
		}
	}

	std::array<std::size_t, kFreeSamplingProposalCount> proposal_counts{};
	double total_weight = 0.;
	for (std::size_t sample_index = 0;
		 sample_index < sample_proposals.size();
		 ++sample_index) {
		const std::size_t proposal_index =
			free_sampling_proposal_index(sample_proposals[sample_index]);
		if (proposal_index >= kFreeSamplingProposalCount) {
			throw std::runtime_error(
				"Invalid overlapping importance-sampling proposal label");
		}
		const double weight = normalized_sample_weights[sample_index];
		if (weight < 0. || !std::isfinite(weight)) {
			throw std::runtime_error(
				"Normalized importance-sampling weights must be non-negative and finite");
		}
		++proposal_counts[proposal_index];
		total_weight += weight;
	}
	if (!(total_weight > 0.) || !std::isfinite(total_weight) ||
		std::abs(total_weight - 1.0) > 1e-10) {
		throw std::runtime_error(
			"Normalized importance-sampling weights must sum to one");
	}

	for (std::size_t proposal_index = 0;
		 proposal_index < kFreeSamplingProposalCount;
		 ++proposal_index) {
		if (proposal_counts[proposal_index] == 0) {
			continue;
		}
		result.proposal_means[proposal_index].assign(
			time_points, std::complex<double>(0., 0.));
		result.proposal_contributions[proposal_index].assign(
			time_points, std::complex<double>(0., 0.));
	}
	result.combined_signal.assign(
		time_points, std::complex<double>(0., 0.));

	const std::complex<double> i_compl(0., 1.);
	for (std::size_t sample_index = 0;
		 sample_index < sample_proposals.size();
		 ++sample_index) {
		const std::size_t proposal_index =
			free_sampling_proposal_index(sample_proposals[sample_index]);
		for (std::size_t time_index = 0;
			 time_index < time_points;
			 ++time_index) {
			const std::complex<double> sample_signal =
				std::exp(i_compl * phase_vector[sample_index][time_index]);
			result.proposal_means[proposal_index][time_index] +=
				sample_signal /
				static_cast<double>(proposal_counts[proposal_index]);
			const std::complex<double> contribution =
				normalized_sample_weights[sample_index] * sample_signal;
			result.proposal_contributions[proposal_index][time_index] +=
				contribution;
			result.combined_signal[time_index] += contribution;
		}
	}
	return result;
}
