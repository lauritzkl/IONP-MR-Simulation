/**
 * @file SimSpace.cpp
 *
 * Geometry construction: optional input selection, rejection placement,
 * cellular/extracellular populations and surrounding magnetic-source copies.
 * The overlap grid accelerates placement checks; it does not alter the exact
 * sphere-overlap predicate. Random draw and insertion order affect seeded runs.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <complex>
#include <omp.h>
#include <filesystem>
#include <limits>
#include <optional>
#include <stdexcept>
#include <unordered_map>

#include "mptmacros.h"

#include "Point3D.hpp"
#include "Vector3D.hpp"
#include "Particle.hpp"
#include "Proton.hpp"
#include "Ionp.hpp"
#include "Aggregate.hpp"
#include "SimMatrix.hpp"
#include "Data.hpp"

#include "SimSpace.hpp"

namespace {

// Lightweight spatial hash grid to avoid O(N^2) overlap checks when placing IONPs.
struct CellKey {
	int x;
	int y;
	int z;

	bool operator==(const CellKey& other) const {
		return x == other.x && y == other.y && z == other.z;
	}
};

struct CellKeyHash {
	std::size_t operator()(const CellKey& key) const {
		// Simple mix that stays deterministic and cheap.
		std::size_t h = static_cast<std::size_t>(key.x);
		// Cast operands to size_t before multiplication to avoid signed overflow UB.
		h = h * 73856093u ^ (static_cast<std::size_t>(key.y) * 19349663u);
		h ^= (static_cast<std::size_t>(key.z) * 83492791u);
		return h;
	}
};

// Placement-only broad phase: bucket accepted spheres and test nearby candidates
// with the same exact overlap predicate used by the legacy exhaustive path.
class OverlapGrid {
public:
	explicit OverlapGrid(double cell_size) : m_cell_size(cell_size) {}

	void insert(const Point3D& p, unsigned int index, double outer_radius) {
		cell(p).push_back(index);
		m_max_outer_radius = std::max(m_max_outer_radius, outer_radius);
	}

	bool overlaps(const Point3D& p, double radius, const std::vector<Ionp>& ionps) const {
		const CellKey c = key(p);
		// Include both the candidate radius and the largest previously inserted
		// radius. Cell-load-derived equivalent particles can be much larger than
		// the configured radii used to choose m_cell_size.
		const int reach = std::max(1, static_cast<int>(std::ceil(
			(radius + m_max_outer_radius) / m_cell_size)));
		for (int dz = -reach; dz <= reach; ++dz) {
			for (int dy = -reach; dy <= reach; ++dy) {
				for (int dx = -reach; dx <= reach; ++dx) {
					const CellKey nkey{c.x + dx, c.y + dy, c.z + dz};
					const auto it = m_cells.find(nkey);
					if (it == m_cells.end()) {
						continue;
					}
					for (unsigned int idx : it->second) {
						const double allowed = radius + ionps[idx].outer_radius();
						if (distance(p, ionps[idx].xt()) < allowed) {
							return true;
						}
					}
				}
			}
		}
		return false;
	}

private:
	CellKey key(const Point3D& p) const {
		return {
			static_cast<int>(std::floor(p.x() / m_cell_size)),
			static_cast<int>(std::floor(p.y() / m_cell_size)),
			static_cast<int>(std::floor(p.z() / m_cell_size))
		};
	}

	std::vector<unsigned int>& cell(const Point3D& p) {
		return m_cells[key(p)];
	}

	double m_cell_size;
	double m_max_outer_radius = 0.;
	std::unordered_map<CellKey, std::vector<unsigned int>, CellKeyHash> m_cells;
};

// Tries to place a single IONP using the provided placement functor until it no longer overlaps neighbors.
template <typename PlaceFunc>
// Retry the supplied placement callback until its sphere is accepted; each
// retry consumes random draws. Exhausting the existing limit reports failure.
void place_ionp_no_overlap(
    Ionp& ionp,
    std::vector<Ionp>& ionpInstances,
    OverlapGrid& grid,
    PlaceFunc&& placer,
    unsigned int current_index,
    unsigned int max_attempts = 5000) {
	const double radius = ionp.outer_radius();
	for (unsigned int attempt = 0; attempt < max_attempts; ++attempt) {
		if (!placer()) {
			continue; // placement failed validation (e.g. inside a cell), try again
		}
		if (!grid.overlaps(ionp.xt(), radius, ionpInstances)) {
			grid.insert(ionp.xt(), current_index, radius);
			return;
		}
	}
	throw std::runtime_error("Unable to place IONP without overlap after many attempts");
}

std::filesystem::path require_ini_subdir(const std::filesystem::path& ini_root_dir, const char* subdir) {
	const std::filesystem::path path = ini_root_dir / subdir;
	if (!std::filesystem::is_directory(path)) {
		throw std::runtime_error("Required input directory does not exist: " + path.string());
	}
	return path;
}

bool is_supported_geometry_file(const std::filesystem::path& path) {
	return path.extension() == ".txt" || path.extension() == ".dat";
}

// An explicit file takes precedence; otherwise collect supported files from
// the legacy directory and sort them before appending their records.
std::vector<std::filesystem::path> geometry_input_files(
		const std::filesystem::path& input_dir,
		const std::optional<std::filesystem::path>& selected_path,
		const char* geometry_name) {
	if (selected_path) {
		const std::filesystem::path absolute_path =
			std::filesystem::absolute(*selected_path).lexically_normal();
		if (!is_supported_geometry_file(absolute_path)) {
			throw std::runtime_error(std::string("Selected ") + geometry_name +
				" file must use the .txt or .dat extension: " + absolute_path.string());
		}

		std::error_code ec;
		if (!std::filesystem::is_regular_file(absolute_path, ec) || ec) {
			throw std::runtime_error(std::string("Selected ") + geometry_name +
				" path is not a regular file: " + absolute_path.string());
		}
		const std::filesystem::path canonical_dir =
			std::filesystem::canonical(std::filesystem::absolute(input_dir), ec);
		if (ec) {
			throw std::runtime_error("Unable to resolve geometry directory " +
				input_dir.string() + ": " + ec.message());
		}
		const std::filesystem::path canonical_path =
			std::filesystem::canonical(absolute_path, ec);
		if (ec || canonical_path.parent_path() != canonical_dir) {
			throw std::runtime_error(std::string("Selected ") + geometry_name +
				" file must be directly inside " + input_dir.string());
		}

		std::cout << "Loading explicit " << geometry_name << " file: "
			<< absolute_path << std::endl;
		return {absolute_path};
	}

	std::vector<std::filesystem::path> files;
	for (const auto& entry : std::filesystem::directory_iterator(input_dir)) {
		if (entry.is_regular_file() && is_supported_geometry_file(entry.path())) {
			files.push_back(std::filesystem::absolute(entry.path()).lexically_normal());
		}
	}
	std::sort(files.begin(), files.end());

	std::cout << "Using legacy " << geometry_name << " directory scan: "
		<< std::filesystem::absolute(input_dir).lexically_normal() << std::endl;
	for (const std::filesystem::path& path : files) {
		std::cout << "Legacy " << geometry_name << " input file: " << path << std::endl;
	}
	if (files.empty()) {
		throw std::runtime_error(std::string("No supported ") + geometry_name +
			" files found in " + input_dir.string());
	}
	return files;
}

Point3D& cell_center_for_placement(std::vector<Point3D>& cell_origins,
								  Point3D& fallback_origin,
								  std::size_t index) {
	if (cell_origins.empty()) {
		return fallback_origin;
	}
	if (index >= cell_origins.size()) {
		throw std::runtime_error("Cell placement requested a center outside the configured cell geometry");
	}
	return cell_origins[index];
}

constexpr double kSurroundingGeometryTolerance = 1e-12;
constexpr double kTwoPi = 6.2831853071795864769;
constexpr std::size_t kSurroundingReplicaCount = 26;
constexpr std::size_t kSuperparticlePlacementAttempts = 100000;

void require_valid_surrounding_dimensions(const SimMatrix& sim_matrix) {
	if (!(sim_matrix.xLength() > 0.0) ||
		!(sim_matrix.yLength() > 0.0) ||
		!(sim_matrix.zLength() > 0.0) ||
		!std::isfinite(sim_matrix.xLength()) ||
		!std::isfinite(sim_matrix.yLength()) ||
		!std::isfinite(sim_matrix.zLength())) {
		throw std::invalid_argument(
			"Surrounding IONP voxels require positive finite voxel dimensions");
	}
}

bool finite_point(const Point3D& point) {
	return std::isfinite(point.x()) &&
		std::isfinite(point.y()) &&
		std::isfinite(point.z());
}

std::vector<Point3D> surrounding_cell_centers(
	const std::vector<Point3D>& cell_origins,
	const Point3D& fallback_cell_origin,
	double cell_radius) {
	if (!std::isfinite(cell_radius) || cell_radius < 0.0) {
		throw std::invalid_argument(
			"Surrounding-voxel superparticles require a finite non-negative cell radius");
	}
	if (cell_radius == 0.0) {
		return {};
	}

	std::vector<Point3D> centers = cell_origins.empty()
		? std::vector<Point3D>{fallback_cell_origin}
		: cell_origins;
	for (const Point3D& center : centers) {
		if (!finite_point(center)) {
			throw std::invalid_argument(
				"Surrounding-voxel superparticles require finite cell origins");
		}
	}
	return centers;
}

// Translate a template into the 26 adjacent voxel offsets. The central
// physical sources remain separate and are added by the sequence orchestrator.
std::vector<Ionp> surrounding_sources_from_template(
	const std::vector<Ionp>& central_ionps,
	const std::vector<Ionp>& surrounding_template,
	const SimMatrix& sim_matrix) {
	require_valid_surrounding_dimensions(sim_matrix);
	if (surrounding_template.size() >
		(std::numeric_limits<std::size_t>::max() - central_ionps.size()) /
			kSurroundingReplicaCount) {
		throw std::length_error(
			"Surrounding IONP source count exceeds the supported collection size");
	}

	std::vector<Ionp> field_sources;
	field_sources.reserve(
		central_ionps.size() +
		kSurroundingReplicaCount * surrounding_template.size());
	field_sources.insert(
		field_sources.end(), central_ionps.begin(), central_ionps.end());
	for (int x_offset = -1; x_offset <= 1; ++x_offset) {
		for (int y_offset = -1; y_offset <= 1; ++y_offset) {
			for (int z_offset = -1; z_offset <= 1; ++z_offset) {
				if (x_offset == 0 && y_offset == 0 && z_offset == 0) {
					continue;
				}
				const Point3D translation(
					x_offset * sim_matrix.xLength(),
					y_offset * sim_matrix.yLength(),
					z_offset * sim_matrix.zLength());
				for (const Ionp& ionp : surrounding_template) {
					Ionp replica = ionp;
					replica.setX0(ionp.x0() + translation);
					replica.setXt(ionp.xt() + translation);
					field_sources.push_back(replica);
				}
			}
		}
	}
	return field_sources;
}

struct SurroundingIonpGroup {
	std::vector<std::size_t> indices;
	std::optional<std::size_t> cell_index;
};

struct SurroundingIonpPartition {
	std::vector<SurroundingIonpGroup> groups;
	std::vector<std::size_t> group_for_ionp;
};

// Partition the central geometry by physical cell or extracellular space
// before any magnetic-source compression; crossing/ambiguous geometry is rejected.
SurroundingIonpPartition partition_surrounding_ionps(
	const std::vector<Ionp>& ionps,
	const std::vector<Point3D>& cell_centers,
	double cell_radius) {
	SurroundingIonpPartition partition;
	partition.groups.resize(cell_centers.size() + 1);
	partition.group_for_ionp.resize(ionps.size());
	for (std::size_t cell_index = 0;
		 cell_index < cell_centers.size(); ++cell_index) {
		partition.groups[cell_index].cell_index = cell_index;
	}
	const std::size_t extracellular_group = cell_centers.size();

	for (std::size_t ionp_index = 0; ionp_index < ionps.size(); ++ionp_index) {
		const Ionp& ionp = ionps[ionp_index];
		if (!finite_point(ionp.xt()) ||
			!std::isfinite(ionp.radius()) || ionp.radius() < 0.0 ||
			!std::isfinite(ionp.coating()) || ionp.coating() < 0.0 ||
			!std::isfinite(ionp.outer_radius())) {
			throw std::invalid_argument(
				"Surrounding-voxel superparticle input IONP " +
				std::to_string(ionp_index) +
				" has invalid position, radius, or coating geometry");
		}

		std::optional<std::size_t> containing_cell;
		bool straddles_membrane = false;
		for (std::size_t cell_index = 0;
			 cell_index < cell_centers.size(); ++cell_index) {
			const double center_distance =
				distance(ionp.xt(), cell_centers[cell_index]);
			if (center_distance + ionp.outer_radius() <=
				cell_radius + kSurroundingGeometryTolerance) {
				if (containing_cell.has_value()) {
					throw std::runtime_error(
						"IONP " + std::to_string(ionp_index) +
						" is fully contained by more than one cell and cannot be "
						"assigned unambiguously for surrounding-voxel superparticles");
				}
				containing_cell = cell_index;
			}
			else if (center_distance <
				cell_radius + ionp.outer_radius() -
					kSurroundingGeometryTolerance) {
				straddles_membrane = true;
			}
		}
		if (straddles_membrane) {
			throw std::runtime_error(
				"IONP " + std::to_string(ionp_index) +
				" straddles a cell membrane when classified by its full outer "
				"radius for surrounding-voxel superparticles");
		}

		const std::size_t group_index = containing_cell.value_or(
			extracellular_group);
		partition.groups[group_index].indices.push_back(ionp_index);
		partition.group_for_ionp[ionp_index] = group_index;
	}
	return partition;
}

bool overlaps_template_source(
	const Point3D& position,
	double radius,
	const std::vector<Ionp>& template_sources) {
	for (const Ionp& source : template_sources) {
		if (distance(position, source.xt()) <
			radius + source.outer_radius()) {
			return true;
		}
	}
	return false;
}

Point3D random_point_in_sphere(
	const Point3D& center,
	double maximum_center_radius,
	RandNumberBetween& random_unit) {
	const double cosine = 1.0 - 2.0 * random_unit();
	const double sine = std::sqrt(std::max(0.0, 1.0 - cosine * cosine));
	const double azimuth = kTwoPi * random_unit();
	const double radial_distance =
		maximum_center_radius * std::cbrt(random_unit());
	return Point3D(
		center.x() + radial_distance * sine * std::cos(azimuth),
		center.y() + radial_distance * sine * std::sin(azimuth),
		center.z() + radial_distance * cosine);
}

Point3D random_point_in_voxel(
	const SimMatrix& sim_matrix,
	double radius,
	RandNumberBetween& random_unit) {
	const double x_low = -sim_matrix.xLength() / 2.0 + radius;
	const double y_low = -sim_matrix.yLength() / 2.0 + radius;
	const double z_low = -sim_matrix.zLength() / 2.0 + radius;
	return Point3D(
		x_low + (sim_matrix.xLength() - 2.0 * radius) * random_unit(),
		y_low + (sim_matrix.yLength() - 2.0 * radius) * random_unit(),
		z_low + (sim_matrix.zLength() - 2.0 * radius) * random_unit());
}

bool proxy_is_outside_cells(
	const Point3D& position,
	double radius,
	const std::vector<Point3D>& cell_centers,
	double cell_radius,
	const std::optional<std::size_t>& containing_cell) {
	for (std::size_t cell_index = 0;
		 cell_index < cell_centers.size(); ++cell_index) {
		if (containing_cell.has_value() &&
			*containing_cell == cell_index) {
			continue;
		}
		if (distance(position, cell_centers[cell_index]) <
			cell_radius + radius) {
			return false;
		}
	}
	return true;
}

std::string superparticle_group_name(
	const std::optional<std::size_t>& cell_index) {
	return cell_index.has_value()
		? "cell " + std::to_string(*cell_index)
		: "the extracellular group";
}

[[noreturn]] void throw_superparticle_placement_failure(
	const std::optional<std::size_t>& cell_index,
	double proxy_radius) {
	throw std::runtime_error(
		"Unable to place non-overlapping surrounding-voxel superparticles for " +
		superparticle_group_name(cell_index) + " with core radius " +
		std::to_string(proxy_radius) +
		". Increase Surrounding_voxel_superparticles to use smaller proxies.");
}

// Represent one group using equal-radius magnetic proxies, conserving its
// sum of core radii cubed. These proxies are not proton collision obstacles.
void append_compressed_surrounding_group(
	std::vector<Ionp>& template_sources,
	const std::vector<Ionp>& physical_ionps,
	const SurroundingIonpGroup& group,
	std::size_t superparticle_cap,
	const std::vector<Point3D>& cell_centers,
	double cell_radius,
	const SimMatrix& sim_matrix,
	std::optional<RandNumberBetween>& random_unit) {
	long double core_radius_cubed_sum = 0.0L;
	for (std::size_t ionp_index : group.indices) {
		const long double radius = physical_ionps[ionp_index].radius();
		core_radius_cubed_sum += radius * radius * radius;
	}
	if (core_radius_cubed_sum == 0.0L) {
		return;
	}
	if (!std::isfinite(core_radius_cubed_sum)) {
		throw std::invalid_argument(
			"Surrounding-voxel superparticle iron volume is not finite for " +
			superparticle_group_name(group.cell_index));
	}

	const long double proxy_radius_long = std::cbrt(
		core_radius_cubed_sum /
		static_cast<long double>(superparticle_cap));
	const double proxy_radius = static_cast<double>(proxy_radius_long);
	if (!std::isfinite(proxy_radius) || !(proxy_radius > 0.0)) {
		throw std::invalid_argument(
			"Surrounding-voxel superparticle radius is not finite and positive for " +
			superparticle_group_name(group.cell_index));
	}

	if (group.cell_index.has_value()) {
		if (proxy_radius > cell_radius + kSurroundingGeometryTolerance) {
			throw_superparticle_placement_failure(group.cell_index, proxy_radius);
		}
	}
	else if (2.0 * proxy_radius > sim_matrix.xLength() ||
			  2.0 * proxy_radius > sim_matrix.yLength() ||
			  2.0 * proxy_radius > sim_matrix.zLength()) {
		throw_superparticle_placement_failure(group.cell_index, proxy_radius);
	}

	if (!random_unit.has_value()) {
		random_unit.emplace(0.0, 1.0);
	}
	for (std::size_t proxy_index = 0;
		 proxy_index < superparticle_cap; ++proxy_index) {
		bool placed = false;
		for (std::size_t attempt = 0;
			 attempt < kSuperparticlePlacementAttempts; ++attempt) {
			Point3D candidate;
			if (group.cell_index.has_value()) {
				candidate = random_point_in_sphere(
					cell_centers[*group.cell_index],
					std::max(0.0, cell_radius - proxy_radius),
					*random_unit);
			}
			else {
				candidate = random_point_in_voxel(
					sim_matrix, proxy_radius, *random_unit);
			}

			if (!proxy_is_outside_cells(
					candidate,
					proxy_radius,
					cell_centers,
					cell_radius,
					group.cell_index) ||
				overlaps_template_source(
					candidate, proxy_radius, template_sources)) {
				continue;
			}
			template_sources.emplace_back(
				candidate, 0.0, proxy_radius, 0.0);
			placed = true;
			break;
		}
		if (!placed) {
			throw_superparticle_placement_failure(
				group.cell_index, proxy_radius);
		}
	}
}

} // namespace

// Build/load cell centres and cell loads, then place each cell's IONPs.
// c_cell_vec and voxel_pos communicate the realized geometry back to main.
void SimSpace::simulationSpace_inCell(
    std::vector<Ionp>& ionpInstances,
    SimMatrix& simMatrix,
    int ionp_n,
    double conc,
    double r_ionp,
    double r_ionp_1,
    double ionp_per,
    double r_ionp_std,
    double r_ionp_std_1,
    double thickness_coating,
    Point3D& c_cell,
    double r_cell,
    int n_cell,
    double conc_cell,
    double conc_cell_std,
    int conc_cell_dist,
    std::vector<Point3D>& c_cell_vec,
    Point3D& voxel_pos,
    bool conc_hist,
    const std::filesystem::path& ini_root_dir,
    bool use_overlap_grid,
    std::optional<std::filesystem::path> cell_position_path,
    std::optional<std::filesystem::path> cell_load_path) {

	std::cout << "Creating simulation space...." << std::endl;
	const bool supplied_cell_origins = !c_cell_vec.empty();
	bool shift_generated_cell_origins = false;

	std::vector<double> rand_conc;
	std::vector<Point3D> cell_positions;

	if (!supplied_cell_origins && n_cell == 0 && r_cell != 0.){
		if (conc_hist && cell_position_path.has_value() != cell_load_path.has_value()) {
			throw std::runtime_error(
				"Cell_hist=Yes with file-based cells requires selected cell position "
				"and cell load files to be specified together, or both omitted");
		}

		const std::filesystem::path path_to_CellPos = require_ini_subdir(ini_root_dir, "Cell_pos");
		Data dat(path_to_CellPos.string());
		const std::vector<std::filesystem::path> position_files =
			geometry_input_files(path_to_CellPos, cell_position_path, "cell position");
		for (const std::filesystem::path& pathCell : position_files) {
			const std::size_t previous_size = cell_positions.size();
			dat.readCell(pathCell.string(), cell_positions, r_cell);
			if (cell_position_path && cell_positions.size() == previous_size) {
				throw std::runtime_error("Selected cell position file produced no records: " +
					pathCell.string());
			}
		}

		c_cell_vec = cell_positions;
		shift_generated_cell_origins = !c_cell_vec.empty();

		const std::filesystem::path path_to_CellLoad = require_ini_subdir(ini_root_dir, "Cell_load");
		Data dat_load(path_to_CellLoad.string());
		const std::vector<std::filesystem::path> load_files =
			geometry_input_files(path_to_CellLoad, cell_load_path, "cell load");
		for (const std::filesystem::path& pathCell : load_files) {
			const std::size_t previous_size = rand_conc.size();
			dat_load.readCell(pathCell.string(), rand_conc);
			if (cell_load_path && rand_conc.size() == previous_size) {
				throw std::runtime_error("Selected cell load file produced no records: " +
					pathCell.string());
			}
		}
		if (conc_hist && rand_conc.size() != c_cell_vec.size()) {
			throw std::runtime_error(
				"Cell load count (" + std::to_string(rand_conc.size()) +
				") does not match cell position count (" +
				std::to_string(c_cell_vec.size()) + ")");
		}
		r_ionp_1 = 0.;
		r_ionp_std = 0.;
		r_ionp_std_1 = 0.;
		thickness_coating = 0.;
		conc = 0.;
	}
	if (conc_hist && (supplied_cell_origins || n_cell > 0)) {
		// Read the histogram independently of how the authoritative centers were supplied.
		const std::filesystem::path path_to_Hist = require_ini_subdir(ini_root_dir, "Histogram");
		Data dat(path_to_Hist.string());
		for (const auto& entryHist : std::filesystem::directory_iterator(path_to_Hist)) {
			const std::string pathHist = entryHist.path();
			if (pathHist.find(".txt") != std::string::npos || pathHist.find(".dat") != std::string::npos) {
				dat.readHist(pathHist, rand_conc);
			}
			else {
				std::cerr << "No Histogram file found!" << std::endl;
				exit(1);
			}
		}
		r_ionp_1 = 0.;
		r_ionp_std = 0.;
		r_ionp_std_1 = 0.;
		thickness_coating = 0.;
		conc = 0.;
	}

	if (!supplied_cell_origins && n_cell > 0){
		for (unsigned int cellNum = 1; cellNum <= SC_UI(n_cell); ++cellNum) {
			c_cell_vec.push_back(c_cell);
			if (n_cell == 1) {
				// A single configured cell is centered at Cell Origin.  Keeping the
				// vector entry identical makes placement, sampling, and membrane
				// propagation use the same authoritative geometry.
				break;
			}

			while (true) {
				bool hit_cell = false;
				c_cell_vec[cellNum - 1].placeRandom(-(simMatrix.xLength() / 2. - r_cell), (simMatrix.xLength() / 2. - r_cell), -(simMatrix.yLength() / 2. - r_cell), (simMatrix.yLength() / 2. - r_cell), -(simMatrix.zLength() / 2. - r_cell), (simMatrix.zLength() / 2. - r_cell));
				for (unsigned int k = 0; k < cellNum - 1; ++k) {
					if (cellNum == 1) {
						break;//skip this validation for the first ionp
					}
					if (distance(c_cell_vec[cellNum - 1], c_cell_vec[k]) < (2 * r_cell)) {
						hit_cell = true;
						break;
					}
				}
				if (!hit_cell) {
					break;
				}
			}

		}
		shift_generated_cell_origins = n_cell > 1;
	}

	if (!c_cell_vec.empty()){
		n_cell = int(c_cell_vec.size());
	}
	else if (r_cell != 0.){
		n_cell = 1;
	}
	else{
		n_cell = 0;
	}

	voxel_pos = (-1.) * voxel_pos;
	for (unsigned int cellNum = 1; cellNum <= SC_UI(n_cell); ++cellNum) {
		if (shift_generated_cell_origins) {
			c_cell_vec[cellNum - 1] = c_cell_vec[cellNum - 1] + voxel_pos;
		}
		const Point3D& center = cell_center_for_placement(c_cell_vec, c_cell, cellNum - 1);
		std::cout << center.x() << " " << center.y() << " " << center.z() << std::endl;
	}

	RandNumberBase* rand_cell_c = nullptr;

	if (conc_cell_dist == 0) {
		//Use lognormal distribution for cell concentration
		rand_cell_c = new RandNumberLognorm(log(conc_cell), conc_cell_std);
	}
	else if (conc_cell_dist == 1) {
		//Use normal distribution for cell concentration
		rand_cell_c = new RandNumberGaussian(conc_cell, conc_cell_std);
	}

	RandNumberLognorm rand(log(r_ionp*2.), r_ionp_std);
	RandNumberLognorm rand1(log(r_ionp_1*2.), r_ionp_std_1);
	RandNumberBetween rand_norm(0., 1.);

	const double max_r_for_grid = std::max({r_ionp, r_ionp_1, thickness_coating});
	const double grid_cell_size = std::max(1e-9, 2.0 * max_r_for_grid);
	OverlapGrid overlap_grid(grid_cell_size);
	if (use_overlap_grid) {
		for (unsigned int idx = 0; idx < SC_UI(ionpInstances.size()); ++idx) {
			overlap_grid.insert(
				ionpInstances[idx].xt(), idx, ionpInstances[idx].outer_radius());
		}
	} else {
		(void) overlap_grid;
	}

	if ( conc != 0. ){
		double concentration;

		if ( r_cell != 0. ){
			concentration = conc * 5.5845e-5 * 4./3. * D_PI * cube(r_cell) * 1000. * n_cell;
		}
		else{
			concentration = conc * 5.5845e-5 * (simMatrix.xLength()*simMatrix.yLength()*simMatrix.zLength()*1000.); //This is for Iron concentration in the TMAOH and PMAO particles specificly (in µmol/l of Fe) (maybe it would be possible to make this more general)
		}

		unsigned int i = 0, small = 0, big = 0;

		for (unsigned int cellNum = 1; cellNum <= SC_UI(n_cell); ++cellNum) {
			double actual_c = 0.;

			while ( actual_c <= concentration / n_cell ){

				//RandNumberBetween r(-simMatrix.xLength() / 2., simMatrix.xLength() / 2.);
				if ( r_ionp_std == 0. ){
					double coat;
					if ( thickness_coating == 0. ) {
						coat = 0.;
					}
					else{
						coat = ionp_shell_thickness(r_ionp, thickness_coating);
					}
					ionpInstances.push_back(Ionp(simMatrix.x0().x(),simMatrix.x0().y(),simMatrix.x0().z(), 0., r_ionp, coat));
				}
				else{
					if (r_ionp_1 != 0.){
						double prop = rand_norm();
						if (prop <= ionp_per){
							double r = rand()/2.;
							double coat;
							if ( thickness_coating == 0. ) {
								coat = 0.;
							}
							else{
								coat = ionp_shell_thickness(r, thickness_coating);
							}
							ionpInstances.push_back(Ionp(simMatrix.x0().x(),simMatrix.x0().y(),simMatrix.x0().z(), 0., r, coat));
							small += 1;
						}
						else{
							double r = rand1()/2.;
							double coat;
							if ( thickness_coating == 0. ) {
								coat = 0.;
							}
							else{
								coat = ionp_shell_thickness(r, thickness_coating);
							}
							ionpInstances.push_back(Ionp(simMatrix.x0().x(),simMatrix.x0().y(),simMatrix.x0().z(), 0., r, coat));
							big += 1;
						}
					}
					else{
						double r = rand()/2.;
						double coat;
						if ( thickness_coating == 0. ) {
							coat = 0.;
						}
						else{
							coat = ionp_shell_thickness(r, thickness_coating);
						}
						ionpInstances.push_back(Ionp(simMatrix.x0().x(),simMatrix.x0().y(),simMatrix.x0().z(), 0., r, coat));
					}
				}

					if (use_overlap_grid) {
						place_ionp_no_overlap(
							ionpInstances[i],
							ionpInstances,
							overlap_grid,
							[&, i]() {
								if (r_cell != 0.){
									ionpInstances[i].placeRandom(cell_center_for_placement(c_cell_vec, c_cell, cellNum - 1), 0., (r_cell - ionpInstances[i].outer_radius()), true);
								}
								else{
									ionpInstances[i].placeRandom(-(simMatrix.xLength() / 2. - ionpInstances[i].outer_radius()), (simMatrix.xLength() / 2. - ionpInstances[i].outer_radius()), -(simMatrix.yLength() / 2. - ionpInstances[i].outer_radius()), (simMatrix.yLength() / 2. - ionpInstances[i].outer_radius()), -(simMatrix.zLength() / 2. - ionpInstances[i].outer_radius()), (simMatrix.zLength() / 2. - ionpInstances[i].outer_radius()), true);
								}
								return true;
							},
							i
						);
					} else {
						while (true) { // random placement and validation of position
							bool hit_ionp = false; //flag if proton enters IONP
							if (r_cell != 0.){
								ionpInstances[i].placeRandom(cell_center_for_placement(c_cell_vec, c_cell, cellNum - 1), 0., (r_cell - ionpInstances[i].outer_radius()), true);
							}
							else{
								ionpInstances[i].placeRandom(-(simMatrix.xLength() / 2. - ionpInstances[i].outer_radius()), (simMatrix.xLength() / 2. - ionpInstances[i].outer_radius()), -(simMatrix.yLength() / 2. - ionpInstances[i].outer_radius()), (simMatrix.yLength() / 2. - ionpInstances[i].outer_radius()), -(simMatrix.zLength() / 2. - ionpInstances[i].outer_radius()), (simMatrix.zLength() / 2. - ionpInstances[i].outer_radius()), true);
							}

							for (unsigned int k = 0; k < SC_UI(ionpInstances.size()) - 1; ++k) { //not inside IONP // don't check final entry(is this ionp)
								if (i == 0) {
									break;//skip this validation for the first ionp
								}
								if (distance(ionpInstances[i].xt(), ionpInstances[k].xt()) < (ionpInstances[i].outer_radius() + ionpInstances[k].outer_radius())) {
									hit_ionp = true;
									break;
								}
							}
							if (!hit_ionp) {
								break;
							}
						}
					}
				actual_c += 4./3. * D_PI * 5.24 * (167.4/231.4) * cube(ionpInstances[i].radius()*100.);
				i++;
			}
		}

		std::cout << "Small: " << small << std::endl;
		std::cout << "Big: " << big << std::endl;
	}
	else{
		const bool cell_load_controls_radius = conc_hist || conc_cell != 0.;
		for (unsigned int cellNum = 1; cellNum <= SC_UI(n_cell); ++cellNum) {
			if ( conc_hist == true){
				double r_ionp_pow3 = (6.299e-20*rand_conc[cellNum-1])/(1000.*ionp_n);
				r_ionp = std::pow(r_ionp_pow3, 1./3.);
				std::cout << rand_conc[cellNum-1] << std::endl;
			}
			if ( conc_cell != 0. && conc_cell_std != 0. ){
				double cell_c_loc = (*rand_cell_c)();
				while (cell_c_loc < 0.){
					cell_c_loc = (*rand_cell_c)();
				}
				double r_ionp_pow3 = (6.299e-20*cell_c_loc)/(1000.*ionp_n);
				r_ionp = std::pow(r_ionp_pow3, 1./3.);
				std::cout << cell_c_loc << std::endl;
			}
			if ( conc_cell != 0. && conc_cell_std == 0. ){
				double cell_c_loc = conc_cell;
				double r_ionp_pow3 = (6.299e-20*cell_c_loc)/(1000.*ionp_n);
				r_ionp = std::pow(r_ionp_pow3, 1./3.);
				std::cout << cell_c_loc << std::endl;
			}

			for (unsigned int i = 0; i < SC_UI(ionp_n); ++i) {
				//RandNumberBetween r(-simMatrix.xLength() / 2., simMatrix.xLength() / 2.);
				// Cell loads are represented by ionp_n equal iron-equivalent
				// particles.  Their derived radius must take precedence over the
				// configured physical-particle distribution, which remains available
				// to the separate extracellular placement path.
				if ( cell_load_controls_radius || r_ionp_std == 0. ){
					double coat;
					if ( thickness_coating == 0. ) {
						coat = 0.;
					}
					else{
						coat = ionp_shell_thickness(r_ionp, thickness_coating);
					}
					ionpInstances.push_back(Ionp(simMatrix.x0().x(),simMatrix.x0().y(),simMatrix.x0().z(), 0., r_ionp, coat));
				}
				else{
					if (r_ionp_1 != 0.){
						double prop = rand_norm();
						if (prop <= ionp_per){
							double r = rand()/2.;
							double coat;
							if ( thickness_coating == 0. ) {
								coat = 0.;
							}
							else{
								coat = ionp_shell_thickness(r, thickness_coating);
							}
							ionpInstances.push_back(Ionp(simMatrix.x0().x(),simMatrix.x0().y(),simMatrix.x0().z(), 0., r, coat));
						}
						else{
							double r = rand1()/2.;
							double coat;
							if ( thickness_coating == 0. ) {
								coat = 0.;
							}
							else{
								coat = ionp_shell_thickness(r, thickness_coating);
							}
							ionpInstances.push_back(Ionp(simMatrix.x0().x(),simMatrix.x0().y(),simMatrix.x0().z(), 0., r, coat));
						}
					}
					else{
						double r = rand()/2.;
						double coat;
						if ( thickness_coating == 0. ) {
							coat = 0.;
						}
						else{
							coat = ionp_shell_thickness(r, thickness_coating);
						}
						ionpInstances.push_back(Ionp(simMatrix.x0().x(),simMatrix.x0().y(),simMatrix.x0().z(), 0., r, coat));
					}
				}

						if (use_overlap_grid) {
							place_ionp_no_overlap(
								ionpInstances[i + SC_UI(ionp_n)*(cellNum - 1)],
								ionpInstances,
								overlap_grid,
								[&, i, cellNum]() {
									if (r_cell != 0.){
										ionpInstances[i + SC_UI(ionp_n)*(cellNum - 1)].placeRandom(cell_center_for_placement(c_cell_vec, c_cell, cellNum - 1), 0., (r_cell - ionpInstances[i + SC_UI(ionp_n)*(cellNum - 1)].outer_radius()), true);
									}
									else{
										const double placed_outer_radius =
											ionpInstances[i + SC_UI(ionp_n)*(cellNum - 1)].outer_radius();
										ionpInstances[i + SC_UI(ionp_n)*(cellNum - 1)].placeRandom(-(simMatrix.xLength() / 2. - placed_outer_radius), (simMatrix.xLength() / 2. - placed_outer_radius), -(simMatrix.yLength() / 2. - placed_outer_radius), (simMatrix.yLength() / 2. - placed_outer_radius), -(simMatrix.zLength() / 2. - placed_outer_radius), (simMatrix.zLength() / 2. - placed_outer_radius), true);
									}
									return true;
								},
								i + SC_UI(ionp_n)*(cellNum - 1)
							);
						} else {
							while (true) { // random placement and validation of position
								bool hit_ionp = false; //flag if proton enters IONP
								if (r_cell != 0.){
									ionpInstances[i + SC_UI(ionp_n)*(cellNum - 1)].placeRandom(cell_center_for_placement(c_cell_vec, c_cell, cellNum - 1), 0., (r_cell - ionpInstances[i + SC_UI(ionp_n)*(cellNum - 1)].outer_radius()), true);
								}
								else{
									const double placed_outer_radius =
										ionpInstances[i + SC_UI(ionp_n)*(cellNum - 1)].outer_radius();
									ionpInstances[i + SC_UI(ionp_n)*(cellNum - 1)].placeRandom(-(simMatrix.xLength() / 2. - placed_outer_radius), (simMatrix.xLength() / 2. - placed_outer_radius), -(simMatrix.yLength() / 2. - placed_outer_radius), (simMatrix.yLength() / 2. - placed_outer_radius), -(simMatrix.zLength() / 2. - placed_outer_radius), (simMatrix.zLength() / 2. - placed_outer_radius), true);
								}

								for (unsigned int k = 0; k < SC_UI(ionpInstances.size()) - 1; ++k) { //not inside IONP // don't check final entry(is this ionp)
									if (i == 0) {
										break;//skip this validation for the first ionp
									}
									if (distance(ionpInstances[i + SC_UI(ionp_n)*(cellNum - 1)].xt(), ionpInstances[k].xt()) < (ionpInstances[i + SC_UI(ionp_n)*(cellNum - 1)].outer_radius() + ionpInstances[k].outer_radius())) {
										hit_ionp = true;
										break;
									}
								}
								if (!hit_ionp) {
									break;
								}
							}
						}
			}
		}

	}

	delete rand_cell_c;
	std::cout << "IONPs created: " << ionpInstances.size() << std::endl;
}

// Add the extracellular IONP population, avoiding cell interiors and
// previously accepted particles. Optional radius histograms supply sampled sizes.
void SimSpace::simulationSpace_outCell(
    std::vector<Ionp>& ionpInstances,
    SimMatrix& simMatrix,
    int ionp_n,
    double conc,
    double r_ionp,
    double r_ionp_1,
    double ionp_per,
    double r_ionp_std,
    double r_ionp_std_1,
    double thickness_coating,
    Point3D& c_cell,
    double r_cell,
    int n_cell,
    std::vector<Point3D>& c_cell_vec,
    Point3D& voxel_pos,
    bool radius_hist,
    const std::filesystem::path& ini_root_dir,
    bool use_overlap_grid) {

	std::cout << "Creating simulation space outside of cells...." << std::endl;

	(void) voxel_pos;

	RandNumberLognorm rand(log(r_ionp*2.), r_ionp_std);
	RandNumberLognorm rand1(log(r_ionp_1*2.), r_ionp_std_1);
	RandNumberBetween rand_norm(0., 1.);

	if (c_cell_vec.size() > 0){
		n_cell = int(c_cell_vec.size());
	}
	else if ( r_cell != 0. ){
		n_cell = 1;
	}
	else{
		n_cell = 0;
	}

	std::vector<double> rand_hist_rad;

	if ( radius_hist == true ){
		//Read the Histogram File for IONP radius distribution
		const std::filesystem::path path_to_Hist = require_ini_subdir(ini_root_dir, "Radius_Histogram");

		Data dat(path_to_Hist.string());

		for (const auto & entryHist : std::filesystem::directory_iterator(path_to_Hist)){
			std::string pathHist = entryHist.path();
			if (pathHist.find(".txt") != std::string::npos || pathHist.find(".dat") != std::string::npos ){
				dat.readHist( pathHist, rand_hist_rad );
			}
			else{
				std::cerr << "No Histogram file found!" << std::endl;
				exit(1);
			}
		}
	}

	const double max_r_for_grid = std::max({r_ionp, r_ionp_1, thickness_coating});
	const double grid_cell_size = std::max(1e-9, 2.0 * max_r_for_grid);
	OverlapGrid overlap_grid(grid_cell_size);
	if (use_overlap_grid) {
		for (unsigned int idx = 0; idx < SC_UI(ionpInstances.size()); ++idx) {
			overlap_grid.insert(
				ionpInstances[idx].xt(), idx, ionpInstances[idx].outer_radius());
		}
	} else {
		(void) overlap_grid;
	}

	if ( conc != 0. ){
		double concentration;

		concentration = conc * 5.5845e-5 * (simMatrix.xLength()*simMatrix.yLength()*simMatrix.zLength()*1000.);

		unsigned int i = uint(ionpInstances.size());
		unsigned int small = 0, big = 0;

		double actual_c = 0.;
		while ( actual_c <= concentration ){
			if ( r_ionp_std == 0. ){
				if (radius_hist == true){
					r_ionp = rand_hist_rad[i];
					std::cout << r_ionp << std::endl;
					if (i >= rand_hist_rad.size()){
						std::cerr << "Not enough entries in the radius histogram for the desired concentration!" << std::endl;
						exit(1);
					}
				}
				double coat;
				if ( thickness_coating == 0. ) {
					coat = 0.;
				}
				else{
					coat = ionp_shell_thickness(r_ionp, thickness_coating);
				}
				ionpInstances.push_back(Ionp(simMatrix.x0().x(),simMatrix.x0().y(),simMatrix.x0().z(), 0., r_ionp, coat));
			}
			else{
				if (r_ionp_1 != 0.){
					double prop = rand_norm();
					if (prop <= ionp_per){
						double r = rand()/2.;
						double coat;
						if ( thickness_coating == 0. ) {
							coat = 0.;
						}
						else{
							coat = ionp_shell_thickness(r, thickness_coating);
						}
						ionpInstances.push_back(Ionp(simMatrix.x0().x(),simMatrix.x0().y(),simMatrix.x0().z(), 0., r, coat));
						small += 1;
					}
					else{
						double r = rand1()/2.;
						double coat;
						if ( thickness_coating == 0. ) {
							coat = 0.;
						}
						else{
							coat = ionp_shell_thickness(r, thickness_coating);
						}
						ionpInstances.push_back(Ionp(simMatrix.x0().x(),simMatrix.x0().y(),simMatrix.x0().z(), 0., r, coat));
						big += 1;
					}
				}
				else{
					double r = rand()/2.;
					double coat;
					if ( thickness_coating == 0. ) {
						coat = 0.;
					}
					else{
						coat = ionp_shell_thickness(r, thickness_coating);
					}
					ionpInstances.push_back(Ionp(simMatrix.x0().x(),simMatrix.x0().y(),simMatrix.x0().z(), 0., r, coat));
				}
			}

			if (use_overlap_grid) {
				place_ionp_no_overlap(
				ionpInstances[i],
				ionpInstances,
				overlap_grid,
				[&, i]() {
					ionpInstances[i].placeRandom(-(simMatrix.xLength() / 2. - ionpInstances[i].outer_radius()), (simMatrix.xLength() / 2. - ionpInstances[i].outer_radius()), -(simMatrix.yLength() / 2. - ionpInstances[i].outer_radius()), (simMatrix.yLength() / 2. - ionpInstances[i].outer_radius()), -(simMatrix.zLength() / 2. - ionpInstances[i].outer_radius()), (simMatrix.zLength() / 2. - ionpInstances[i].outer_radius()), true);
					// Reject placements that would intersect existing cells
					for (unsigned int cellNum = 1; cellNum <= SC_UI(n_cell); ++cellNum) {
						if (distance(cell_center_for_placement(c_cell_vec, c_cell, cellNum - 1), ionpInstances[i].xt()) < (r_cell + ionpInstances[i].outer_radius())) {
							return false;
						}
					}
					return true;
				},
							i
				);
			} else {
				while (true) { // random placement and validation of position
					bool hit_ionp = false; //flag if proton enters IONP
					ionpInstances[i].placeRandom(-(simMatrix.xLength() / 2. - ionpInstances[i].outer_radius()), (simMatrix.xLength() / 2. - ionpInstances[i].outer_radius()), -(simMatrix.yLength() / 2. - ionpInstances[i].outer_radius()), (simMatrix.yLength() / 2. - ionpInstances[i].outer_radius()), -(simMatrix.zLength() / 2. - ionpInstances[i].outer_radius()), (simMatrix.zLength() / 2. - ionpInstances[i].outer_radius()), true);

					for (unsigned int cellNum = 1; cellNum <= SC_UI(n_cell); ++cellNum) {
						if (distance(cell_center_for_placement(c_cell_vec, c_cell, cellNum - 1), ionpInstances[i].xt()) < (r_cell + ionpInstances[i].outer_radius())) {
							hit_ionp = true;
							break;
						}
					}
					if (!hit_ionp) {
						for (unsigned int k = 0; k < SC_UI(ionpInstances.size()) - 1; ++k) { //not inside IONP // don't check final entry(is this ionp)
							if (i == 0) {
								break;//skip this validation for the first ionp
							}
							if (distance(ionpInstances[i].xt(), ionpInstances[k].xt()) < (ionpInstances[i].outer_radius() + ionpInstances[k].outer_radius())) {
								hit_ionp = true;
								break;
							}
						}
					}
					if (!hit_ionp) {
						break;
					}
				}
			}
			actual_c += 4./3. * D_PI * 5.24 * (167.4/231.4) * cube(ionpInstances[i].radius()*100.);
			i++;
		}

		std::cout << "Small: " << small << std::endl;
		std::cout << "Big: " << big << std::endl;
	}
	else{
		unsigned int ionpInstances_len = uint(ionpInstances.size());
		unsigned int initial_ionp_n = SC_UI(ionp_n) + ionpInstances_len;
		for (unsigned int i = ionpInstances_len; i < initial_ionp_n; ++i) {
			if ( r_ionp_std == 0. ){
				if (radius_hist == true){
					r_ionp = rand_hist_rad[i];
					std::cout << r_ionp << std::endl;
					if (i >= rand_hist_rad.size()){
						std::cerr << "Not enough entries in the radius histogram for the desired IONP number!" << std::endl;
						exit(1);
					}
				}
				double coat;
				if ( thickness_coating == 0. ) {
					coat = 0.;
				}
				else{
					coat = ionp_shell_thickness(r_ionp, thickness_coating);
				}
				ionpInstances.push_back(Ionp(simMatrix.x0().x(),simMatrix.x0().y(),simMatrix.x0().z(), 0., r_ionp, coat));
			}
			else{
				if (r_ionp_1 != 0.){
					double prop = rand_norm();
					if (prop <= ionp_per){
						double r = rand()/2.;
						double coat;
						if ( thickness_coating == 0. ) {
							coat = 0.;
						}
						else{
							coat = ionp_shell_thickness(r, thickness_coating);
						}
						ionpInstances.push_back(Ionp(simMatrix.x0().x(),simMatrix.x0().y(),simMatrix.x0().z(), 0., r, coat));
					}
					else{
						double r = rand1()/2.;
						double coat;
						if ( thickness_coating == 0. ) {
							coat = 0.;
						}
						else{
							coat = ionp_shell_thickness(r, thickness_coating);
						}
						ionpInstances.push_back(Ionp(simMatrix.x0().x(),simMatrix.x0().y(),simMatrix.x0().z(), 0., r, coat));
					}
				}
				else{
					double r = rand()/2.;
					double coat;
					if ( thickness_coating == 0. ) {
						coat = 0.;
					}
					else{
						coat = ionp_shell_thickness(r, thickness_coating);
					}
					ionpInstances.push_back(Ionp(simMatrix.x0().x(),simMatrix.x0().y(),simMatrix.x0().z(), 0., r, coat));
				}
			}

					if (use_overlap_grid) {
						place_ionp_no_overlap(
							ionpInstances[i],
							ionpInstances,
							overlap_grid,
							[&, i]() {
								ionpInstances[i].placeRandom(-(simMatrix.xLength() / 2. - ionpInstances[i].outer_radius()), (simMatrix.xLength() / 2. - ionpInstances[i].outer_radius()), -(simMatrix.yLength() / 2. - ionpInstances[i].outer_radius()), (simMatrix.yLength() / 2. - ionpInstances[i].outer_radius()), -(simMatrix.zLength() / 2. - ionpInstances[i].outer_radius()), (simMatrix.zLength() / 2. - ionpInstances[i].outer_radius()), true);
								for (unsigned int cellNum = 1; cellNum <= SC_UI(n_cell); ++cellNum) {
									if (distance(cell_center_for_placement(c_cell_vec, c_cell, cellNum - 1), ionpInstances[i].xt()) < (r_cell + ionpInstances[i].outer_radius())) {
										return false;
									}
								}
								return true;
							},
							i
						);
					} else {
						while (true) { // random placement and validation of position
							bool hit_ionp = false; //flag if proton enters IONP
							ionpInstances[i].placeRandom(-(simMatrix.xLength() / 2. - ionpInstances[i].outer_radius()), (simMatrix.xLength() / 2. - ionpInstances[i].outer_radius()), -(simMatrix.yLength() / 2. - ionpInstances[i].outer_radius()), (simMatrix.yLength() / 2. - ionpInstances[i].outer_radius()), -(simMatrix.zLength() / 2. - ionpInstances[i].outer_radius()), (simMatrix.zLength() / 2. - ionpInstances[i].outer_radius()), true);
							for (unsigned int cellNum = 1; cellNum <= SC_UI(n_cell); ++cellNum) {
								if (distance(cell_center_for_placement(c_cell_vec, c_cell, cellNum - 1), ionpInstances[i].xt()) < (r_cell + ionpInstances[i].outer_radius())) {
									hit_ionp = true;
									break;
								}
							}
							if (!hit_ionp) {
								for (unsigned int k = 0; k < SC_UI(ionpInstances.size()) - 1; ++k) {
									if (i == 0) {
										break;//skip this validation for the first ionp
									}
									if (distance(ionpInstances[i].xt(), ionpInstances[k].xt()) < (ionpInstances[i].outer_radius() + ionpInstances[k].outer_radius())) {
										hit_ionp = true;
										break;
									}
								}
							}
							if (!hit_ionp) {
								break;
							}
						}
					}
		}
	}

	std::cout << "Total IONPs created: " << ionpInstances.size() << std::endl;
}

// Place extracellular aggregates first, then distribute their IONPs while
// retaining the configured radius distribution and coating conversion.
void SimSpace::simulationSpace_outCellAggregates(
	std::vector<Aggregate>& aggregateInstances,
	std::vector<Ionp>& ionpInstances,
	SimMatrix& simMatrix,
	int agg_n,
	double conc_ext,
	double r_agg,
	double r_agg_std,
	double r_ionp,
	double r_ionp_1,
	double ionp_per,
	double r_ionp_std,
	double r_ionp_std_1,
	double thickness_coating,
	const Point3D& c_cell,
	double r_cell,
	const std::vector<Point3D>& c_cell_vec) {
	std::cout << "Creating extracellular spherical aggregates...." << std::endl;

	if (agg_n <= 0 || !(conc_ext > 0.) || !(r_agg > 0.) ||
		!std::isfinite(conc_ext) || !std::isfinite(r_agg) ||
		!std::isfinite(r_agg_std) || r_agg_std < 0.) {
		throw std::invalid_argument(
			"Extracellular aggregates require positive finite Agg_num, Conc_ext, "
			"and Agg_radius, plus a finite non-negative Agg_radius_std");
	}
	if (!(simMatrix.xLength() > 2. * r_agg) ||
		!(simMatrix.yLength() > 2. * r_agg) ||
		!(simMatrix.zLength() > 2. * r_agg)) {
		throw std::runtime_error(
			"The configured extracellular aggregate radius does not fit inside the simulation voxel");
	}

	std::vector<Point3D> cell_centers = c_cell_vec;
	if (cell_centers.empty() && r_cell > 0.) {
		cell_centers.push_back(c_cell);
	}

	RandNumberBetween random_unit(0., 1.);
	RandNumberLognorm random_aggregate_radius(log(r_agg * 2.), r_agg_std);
	constexpr unsigned int aggregate_placement_attempts = 100000;
	const std::size_t first_aggregate = aggregateInstances.size();

	for (int aggregate_number = 0; aggregate_number < agg_n; ++aggregate_number) {
		const double aggregate_radius = r_agg_std == 0.
			? r_agg
			: random_aggregate_radius() / 2.;
		if (!(aggregate_radius > 0.) || !std::isfinite(aggregate_radius) ||
			simMatrix.xLength() <= 2. * aggregate_radius ||
			simMatrix.yLength() <= 2. * aggregate_radius ||
			simMatrix.zLength() <= 2. * aggregate_radius) {
			throw std::runtime_error(
				"A sampled extracellular aggregate radius does not fit inside the simulation voxel");
		}

		bool placed = false;
		for (unsigned int attempt = 0;
			 attempt < aggregate_placement_attempts;
			 ++attempt) {
			const Point3D candidate = random_point_in_voxel(
				simMatrix, aggregate_radius, random_unit);
			bool valid = true;
			for (const Point3D& cell_center : cell_centers) {
				if (distance(candidate, cell_center) <
					r_cell + aggregate_radius) {
					valid = false;
					break;
				}
			}
			if (!valid) {
				continue;
			}
			for (const Aggregate& aggregate : aggregateInstances) {
				if (distance(candidate, aggregate.xt()) <
					aggregate_radius + aggregate.radius()) {
					valid = false;
					break;
				}
			}
			if (!valid) {
				continue;
			}

			aggregateInstances.emplace_back(
				candidate, 0., aggregate_radius, 0);
			placed = true;
			break;
		}
		if (!placed) {
			throw std::runtime_error(
				"Unable to place a non-overlapping extracellular aggregate outside all cells");
		}
	}

	RandNumberLognorm random_primary(log(r_ionp * 2.), r_ionp_std);
	RandNumberLognorm random_secondary(log(r_ionp_1 * 2.), r_ionp_std_1);
	RandNumberBetween random_mode(0., 1.);
	const double total_concentration_mass =
		conc_ext * 5.5845e-5 *
		(simMatrix.xLength() * simMatrix.yLength() *
		 simMatrix.zLength() * 1000.);
	const double target_mass_per_aggregate =
		total_concentration_mass / static_cast<double>(agg_n);
	const double configured_outer_radius = std::max({
		r_ionp, r_ionp_1, thickness_coating});
	OverlapGrid overlap_grid(std::max(1e-9, 2. * configured_outer_radius));
	for (unsigned int index = 0; index < SC_UI(ionpInstances.size()); ++index) {
		overlap_grid.insert(
			ionpInstances[index].xt(), index,
			ionpInstances[index].outer_radius());
	}

	std::size_t total_primary = 0;
	std::size_t total_secondary = 0;
	std::size_t total_aggregate_ionps = 0;
	for (std::size_t aggregate_index = first_aggregate;
		 aggregate_index < aggregateInstances.size();
		 ++aggregate_index) {
		const Point3D aggregate_center =
			aggregateInstances[aggregate_index].x0();
		std::size_t aggregate_ionps = 0;
		double sampled_mass = 0.;
		while (sampled_mass < target_mass_per_aggregate) {
			bool primary_mode = true;
			double core_radius = r_ionp;
			if (r_ionp_std != 0.) {
				if (r_ionp_1 != 0. && random_mode() > ionp_per) {
					primary_mode = false;
					core_radius = random_secondary() / 2.;
				}
				else {
					core_radius = random_primary() / 2.;
				}
			}
			if (!(core_radius > 0.) || !std::isfinite(core_radius)) {
				throw std::runtime_error(
					"The extracellular aggregate radius distribution produced an invalid IONP radius");
			}
			const double coating = ionp_shell_thickness(
				core_radius, thickness_coating);
			const double outer_radius = core_radius + coating;
			if (outer_radius > aggregateInstances[aggregate_index].radius()) {
				throw std::runtime_error(
					"A sampled coated IONP is larger than its extracellular aggregate");
			}
			ionpInstances.emplace_back(
				aggregate_center, 0., core_radius, coating);
			const unsigned int ionp_index =
				SC_UI(ionpInstances.size() - 1);
			const double placement_radius =
				aggregateInstances[aggregate_index].radius() -
				outer_radius;
			try {
				place_ionp_no_overlap(
					ionpInstances[ionp_index],
					ionpInstances,
					overlap_grid,
					[&, ionp_index, placement_radius]() {
						ionpInstances[ionp_index].placeRandomSphere(
							aggregate_center.x(),
							aggregate_center.y(),
							aggregate_center.z(),
							0.,
							placement_radius,
							true);
						return true;
					},
					ionp_index);
			}
			catch (const std::runtime_error&) {
				throw std::runtime_error(
					"Unable to place all coated IONPs randomly inside extracellular aggregate " +
					std::to_string(aggregate_index - first_aggregate) +
					" without overlap. Increase Agg_radius or Agg_num.");
			}
			if (primary_mode) {
				++total_primary;
			}
			else {
				++total_secondary;
			}
			sampled_mass += 4. / 3. * D_PI * 5.24 *
				(167.4 / 231.4) * cube(core_radius * 100.);
			++aggregate_ionps;
		}

		if (aggregate_ionps >
			static_cast<std::size_t>(std::numeric_limits<int>::max())) {
			throw std::length_error(
				"An extracellular aggregate contains too many IONPs");
		}
		aggregateInstances[aggregate_index] = Aggregate(
			aggregate_center, 0.,
			aggregateInstances[aggregate_index].radius(),
			static_cast<int>(aggregate_ionps));
		total_aggregate_ionps += aggregate_ionps;
	}

	std::cout << "Extracellular aggregates created: "
		<< agg_n << std::endl;
	std::cout << "Extracellular aggregate IONPs created: "
		<< total_aggregate_ionps << std::endl;
	std::cout << "Extracellular aggregate bimodal counts (primary / secondary): "
		<< total_primary << " / " << total_secondary << std::endl;
}

std::vector<Ionp> SimSpace::surroundingSpace(
	const std::vector<Ionp>& ionpInstances,
	const SimMatrix& simMatrix) {
	return surrounding_sources_from_template(
		ionpInstances, ionpInstances, simMatrix);
}

std::vector<Ionp> SimSpace::surroundingSpace(
	const std::vector<Ionp>& ionpInstances,
	const SimMatrix& simMatrix,
	std::size_t superparticle_cap,
	const std::vector<Point3D>& cell_origins,
	const Point3D& fallback_cell_origin,
	double cell_radius) {
	if (superparticle_cap == 0) {
		return surroundingSpace(ionpInstances, simMatrix);
	}
	require_valid_surrounding_dimensions(simMatrix);
	const std::vector<Point3D> cell_centers = surrounding_cell_centers(
		cell_origins, fallback_cell_origin, cell_radius);
	const SurroundingIonpPartition partition = partition_surrounding_ionps(
		ionpInstances, cell_centers, cell_radius);

	std::vector<bool> compressed_group(partition.groups.size(), false);
	for (std::size_t group_index = 0;
		 group_index < partition.groups.size(); ++group_index) {
		compressed_group[group_index] =
			partition.groups[group_index].indices.size() > superparticle_cap;
	}

	std::vector<Ionp> surrounding_template;
	surrounding_template.reserve(ionpInstances.size());
	for (std::size_t ionp_index = 0;
		 ionp_index < ionpInstances.size(); ++ionp_index) {
		if (!compressed_group[partition.group_for_ionp[ionp_index]]) {
			surrounding_template.push_back(ionpInstances[ionp_index]);
		}
	}

	std::optional<RandNumberBetween> random_unit;
	for (std::size_t group_index = 0;
		 group_index < partition.groups.size(); ++group_index) {
		if (!compressed_group[group_index]) {
			continue;
		}
		append_compressed_surrounding_group(
			surrounding_template,
			ionpInstances,
			partition.groups[group_index],
			superparticle_cap,
			cell_centers,
			cell_radius,
			simMatrix,
			random_unit);
	}

	return surrounding_sources_from_template(
		ionpInstances, surrounding_template, simMatrix);
}

void SimSpace::simulationSpace(
    std::vector<Aggregate>& aggregateInstances,
    std::vector<Ionp>& ionpInstances,
    SimMatrix& simMatrix,
    int agg_n,
    int ionp_n,
    double conc,
    double r_agg,
    double r_ionp,
    double r_ionp_1,
    double ionp_per,
    double r_ionp_std,
    double r_ionp_std_1,
    double r_agg_std,
    double thickness_coating,
    Point3D& c_cell,
    double r_cell,
    int n_cell,
    std::vector<Point3D>& c_cell_vec,
    bool conc_hist,
    const std::filesystem::path& ini_root_dir,
    std::optional<std::filesystem::path> cell_position_path) {
	std::cout << "Creating simulation space containing spherical aggregates...." << std::endl;

	std::vector<double> rand_conc;

	if ( conc_hist == true ){
		//Read the Histogram File from single cell ICP-MS
		const std::filesystem::path path_to_Hist = require_ini_subdir(ini_root_dir, "Histogram");

		Data dat(path_to_Hist.string());

		for (const auto & entryHist : std::filesystem::directory_iterator(path_to_Hist)){
			std::string pathHist = entryHist.path();
			if (pathHist.find(".txt") != std::string::npos || pathHist.find(".dat") != std::string::npos ){
				dat.readHist( pathHist, rand_conc );
			}
			else{
				std::cerr << "No Histogram file found!" << std::endl;
				exit(1);
			}
		}
	}

	std::vector<Point3D> cell_positions;
	const bool supplied_cell_origins = !c_cell_vec.empty();

	if (supplied_cell_origins) {
		n_cell = int(c_cell_vec.size());
	}
	else if (n_cell == 0 && r_cell != 0.){
		const std::filesystem::path path_to_CellPos = require_ini_subdir(ini_root_dir, "Cell_pos");
		Data dat(path_to_CellPos.string());
		const std::vector<std::filesystem::path> position_files =
			geometry_input_files(path_to_CellPos, cell_position_path, "cell position");
		for (const std::filesystem::path& pathCell : position_files) {
			const std::size_t previous_size = cell_positions.size();
			dat.readCell(pathCell.string(), cell_positions, r_cell);
			if (cell_position_path && cell_positions.size() == previous_size) {
				throw std::runtime_error("Selected cell position file produced no records: " +
					pathCell.string());
			}
		}

		c_cell_vec = cell_positions;
		n_cell = int(cell_positions.size());
	}
	else{
		for (unsigned int cellNum = 1; cellNum <= SC_UI(n_cell); ++cellNum) {
			c_cell_vec.push_back(c_cell);
			if (n_cell == 1) {
				// Match the single-cell geometry used by the non-aggregate path.
				break;
			}

			while (true) {
				bool hit_cell = false;
				c_cell_vec[cellNum - 1].placeRandom(-(simMatrix.xLength() / 2. - r_cell), (simMatrix.xLength() / 2. - r_cell), -(simMatrix.yLength() / 2. - r_cell), (simMatrix.yLength() / 2. - r_cell), -(simMatrix.zLength() / 2. - r_cell), (simMatrix.zLength() / 2. - r_cell));
				for (unsigned int k = 0; k < cellNum - 1; ++k) {
					if (cellNum == 1) {
						break;//skip this validation for the first ionp
					}
					if (distance(c_cell_vec[cellNum - 1], c_cell_vec[k]) < (2 * r_cell)) {
						hit_cell = true;
						break;
					}
				}
				if (!hit_cell) {
					break;
				}
			}
			std::cout << c_cell_vec[cellNum - 1].x() << " " << c_cell_vec[cellNum - 1].y() << " " << c_cell_vec[cellNum - 1].z() << std::endl;
		}
	}

	if (!c_cell_vec.empty()) {
		n_cell = int(c_cell_vec.size());
	}
	else if (r_cell != 0.) {
		n_cell = 1;
	}
	else {
		n_cell = 0;
	}

	RandNumberLognorm rand(log(r_ionp*2.), r_ionp_std);
	RandNumberLognorm rand1(log(r_ionp_1*2.), r_ionp_std_1);
	RandNumberBetween rand_norm(0., 1.);
	RandNumberLognorm randAgg(log(r_agg*2.), r_agg_std);

	double r_agg_rand;

	if (aggregateInstances.empty()){
		for (unsigned int cellNum = 1; cellNum <= SC_UI(n_cell); ++cellNum) {

			if ( conc_hist == true){
				double r_ionp_pow3 = (6.299e-20*rand_conc[cellNum-1])/(1000.*ionp_n);
				r_ionp = std::pow(r_ionp_pow3, 1./3.);
				std::cout << rand_conc[cellNum-1] << std::endl;
			}

			for (unsigned int i = 0; i < SC_UI(agg_n); ++i) {

				if ( r_agg_std == 0. ){
					aggregateInstances.push_back(Aggregate(r_agg, ionp_n));
					if (r_cell != 0.){
						aggregateInstances[i + SC_UI(agg_n)*(cellNum - 1)].placeRandom(cell_center_for_placement(c_cell_vec, c_cell, cellNum - 1), 0., (r_cell - r_agg), true);
					}
					else{

						aggregateInstances[i].placeRandom(-(simMatrix.xLength() / 2. - r_agg), (simMatrix.xLength() / 2. - r_agg), -(simMatrix.yLength() / 2. - r_agg), (simMatrix.yLength() / 2. - r_agg), -(simMatrix.zLength() / 2. - r_agg), (simMatrix.zLength() / 2. - r_agg), true);
					}
				}
				else{
					r_agg_rand = randAgg()/2.;
					aggregateInstances.push_back(Aggregate(r_agg_rand, ionp_n));
					if (r_cell != 0.){
						aggregateInstances[i + SC_UI(agg_n)*(cellNum - 1)].placeRandom(cell_center_for_placement(c_cell_vec, c_cell, cellNum - 1), 0., (r_cell - r_agg_rand), true);
					}
					else{
						aggregateInstances[i].placeRandom(-(simMatrix.xLength() / 2. - r_agg_rand), (simMatrix.xLength() / 2. - r_agg_rand), -(simMatrix.yLength() / 2. - r_agg_rand), (simMatrix.yLength() / 2. - r_agg_rand), -(simMatrix.zLength() / 2. - r_agg_rand), (simMatrix.zLength() / 2. - r_agg_rand), true);
					}
				}
			}
		}
	}

	const double max_outer_for_grid = std::max({r_ionp, r_ionp_1, thickness_coating});
	OverlapGrid overlap_grid(std::max(1e-9, 2.0 * max_outer_for_grid));
	for (unsigned int idx = 0; idx < SC_UI(ionpInstances.size()); ++idx) {
		overlap_grid.insert(
			ionpInstances[idx].xt(), idx, ionpInstances[idx].outer_radius());
	}

	auto validate_aggregate_inside_cell = [&](const Aggregate& aggregate, unsigned int aggregate_index) {
		if (r_cell == 0.) {
			return;
		}
		for (unsigned int cellNum = 1; cellNum <= SC_UI(n_cell); ++cellNum) {
			const Point3D& center = cell_center_for_placement(c_cell_vec, c_cell, cellNum - 1);
			if (distance(center, aggregate.x0()) <= r_cell - aggregate.radius()) {
				return;
			}
		}
		throw std::runtime_error("Aggregate " + std::to_string(aggregate_index) + " is not fully inside any configured cell");
	};

	auto append_ionp_at = [&](const Point3D& center, unsigned int& small, unsigned int& big) {
		double r = r_ionp;
		if (r_ionp_std != 0.) {
			if (r_ionp_1 != 0.) {
				const double prop = rand_norm();
				if (prop <= ionp_per) {
					r = rand() / 2.;
					small += 1;
				}
				else {
					r = rand1() / 2.;
					big += 1;
				}
			}
			else {
				r = rand() / 2.;
			}
		}
		const double coat = ionp_shell_thickness(r, thickness_coating);
		ionpInstances.push_back(Ionp(center.x(), center.y(), center.z(), 0., r, coat));
		return SC_UI(ionpInstances.size() - 1);
	};

	auto place_ionp_in_aggregate = [&](unsigned int ionp_index, const Aggregate& aggregate) {
		Ionp& ionp = ionpInstances[ionp_index];
		const double outer_radius = ionp.outer_radius();
		const double placement_radius = aggregate.radius() - outer_radius;
		if (placement_radius < 0.) {
			throw std::runtime_error("Aggregate radius is smaller than an IONP outer radius");
		}
		place_ionp_no_overlap(
			ionp,
			ionpInstances,
			overlap_grid,
			[&, ionp_index, placement_radius]() {
				ionpInstances[ionp_index].placeRandomSphere(
					aggregate.x0().x(),
					aggregate.x0().y(),
					aggregate.x0().z(),
					0.,
					placement_radius,
					true);
				return true;
			},
			ionp_index
		);
	};

	const double total_concentration_mass =
		(r_cell != 0.)
			? conc * 5.5845e-5 * (4. / 3. * D_PI * cube(r_cell * 10.)) * std::max(1, n_cell)
			: conc * 5.5845e-5 * (simMatrix.xLength() * simMatrix.yLength() * simMatrix.zLength() * 1000.);
	const double target_mass_per_aggregate =
		(conc != 0. && !aggregateInstances.empty())
			? total_concentration_mass / static_cast<double>(aggregateInstances.size())
			: 0.;

	for (unsigned int i = 0; i < SC_UI(aggregateInstances.size()); ++i) {
		validate_aggregate_inside_cell(aggregateInstances[i], i);
		unsigned int small = 0;
		unsigned int big = 0;

		if (conc != 0.) {
			double actual_c = 0.;
			while (actual_c <= target_mass_per_aggregate) {
				const unsigned int ionp_index = append_ionp_at(aggregateInstances[i].x0(), small, big);
				place_ionp_in_aggregate(ionp_index, aggregateInstances[i]);
				actual_c += 4. / 3. * D_PI * 5.24 * (167.4 / 231.4) * cube(ionpInstances[ionp_index].radius() * 100.);
			}
			std::cout << "Aggregate " << i << ":" << std::endl;
			std::cout << "Small: " << small << std::endl;
			std::cout << "Big: " << big << std::endl;
		}
		else {
			for (unsigned int j = 0; j < SC_UI(ionp_n); ++j) {
				const unsigned int ionp_index = append_ionp_at(aggregateInstances[i].x0(), small, big);
				place_ionp_in_aggregate(ionp_index, aggregateInstances[i]);
			}
		}
	}

	std::cout << "Aggregates created: " << aggregateInstances.size() << std::endl;
	std::cout << "IONPs created: " << ionpInstances.size() << std::endl;
}

void SimSpace::simulationSpace(
    std::vector<Aggregate>& aggregateInstances,
    std::vector<Ionp>& ionpInstances,
    SimMatrix& simMatrix,
    int agg_n,
    int ionp_n,
    double r_ionp) {
	std::cout << "Creating simulation space spherical dense aggregates...." << std::endl;

	double layers = ceil(ionp_n/26);

	double r_agg = (layers + 0.5) * r_ionp;

	std::cout << "Approximated Aggregate Radius = " << r_agg << std::endl;

	for (unsigned int i = 0; i < SC_UI(agg_n); ++i) {
		aggregateInstances.push_back(Aggregate(ionp_n));

		while (true) { // random placement and validation of position
				bool hit_ionp = false; //flag if proton enters IONP

				aggregateInstances[i].placeRandom(-simMatrix.xLength() / 2., simMatrix.xLength() / 2.,-simMatrix.yLength() / 2., simMatrix.yLength() / 2.,-simMatrix.zLength() / 2., simMatrix.zLength() / 2., true);

				for (unsigned int k = 0; k < SC_UI(aggregateInstances.size()) - 1; ++k) { //not inside IONP // don't check final entry(is this ionp)
					if (i == 0) {
						break;//skip this validation for the first ionp
					}
					if (distance(aggregateInstances[i].xt(), aggregateInstances[k].xt()) < (2. * r_agg)) {
						hit_ionp = true;
						break;
					}
				}
				if (!hit_ionp) {
					break;
				}
		}

		unsigned int k = 1;

		for (unsigned int j = 0; j < SC_UI(ionp_n); ++j) {

			if (!j){
				ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z(), 0., r_ionp));
			}
			else if ( j == 1 + (k - 1)*26 ){
				ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x() + 2*k*r_ionp, aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z(), 0., r_ionp));
			}
			else if ( j == 2 + (k - 1)*26 ){
				ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x() - 2*k*r_ionp, aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z(), 0., r_ionp));
			}
			else if ( j == 3 + (k - 1)*26 ){
				ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y() + 2*k*r_ionp, aggregateInstances[i].x0().z(), 0., r_ionp));
			}
			else if ( j == 4 + (k - 1)*26 ){
				ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y() - 2*k*r_ionp, aggregateInstances[i].x0().z(), 0., r_ionp));
			}
			else if ( j == 5 + (k - 1)*26 ){
				ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z() + 2*k*r_ionp, 0., r_ionp));
			}
			else if ( j == 6 + (k - 1)*26 ){
				ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z() - 2*k*r_ionp, 0., r_ionp));
			}
			else if ( j == 7 + (k - 1)*26 ){
				ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x() + 2*k*r_ionp, aggregateInstances[i].x0().y() + 2*k*r_ionp, aggregateInstances[i].x0().z(), 0., r_ionp));
			}
			else if ( j == 8 + (k - 1)*26 ){
				ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x() + 2*k*r_ionp, aggregateInstances[i].x0().y() - 2*k*r_ionp, aggregateInstances[i].x0().z(), 0., r_ionp));
			}
			else if ( j == 9 + (k - 1)*26 ){
				ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x() - 2*k*r_ionp, aggregateInstances[i].x0().y() + 2*k*r_ionp, aggregateInstances[i].x0().z(), 0., r_ionp));
			}
			else if ( j == 10 + (k - 1)*26 ){
				ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x() - 2*k*r_ionp, aggregateInstances[i].x0().y() - 2*k*r_ionp, aggregateInstances[i].x0().z(), 0., r_ionp));
			}
			else if ( j == 11 + (k - 1)*26 ){
				ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x() + 2*k*r_ionp, aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z() + 2*k*r_ionp, 0., r_ionp));
			}
			else if ( j == 12 + (k - 1)*26 ){
				ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x() + 2*k*r_ionp, aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z() - 2*k*r_ionp, 0., r_ionp));
			}
			else if ( j == 13 + (k - 1)*26 ){
				ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x() - 2*k*r_ionp, aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z() + 2*k*r_ionp, 0., r_ionp));
			}
			else if ( j == 14 + (k - 1)*26 ){
				ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x() - 2*k*r_ionp, aggregateInstances[i].x0().y(), aggregateInstances[i].x0().z() - 2*k*r_ionp, 0., r_ionp));
			}
			else if ( j == 15 + (k - 1)*26 ){
				ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y() + 2*k*r_ionp, aggregateInstances[i].x0().z() + 2*k*r_ionp, 0., r_ionp));
			}
			else if ( j == 16 + (k - 1)*26 ){
				ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y() + 2*k*r_ionp, aggregateInstances[i].x0().z() - 2*k*r_ionp, 0., r_ionp));
			}
			else if ( j == 17 + (k - 1)*26 ){
				ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y() - 2*k*r_ionp, aggregateInstances[i].x0().z() + 2*k*r_ionp, 0., r_ionp));
			}
			else if ( j == 18 + (k - 1)*26 ){
				ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x(), aggregateInstances[i].x0().y() - 2*k*r_ionp, aggregateInstances[i].x0().z() - 2*k*r_ionp, 0., r_ionp));
			}
			else if ( j == 19 + (k - 1)*26 ){
				ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x() + 2*k*r_ionp, aggregateInstances[i].x0().y() + 2*k*r_ionp, aggregateInstances[i].x0().z() + 2*k*r_ionp, 0., r_ionp));
			}
			else if ( j == 20 + (k - 1)*26 ){
				ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x() + 2*k*r_ionp, aggregateInstances[i].x0().y() + 2*k*r_ionp, aggregateInstances[i].x0().z() - 2*k*r_ionp, 0., r_ionp));
			}
			else if ( j == 21 + (k - 1)*26 ){
				ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x() + 2*k*r_ionp, aggregateInstances[i].x0().y() - 2*k*r_ionp, aggregateInstances[i].x0().z() + 2*k*r_ionp, 0., r_ionp));
			}
			else if ( j == 22 + (k - 1)*26 ){
				ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x() + 2*k*r_ionp, aggregateInstances[i].x0().y() - 2*k*r_ionp, aggregateInstances[i].x0().z() - 2*k*r_ionp, 0., r_ionp));
			}
			else if ( j == 23 + (k - 1)*26 ){
				ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x() - 2*k*r_ionp, aggregateInstances[i].x0().y() + 2*k*r_ionp, aggregateInstances[i].x0().z() + 2*k*r_ionp, 0., r_ionp));
			}
			else if ( j == 24 + (k - 1)*26 ){
				ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x() - 2*k*r_ionp, aggregateInstances[i].x0().y() + 2*k*r_ionp, aggregateInstances[i].x0().z() - 2*k*r_ionp, 0., r_ionp));
			}
			else if ( j == 25 + (k - 1)*26 ){
				ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x() - 2*k*r_ionp, aggregateInstances[i].x0().y() - 2*k*r_ionp, aggregateInstances[i].x0().z() + 2*k*r_ionp, 0., r_ionp));
			}
			else if ( j == 26 + (k - 1)*26 ){
				ionpInstances.push_back(Ionp(aggregateInstances[i].x0().x() - 2*k*r_ionp, aggregateInstances[i].x0().y() - 2*k*r_ionp, aggregateInstances[i].x0().z() - 2*k*r_ionp, 0., r_ionp));
			}

			if (j == k*26){
				k += 1;
			}

			for (unsigned int v = 0; v < SC_UI(ionpInstances.size()) - 1; ++v) { //not inside IONP // don't check final entry(is this ionp)
				if (j == 0 && i == 0) {
					break;//skip this validation for the first ionp
				}
				if (distance(ionpInstances[i * SC_UI(ionp_n) + j].xt(), ionpInstances[v].xt()) - (2. * r_ionp) < -1e-10) {
					std::cerr << "Some IONP overlap!!! " << i * SC_UI(ionp_n) + j << " " << v << std::endl;
					exit(1);
				}
			}

		}

	}

	std::cout << "Dense spherical Aggregates created: " << aggregateInstances.size() << std::endl;
	std::cout << "IONPs created: " << ionpInstances.size() << std::endl;
}
