/**
 * @file main.cpp
 *
 * Executable entry point: resolve paths and CLI options, select configuration
 * files, load/build geometry, choose the time step, run MRsequence and write
 * outputs. Each configuration is handled independently; inspect completion
 * diagnostics because a caught simulation failure does not make main return failure.
 */

/*****************************************************************************
*16/05/2022 Lauritz Klünder
******************************************************************************/

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "mptmacros.h"

#include "Aggregate.hpp"
#include "Data.hpp"
#include "Ionp.hpp"
#include "MRsequence.hpp"
#include "Point3D.hpp"
#include "Proton.hpp"
#include "RandomList.hpp"
#include "RandNumberBetween.hpp"
#include "SimMatrix.hpp"
#include "SimSpace.hpp"
#include "OutputUtils.hpp"

namespace fs = std::filesystem;

struct SimulationPaths {
	fs::path ini_root_dir;
	fs::path ini_dir;
	fs::path tmp_output_dir;
	fs::path log_output_dir;
};

struct SimulationConfig {
	SimulationPaths paths;
	bool redirect_logs = true;
	std::optional<fs::path> single_config;
	std::optional<std::uint32_t> seed;
	bool example_mode = false;
	bool record_positions = false;
	bool save_region_signals = false;
	bool use_legacy_bfield_grid = false;
	// Preserve the historical placement by default; the CLI flag remains compatible.
	bool legacy_placement = true;
	bool legacy_ionp_check = false;
	std::optional<fs::path> save_bfield_grid;
};

double automatic_base_time_step(double radius_ionp,
								double realized_outer_radius,
								double max_diffusion_coefficient) {
	const double diffusion_length =
		(realized_outer_radius > 0. ? realized_outer_radius : radius_ionp) * 2.;
	if (diffusion_length <= 0.) {
		throw std::runtime_error(
			"Cannot derive automatic time_step from a non-positive realized IONP radius");
	}
	if (max_diffusion_coefficient <= 0.) {
		throw std::runtime_error("Cannot derive automatic time_step from non-positive diffusion coefficient");
	}
	return square(diffusion_length) / (6 * max_diffusion_coefficient);
}

// The ini marker is also needed in a synthetic-only public source export.
bool is_repository_root(const fs::path& candidate) {
	return fs::is_directory(candidate / "ini") &&
		fs::is_directory(candidate / "src" / "ionp");
}

std::optional<fs::path> find_repository_root_from(fs::path start) {
	std::error_code ec;
	start = fs::weakly_canonical(start, ec);
	if (ec) {
		start = fs::absolute(start);
	}

	if (fs::is_regular_file(start)) {
		start = start.parent_path();
	}

	for (fs::path current = start; !current.empty(); current = current.parent_path()) {
		if (is_repository_root(current)) {
			return current;
		}
		if (current == current.root_path()) {
			break;
		}
	}
	return std::nullopt;
}

fs::path find_repository_root(char** argv) {
	if (auto root = find_repository_root_from(fs::current_path())) {
		return *root;
	}
	if (argv != nullptr && argv[0] != nullptr) {
		if (auto root = find_repository_root_from(fs::absolute(argv[0]))) {
			return *root;
		}
	}
	return fs::current_path();
}

SimulationPaths default_paths(const fs::path& repo_root) {
	return {
		repo_root / "ini",
		repo_root / "ini",
		repo_root / "data" / "tmp",
		repo_root / "data" / "log"
	};
}

bool contains_auxiliary_ini_subdir(const fs::path& path) {
	static const char* subdirs[] = {
		"Cell_pos",
		"Cell_load",
		"Histogram",
		"Radius_Histogram",
		"IONP_pos",
		"Agg_pos"
	};

	for (const char* subdir : subdirs) {
		if (fs::is_directory(path / subdir)) {
			return true;
		}
	}
	return false;
}

fs::path infer_ini_root_dir(const fs::path& ini_dir) {
	if (ini_dir.filename() == "ini" || contains_auxiliary_ini_subdir(ini_dir)) {
		return ini_dir;
	}

	const fs::path parent = ini_dir.parent_path();
	if (!parent.empty() && (parent.filename() == "ini" || contains_auxiliary_ini_subdir(parent))) {
		return parent;
	}

	return ini_dir;
}

bool is_supported_geometry_file(const fs::path& path) {
	return path.extension() == ".txt" || path.extension() == ".dat";
}

// Resolve a configured geometry selector against its auxiliary-input directory
// and reject unsupported or missing selections before loading simulation state.
fs::path resolve_selected_geometry_file(const fs::path& base_dir,
										const std::string& filename,
										const char* selector_name,
										const char* geometry_name) {
	if (filename.empty()) {
		throw std::runtime_error(std::string(selector_name) + " must not be empty");
	}

	const fs::path configured_path(filename);
	if (configured_path.is_absolute() || configured_path.has_root_path() ||
		configured_path != configured_path.filename() ||
		filename.find('/') != std::string::npos ||
		filename.find('\\') != std::string::npos) {
		throw std::runtime_error(std::string(selector_name) +
			" must be a leaf filename without directory components: " + filename);
	}
	if (!is_supported_geometry_file(configured_path)) {
		throw std::runtime_error(std::string(selector_name) +
			" must name a .txt or .dat file: " + filename);
	}

	std::error_code ec;
	const fs::path absolute_base = fs::absolute(base_dir).lexically_normal();
	const fs::path selected_path = (absolute_base / configured_path).lexically_normal();
	if (!fs::is_regular_file(selected_path, ec) || ec) {
		throw std::runtime_error(std::string(selector_name) +
			" does not resolve to a regular file: " + selected_path.string());
	}

	const fs::path canonical_base = fs::canonical(absolute_base, ec);
	if (ec) {
		throw std::runtime_error("Unable to resolve geometry directory " +
			absolute_base.string() + ": " + ec.message());
	}
	const fs::path canonical_selected = fs::canonical(selected_path, ec);
	if (ec || canonical_selected.parent_path() != canonical_base) {
		throw std::runtime_error(std::string(selector_name) +
			" must resolve to a file directly inside " + absolute_base.string());
	}

	std::cout << "Using explicit " << geometry_name << " file: "
		<< selected_path << std::endl;
	return selected_path;
}

std::vector<fs::path> collect_legacy_geometry_files(const fs::path& base_dir,
											 const char* geometry_name) {
	if (!fs::is_directory(base_dir)) {
		throw std::runtime_error("Required input directory does not exist: " +
			base_dir.string());
	}

	std::vector<fs::path> files;
	for (const auto& entry : fs::directory_iterator(base_dir)) {
		if (entry.is_regular_file() && is_supported_geometry_file(entry.path())) {
			files.push_back(fs::absolute(entry.path()).lexically_normal());
		}
	}
	std::sort(files.begin(), files.end());

	std::cout << "Using legacy " << geometry_name << " directory scan: "
		<< fs::absolute(base_dir).lexically_normal() << std::endl;
	for (const fs::path& path : files) {
		std::cout << "Legacy " << geometry_name << " input file: " << path << std::endl;
	}
	return files;
}

void warn_ignored_geometry_selector(const std::optional<std::string>& selector,
									const char* selector_name,
									const char* reason) {
	if (selector) {
		std::cout << "Warning: ignoring " << selector_name << "='" << *selector
			<< "' because " << reason << std::endl;
	}
}

class ScopedLogRedirect {
public:
	ScopedLogRedirect() = default;
	explicit ScopedLogRedirect(const fs::path& log_path) {
		open(log_path);
	}

	void open(const fs::path& log_path) {
		if (log_path.empty()) {
			return;
		}
		fs::create_directories(log_path.parent_path());
		log_stream.open(log_path);
		if (!log_stream.is_open()) {
			throw std::runtime_error("Unable to open log file: " + log_path.string());
		}
		old_buf = std::cout.rdbuf(log_stream.rdbuf());
	}

	~ScopedLogRedirect() {
		if (old_buf) {
			std::cout.rdbuf(old_buf);
		}
	}

private:
	std::ofstream log_stream;
	std::streambuf* old_buf = nullptr;
};

// Options are processed in argument order. Unknown options and missing values
// retain the historical permissive handling; --example also resets path defaults.
SimulationConfig parse_cli(int argc, char** argv, const fs::path& repo_root) {
	SimulationConfig config;
	config.paths = default_paths(repo_root);
	bool ini_dir_was_set = false;
	bool ini_root_was_set = false;
	for (int i = 1; i < argc; ++i) {
		std::string arg = argv[i];
		if ((arg == "--ini-dir" || arg == "-i") && i + 1 < argc) {
			config.paths.ini_dir = argv[++i];
			ini_dir_was_set = true;
		} else if (arg == "--ini-root" && i + 1 < argc) {
			config.paths.ini_root_dir = argv[++i];
			ini_root_was_set = true;
		} else if (arg == "--output-dir" && i + 1 < argc) {
			config.paths.tmp_output_dir = argv[++i];
		} else if (arg == "--log-dir" && i + 1 < argc) {
			config.paths.log_output_dir = argv[++i];
		} else if (arg == "--config" && i + 1 < argc) {
			config.single_config = fs::path(argv[++i]);
		} else if (arg == "--no-log") {
			config.redirect_logs = false;
		} else if (arg == "--seed" && i + 1 < argc) {
			config.seed = static_cast<std::uint32_t>(std::stoul(argv[++i]));
		} else if (arg == "--example") {
			config.example_mode = true;
			config.paths.ini_root_dir = repo_root / "Testing" / "example" / "ini";
			config.paths.ini_dir = config.paths.ini_root_dir;
			config.paths.tmp_output_dir = repo_root / "Testing" / "example" / "data" / "tmp";
			config.paths.log_output_dir = repo_root / "Testing" / "example" / "data" / "log";
			fs::create_directories(config.paths.ini_dir);
			fs::create_directories(config.paths.tmp_output_dir);
			fs::create_directories(config.paths.log_output_dir);
		} else if (arg == "--record-positions") {
			config.record_positions = true;
		} else if (arg == "--save-region-signals") {
			config.save_region_signals = true;
		} else if (arg == "--legacy-bgrid") {
			config.use_legacy_bfield_grid = true;
		} else if (arg == "--legacy-placement") {
			config.legacy_placement = true;
		} else if (arg == "--legacy-ionp-check") {
			config.legacy_ionp_check = true;
		} else if (arg == "--save-bgrid" && i + 1 < argc) {
			config.save_bfield_grid = fs::path(argv[++i]);
		}
	}
	if (ini_dir_was_set && !ini_root_was_set && !config.example_mode) {
		config.paths.ini_root_dir = infer_ini_root_dir(config.paths.ini_dir);
	}
	return config;
}

fs::path write_example_config(const SimulationConfig& config) {
	const fs::path example_path = config.paths.ini_dir / "example_config.txt";
	std::ofstream out(example_path);
	out << "IONP_radius=1\n";
	out << "Matrix_xlength=10\n";
	out << "Matrix_ylength=10\n";
	out << "Matrix_zlength=10\n";
	out << "Matrix_x0=0\n";
	out << "Matrix_y0=0\n";
	out << "Matrix_z0=0\n";
	out << "Matrix_xnodes=5\n";
	out << "Matrix_ynodes=5\n";
	out << "Matrix_znodes=5\n";
	out << "IONP_number=10\n";
	out << "Agg_number=0\n";
	out << "Prot_number=50\n";
	out << "echo_time=0.001\n";
	out << "time_step=0.0005\n";
	out << "variable_Time=No\n";
	out << "Surrounding_voxels=No\n";
	out << "Surrounding_voxel_superparticles=0\n";
	out << "echo_spacing=0.001\n";
	out << "Threads=1\n";
	out << "Sequence=MSE\n";
	out << "Cell_hist=No\n";
	out << "Coating=0\n";
	out << "Bound_Prot_n=0\n";
	out << "T2out=0\n";
	out << "T2in=0\n";
	out << "IONP_r_std=0\n";
	out << "Agg_r_std=0\n";
	out << "Concentration=0\n";
	out << "Concentration_extern=0\n";
	out << "Magnetic Moment=0\n";
	out << "Cell Origin= 0 0 0\n";
	out << "Voxel position= 0 0 0\n";
	out << "Cube Origin= 0 0 0\n";
	out << "Cube Size= 0 0 0\n";
	out << "Cell Radius=0\n";
	out << "Cell_permeability=0\n";
	out << "Cell Concentration=0\n";
	out << "Cell Concentration Std=0\n";
	out << "Cell Concentration Distribution=Normal\n";
	out << "Cell Number=0\n";
	return example_path;
}

// Batch configurations are sorted for a stable run order. The run-wide seed
// counter is shared across the batch, so earlier inputs affect later streams.
std::vector<fs::path> collect_config_files(const SimulationConfig& config) {
	std::vector<fs::path> files;
	if (config.single_config) {
		files.push_back(*config.single_config);
		return files;
	}

	if (!fs::exists(config.paths.ini_dir)) {
		throw std::runtime_error("INI directory does not exist: " + config.paths.ini_dir.string());
	}

	for (const auto& entry : fs::directory_iterator(config.paths.ini_dir)) {
		if (entry.is_regular_file() && entry.path().extension() == ".txt") {
			files.push_back(entry.path());
		}
	}
	std::sort(files.begin(), files.end());
	return files;
}

// One configuration owns its Data object, geometry, sequence and output vectors.
// Resolve input paths and material parameters before constructing MRsequence.
int run_single_simulation(const fs::path& config_path, const SimulationConfig& config) {
	const std::size_t next_tmp_index = output::next_file_index(config.paths.tmp_output_dir);
	const std::size_t next_log_index = output::next_file_index(config.paths.log_output_dir);

	ScopedLogRedirect log_redirect;
	if (config.redirect_logs) {
		const std::string log_name = std::to_string(next_log_index) + "_" + std::to_string(next_tmp_index) + "_log.txt";
		log_redirect.open(config.paths.log_output_dir / log_name);
	}

	auto start = std::chrono::system_clock::now();
	std::time_t start_time = std::chrono::system_clock::to_time_t(start);
	std::cout << "Started computation at: " << std::ctime(&start_time) << "\n";
	const IonpCheckMode ionp_check_mode = config.legacy_ionp_check ? IonpCheckMode::All : IonpCheckMode::Spatial;
	std::cout << "IONP check mode: " << (ionp_check_mode == IonpCheckMode::All ? "all IONPs" : "spatial index") << std::endl;

	std::cout << config_path << std::endl;

	Data dat(config_path.string());
	dat.set_output_context(config.paths.tmp_output_dir, config.paths.log_output_dir);

	dat.readData();

	double matrix_x0 = dat.get_matrix_x0() ;
	double matrix_y0 = dat.get_matrix_y0() ;
	double matrix_z0 = dat.get_matrix_z0() ;
	double matrix_xlength = dat.get_matrix_xlength() ;
	double matrix_ylength = dat.get_matrix_ylength() ;
	double matrix_zlength = dat.get_matrix_zlength() ;
	int matrix_nxnodes = dat.get_matrix_nxnodes() ;
	int matrix_nynodes = dat.get_matrix_nynodes() ;
	int matrix_nznodes = dat.get_matrix_nznodes() ;
	double radius_ionp = dat.get_radius_ionp() ;
	double radius_ionp_1 = dat.get_radius_ionp_1() ;
	double ionp_per = dat.get_ionp_per() ;
	double radius_agg = dat.get_radius_agg() ;
	double radius_ionp_std = dat.get_radius_ionp_std() ;
	double radius_ionp_std_1 = dat.get_radius_ionp_std_1() ;
	double radius_agg_std = dat.get_radius_agg_std() ;
	int ionp_number = dat.get_ionp_number() ;
	int agg_number = dat.get_agg_number() ;
	int proton_number = dat.get_proton_number();
	double echo_time = dat.get_echo_time() ;
	double time_step = dat.get_time_step() ;
	bool var_TimeStep = dat.get_variable_time();
	const bool surrounding_voxels = dat.get_surrounding_voxels();
	const std::size_t surrounding_voxel_superparticles =
		dat.get_surrounding_voxel_superparticles();
	double echo_spacing = dat.get_echo_spacing() ;
	int numThreads = dat.get_numThreads() ;
	int sequence_num = dat.get_sequence_num() ;
	double velocityIONP = dat.get_velocity_IONP() ;
	Point3D directionIONP = dat.get_direction_IONP() ;
	double thickness_ionp = dat.get_thickness_IONPcoating() ;
	double coating_diffusion_coefficient = dat.has_coating_diffusion_coefficient() ? dat.get_coating_diffusion_coefficient() : D_D_CONST;
	double cell_diffusion_coefficient = dat.has_cell_diffusion_coefficient() ? dat.get_cell_diffusion_coefficient() : D_D_CONST;
	const double max_diffusion_coefficient = std::max({D_D_CONST, coating_diffusion_coefficient, cell_diffusion_coefficient});
	bool check_coating_permeability = dat.has_coating_permeability() && dat.get_coating_permeability() < 1.;
	double coating_permeability = dat.get_coating_permeability();
	int bound_proton_number = dat.get_bound_proton_number() ;
	double background_T2 = dat.get_background_T2() ;
	double cell_T2 = dat.get_cell_T2() ;
	double concentration = dat.get_concentration() ;
	double concentration_ext = dat.get_concentration_ext() ;
	double magnetic_moment = dat.get_mag_moment() ;
	Point3D Cell_orig = dat.get_Cell_Origin() ;
	double Cell_rad = dat.get_Cell_radius() ;
	double Cell_permeability = dat.get_cell_permeability() ;
	int Cell_num = dat.get_Cell_number() ;
	double Cell_concentration = dat.get_Cell_concentration() ;
	double Cell_concentration_std = dat.get_Cell_concentration_std() ;
	Point3D Cube_orig = dat.get_Cube_Origin() ;
	Point3D voxel_pos = dat.get_Voxel_position() ;
	Point3D Cube_size = dat.get_Cube_size() ;
	bool hist_bool = dat.get_Cell_hist() ;
	bool radius_hist_bool = dat.get_Radius_hist();
	int Cell_concentration_dist = dat.get_Cell_concentration_dist();
	ImportanceSamplingConfig importance_sampling = dat.get_importance_sampling();
	const BFieldSolver bfield_solver = dat.get_bfield_solver();
	const std::optional<std::string> ionp_position_selector = dat.get_ionp_position_file();
	const std::optional<std::string> agg_position_selector = dat.get_agg_position_file();
	const std::optional<std::string> cell_position_selector = dat.get_cell_position_file();
	const std::optional<std::string> cell_load_selector = dat.get_cell_load_file();

	const bool uses_explicit_ionp_geometry = ionp_number == 0;
	const bool uses_explicit_aggregate_geometry =
		ionp_number > 0 && agg_number == 0 && radius_agg > 0.;
	const bool uses_nonaggregate_cell_geometry =
		ionp_number > 0 && agg_number == 0 && radius_agg == 0. &&
		Cell_num == 0 && Cell_rad > 0.;
	const bool uses_aggregate_cell_geometry =
		ionp_number > 0 && radius_agg > 0. && Cell_num == 0 && Cell_rad > 0.;
	const bool uses_cell_position_geometry =
		uses_nonaggregate_cell_geometry || uses_aggregate_cell_geometry;
	const bool uses_cell_load_geometry =
		uses_nonaggregate_cell_geometry;

	if (!uses_explicit_ionp_geometry) {
		warn_ignored_geometry_selector(
			ionp_position_selector, "IONP_position_file",
			"the selected placement branch does not load manual IONP positions");
	}
	if (!uses_explicit_aggregate_geometry) {
		warn_ignored_geometry_selector(
			agg_position_selector, "Agg_position_file",
			"the selected placement branch does not load manual aggregate positions");
	}
	if (!uses_cell_position_geometry) {
		warn_ignored_geometry_selector(
			cell_position_selector, "Cell_position_file",
			"the selected placement branch does not load cell positions from a file");
	}
	if (!uses_cell_load_geometry) {
		warn_ignored_geometry_selector(
			cell_load_selector, "Cell_load_file",
			uses_aggregate_cell_geometry
				? "aggregate placement uses Histogram data for cell loading"
				: "the selected placement branch does not load per-cell iron loads");
	}

	if (uses_nonaggregate_cell_geometry && hist_bool &&
		cell_position_selector.has_value() != cell_load_selector.has_value()) {
		throw std::runtime_error(
			"Cell_hist=Yes with file-based cells requires Cell_position_file and "
			"Cell_load_file to be specified together, or both omitted");
	}

	std::optional<fs::path> selected_ionp_position_path;
	std::optional<fs::path> selected_agg_position_path;
	std::optional<fs::path> selected_cell_position_path;
	std::optional<fs::path> selected_cell_load_path;
	if (uses_explicit_ionp_geometry && ionp_position_selector) {
		selected_ionp_position_path = resolve_selected_geometry_file(
			config.paths.ini_dir / "IONP_pos", *ionp_position_selector,
			"IONP_position_file", "IONP position");
	}
	if (uses_explicit_aggregate_geometry && agg_position_selector) {
		selected_agg_position_path = resolve_selected_geometry_file(
			config.paths.ini_dir / "Agg_pos", *agg_position_selector,
			"Agg_position_file", "aggregate position");
	}
	if (uses_cell_position_geometry && cell_position_selector) {
		selected_cell_position_path = resolve_selected_geometry_file(
			config.paths.ini_root_dir / "Cell_pos", *cell_position_selector,
			"Cell_position_file", "cell position");
	}
	if (uses_cell_load_geometry && cell_load_selector) {
		selected_cell_load_path = resolve_selected_geometry_file(
			config.paths.ini_root_dir / "Cell_load", *cell_load_selector,
			"Cell_load_file", "cell load");
	}

	if (importance_sampling.enabled) {
		proton_number = importance_sampling.total_free_protons();
		dat.set_effective_proton_number(proton_number);
		std::cout << "Importance-sampled free proton count: " << proton_number << std::endl;
		if (config.save_region_signals) {
			std::cout << "Importance-sampling diagnostic signal output enabled" << std::endl;
		}
	} else {
		dat.set_effective_proton_number(proton_number);
		if (config.save_region_signals) {
			std::cout << "--save-region-signals has no effect because importance sampling is disabled" << std::endl;
		}
	}

	(void) radius_agg;
	(void) agg_number;
	(void) echo_spacing;

	SimMatrix sMatrix = SimMatrix(
	    matrix_x0,
	    matrix_y0,
	    matrix_z0,
	    matrix_xlength,
	    matrix_ylength,
	    matrix_zlength,
	    matrix_nxnodes,
	    matrix_nynodes,
	    matrix_nznodes);
	if (config.use_legacy_bfield_grid) {
		std::cout << "--legacy-bgrid is retained for compatibility; the full all-source magnetic grid is already the default" << std::endl;
	}
	sMatrix.set_save_grid_path(config.save_bfield_grid);

	std::vector<Ionp> ionps;
	ionps.reserve(SC_UI(ionp_number));

	std::vector<Point3D> cell_origins;
	cell_origins.reserve(SC_UI(Cell_num));
	if (Cell_num == 1) {
		cell_origins.push_back(Cell_orig);
	}

	std::vector<Ionp> ionpMove;
	ionpMove.reserve(SC_UI(ionp_number));

	std::vector<Aggregate> aggregates;
	aggregates.reserve(SC_UI(agg_number));

	if (ionp_number == 0){
		if (matrix_xlength == 0.){
			throw std::runtime_error("Size of Simulation volume has to be given, if IONPs are placed by hand!");
		}
		//This does not work for dense aggregates yet ...

		const fs::path path_to_ionp = config.paths.ini_dir / "IONP_pos";
		if (selected_ionp_position_path) {
			const std::size_t previous_size = ionps.size();
			dat.readIONP(selected_ionp_position_path->string(), ionps, sMatrix,
				radius_ionp, thickness_ionp);
			if (ionps.size() == previous_size) {
				throw std::runtime_error("Selected IONP position file produced no records: " +
					selected_ionp_position_path->string());
			}
		}
		else {
			const std::vector<fs::path> position_files =
				collect_legacy_geometry_files(path_to_ionp, "IONP position");
			for (const fs::path& pathIONP : position_files) {
				dat.readIONP(pathIONP.string(), ionps, sMatrix, radius_ionp, thickness_ionp);
			}
			if (position_files.empty()) {
				throw std::runtime_error("No IONP number given and no IONP position file found in " + path_to_ionp.string());
			}
		}

		if (ionps.empty()) {
			throw std::runtime_error("No IONP number given and no IONP position file found in " + path_to_ionp.string());
		}

	}
	else{
		if (agg_number > 0 && radius_agg > 0. && concentration_ext > 0.) {
			// Preserve the normal intracellular population, but group only the
			// extracellular iron into the requested aggregates.
			if (Cell_rad != 0. || Cell_num != 0) {
				SimSpace::simulationSpace_inCell(ionps, sMatrix, ionp_number, concentration, radius_ionp, radius_ionp_1, ionp_per, radius_ionp_std, radius_ionp_std_1, thickness_ionp, Cell_orig, Cell_rad, Cell_num, Cell_concentration, Cell_concentration_std, Cell_concentration_dist, cell_origins, voxel_pos, hist_bool, config.paths.ini_root_dir, !config.legacy_placement, selected_cell_position_path, selected_cell_load_path);
			}
			SimSpace::simulationSpace_outCellAggregates(
				aggregates, ionps, sMatrix, agg_number, concentration_ext,
				radius_agg, radius_agg_std, radius_ionp, radius_ionp_1,
				ionp_per, radius_ionp_std, radius_ionp_std_1,
				thickness_ionp, Cell_orig, Cell_rad, cell_origins);
		}
		else if ( agg_number == 0 && radius_agg == 0){
			if ( Cell_rad != 0. || Cell_num != 0 ){
				SimSpace::simulationSpace_inCell(ionps, sMatrix, ionp_number, concentration, radius_ionp, radius_ionp_1, ionp_per, radius_ionp_std, radius_ionp_std_1, thickness_ionp, Cell_orig, Cell_rad, Cell_num, Cell_concentration, Cell_concentration_std, Cell_concentration_dist, cell_origins, voxel_pos, hist_bool, config.paths.ini_root_dir, !config.legacy_placement, selected_cell_position_path, selected_cell_load_path);
				if ( concentration_ext != 0 ){
					SimSpace::simulationSpace_outCell(ionps, sMatrix, ionp_number, concentration_ext, radius_ionp, radius_ionp_1, ionp_per, radius_ionp_std, radius_ionp_std_1, thickness_ionp, Cell_orig, Cell_rad, Cell_num, cell_origins, voxel_pos, radius_hist_bool, config.paths.ini_root_dir, !config.legacy_placement);
				}
			}
			else{
				if ( concentration_ext != 0. ){
					concentration = concentration_ext;
				}
				SimSpace::simulationSpace_outCell(ionps, sMatrix, ionp_number, concentration, radius_ionp, radius_ionp_1, ionp_per, radius_ionp_std, radius_ionp_std_1, thickness_ionp, Cell_orig, Cell_rad, Cell_num, cell_origins, voxel_pos, radius_hist_bool, config.paths.ini_root_dir, !config.legacy_placement);
			}
		}
		else if ( agg_number != 0 && radius_agg == 0 ){
			SimSpace::simulationSpace(aggregates, ionps, sMatrix, agg_number, ionp_number, radius_ionp);
		}
		else if ( agg_number == 0 && radius_agg != 0 ){
			const fs::path path_to_agg = config.paths.ini_dir / "Agg_pos";

			if (selected_agg_position_path) {
				const std::size_t previous_size = aggregates.size();
				dat.readAgg(selected_agg_position_path->string(), aggregates, sMatrix, radius_agg);
				if (aggregates.size() == previous_size) {
					throw std::runtime_error("Selected aggregate position file produced no records: " +
						selected_agg_position_path->string());
				}
			}
			else {
				const std::vector<fs::path> aggregate_files =
					collect_legacy_geometry_files(path_to_agg, "aggregate position");
				for (const fs::path& pathAgg : aggregate_files) {
					dat.readAgg(pathAgg.string(), aggregates, sMatrix, radius_agg);
				}
				if (aggregate_files.empty()) {
					throw std::runtime_error("No Aggregates number given and no Aggregates position file found in " + path_to_agg.string());
				}
			}

			if (aggregates.empty()) {
				throw std::runtime_error("No Aggregates number given and no Aggregates position file found in " + path_to_agg.string());
			}

			SimSpace::simulationSpace(aggregates, ionps, sMatrix, agg_number, ionp_number, concentration, radius_agg, radius_ionp, radius_ionp_1, ionp_per, radius_ionp_std, radius_ionp_std_1, radius_agg_std, thickness_ionp, Cell_orig, Cell_rad, Cell_num, cell_origins, hist_bool, config.paths.ini_root_dir, selected_cell_position_path);
		}
		else{
			SimSpace::simulationSpace(aggregates, ionps, sMatrix, agg_number, ionp_number, concentration, radius_agg, radius_ionp, radius_ionp_1, ionp_per, radius_ionp_std, radius_ionp_std_1, radius_agg_std, thickness_ionp, Cell_orig, Cell_rad, Cell_num, cell_origins, hist_bool, config.paths.ini_root_dir, selected_cell_position_path);
		}
	}

	const IonpPopulationMap ionp_populations =
		classify_ionp_populations(
			ionps, cell_origins, Cell_orig, Cell_rad);

	double minimum_outer_radius = std::numeric_limits<double>::infinity();
	double maximum_outer_radius = 0.;
	double minimum_core_radius = std::numeric_limits<double>::infinity();
	double maximum_core_radius = 0.;
	double minimum_shell_thickness = std::numeric_limits<double>::infinity();
	double maximum_shell_thickness = 0.;
	std::size_t zero_shell_particles = 0;
	for (const Ionp& ionp : ionps) {
		minimum_outer_radius =
			std::min(minimum_outer_radius, ionp.outer_radius());
		maximum_outer_radius =
			std::max(maximum_outer_radius, ionp.outer_radius());
		minimum_core_radius =
			std::min(minimum_core_radius, ionp.radius());
		maximum_core_radius =
			std::max(maximum_core_radius, ionp.radius());
		minimum_shell_thickness =
			std::min(minimum_shell_thickness, ionp.coating());
		maximum_shell_thickness =
			std::max(maximum_shell_thickness, ionp.coating());
		if (ionp.coating() <= 0.) {
			++zero_shell_particles;
		}
	}
	if (!ionps.empty()) {
		std::cout << "Realized IONP core-radius range: "
			<< minimum_core_radius << " - " << maximum_core_radius
			<< std::endl;
		std::cout << "Stored IONP shell-thickness range: "
			<< minimum_shell_thickness << " - "
			<< maximum_shell_thickness << std::endl;
		std::cout << "Realized IONP outer-radius range: "
			<< minimum_outer_radius << " - " << maximum_outer_radius
			<< std::endl;
		std::cout << "IONPs without a positive shell: "
			<< zero_shell_particles << std::endl;
	}
	std::cout << "IONP populations (intracellular / extracellular): "
		<< ionp_populations.intracellular_indices.size() << " / "
		<< ionp_populations.extracellular_indices.size() << std::endl;
	if (bound_proton_number > 0 &&
		ionp_populations.coated_intracellular_count +
			ionp_populations.coated_extracellular_count == 0) {
		throw std::runtime_error(
			"Bound_Prot_n was requested but no IONP has a positive shell");
	}

	if (importance_sampling.enabled) {
		const bool overlapping_proposals =
			importance_sampling.overlapping_proposals_enabled();
		if (importance_sampling.split_bound_configured) {
			if (importance_sampling.bound_protons_intracellular > 0 &&
				ionp_populations.coated_intracellular_count == 0) {
				throw std::runtime_error(
					"Intracellular bound samples were requested but no intracellular IONP has a positive shell");
			}
			if (importance_sampling.bound_protons_extracellular > 0 &&
				ionp_populations.coated_extracellular_count == 0) {
				throw std::runtime_error(
					"Extracellular bound samples were requested but no extracellular IONP has a positive shell");
			}
		}
		if (overlapping_proposals) {
			const std::array<int, kFreeSamplingProposalCount> requested_counts =
				importance_sampling.proposal_protons;
			const auto proposal_pilot_start =
				std::chrono::steady_clock::now();
			const OverlappingProposalVolumeEstimate estimate =
				estimate_overlapping_proposal_volumes(
					importance_sampling,
					cell_origins,
					Cell_orig,
					Cell_rad,
					ionps,
					ionp_populations,
					sMatrix.xLength(),
					sMatrix.yLength(),
					sMatrix.zLength());
			const std::chrono::duration<double> proposal_pilot_elapsed =
				std::chrono::steady_clock::now() - proposal_pilot_start;
			if (!(estimate.accessible_volume > 0.) ||
				!(estimate.box_volume > 0.)) {
				throw std::runtime_error(
					"Overlapping importance sampling found no accessible free-water volume");
			}
			const EffectiveProposalCounts effective =
				resolve_effective_proposal_counts(
					requested_counts, estimate.region_fractions);
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
			}
			importance_sampling.proposal_volume_fractions =
				estimate.region_fractions;
			importance_sampling.free_water_volume_fraction =
				std::clamp(
					estimate.accessible_volume / estimate.box_volume,
					0.0,
					1.0);
			importance_sampling.proposal_volumes_estimated = true;

			std::cout << "Overlapping importance-sampling proposal volumes and counts:"
				<< std::endl;
			std::cout << "  volume_pilot_seconds: "
				<< proposal_pilot_elapsed.count() << std::endl;
			std::cout << "  accessible_free_water_fraction: "
				<< importance_sampling.free_water_volume_fraction << std::endl;
			for (std::size_t proposal_index = 0;
				 proposal_index < kFreeSamplingProposalCount;
				 ++proposal_index) {
				std::cout << "  "
					<< free_sampling_proposal_name(
						static_cast<FreeSamplingProposal>(proposal_index))
					<< ": volume_fraction="
					<< estimate.region_fractions[proposal_index]
					<< ", requested=" << requested_counts[proposal_index]
					<< ", effective=" << effective.counts[proposal_index];
				if (effective.redistributed[proposal_index]) {
					std::cout << " (moved to uniform: zero estimated support)";
				}
				std::cout << std::endl;
			}
			std::cout << "  transferred_to_uniform: "
				<< effective.redistributed_to_uniform << std::endl;
		}
		else {
			bool requests_intracellular_ionp_region = false;
			bool requests_extracellular_ionp_region = false;
			for (std::size_t region_index = 0;
				 region_index < kFreeSampleRegionCount;
				 ++region_index) {
				if (importance_sampling.joint_protons[region_index] <= 0) {
					continue;
				}
				const JointSampleRegion region =
					decode_joint_sample_region(
						static_cast<ProtonSampleRegion>(region_index));
				requests_intracellular_ionp_region =
					requests_intracellular_ionp_region ||
					region.intracellular_ionp != IonpDistanceRegion::Far;
				requests_extracellular_ionp_region =
					requests_extracellular_ionp_region ||
					region.extracellular_ionp != IonpDistanceRegion::Far;
			}
			if (requests_intracellular_ionp_region &&
				ionp_populations.intracellular_indices.empty()) {
				throw std::runtime_error(
					"Intracellular IONP importance sampling is configured but no intracellular IONPs exist");
			}
			if (requests_extracellular_ionp_region &&
				ionp_populations.extracellular_indices.empty()) {
				throw std::runtime_error(
					"Extracellular IONP importance sampling is configured but no extracellular IONPs exist");
			}

			const bool needs_automatic_weights =
				!has_explicit_free_region_weights(importance_sampling);
			std::array<double, kFreeSampleRegionCount> estimated_weights{};
			if (needs_automatic_weights ||
				importance_sampling.joint_regions_configured) {
				estimated_weights =
					estimate_active_joint_region_volume_weights(
						importance_sampling,
						cell_origins,
						Cell_orig,
						Cell_rad,
						ionps,
						ionp_populations,
						sMatrix.xLength(),
						sMatrix.yLength(),
						sMatrix.zLength());
			}
			if (needs_automatic_weights) {
				importance_sampling.joint_weights = estimated_weights;
				std::cout << "Automatic active joint importance-sampling weights:"
					<< std::endl;
				for (std::size_t region_index = 0;
					 region_index < kFreeSampleRegionCount;
					 ++region_index) {
					if (importance_sampling.joint_protons[region_index] == 0) {
						continue;
					}
					std::cout << "  "
						<< sample_region_name(
							static_cast<ProtonSampleRegion>(region_index))
						<< ": "
						<< importance_sampling.joint_weights[region_index]
						<< std::endl;
				}
			}
		}

		const std::size_t bound_samples =
			importance_sampling.split_bound_configured
				? SC_UI(importance_sampling.bound_protons_intracellular) *
						ionp_populations.coated_intracellular_count +
					SC_UI(importance_sampling.bound_protons_extracellular) *
						ionp_populations.coated_extracellular_count
				: SC_UI(bound_proton_number) *
					(ionp_populations.coated_intracellular_count +
					 ionp_populations.coated_extracellular_count);
		dat.set_effective_proton_number(
			proton_number + SC_I(bound_samples));
	}

	double radius_ionp_mean = 0.;

	for (std::size_t i = 0; i < ionps.size(); ++i){
		radius_ionp_mean += ionps[i].radius();
	}
	radius_ionp_mean = radius_ionp_mean / SC_D(ionps.size());

	radius_ionp = radius_ionp_mean;

	int steps_number;
	if (time_step == 0){
		steps_number = SC_I(round(echo_time / automatic_base_time_step(radius_ionp, maximum_outer_radius, max_diffusion_coefficient))); // number of steps
		std::cout << "Steps_n: " << steps_number <<"\n" << std::endl;
	}
	else{
		steps_number = SC_I(round(echo_time / time_step)); // number of steps
		std::cout << "Steps_n: " << steps_number <<"\n" << std::endl;
	}

	std::vector<Proton> protons;
	protons.reserve(SC_UI(agg_number * ionp_number));

	std::vector<Proton> protonPos;
	protonPos.reserve(SC_UI(agg_number * ionp_number));

	std::vector<std::complex<double>> signalV;
	signalV.reserve(SC_UI(steps_number + 1));

	std::vector<std::complex<double>> signalV_MSE;
	signalV_MSE.reserve(SC_UI(steps_number + 1));

	std::vector<double> timeV;
	timeV.reserve(SC_UI(steps_number + 1));

	std::vector<double> timeV_MSE;
	timeV_MSE.reserve(SC_UI(steps_number + 1));

	std::vector<double> phaseV;
	phaseV.reserve(SC_UI(steps_number + 1));

	bool Grid, var_Steps = false;

	if (matrix_nxnodes == 0 && matrix_nynodes == 0 && matrix_nznodes == 0){
		if ( time_step == 0. || var_TimeStep) {
			std::cout << "Variable diffusion step length (R or R/8) used" << std::endl;
			if ( time_step == 0. ){
				time_step = automatic_base_time_step(radius_ionp, maximum_outer_radius, max_diffusion_coefficient);
			}
			dat.set_time_step(time_step);
			Grid = false;
			var_Steps = true;
		}
		else{
			Grid = false;
			var_Steps = false;
		}
	}
	else {
		if ( time_step == 0. || var_TimeStep ) {
			std::cout << "Variable diffusion step length (R or R/8) AND B-Field Grid with interpolation used" << std::endl;
			if ( time_step == 0. ){
				time_step = automatic_base_time_step(radius_ionp, maximum_outer_radius, max_diffusion_coefficient);
			}
			dat.set_time_step(time_step);
			Grid = true;
			var_Steps = true;
		}
		else{
			std::cout << "B-Field Grid with interpolation used" << std::endl;
			Grid = true;
			var_Steps = false;
		}
	}

	MRsequence MRI(
	    sMatrix,
	    sequence_num,
	    proton_number,
	    bound_proton_number,
	    radius_ionp,
	    echo_time,
	    time_step,
	    echo_spacing,
	    numThreads,
	    matrix_nxnodes,
	    matrix_xlength,
	    directionIONP,
	    velocityIONP,
	    thickness_ionp,
	    coating_diffusion_coefficient,
	    cell_diffusion_coefficient,
	    check_coating_permeability,
	    coating_permeability,
	    background_T2,
	    cell_T2,
	    magnetic_moment,
	    Cell_orig,
	    Cell_rad,
	    Cell_permeability,
	    Cube_orig,
	    Cube_size,
	    Grid,
	    var_Steps,
	    config.record_positions,
	    importance_sampling,
	    ionp_check_mode);
	MRI.set_bfield_solver(bfield_solver);
	MRI.set_surrounding_voxels(surrounding_voxels);
	MRI.set_surrounding_voxel_superparticles(surrounding_voxel_superparticles);
	MRI.set_output_context(config.paths.tmp_output_dir);

	ImportanceRegionSignals regionSignalVectors;
	ImportanceRegionSignals regionSignalVectorsMSE;
	ImportanceProposalSignalDiagnostics proposalSignalDiagnostics;
	ImportanceProposalSignalDiagnostics proposalSignalDiagnosticsMSE;
	const bool save_overlapping_proposal_signals =
		config.save_region_signals &&
		importance_sampling.enabled &&
		importance_sampling.overlapping_proposals_enabled();
	const bool save_legacy_region_signals =
		config.save_region_signals &&
		importance_sampling.enabled &&
		!importance_sampling.overlapping_proposals_enabled();
	MRI.Sequence(
		protons,
		ionps,
		signalV,
		timeV,
		signalV_MSE,
		timeV_MSE,
		cell_origins,
		save_legacy_region_signals
			? &regionSignalVectors
			: nullptr,
		save_legacy_region_signals
			? &regionSignalVectorsMSE
			: nullptr,
		save_overlapping_proposal_signals
			? &proposalSignalDiagnostics
			: nullptr,
		save_overlapping_proposal_signals
			? &proposalSignalDiagnosticsMSE
			: nullptr);

	//============ Outfile magnetization =================================//

	dat.writeData(signalV, timeV);
	if (config.save_region_signals && importance_sampling.enabled) {
		if (save_legacy_region_signals) {
			for (std::size_t region_index = 0;
				 region_index < kSampleRegionCount;
				 ++region_index) {
				if (regionSignalVectors[region_index].empty()) {
					continue;
				}
				dat.writeImportanceRegionData(
					regionSignalVectors[region_index],
					timeV,
					static_cast<ProtonSampleRegion>(region_index));
			}
		}
		else if (save_overlapping_proposal_signals) {
			for (std::size_t proposal_index = 0;
				 proposal_index < kFreeSamplingProposalCount;
				 ++proposal_index) {
				if (proposalSignalDiagnostics.means[proposal_index].empty()) {
					continue;
				}
				const FreeSamplingProposal proposal =
					static_cast<FreeSamplingProposal>(proposal_index);
				dat.writeImportanceProposalData(
					proposalSignalDiagnostics.means[proposal_index],
					timeV,
					proposal,
					false);
				dat.writeImportanceProposalData(
					proposalSignalDiagnostics.contributions[proposal_index],
					timeV,
					proposal,
					true);
			}
		}
	}

	if (sequence_num == 2){
		dat.writeData(signalV_MSE, timeV_MSE);
		if (config.save_region_signals && importance_sampling.enabled) {
			if (save_legacy_region_signals) {
				for (std::size_t region_index = 0;
					 region_index < kSampleRegionCount;
					 ++region_index) {
					if (regionSignalVectorsMSE[region_index].empty()) {
						continue;
					}
					dat.writeImportanceRegionData(
						regionSignalVectorsMSE[region_index],
						timeV_MSE,
						static_cast<ProtonSampleRegion>(region_index));
				}
			}
			else if (save_overlapping_proposal_signals) {
				for (std::size_t proposal_index = 0;
					 proposal_index < kFreeSamplingProposalCount;
					 ++proposal_index) {
					if (proposalSignalDiagnosticsMSE.means[proposal_index].empty()) {
						continue;
					}
					const FreeSamplingProposal proposal =
						static_cast<FreeSamplingProposal>(proposal_index);
					dat.writeImportanceProposalData(
						proposalSignalDiagnosticsMSE.means[proposal_index],
						timeV_MSE,
						proposal,
						false);
					dat.writeImportanceProposalData(
						proposalSignalDiagnosticsMSE.contributions[proposal_index],
						timeV_MSE,
						proposal,
						true);
				}
			}
		}
	}

	dat.writeData(ionps);

	if ( agg_number != 0 || radius_agg != 0 ){
		dat.writeData(aggregates);
	}

	auto end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = end - start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);

	std::cout << "finished computation at " << std::ctime(&end_time)
		<< "elapsed time: " << elapsed_seconds.count() << "s\n";

	return 0;
}

// Set a requested seed once for the entire invocation. Per-configuration
// failures are reported and the loop continues; only outer failures return 1.
int main(int argc, char** argv ) {
	try {
		const fs::path repo_root = find_repository_root(argv);
		auto config = parse_cli(argc, argv, repo_root);
		if (config.seed.has_value()) {
			rng_config::set_seed(*config.seed);
		}

		std::vector<fs::path> configs_to_run;
		if (config.example_mode) {
			configs_to_run.push_back(write_example_config(config));
		} else {
			configs_to_run = collect_config_files(config);
		}

		for (const auto& config_file : configs_to_run) {
			try {
				run_single_simulation(config_file, config);
			} catch (const std::exception& e) {
				std::cerr << "Simulation failed for " << config_file << ": " << e.what() << std::endl;
			}
		}
	} catch (const std::exception& e) {
		std::cerr << "Fatal error: " << e.what() << std::endl;
		return 1;
	}

	return 0;
}
