#include <cassert>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Aggregate.hpp"
#include "Data.hpp"
#include "Ionp.hpp"
#include "RandNumberBetween.hpp"
#include "SimMatrix.hpp"
#include "SimSpace.hpp"

namespace {

void write_common_config(std::ofstream& cfg,
						 double cell_radius,
						 int bound_proton_number,
						 bool variable_time = false) {
	cfg << "IONP_radius=1\n";
	cfg << "Matrix_xlength=10\n";
	cfg << "Matrix_ylength=10\n";
	cfg << "Matrix_zlength=10\n";
	cfg << "Matrix_x0=0\n";
	cfg << "Matrix_y0=0\n";
	cfg << "Matrix_z0=0\n";
	cfg << "Matrix_xnodes=10\n";
	cfg << "Matrix_ynodes=10\n";
	cfg << "Matrix_znodes=10\n";
	cfg << "IONP_number=1\n";
	cfg << "Agg_number=0\n";
	cfg << "Prot_number=1\n";
	cfg << "echo_time=0.01\n";
	cfg << "time_step=0.01\n";
	cfg << "variable_Time=" << (variable_time ? "Yes" : "No") << "\n";
	cfg << "echo_spacing=0.01\n";
	cfg << "Threads=1\n";
	cfg << "Sequence=FID\n";
	cfg << "Cell_hist=No\n";
	cfg << "Coating=0\n";
	cfg << "Bound_Prot_n=" << bound_proton_number << "\n";
	cfg << "T2out=0\n";
	cfg << "T2in=0\n";
	cfg << "IONP_r_std=0\n";
	cfg << "Agg_r_std=0\n";
	cfg << "Concentration=0\n";
	cfg << "Concentration_extern=0\n";
	cfg << "Magnetic Moment=0\n";
	cfg << "Cell Origin= 0 0 0\n";
	cfg << "Voxel position= 0 0 0\n";
	cfg << "Cube Origin= 0 0 0\n";
	cfg << "Cube Size= 0 0 0\n";
	cfg << "Cell Radius=" << cell_radius << "\n";
	cfg << "Cell Concentration=0\n";
	cfg << "Cell Concentration Std=0\n";
	cfg << "Cell Concentration Distribution=Normal\n";
	cfg << "Cell Number=1\n";
}

} // namespace

int main() {
	namespace fs = std::filesystem;

	const fs::path tmp_dir = fs::current_path() / "basic_parsing_tmp";
	fs::remove_all(tmp_dir);
	fs::create_directories(tmp_dir);
	const fs::path config_path = tmp_dir / "ionp_basic_config.txt";
	const fs::path is_config_path = tmp_dir / "ionp_importance_sampling_config.txt";
	const fs::path auto_weight_config_path = tmp_dir / "ionp_auto_weight_config.txt";
	const fs::path joint_is_config_path = tmp_dir / "ionp_joint_importance_sampling_config.txt";
	const fs::path overlapping_is_config_path =
		tmp_dir / "ionp_overlapping_importance_sampling_config.txt";
	const fs::path overlapping_without_uniform_config_path =
		tmp_dir / "ionp_overlapping_without_uniform_config.txt";
	const fs::path overlapping_negative_count_config_path =
		tmp_dir / "ionp_overlapping_negative_count_config.txt";
	const fs::path negative_bound_count_config_path =
		tmp_dir / "ionp_negative_bound_count_config.txt";
	const fs::path overlapping_far_without_cutoffs_config_path =
		tmp_dir / "ionp_overlapping_far_without_cutoffs_config.txt";
	const fs::path mixed_overlapping_legacy_config_path =
		tmp_dir / "ionp_mixed_overlapping_legacy_config.txt";
	const fs::path mixed_overlapping_joint_config_path =
		tmp_dir / "ionp_mixed_overlapping_joint_config.txt";
	const fs::path mixed_overlapping_weight_config_path =
		tmp_dir / "ionp_mixed_overlapping_weight_config.txt";
	const fs::path ionp_only_is_with_cell_config_path = tmp_dir / "ionp_only_importance_sampling_with_cell_config.txt";
	const fs::path near_cell_without_cutoff_config_path = tmp_dir / "ionp_near_cell_without_cutoff_config.txt";
	const fs::path invalid_is_config_path = tmp_dir / "ionp_invalid_importance_sampling_config.txt";
	const fs::path mixed_is_config_path = tmp_dir / "ionp_mixed_importance_sampling_config.txt";
	const fs::path incomplete_cutoff_config_path = tmp_dir / "ionp_incomplete_ionp_cutoff_config.txt";
	const fs::path coating_diffusion_config_path = tmp_dir / "ionp_coating_diffusion_config.txt";
	const fs::path cell_diffusion_config_path = tmp_dir / "ionp_cell_diffusion_config.txt";
	const std::vector<std::pair<std::string, fs::path>> cell_diffusion_alias_configs = {
		{"D_cell", tmp_dir / "ionp_d_cell_config.txt"},
		{"cell_diffusion", tmp_dir / "ionp_cell_diffusion_alias_config.txt"},
		{"cell_diffusion_coefficient", tmp_dir / "ionp_cell_diffusion_coefficient_config.txt"}
	};
	const fs::path invalid_cell_diffusion_config_path = tmp_dir / "ionp_invalid_cell_diffusion_config.txt";
	const fs::path coating_permeability_config_path = tmp_dir / "ionp_coating_permeability_config.txt";
	const fs::path invalid_low_coating_permeability_config_path = tmp_dir / "ionp_invalid_low_coating_permeability_config.txt";
	const fs::path invalid_high_coating_permeability_config_path = tmp_dir / "ionp_invalid_high_coating_permeability_config.txt";
	const fs::path cell_permeability_config_path = tmp_dir / "ionp_cell_permeability_config.txt";
	const std::vector<std::pair<std::string, fs::path>> cell_permeability_alias_configs = {
		{"Cell_P", tmp_dir / "ionp_cell_p_config.txt"},
		{"cell_probability", tmp_dir / "ionp_cell_probability_config.txt"},
		{"cell_permeation_probability", tmp_dir / "ionp_cell_permeation_probability_config.txt"},
		{"Cellular_permeability", tmp_dir / "ionp_cellular_permeability_config.txt"}
	};
	const fs::path invalid_low_cell_permeability_config_path = tmp_dir / "ionp_invalid_low_cell_permeability_config.txt";
	const fs::path invalid_high_cell_permeability_config_path = tmp_dir / "ionp_invalid_high_cell_permeability_config.txt";
	const fs::path invalid_percent_cell_permeability_config_path = tmp_dir / "ionp_invalid_percent_cell_permeability_config.txt";
	const fs::path ionp_per_config_path = tmp_dir / "ionp_per_config.txt";
	const fs::path ionp_per_percent_config_path = tmp_dir / "ionp_per_percent_config.txt";
	const fs::path invalid_low_ionp_per_config_path = tmp_dir / "ionp_invalid_low_per_config.txt";
	const fs::path invalid_high_ionp_per_config_path = tmp_dir / "ionp_invalid_high_per_config.txt";
	const fs::path hybrid_solver_config_path = tmp_dir / "ionp_hybrid_solver_config.txt";
	const fs::path invalid_solver_config_path = tmp_dir / "ionp_invalid_solver_config.txt";
	const fs::path surrounding_voxels_yes_config_path =
		tmp_dir / "ionp_surrounding_voxels_yes_config.txt";
	const fs::path surrounding_voxels_no_config_path =
		tmp_dir / "ionp_surrounding_voxels_no_config.txt";
	const fs::path invalid_surrounding_voxels_config_path =
		tmp_dir / "ionp_invalid_surrounding_voxels_config.txt";
	const fs::path surrounding_voxel_superparticles_config_path =
		tmp_dir / "ionp_surrounding_voxel_superparticles_config.txt";
	const fs::path negative_surrounding_voxel_superparticles_config_path =
		tmp_dir / "ionp_negative_surrounding_voxel_superparticles_config.txt";
	const fs::path superparticles_without_surrounding_voxels_config_path =
		tmp_dir / "ionp_superparticles_without_surrounding_voxels_config.txt";
	const fs::path superparticles_with_magnetic_moment_config_path =
		tmp_dir / "ionp_superparticles_with_magnetic_moment_config.txt";
	const fs::path custom_input_selector_config_path =
		tmp_dir / "custom_input_selector_config.txt";
	const fs::path ionp_path = tmp_dir / "ionp_basic_positions.dat";

	{
		std::ofstream cfg(config_path);
		write_common_config(cfg, 0.0, 0);
	}

	{
		std::ofstream cfg(is_config_path);
		write_common_config(cfg, 1.0, 2, true);
		cfg << "ImportanceSampling=Yes\n";
		cfg << "IS_near_cutoff=2.5\n";
		cfg << "IS_protons_inside=3\n";
		cfg << "IS_protons_near=4\n";
		cfg << "IS_protons_far=5\n";
		cfg << "IS_weight_inside=2\n";
		cfg << "IS_weight_near=3\n";
		cfg << "IS_weight_far=4\n";
		cfg << "IS_weight_bound=5\n";
	}

	{
		std::ofstream cfg(auto_weight_config_path);
		write_common_config(cfg, 1.0, 0);
		cfg << "ImportanceSampling=Yes\n";
		cfg << "IS_near_cutoff=2.5\n";
		cfg << "IS_protons_inside=3\n";
		cfg << "IS_protons_near=4\n";
		cfg << "IS_protons_far=5\n";
	}

	{
		std::ofstream cfg(joint_is_config_path);
		write_common_config(cfg, 5.0, 0);
		cfg << "ImportanceSampling=Yes\n";
		cfg << "IS_cell_near_cutoff=6\n";
		cfg << "IS_intracellular_ionp_near_cutoff=1\n";
		cfg << "IS_intracellular_ionp_far_cutoff=2\n";
		cfg << "IS_extracellular_ionp_near_cutoff=3\n";
		cfg << "IS_extracellular_ionp_far_cutoff=4\n";
		cfg << "IS_protons_cell_inside__intracellular_ionp_near__extracellular_ionp_far=7\n";
		cfg << "IS_weight_cell_inside__intracellular_ionp_near__extracellular_ionp_far=0.8\n";
		cfg << "IS_protons_cell_near__intracellular_ionp_intermediate__extracellular_ionp_near=3\n";
		cfg << "IS_bound_protons_intracellular_ionp=2\n";
		cfg << "IS_weight_bound_intracellular_ionp=0.2\n";
	}

	{
		std::ofstream cfg(overlapping_is_config_path);
		write_common_config(cfg, 5.0, 0);
		cfg << "ImportanceSampling=Yes\n";
		cfg << "IS_cell_near_cutoff=6\n";
		cfg << "IS_intracellular_ionp_near_cutoff=1\n";
		cfg << "IS_intracellular_ionp_far_cutoff=2\n";
		cfg << "IS_extracellular_ionp_near_cutoff=3\n";
		cfg << "IS_extracellular_ionp_far_cutoff=4\n";
		cfg << "IS_protons_uniform=11\n";
		cfg << "IS_protons_cell_inside=12\n";
		cfg << "IS_protons_cell_near=13\n";
		cfg << "IS_protons_cell_far=14\n";
		cfg << "IS_protons_intracellular_ionp_near=15\n";
		cfg << "IS_protons_intracellular_ionp_intermediate=16\n";
		cfg << "IS_protons_intracellular_ionp_far=17\n";
		cfg << "IS_protons_extracellular_ionp_near=18\n";
		cfg << "IS_protons_extracellular_ionp_intermediate=19\n";
		cfg << "IS_protons_extracellular_ionp_far=20\n";
	}

	{
		std::ofstream cfg(overlapping_without_uniform_config_path);
		write_common_config(cfg, 0.0, 0);
		cfg << "ImportanceSampling=Yes\n";
		cfg << "IS_protons_cell_far=1\n";
	}

	{
		std::ofstream cfg(overlapping_negative_count_config_path);
		write_common_config(cfg, 0.0, 0);
		cfg << "ImportanceSampling=Yes\n";
		cfg << "IS_protons_uniform=1\n";
		cfg << "IS_protons_cell_far=-1\n";
	}

	{
		std::ofstream cfg(negative_bound_count_config_path);
		write_common_config(cfg, 0.0, -1);
	}

	{
		std::ofstream cfg(overlapping_far_without_cutoffs_config_path);
		write_common_config(cfg, 0.0, 0);
		cfg << "ImportanceSampling=Yes\n";
		cfg << "IS_protons_uniform=1\n";
		cfg << "IS_protons_extracellular_ionp_far=1\n";
	}

	{
		std::ofstream cfg(mixed_overlapping_legacy_config_path);
		write_common_config(cfg, 1.0, 0);
		cfg << "ImportanceSampling=Yes\n";
		cfg << "IS_protons_uniform=1\n";
		// Presence is invalid even when the legacy value is zero.
		cfg << "IS_protons_inside=0\n";
	}

	{
		std::ofstream cfg(mixed_overlapping_joint_config_path);
		write_common_config(cfg, 0.0, 0);
		cfg << "ImportanceSampling=Yes\n";
		cfg << "IS_protons_uniform=1\n";
		// Presence is invalid even when the joint value is zero.
		cfg << "IS_weight_cell_far__intracellular_ionp_far__extracellular_ionp_far=0\n";
	}

	{
		std::ofstream cfg(mixed_overlapping_weight_config_path);
		write_common_config(cfg, 0.0, 0);
		cfg << "ImportanceSampling=Yes\n";
		cfg << "IS_protons_uniform=1\n";
		// Proposal weights are derived; manual legacy weights are incompatible,
		// including an explicitly configured zero.
		cfg << "IS_weight_far=0\n";
	}

	{
		std::ofstream cfg(ionp_only_is_with_cell_config_path);
		write_common_config(cfg, 5.0, 0);
		cfg << "ImportanceSampling=Yes\n";
		cfg << "IS_ionp_inside_near_cutoff=1\n";
		cfg << "IS_ionp_inside_far_cutoff=2\n";
		cfg << "IS_ionp_outside_near_cutoff=1\n";
		cfg << "IS_ionp_outside_far_cutoff=2\n";
		cfg << "IS_protons_cell_inside_ionp_inside_near_ionp_outside_far=7\n";
		cfg << "IS_protons_cell_far_ionp_inside_far_ionp_outside_near=5\n";
	}

	{
		std::ofstream cfg(near_cell_without_cutoff_config_path);
		write_common_config(cfg, 5.0, 0);
		cfg << "ImportanceSampling=Yes\n";
		cfg << "IS_protons_cell_near_ionp_inside_far_ionp_outside_far=1\n";
	}

	{
		std::ofstream cfg(invalid_is_config_path);
		write_common_config(cfg, 1.0, 1);
		cfg << "ImportanceSampling=Yes\n";
		cfg << "IS_near_cutoff=2.5\n";
		cfg << "IS_protons_inside=1\n";
		cfg << "IS_protons_near=1\n";
		cfg << "IS_protons_far=1\n";
		cfg << "IS_weight_inside=1\n";
		cfg << "IS_weight_near=1\n";
		cfg << "IS_weight_far=1\n";
	}

	{
		std::ofstream cfg(mixed_is_config_path);
		write_common_config(cfg, 5.0, 0);
		cfg << "ImportanceSampling=Yes\n";
		cfg << "IS_cell_near_cutoff=6\n";
		cfg << "IS_protons_inside=1\n";
		cfg << "IS_protons_cell_inside_ionp_inside_far_ionp_outside_far=1\n";
	}

	{
		std::ofstream cfg(incomplete_cutoff_config_path);
		write_common_config(cfg, 0.0, 0);
		cfg << "ImportanceSampling=Yes\n";
		cfg << "IS_ionp_inside_near_cutoff=1\n";
		cfg << "IS_protons_cell_far_ionp_inside_near_ionp_outside_far=1\n";
	}

	{
		std::ofstream cfg(coating_diffusion_config_path);
		write_common_config(cfg, 0.0, 0);
		cfg << "Coating_D=1e-10\n";
	}

	{
		std::ofstream cfg(cell_diffusion_config_path);
		write_common_config(cfg, 1.0, 0);
		cfg << "Cell_D=1e-9\n";
	}

	for (std::size_t alias_index = 0; alias_index < cell_diffusion_alias_configs.size(); ++alias_index) {
		std::ofstream cfg(cell_diffusion_alias_configs[alias_index].second);
		write_common_config(cfg, 1.0, 0);
		cfg << cell_diffusion_alias_configs[alias_index].first << "="
			<< (1e-9 * static_cast<double>(alias_index + 2)) << "\n";
	}

	{
		std::ofstream cfg(invalid_cell_diffusion_config_path);
		write_common_config(cfg, 1.0, 0);
		cfg << "Cell_D=-1e-9\n";
	}

	{
		std::ofstream cfg(coating_permeability_config_path);
		write_common_config(cfg, 0.0, 0);
		cfg << "Coating_permeability=0.25\n";
	}

	{
		std::ofstream cfg(invalid_low_coating_permeability_config_path);
		write_common_config(cfg, 0.0, 0);
		cfg << "Coating_permeability=-0.1\n";
	}

	{
		std::ofstream cfg(invalid_high_coating_permeability_config_path);
		write_common_config(cfg, 0.0, 0);
		cfg << "Coating_permeability=1.1\n";
	}

	{
		std::ofstream cfg(cell_permeability_config_path);
		write_common_config(cfg, 1.0, 0);
		cfg << "Cell_permeability=0.25\n";
	}

	for (std::size_t alias_index = 0; alias_index < cell_permeability_alias_configs.size(); ++alias_index) {
		std::ofstream cfg(cell_permeability_alias_configs[alias_index].second);
		write_common_config(cfg, 1.0, 0);
		cfg << cell_permeability_alias_configs[alias_index].first << "="
			<< (0.1 * static_cast<double>(alias_index + 1)) << "\n";
	}

	{
		std::ofstream cfg(invalid_low_cell_permeability_config_path);
		write_common_config(cfg, 1.0, 0);
		cfg << "Cell_permeability=-0.1\n";
	}

	{
		std::ofstream cfg(invalid_high_cell_permeability_config_path);
		write_common_config(cfg, 1.0, 0);
		cfg << "Cell_permeability=1.1\n";
	}

	{
		std::ofstream cfg(invalid_percent_cell_permeability_config_path);
		write_common_config(cfg, 1.0, 0);
		cfg << "Cell_permeability=0.5%\n";
	}

	{
		std::ofstream cfg(ionp_per_config_path);
		write_common_config(cfg, 0.0, 0);
		cfg << "IONP_per=0.73\n";
	}

	{
		std::ofstream cfg(ionp_per_percent_config_path);
		write_common_config(cfg, 0.0, 0);
		cfg << "IONP_per=21\n";
	}

	{
		std::ofstream cfg(invalid_low_ionp_per_config_path);
		write_common_config(cfg, 0.0, 0);
		cfg << "IONP_per=-0.1\n";
	}

	{
		std::ofstream cfg(invalid_high_ionp_per_config_path);
		write_common_config(cfg, 0.0, 0);
		cfg << "IONP_per=101\n";
	}

	{
		std::ofstream cfg(hybrid_solver_config_path);
		write_common_config(cfg, 0.0, 0);
		cfg << "BField_solver=corrected_hybrid\n";
	}

	{
		std::ofstream cfg(invalid_solver_config_path);
		write_common_config(cfg, 0.0, 0);
		cfg << "BField_solver=tree\n";
	}

	{
		std::ofstream cfg(surrounding_voxels_yes_config_path);
		write_common_config(cfg, 0.0, 0);
		cfg << "Surrounding_voxels=Yes\n";
	}

	{
		std::ofstream cfg(surrounding_voxels_no_config_path);
		write_common_config(cfg, 0.0, 0);
		cfg << "Surrounding_voxels=No\n";
		cfg << "Surrounding_voxel_superparticles=0\n";
	}

	{
		std::ofstream cfg(invalid_surrounding_voxels_config_path);
		write_common_config(cfg, 0.0, 0);
		cfg << "Surrounding_voxels=Sometimes\n";
	}

	{
		std::ofstream cfg(surrounding_voxel_superparticles_config_path);
		write_common_config(cfg, 0.0, 0);
		cfg << "Surrounding_voxels=Yes\n";
		cfg << "Surrounding_voxel_superparticles=10\n";
	}

	{
		std::ofstream cfg(negative_surrounding_voxel_superparticles_config_path);
		write_common_config(cfg, 0.0, 0);
		cfg << "Surrounding_voxels=Yes\n";
		cfg << "Surrounding_voxel_superparticles=-1\n";
	}

	{
		std::ofstream cfg(superparticles_without_surrounding_voxels_config_path);
		write_common_config(cfg, 0.0, 0);
		cfg << "Surrounding_voxel_superparticles=10\n";
	}

	{
		std::ofstream cfg(superparticles_with_magnetic_moment_config_path);
		write_common_config(cfg, 0.0, 0);
		cfg << "Surrounding_voxels=Yes\n";
		cfg << "Surrounding_voxel_superparticles=10\n";
		cfg << "Magnetic Moment=1\n";
	}

	{
		std::ofstream cfg(custom_input_selector_config_path);
		write_common_config(cfg, 0.0, 0);
		cfg << "IONP position file=ionp positions.dat\n";
		cfg << "AGG_POSITION_FILE=aggregate positions.txt\n";
		cfg << "Cell_position_file=cell positions.txt\n";
		cfg << "Cell load file=cell loads.dat\n";
	}

	{
		std::ofstream ionp(ionp_path);
		ionp << "0 0 0\n";
		ionp << "3 0 0\n";
	}

	Data dat(config_path.string());
	dat.readData();

	assert(std::fabs(dat.get_radius_ionp() - 1.0) < 1e-12);
	assert(dat.get_ionp_number() == 1);
	assert(std::fabs(dat.get_matrix_xlength() - 10.0) < 1e-12);
	assert(dat.get_Cell_hist() == false);
	assert(dat.get_variable_time() == false);
	if (dat.get_surrounding_voxels()) {
		return 1;
	}
	if (dat.get_surrounding_voxel_superparticles() != 0) {
		return 1;
	}
	assert(dat.get_importance_sampling().enabled == false);
	assert(dat.has_coating_diffusion_coefficient() == false);
	assert(dat.has_cell_diffusion_coefficient() == false);
	assert(dat.has_coating_permeability() == false);
	assert(std::fabs(dat.get_coating_permeability() - 1.0) < 1e-12);
	assert(dat.has_cell_permeability() == false);
	assert(std::fabs(dat.get_cell_permeability()) < 1e-12);
	assert(std::fabs(dat.get_ionp_per() - 0.21) < 1e-12);
	assert(dat.get_bfield_solver() == BFieldSolver::Current);
	assert(!dat.get_ionp_position_file().has_value());
	assert(!dat.get_agg_position_file().has_value());
	assert(!dat.get_cell_position_file().has_value());
	assert(!dat.get_cell_load_file().has_value());

	Data selector_dat(custom_input_selector_config_path.string());
	selector_dat.readData();
	assert(selector_dat.get_ionp_position_file() ==
		std::optional<std::string>("ionp positions.dat"));
	assert(selector_dat.get_agg_position_file() ==
		std::optional<std::string>("aggregate positions.txt"));
	assert(selector_dat.get_cell_position_file() ==
		std::optional<std::string>("cell positions.txt"));
	assert(selector_dat.get_cell_load_file() ==
		std::optional<std::string>("cell loads.dat"));

	Data hybrid_solver_dat(hybrid_solver_config_path.string());
	hybrid_solver_dat.readData();
	assert(hybrid_solver_dat.get_bfield_solver() ==
		   BFieldSolver::CorrectedHybrid);

	bool invalid_solver_threw = false;
	try {
		Data invalid_solver_dat(invalid_solver_config_path.string());
		invalid_solver_dat.readData();
	} catch (const std::runtime_error&) {
		invalid_solver_threw = true;
	}
	if (!invalid_solver_threw) {
		return 1;
	}

	Data surrounding_voxels_yes_dat(
		surrounding_voxels_yes_config_path.string());
	surrounding_voxels_yes_dat.readData();
	if (!surrounding_voxels_yes_dat.get_surrounding_voxels()) {
		return 1;
	}
	if (surrounding_voxels_yes_dat.get_surrounding_voxel_superparticles() != 0) {
		return 1;
	}

	Data surrounding_voxels_no_dat(
		surrounding_voxels_no_config_path.string());
	surrounding_voxels_no_dat.readData();
	if (surrounding_voxels_no_dat.get_surrounding_voxels()) {
		return 1;
	}
	if (surrounding_voxels_no_dat.get_surrounding_voxel_superparticles() != 0) {
		return 1;
	}

	bool invalid_surrounding_voxels_threw = false;
	try {
		Data invalid_surrounding_voxels_dat(
			invalid_surrounding_voxels_config_path.string());
		invalid_surrounding_voxels_dat.readData();
	} catch (const std::runtime_error& error) {
		const std::string message = error.what();
		invalid_surrounding_voxels_threw =
			message.find("Expected Yes/No value") != std::string::npos &&
			message.find("Surrounding_voxels") != std::string::npos &&
			message.find("Sometimes") != std::string::npos;
	}
	if (!invalid_surrounding_voxels_threw) {
		return 1;
	}

	Data surrounding_voxel_superparticles_dat(
		surrounding_voxel_superparticles_config_path.string());
	surrounding_voxel_superparticles_dat.readData();
	if (!surrounding_voxel_superparticles_dat.get_surrounding_voxels() ||
		surrounding_voxel_superparticles_dat.get_surrounding_voxel_superparticles() != 10) {
		return 1;
	}

	bool negative_superparticles_threw = false;
	try {
		Data negative_superparticles_dat(
			negative_surrounding_voxel_superparticles_config_path.string());
		negative_superparticles_dat.readData();
	} catch (const std::runtime_error& error) {
		const std::string message = error.what();
		negative_superparticles_threw =
			message.find("Surrounding_voxel_superparticles") != std::string::npos &&
			message.find("non-negative integer") != std::string::npos;
	}
	if (!negative_superparticles_threw) {
		return 1;
	}

	bool missing_surrounding_voxels_threw = false;
	try {
		Data missing_surrounding_voxels_dat(
			superparticles_without_surrounding_voxels_config_path.string());
		missing_surrounding_voxels_dat.readData();
	} catch (const std::runtime_error& error) {
		const std::string message = error.what();
		missing_surrounding_voxels_threw =
			message.find("Surrounding_voxel_superparticles") != std::string::npos &&
			message.find("Surrounding_voxels = Yes") != std::string::npos;
	}
	if (!missing_surrounding_voxels_threw) {
		return 1;
	}

	bool magnetic_moment_superparticles_threw = false;
	try {
		Data magnetic_moment_superparticles_dat(
			superparticles_with_magnetic_moment_config_path.string());
		magnetic_moment_superparticles_dat.readData();
	} catch (const std::runtime_error& error) {
		const std::string message = error.what();
		magnetic_moment_superparticles_threw =
			message.find("Surrounding_voxel_superparticles") != std::string::npos &&
			message.find("Magnetic Moment = 0") != std::string::npos;
	}
	if (!magnetic_moment_superparticles_threw) {
		return 1;
	}

	SimMatrix sim(
		0.0,
		0.0,
		0.0,
		dat.get_matrix_xlength(),
		dat.get_matrix_ylength(),
		dat.get_matrix_zlength(),
		dat.get_matrix_nxnodes(),
		dat.get_matrix_nynodes(),
		dat.get_matrix_nznodes());

	std::vector<Ionp> ionps;
	ionps.reserve(2);
	dat.readIONP(ionp_path.string(), ionps, sim, dat.get_radius_ionp(), dat.get_thickness_IONPcoating());

	assert(ionps.size() == 2);
	assert(std::fabs(ionps.back().xt().x() - 3.0) < 1e-12);

	Data is_dat(is_config_path.string());
	is_dat.readData();
	const ImportanceSamplingConfig importance_sampling = is_dat.get_importance_sampling();
	assert(is_dat.get_variable_time());
	assert(importance_sampling.enabled);
	assert(std::fabs(importance_sampling.cell_near_cutoff - 2.5) < 1e-12);
	assert(importance_sampling.protons_inside == 3);
	assert(importance_sampling.protons_near == 4);
	assert(importance_sampling.protons_far == 5);
	assert(std::fabs(importance_sampling.weight_inside - 2.0) < 1e-12);
	assert(std::fabs(importance_sampling.weight_near - 3.0) < 1e-12);
	assert(std::fabs(importance_sampling.weight_far - 4.0) < 1e-12);
	assert(std::fabs(importance_sampling.weight_bound - 5.0) < 1e-12);
	is_dat.set_effective_proton_number(importance_sampling.total_free_protons());

	Data auto_weight_dat(auto_weight_config_path.string());
	auto_weight_dat.readData();
	const ImportanceSamplingConfig auto_weight_importance_sampling = auto_weight_dat.get_importance_sampling();
	if (!auto_weight_importance_sampling.enabled ||
		auto_weight_importance_sampling.weight_inside != 0.0 ||
		auto_weight_importance_sampling.weight_near != 0.0 ||
		auto_weight_importance_sampling.weight_far != 0.0) {
		return 1;
	}

	Data joint_is_dat(joint_is_config_path.string());
	joint_is_dat.readData();
	const ImportanceSamplingConfig joint_importance =
		joint_is_dat.get_importance_sampling();
	const std::size_t joint_region_index = joint_sample_region_index(
		CellSampleRegion::Inside,
		IonpDistanceRegion::Near,
		IonpDistanceRegion::Far);
	const std::size_t intermediate_region_index = joint_sample_region_index(
		CellSampleRegion::Near,
		IonpDistanceRegion::Intermediate,
		IonpDistanceRegion::Near);
	if (!joint_importance.joint_regions_configured ||
		joint_importance.joint_protons[joint_region_index] != 7 ||
		joint_importance.joint_protons[intermediate_region_index] != 3 ||
		std::fabs(joint_importance.joint_weights[joint_region_index] - 0.8) >= 1e-12 ||
		std::fabs(joint_importance.intracellular_ionp_near_cutoff - 1.0) >= 1e-12 ||
		std::fabs(joint_importance.intracellular_ionp_far_cutoff - 2.0) >= 1e-12 ||
		std::fabs(joint_importance.extracellular_ionp_near_cutoff - 3.0) >= 1e-12 ||
		std::fabs(joint_importance.extracellular_ionp_far_cutoff - 4.0) >= 1e-12 ||
		!joint_importance.split_bound_configured ||
		joint_importance.bound_protons_intracellular != 2 ||
		std::fabs(joint_importance.weight_bound_intracellular - 0.2) >= 1e-12) {
		return 1;
	}

	Data overlapping_is_dat(overlapping_is_config_path.string());
	overlapping_is_dat.readData();
	const ImportanceSamplingConfig overlapping_importance =
		overlapping_is_dat.get_importance_sampling();
	if (!overlapping_importance.enabled ||
		overlapping_importance.mode !=
			ImportanceSamplingMode::OverlappingProposal ||
		!overlapping_importance.overlapping_proposals_enabled()) {
		return 1;
	}
	int expected_overlapping_total = 0;
	for (std::size_t index = 0;
		 index < kFreeSamplingProposalCount;
		 ++index) {
		const int expected_count = static_cast<int>(index) + 11;
		if (overlapping_importance.proposal_protons[index] != expected_count) {
			return 1;
		}
		expected_overlapping_total += expected_count;
	}
	if (overlapping_importance.total_proposal_protons() !=
			expected_overlapping_total ||
		overlapping_importance.total_free_protons() !=
			expected_overlapping_total) {
		return 1;
	}
	const fs::path proposal_output_dir = tmp_dir / "proposal_output";
	overlapping_is_dat.set_output_context(
		proposal_output_dir, tmp_dir / "proposal_logs");
	overlapping_is_dat.set_effective_proton_number(expected_overlapping_total);
	const std::vector<std::complex<double>> proposal_signal = {
		std::complex<double>(1.0, 0.25)
	};
	const std::vector<double> proposal_time = {0.0};
	overlapping_is_dat.writeImportanceProposalData(
		proposal_signal,
		proposal_time,
		FreeSamplingProposal::CellInside,
		false);
	overlapping_is_dat.writeImportanceProposalData(
		proposal_signal,
		proposal_time,
		FreeSamplingProposal::CellInside,
		true);
	bool found_proposal_mean = false;
	bool found_proposal_contribution = false;
	for (const fs::directory_entry& entry :
		 fs::directory_iterator(proposal_output_dir)) {
		const std::string name = entry.path().filename().string();
		found_proposal_mean = found_proposal_mean ||
			name.find("_proposal_cell_inside_mean_") != std::string::npos;
		found_proposal_contribution = found_proposal_contribution ||
			name.find("_proposal_cell_inside_mis_contribution_") !=
				std::string::npos;
	}
	if (!found_proposal_mean || !found_proposal_contribution) {
		return 1;
	}

	bool overlapping_without_uniform_threw = false;
	try {
		Data invalid_dat(overlapping_without_uniform_config_path.string());
		invalid_dat.readData();
	}
	catch (const std::runtime_error&) {
		overlapping_without_uniform_threw = true;
	}
	if (!overlapping_without_uniform_threw) {
		return 1;
	}

	bool overlapping_negative_count_threw = false;
	try {
		Data invalid_dat(overlapping_negative_count_config_path.string());
		invalid_dat.readData();
	}
	catch (const std::runtime_error&) {
		overlapping_negative_count_threw = true;
	}
	if (!overlapping_negative_count_threw) {
		return 1;
	}

	bool negative_bound_count_threw = false;
	try {
		Data invalid_dat(negative_bound_count_config_path.string());
		invalid_dat.readData();
	}
	catch (const std::runtime_error&) {
		negative_bound_count_threw = true;
	}
	if (!negative_bound_count_threw) {
		return 1;
	}

	bool overlapping_far_without_cutoffs_threw = false;
	try {
		Data invalid_dat(overlapping_far_without_cutoffs_config_path.string());
		invalid_dat.readData();
	}
	catch (const std::runtime_error&) {
		overlapping_far_without_cutoffs_threw = true;
	}
	if (!overlapping_far_without_cutoffs_threw) {
		return 1;
	}

	bool mixed_overlapping_legacy_threw = false;
	try {
		Data invalid_dat(mixed_overlapping_legacy_config_path.string());
		invalid_dat.readData();
	}
	catch (const std::runtime_error&) {
		mixed_overlapping_legacy_threw = true;
	}
	if (!mixed_overlapping_legacy_threw) {
		return 1;
	}

	bool mixed_overlapping_joint_threw = false;
	try {
		Data invalid_dat(mixed_overlapping_joint_config_path.string());
		invalid_dat.readData();
	}
	catch (const std::runtime_error&) {
		mixed_overlapping_joint_threw = true;
	}
	if (!mixed_overlapping_joint_threw) {
		return 1;
	}

	bool mixed_overlapping_weight_threw = false;
	try {
		Data invalid_dat(mixed_overlapping_weight_config_path.string());
		invalid_dat.readData();
	}
	catch (const std::runtime_error&) {
		mixed_overlapping_weight_threw = true;
	}
	if (!mixed_overlapping_weight_threw) {
		return 1;
	}

	Data ionp_only_is_with_cell_dat(
		ionp_only_is_with_cell_config_path.string());
	ionp_only_is_with_cell_dat.readData();
	const ImportanceSamplingConfig ionp_only_is_with_cell =
		ionp_only_is_with_cell_dat.get_importance_sampling();
	const std::size_t ionp_only_region_index = joint_sample_region_index(
		CellSampleRegion::Inside,
		IonpDistanceRegion::Near,
		IonpDistanceRegion::Far);
	const std::size_t ionp_only_outside_region_index = joint_sample_region_index(
		CellSampleRegion::Far,
		IonpDistanceRegion::Far,
		IonpDistanceRegion::Near);
	if (!ionp_only_is_with_cell.enabled ||
		ionp_only_is_with_cell.cell_near_cutoff != 0. ||
		ionp_only_is_with_cell.joint_protons[ionp_only_region_index] != 7 ||
		ionp_only_is_with_cell.joint_protons[ionp_only_outside_region_index] != 5) {
		return 1;
	}

	bool near_cell_without_cutoff_threw = false;
	try {
		Data near_cell_without_cutoff_dat(
			near_cell_without_cutoff_config_path.string());
		near_cell_without_cutoff_dat.readData();
	} catch (const std::runtime_error&) {
		near_cell_without_cutoff_threw = true;
	}
	if (!near_cell_without_cutoff_threw) {
		return 1;
	}

	bool threw = false;
	try {
		Data invalid_dat(invalid_is_config_path.string());
		invalid_dat.readData();
	} catch (const std::runtime_error&) {
		threw = true;
	}
	if (!threw) {
		return 1;
	}

	bool mixed_is_threw = false;
	try {
		Data mixed_is_dat(mixed_is_config_path.string());
		mixed_is_dat.readData();
	} catch (const std::runtime_error&) {
		mixed_is_threw = true;
	}
	if (!mixed_is_threw) {
		return 1;
	}

	bool incomplete_cutoff_threw = false;
	try {
		Data incomplete_cutoff_dat(incomplete_cutoff_config_path.string());
		incomplete_cutoff_dat.readData();
	} catch (const std::runtime_error&) {
		incomplete_cutoff_threw = true;
	}
	if (!incomplete_cutoff_threw) {
		return 1;
	}

	Data coating_diffusion_dat(coating_diffusion_config_path.string());
	coating_diffusion_dat.readData();
	assert(coating_diffusion_dat.has_coating_diffusion_coefficient());
	assert(std::fabs(coating_diffusion_dat.get_coating_diffusion_coefficient() - 1e-10) < 1e-20);
	assert(coating_diffusion_dat.has_coating_permeability() == false);

	Data cell_diffusion_dat(cell_diffusion_config_path.string());
	cell_diffusion_dat.readData();
	assert(cell_diffusion_dat.has_cell_diffusion_coefficient());
	assert(std::fabs(cell_diffusion_dat.get_cell_diffusion_coefficient() - 1e-9) < 1e-20);

	for (std::size_t alias_index = 0; alias_index < cell_diffusion_alias_configs.size(); ++alias_index) {
		Data alias_dat(cell_diffusion_alias_configs[alias_index].second.string());
		alias_dat.readData();
		assert(alias_dat.has_cell_diffusion_coefficient());
		assert(std::fabs(alias_dat.get_cell_diffusion_coefficient() -
			(1e-9 * static_cast<double>(alias_index + 2))) < 1e-20);
	}

	bool invalid_cell_diffusion_threw = false;
	try {
		Data invalid_cell_diffusion_dat(invalid_cell_diffusion_config_path.string());
		invalid_cell_diffusion_dat.readData();
	} catch (const std::runtime_error&) {
		invalid_cell_diffusion_threw = true;
	}
	if (!invalid_cell_diffusion_threw) {
		return 1;
	}

	Data coating_permeability_dat(coating_permeability_config_path.string());
	coating_permeability_dat.readData();
	assert(coating_permeability_dat.has_coating_permeability());
	assert(std::fabs(coating_permeability_dat.get_coating_permeability() - 0.25) < 1e-12);
	assert(coating_permeability_dat.has_coating_diffusion_coefficient() == false);

	bool invalid_low_permeability_threw = false;
	try {
		Data invalid_low_permeability_dat(invalid_low_coating_permeability_config_path.string());
		invalid_low_permeability_dat.readData();
	} catch (const std::runtime_error&) {
		invalid_low_permeability_threw = true;
	}
	if (!invalid_low_permeability_threw) {
		return 1;
	}

	bool invalid_high_permeability_threw = false;
	try {
		Data invalid_high_permeability_dat(invalid_high_coating_permeability_config_path.string());
		invalid_high_permeability_dat.readData();
	} catch (const std::runtime_error&) {
		invalid_high_permeability_threw = true;
	}
	if (!invalid_high_permeability_threw) {
		return 1;
	}

	Data cell_permeability_dat(cell_permeability_config_path.string());
	cell_permeability_dat.readData();
	assert(cell_permeability_dat.has_cell_permeability());
	assert(std::fabs(cell_permeability_dat.get_cell_permeability() - 0.25) < 1e-12);

	for (std::size_t alias_index = 0; alias_index < cell_permeability_alias_configs.size(); ++alias_index) {
		Data alias_dat(cell_permeability_alias_configs[alias_index].second.string());
		alias_dat.readData();
		assert(alias_dat.has_cell_permeability());
		assert(std::fabs(alias_dat.get_cell_permeability() -
			(0.1 * static_cast<double>(alias_index + 1))) < 1e-12);
	}

	bool invalid_low_cell_permeability_threw = false;
	try {
		Data invalid_low_cell_permeability_dat(invalid_low_cell_permeability_config_path.string());
		invalid_low_cell_permeability_dat.readData();
	} catch (const std::runtime_error&) {
		invalid_low_cell_permeability_threw = true;
	}
	if (!invalid_low_cell_permeability_threw) {
		return 1;
	}

	bool invalid_high_cell_permeability_threw = false;
	try {
		Data invalid_high_cell_permeability_dat(invalid_high_cell_permeability_config_path.string());
		invalid_high_cell_permeability_dat.readData();
	} catch (const std::runtime_error&) {
		invalid_high_cell_permeability_threw = true;
	}
	if (!invalid_high_cell_permeability_threw) {
		return 1;
	}

	bool invalid_percent_cell_permeability_threw = false;
	try {
		Data invalid_percent_cell_permeability_dat(invalid_percent_cell_permeability_config_path.string());
		invalid_percent_cell_permeability_dat.readData();
	} catch (const std::runtime_error&) {
		invalid_percent_cell_permeability_threw = true;
	}
	if (!invalid_percent_cell_permeability_threw) {
		return 1;
	}

	Data ionp_per_dat(ionp_per_config_path.string());
	ionp_per_dat.readData();
	assert(std::fabs(ionp_per_dat.get_ionp_per() - 0.73) < 1e-12);

	Data ionp_per_percent_dat(ionp_per_percent_config_path.string());
	ionp_per_percent_dat.readData();
	assert(std::fabs(ionp_per_percent_dat.get_ionp_per() - 0.21) < 1e-12);

	bool invalid_low_ionp_per_threw = false;
	try {
		Data invalid_low_ionp_per_dat(invalid_low_ionp_per_config_path.string());
		invalid_low_ionp_per_dat.readData();
	} catch (const std::runtime_error&) {
		invalid_low_ionp_per_threw = true;
	}
	if (!invalid_low_ionp_per_threw) {
		return 1;
	}

	bool invalid_high_ionp_per_threw = false;
	try {
		Data invalid_high_ionp_per_dat(invalid_high_ionp_per_config_path.string());
		invalid_high_ionp_per_dat.readData();
	} catch (const std::runtime_error&) {
		invalid_high_ionp_per_threw = true;
	}
	if (!invalid_high_ionp_per_threw) {
		return 1;
	}

	{
		SimMatrix selection_sim(0.0, 0.0, 0.0, 1e-4, 1e-4, 1e-4, 0, 0, 0);
		Point3D origin;
		Point3D voxel_pos;
		std::vector<Point3D> cells;
		std::vector<Ionp> primary_ionps;
		const double primary_radius = 1e-7;
		const double secondary_radius = 5e-7;
		SimSpace::simulationSpace_outCell(primary_ionps,
										   selection_sim,
										   4,
										   0.0,
										   primary_radius,
										   secondary_radius,
										   1.0,
										   1e-12,
										   1e-12,
										   0.0,
										   origin,
										   0.0,
										   0,
										   cells,
										   voxel_pos,
										   false,
										   tmp_dir);
		for (const Ionp& ionp : primary_ionps) {
			if (std::fabs(ionp.radius() - primary_radius) >= 1e-8) {
				return 1;
			}
		}

		std::vector<Ionp> secondary_ionps;
		SimSpace::simulationSpace_outCell(secondary_ionps,
										   selection_sim,
										   4,
										   0.0,
										   primary_radius,
										   secondary_radius,
										   0.0,
										   1e-12,
										   1e-12,
										   0.0,
										   origin,
										   0.0,
										   0,
										   cells,
										   voxel_pos,
										   false,
										   tmp_dir);
		for (const Ionp& ionp : secondary_ionps) {
			if (std::fabs(ionp.radius() - secondary_radius) >= 1e-8) {
				return 1;
			}
		}

		std::vector<Point3D> fallback_cells;
		std::vector<Ionp> outside_fallback_cell;
		Point3D fallback_voxel_position;
		const double fallback_cell_radius = 1e-5;
		SimSpace::simulationSpace_outCell(outside_fallback_cell,
										  selection_sim,
										  1,
										  0.0,
										  primary_radius,
										  secondary_radius,
										  1.0,
										  0.0,
										  0.0,
										  0.0,
										  origin,
										  fallback_cell_radius,
										  1,
										  fallback_cells,
										  fallback_voxel_position,
										  false,
										  tmp_dir,
										  false);
		if (outside_fallback_cell.size() != 1 ||
			distance(outside_fallback_cell.front().xt(), origin) < fallback_cell_radius + primary_radius) {
			return 1;
		}

		Point3D configured_cell_origin(1e-5, -2e-5, 3e-5);
		Point3D nonzero_voxel_position(4e-6, 5e-6, -6e-6);
		std::vector<Point3D> configured_cells;
		std::vector<Ionp> no_ionps;
		SimSpace::simulationSpace_inCell(no_ionps,
										 selection_sim,
										 0,
										 0.0,
										 primary_radius,
										 secondary_radius,
										 1.0,
										 0.0,
										 0.0,
										 0.0,
										 configured_cell_origin,
										 1e-5,
										 1,
										 0.0,
										 0.0,
										 1,
										 configured_cells,
										 nonzero_voxel_position,
										 false,
										 tmp_dir,
										 false);
		if (configured_cells.size() != 1 || distance(configured_cells.front(), configured_cell_origin) > 1e-15) {
			return 1;
		}

		Point3D loaded_center(-2e-5, 1e-5, 4e-5);
		Point3D loaded_voxel_position(7e-6, -3e-6, 2e-6);
		std::vector<Point3D> loaded_cells = {loaded_center};
		SimSpace::simulationSpace_inCell(no_ionps,
										 selection_sim,
										 0,
										 0.0,
										 primary_radius,
										 secondary_radius,
										 1.0,
										 0.0,
										 0.0,
										 0.0,
										 configured_cell_origin,
										 1e-5,
										 1,
										 0.0,
										 0.0,
										 1,
										 loaded_cells,
										 loaded_voxel_position,
										 false,
										 tmp_dir,
										 false);
		if (loaded_cells.size() != 1 || distance(loaded_cells.front(), loaded_center) > 1e-15) {
			return 1;
		}
	}

	{
		rng_config::set_seed(20260817u);
		const double voxel_length = 12e-6;
		const double aggregate_radius = 0.5e-6;
		const double solid_equivalent_radius = 50e-9;
		const int aggregate_count = 2;
		const double pi = std::acos(-1.);
		const double target_mass = aggregate_count * 4. / 3. * pi * 5.24 *
			(167.4 / 231.4) * std::pow(solid_equivalent_radius * 100., 3.);
		const double extracellular_concentration = target_mass /
			(5.5845e-5 * std::pow(voxel_length, 3.) * 1000.);
		SimMatrix aggregate_sim(
			0., 0., 0., voxel_length, voxel_length, voxel_length, 0, 0, 0);
		Point3D cell_origin;
		const double cell_radius = 1e-6;
		std::vector<Point3D> aggregate_cells = {cell_origin};
		std::vector<Aggregate> extracellular_aggregates;
		std::vector<Ionp> mixed_ionps = {
			Ionp(cell_origin, 0., 100e-9, 0.)
		};
		const std::size_t intracellular_count = mixed_ionps.size();
		SimSpace::simulationSpace_outCellAggregates(
			extracellular_aggregates, mixed_ionps, aggregate_sim,
			aggregate_count, extracellular_concentration,
			aggregate_radius, 0., 4e-9, 12.5e-9, 0.21,
			0.48, 0.3, 50e-9, cell_origin, cell_radius,
			aggregate_cells);
		if (extracellular_aggregates.size() !=
			static_cast<std::size_t>(aggregate_count) ||
			mixed_ionps.size() <= intracellular_count) {
			return 1;
		}
		std::size_t recorded_aggregate_particles = 0;
		for (std::size_t left = 0;
			 left < extracellular_aggregates.size();
			 ++left) {
			const Aggregate& aggregate = extracellular_aggregates[left];
			recorded_aggregate_particles +=
				static_cast<std::size_t>(aggregate.nparticles());
			if (distance(aggregate.xt(), cell_origin) <
				cell_radius + aggregate.radius()) {
				return 1;
			}
			for (std::size_t right = left + 1;
				 right < extracellular_aggregates.size();
				 ++right) {
				if (distance(aggregate.xt(), extracellular_aggregates[right].xt()) <
					aggregate.radius() + extracellular_aggregates[right].radius()) {
					return 1;
				}
			}
		}
		if (recorded_aggregate_particles !=
			mixed_ionps.size() - intracellular_count) {
			return 1;
		}
		for (std::size_t index = intracellular_count;
			 index < mixed_ionps.size();
			 ++index) {
			const Ionp& ionp = mixed_ionps[index];
			if (distance(ionp.xt(), cell_origin) <
				cell_radius + ionp.outer_radius()) {
				return 1;
			}
			bool contained = false;
			for (const Aggregate& aggregate : extracellular_aggregates) {
				if (distance(ionp.xt(), aggregate.xt()) + ionp.outer_radius() <=
					aggregate.radius() + 1e-15) {
					contained = true;
					break;
				}
			}
			if (!contained) {
				return 1;
			}
		}
	}

	{
		// Exercise the no-cell branch and verify that extracellular aggregate
		// particles occupy the full spherical volume rather than a compact
		// center-outward lattice.
		rng_config::set_seed(20260818u);
		const double voxel_length = 12e-6;
		const double aggregate_radius = 1.9e-6;
		const double solid_equivalent_radius = 100e-9;
		const double pi = std::acos(-1.);
		const double target_mass = 4. / 3. * pi * 5.24 *
			(167.4 / 231.4) * std::pow(solid_equivalent_radius * 100., 3.);
		const double extracellular_concentration = target_mass /
			(5.5845e-5 * std::pow(voxel_length, 3.) * 1000.);
		SimMatrix aggregate_sim(
			0., 0., 0., voxel_length, voxel_length, voxel_length, 0, 0, 0);
		Point3D no_cell_origin;
		std::vector<Point3D> no_cells;
		std::vector<Aggregate> extracellular_aggregates;
		std::vector<Ionp> extracellular_ionps;
		SimSpace::simulationSpace_outCellAggregates(
			extracellular_aggregates, extracellular_ionps, aggregate_sim,
			1, extracellular_concentration, aggregate_radius, 0.,
			4e-9, 12.5e-9, 0.21, 0.48, 0.3, 50e-9,
			no_cell_origin, 0., no_cells);
		if (extracellular_aggregates.size() != 1 ||
			extracellular_ionps.size() < 100 ||
			extracellular_aggregates.front().nparticles() !=
				static_cast<int>(extracellular_ionps.size())) {
			return 1;
		}
		bool outer_volume_occupied = false;
		double normalized_volume_sum = 0.;
		for (std::size_t left = 0; left < extracellular_ionps.size(); ++left) {
			const Ionp& ionp = extracellular_ionps[left];
			const double center_distance =
				distance(ionp.xt(), extracellular_aggregates.front().xt());
			if (center_distance +
				ionp.outer_radius() > aggregate_radius + 1e-15) {
				return 1;
			}
			const double available_radius =
				aggregate_radius - ionp.outer_radius();
			if (center_distance > 0.8 * available_radius) {
				outer_volume_occupied = true;
			}
			const double normalized_radius = center_distance / available_radius;
			normalized_volume_sum += normalized_radius * normalized_radius *
				normalized_radius;
			for (std::size_t right = left + 1;
				 right < extracellular_ionps.size();
				 ++right) {
				if (distance(ionp.xt(), extracellular_ionps[right].xt()) <
					ionp.outer_radius() +
						extracellular_ionps[right].outer_radius()) {
					return 1;
				}
			}
		}
		const double mean_normalized_volume =
			normalized_volume_sum /
			static_cast<double>(extracellular_ionps.size());
		if (!outer_volume_occupied ||
			mean_normalized_volume < 0.35 ||
			mean_normalized_volume > 0.65) {
			return 1;
		}
	}

	{
		const fs::path selected_ini_root = tmp_dir / "selected_cell_inputs";
		const fs::path selected_cell_dir = selected_ini_root / "Cell_pos";
		const fs::path selected_load_dir = selected_ini_root / "Cell_load";
		fs::create_directories(selected_cell_dir);
		fs::create_directories(selected_load_dir);
		const fs::path selected_cell_path = selected_cell_dir / "selected cells.txt";
		const fs::path unselected_cell_path = selected_cell_dir / "other cells.dat";
		const fs::path selected_load_path = selected_load_dir / "selected loads.txt";
		const fs::path unselected_load_path = selected_load_dir / "other loads.dat";
		const fs::path mismatch_cell_path = selected_cell_dir / "mismatch cells.txt";
		const fs::path mismatch_load_path = selected_load_dir / "mismatch loads.txt";
		{
			std::ofstream out(selected_cell_path);
			out << "0 0 0\n";
		}
		{
			std::ofstream out(unselected_cell_path);
			out << "2 0 0\n";
		}
		{
			std::ofstream out(selected_load_path);
			out << "1\n";
		}
		{
			std::ofstream out(unselected_load_path);
			out << "2\n";
		}
		{
			std::ofstream out(mismatch_cell_path);
			out << "-2 0 0\n2 0 0\n";
		}
		{
			std::ofstream out(mismatch_load_path);
			out << "1\n";
		}

		SimMatrix cell_file_sim(0.0, 0.0, 0.0, 10.0, 10.0, 10.0, 0, 0, 0);
		Point3D cell_origin;
		Point3D voxel_position;
		std::vector<Point3D> selected_cells;
		std::vector<Ionp> no_selected_ionps;
		SimSpace::simulationSpace_inCell(no_selected_ionps,
			cell_file_sim, 0, 0.0, 0.01, 0.02, 1.0, 0.0, 0.0, 0.0,
			cell_origin, 0.5, 0, 0.0, 0.0, 1, selected_cells,
			voxel_position, false, selected_ini_root, false,
			selected_cell_path, selected_load_path);
		if (selected_cells.size() != 1 ||
			distance(selected_cells.front(), Point3D(0.0, 0.0, 0.0)) > 1e-15) {
			return 1;
		}

		bool missing_cell_pair_threw = false;
		try {
			std::vector<Point3D> cells;
			std::vector<Ionp> ionps_for_cells;
			Point3D local_voxel;
			SimSpace::simulationSpace_inCell(ionps_for_cells,
				cell_file_sim, 0, 0.0, 0.01, 0.02, 1.0, 0.0, 0.0, 0.0,
				cell_origin, 0.5, 0, 0.0, 0.0, 1, cells, local_voxel,
				true, selected_ini_root, false, selected_cell_path);
		} catch (const std::runtime_error&) {
			missing_cell_pair_threw = true;
		}
		if (!missing_cell_pair_threw) {
			return 1;
		}

		bool mismatched_cell_loads_threw = false;
		try {
			std::vector<Point3D> cells;
			std::vector<Ionp> ionps_for_cells;
			Point3D local_voxel;
			SimSpace::simulationSpace_inCell(ionps_for_cells,
				cell_file_sim, 0, 0.0, 0.01, 0.02, 1.0, 0.0, 0.0, 0.0,
				cell_origin, 0.5, 0, 0.0, 0.0, 1, cells, local_voxel,
				true, selected_ini_root, false,
				mismatch_cell_path, mismatch_load_path);
		} catch (const std::runtime_error&) {
			mismatched_cell_loads_threw = true;
		}
		if (!mismatched_cell_loads_threw) {
			return 1;
		}

		std::vector<Point3D> aggregate_cells;
		std::vector<Aggregate> no_aggregates;
		std::vector<Ionp> no_aggregate_ionps;
		SimSpace::simulationSpace(no_aggregates, no_aggregate_ionps,
			cell_file_sim, 0, 0, 0.0, 0.5, 0.01, 0.02, 1.0, 0.0, 0.0,
			0.0, 0.0, cell_origin, 0.5, 0, aggregate_cells, false,
			selected_ini_root, selected_cell_path);
		if (aggregate_cells.size() != 1 ||
			distance(aggregate_cells.front(), Point3D(0.0, 0.0, 0.0)) > 1e-15) {
			return 1;
		}

		const fs::path legacy_ini_root = tmp_dir / "legacy_cell_inputs";
		fs::create_directories(legacy_ini_root / "Cell_pos");
		fs::create_directories(legacy_ini_root / "Cell_load");
		{
			std::ofstream out(legacy_ini_root / "Cell_pos" / "a.txt");
			out << "-2 0 0\n";
		}
		{
			std::ofstream out(legacy_ini_root / "Cell_pos" / "b.dat");
			out << "2 0 0\n";
		}
		{
			std::ofstream out(legacy_ini_root / "Cell_load" / "a.txt");
			out << "1\n";
		}
		{
			std::ofstream out(legacy_ini_root / "Cell_load" / "b.dat");
			out << "2\n";
		}
		std::vector<Point3D> legacy_cells;
		std::vector<Ionp> no_legacy_ionps;
		Point3D legacy_voxel;
		SimSpace::simulationSpace_inCell(no_legacy_ionps,
			cell_file_sim, 0, 0.0, 0.01, 0.02, 1.0, 0.0, 0.0, 0.0,
			cell_origin, 0.5, 0, 0.0, 0.0, 1, legacy_cells, legacy_voxel,
			false, legacy_ini_root, false);
		if (legacy_cells.size() != 2) {
			return 1;
		}
	}

	fs::remove(config_path);
	fs::remove(is_config_path);
	fs::remove(auto_weight_config_path);
	fs::remove(joint_is_config_path);
	fs::remove(ionp_only_is_with_cell_config_path);
	fs::remove(near_cell_without_cutoff_config_path);
	fs::remove(invalid_is_config_path);
	fs::remove(coating_diffusion_config_path);
	fs::remove(cell_diffusion_config_path);
	for (const auto& alias_config : cell_diffusion_alias_configs) {
		fs::remove(alias_config.second);
	}
	fs::remove(invalid_cell_diffusion_config_path);
	fs::remove(coating_permeability_config_path);
	fs::remove(invalid_low_coating_permeability_config_path);
	fs::remove(invalid_high_coating_permeability_config_path);
	fs::remove(cell_permeability_config_path);
	for (const auto& alias_config : cell_permeability_alias_configs) {
		fs::remove(alias_config.second);
	}
	fs::remove(invalid_low_cell_permeability_config_path);
	fs::remove(invalid_high_cell_permeability_config_path);
	fs::remove(invalid_percent_cell_permeability_config_path);
	fs::remove(ionp_per_config_path);
	fs::remove(ionp_per_percent_config_path);
	fs::remove(invalid_low_ionp_per_config_path);
	fs::remove(invalid_high_ionp_per_config_path);
	fs::remove(hybrid_solver_config_path);
	fs::remove(invalid_solver_config_path);
	fs::remove(surrounding_voxels_yes_config_path);
	fs::remove(surrounding_voxels_no_config_path);
	fs::remove(invalid_surrounding_voxels_config_path);
	fs::remove(surrounding_voxel_superparticles_config_path);
	fs::remove(negative_surrounding_voxel_superparticles_config_path);
	fs::remove(superparticles_without_surrounding_voxels_config_path);
	fs::remove(superparticles_with_magnetic_moment_config_path);
	fs::remove(custom_input_selector_config_path);
	fs::remove(ionp_path);
	fs::remove_all(tmp_dir);
	return 0;
}
