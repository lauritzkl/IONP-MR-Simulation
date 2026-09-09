/**
 * @file Data.hpp
 *
 * Configuration state, auxiliary input readers and result writers. Call
 * readData() before consuming configuration getters. Geometry readers fill
 * caller-owned vectors; output paths are supplied by main.cpp. Units, keys
 * and file columns are described in docs/simulations.md.
 */

/*****************************************************************************
*Class to read configuration files

*10-05-2022 (Lauritz Klünder)
******************************************************************************/

#ifndef DATA_HPP_
#define DATA_HPP_

#include <complex>
#include <filesystem>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

#include "Aggregate.hpp"
#include "BFieldSolver.hpp"
#include "ImportanceSampling.hpp"
#include "Ionp.hpp"
#include "Point3D.hpp"
#include "Proton.hpp"
#include "SimMatrix.hpp"

class Data {

public:

    Data(std::string filename_in);

    // Parse configuration, validate supported combinations and resolve defaults.
    void readData(  );
    void set_output_context( std::filesystem::path tmp_dir, std::filesystem::path log_dir );
    // Geometry readers append records; readHist instead replaces its sample vector.
    // Numeric stream extraction stops at malformed input; no header row is skipped.
    void readIONP(
        const std::string& path_to_ionp_data,
        std::vector<Ionp>& ionpInstances,
        SimMatrix& simMatrix,
        double r_ionp,
        double thickness_ionp_coating = 0.);
    void readCell( const std::string& path_to_cell_data, std::vector<Point3D>& cell_origins, double r_cell );
    void readCell( const std::string& path_to_cell_data, std::vector<double>& cell_loads );
    void readAgg( const std::string& path_to_ionp_data, std::vector<Aggregate>& ionpInstances, SimMatrix& simMatrix, double r_ionp );
    void readHist( const std::string& path_to_ionp_data, std::vector<double>& randomHist );

    // Serialize signals/geometry using the existing filenames and column formats.
    // Proposal means and weighted contributions have different interpretations.
    void writeData( const std::vector<std::complex<double>>& Signal_vec, const std::vector<double>& timeV ) const;
    void writeImportanceRegionData( const std::vector<std::complex<double>>& Signal_vec, const std::vector<double>& timeV, ProtonSampleRegion region ) const;
    void writeImportanceProposalData(
        const std::vector<std::complex<double>>& Signal_vec,
        const std::vector<double>& timeV,
        FreeSamplingProposal proposal,
        bool weighted_contribution) const;
    void writeData( const std::vector<double>& Phase_vec ) const;
    void writeData( const std::vector<Ionp>& IONP_vec ) const;
    void writeData( const std::vector<Aggregate>& Agg_vec ) const;
    void writeData( const std::vector<Proton>& Proton_vec ) const;
    void log( bool reset );

    // Geometry: origins/lengths/radii in metres, node counts dimensionless.
    double get_matrix_x0() const;
    double get_matrix_y0() const;
    double get_matrix_z0() const;
    double get_matrix_xlength() const;
    double get_matrix_ylength() const;
    double get_matrix_zlength() const;
    int get_matrix_nxnodes() const;
    int get_matrix_nynodes() const;
    int get_matrix_nznodes() const;
    double get_radius_ionp() const;
    double get_radius_agg() const;
    int get_ionp_number() const;
    int get_agg_number() const;
    std::optional<std::string> get_ionp_position_file() const;
    std::optional<std::string> get_agg_position_file() const;
    std::optional<std::string> get_cell_position_file() const;
    std::optional<std::string> get_cell_load_file() const;
    int get_proton_number() const;
    // Timing in seconds; sequence 0=FID, 1=SE, 2=MSE.
    double get_echo_time() const;
    double get_time_step() const;
    bool get_variable_time() const;
    bool get_surrounding_voxels() const;
    std::size_t get_surrounding_voxel_superparticles() const;
    double get_echo_spacing() const;
    int get_numThreads() const;
    int get_sequence_num() const;
    double get_velocity_IONP() const;
    Point3D get_direction_IONP() const;
    double get_thickness_IONPcoating() const;
    // has_* reports whether a value was explicitly configured, including zero.
    // Diffusion uses m^2/s; permeability is a probability in [0,1].
    bool has_coating_diffusion_coefficient() const;
    double get_coating_diffusion_coefficient() const;
    bool has_cell_diffusion_coefficient() const;
    double get_cell_diffusion_coefficient() const;
    bool has_coating_permeability() const;
    double get_coating_permeability() const;
    bool has_cell_permeability() const;
    double get_cell_permeability() const;
    int get_bound_proton_number() const;
    double get_background_T2() const;
    double get_cell_T2() const;
    double get_radius_ionp_1() const;
    double get_ionp_per() const;
    double get_radius_ionp_std() const;
    double get_radius_ionp_std_1() const;
    double get_radius_agg_std() const;
    double get_concentration() const;
    double get_concentration_ext() const;
    double get_mag_moment() const;
    Point3D get_Cell_Origin() const;
    Point3D get_Cube_Origin() const;
    Point3D get_Cube_size() const;
    Point3D get_Voxel_position() const;
    double get_Cell_radius() const;
    double get_Cell_concentration() const;
    double get_Cell_concentration_std() const;
    int get_Cell_concentration_dist() const;
    int get_Cell_number() const;
    bool get_Cell_hist() const;
    bool get_Radius_hist() const;
    BFieldSolver get_bfield_solver() const;
    ImportanceSamplingConfig get_importance_sampling() const;

    void set_time_step(double dt);
    // Update the sample count used in output naming after allocations resolve.
    void set_effective_proton_number(int proton_count);

private:
    void writeSignalData(const std::vector<std::complex<double>>& Signal_vec,
        const std::vector<double>& timeV,
        const std::string& region_suffix,
        bool complex_output) const;

    std::string filename;
    double matrix_x0, matrix_y0, matrix_z0; //matrix
    double matrix_xlength, matrix_ylength, matrix_zlength;
    int matrix_nxnodes, matrix_nynodes, matrix_nznodes;
    double radius_ionp, radius_agg; //particles
    int ionp_number, agg_number;
    std::optional<std::string> ionp_position_file;
    std::optional<std::string> agg_position_file;
    std::optional<std::string> cell_position_file;
    std::optional<std::string> cell_load_file;
    int proton_number, logID, numThreads;
    double echo_time, time_step, echo_spacing;
    // Presence flags distinguish omitted options from explicitly supplied values.
    bool variable_time, variable_time_explicit;
    bool surrounding_voxels;
    std::size_t surrounding_voxel_superparticles;
    int sequence_num;
    double velocityIONP;
    Point3D directionIONP;
    double thickness_IONPcoating;
    bool coating_diffusion_coefficient_explicit;
    double coating_diffusion_coefficient;
    bool cell_diffusion_coefficient_explicit;
    double cell_diffusion_coefficient;
    bool coating_permeability_explicit;
    double coating_permeability;
    bool cell_permeability_explicit;
    double cell_permeability;
    int bound_proton_number;
    double background_T2, cell_T2;
    double radius_ionp_1, ionp_per, radius_ionp_std, radius_ionp_std_1, radius_agg_std, concentration, concentration_ext, mag_moment;
    Point3D Cell_Origin, Cube_Origin, Cube_size, Voxel_Pos;
    double Cell_radius, Cell_concentration, Cell_concentration_std;
    int Cell_concentration_dist;
    int Cell_number;
    bool Cell_hist, Radius_hist;
    BFieldSolver bfield_solver;
    ImportanceSamplingConfig importance_sampling;
    int effective_proton_number;
    std::filesystem::path output_tmp_dir;
    std::filesystem::path output_log_dir;
    std::size_t first_tmp_output_index;

};

inline Data::Data(std::string filename_in) :
    filename(std::move(filename_in)),
    matrix_x0(),
    matrix_y0(),
    matrix_z0(),
    matrix_xlength(),
    matrix_ylength(),
    matrix_zlength(),
    matrix_nxnodes(),
    matrix_nynodes(),
    matrix_nznodes(),
    radius_ionp(),
    radius_agg(),
    ionp_number(),
    agg_number(),
    ionp_position_file(),
    agg_position_file(),
    cell_position_file(),
    cell_load_file(),
    proton_number(),
    logID(0),
    numThreads(),
    echo_time(),
    time_step(),
    echo_spacing(),
    variable_time(),
    variable_time_explicit(),
    surrounding_voxels(false),
    surrounding_voxel_superparticles(0),
    sequence_num(0),
    velocityIONP(),
    directionIONP(),
    thickness_IONPcoating(),
    coating_diffusion_coefficient_explicit(false),
    coating_diffusion_coefficient(),
    cell_diffusion_coefficient_explicit(false),
    cell_diffusion_coefficient(),
    coating_permeability_explicit(false),
    coating_permeability(1.0),
    cell_permeability_explicit(false),
    cell_permeability(0.0),
    bound_proton_number(),
    background_T2(),
    cell_T2(),
    radius_ionp_1(),
    ionp_per(0.21),
    radius_ionp_std(),
    radius_ionp_std_1(),
    radius_agg_std(),
    concentration(),
    concentration_ext(),
    mag_moment(),
    Cell_Origin(),
    Cube_Origin(),
    Cube_size(),
    Voxel_Pos(),
    Cell_radius(),
    Cell_concentration(),
    Cell_concentration_std(),
    Cell_concentration_dist(0),
    Cell_number(),
    Cell_hist(),
    Radius_hist(),
    bfield_solver(BFieldSolver::Current),
    importance_sampling(),
    effective_proton_number(),
    output_tmp_dir(),
    output_log_dir(),
    first_tmp_output_index(0) {}

#endif
