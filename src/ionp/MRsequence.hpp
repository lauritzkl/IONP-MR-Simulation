/**
 * @file MRsequence.hpp
 *
 * Sequence execution and its collision-query helpers. Lengths are metres,
 * times seconds, phase radians, diffusion coefficients m^2/s, and membrane
 * permeabilities are crossing probabilities. MRsequence stores configuration
 * and a SimMatrix copy; Sequence receives geometry and output buffers by reference.
 */

/*****************************************************************************
*Class to apply different MR sequences

*11-05-2022 (Lauritz Klünder)
******************************************************************************/

#ifndef MRSEQUENCE_HPP_
#define MRSEQUENCE_HPP_

#include <complex>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

#include "ImportanceSampling.hpp"
#include "BFieldSolver.hpp"
#include "mptmacros.h"
#include "SimMatrix.hpp"

class Aggregate;
class Ionp;
class Point3D;
class Proton;
class RandomList;

// Optional path diagnostics: time in seconds, position in metres.
struct TrajectoryEvent {
    double time;
    Point3D position;
    int kind; // 0 = regular sample, 1 = hit point, 2 = post-reflection
};

enum class IonpCheckMode {
    Spatial,
    All
};

struct ImportanceProposalSignalDiagnostics {
    ImportanceProposalSignals means;
    ImportanceProposalSignals contributions;
};

// Conservative collision/near-field candidate lookup. Candidates still need
// exact geometric tests. Scratch-buffer overloads avoid per-query allocation;
// each thread must supply its own scratch vector. Translation offsets let
// moving-source queries reuse an index built at the original coordinates.
class IonpSpatialIndex {
public:
    IonpSpatialIndex() = default;
    explicit IonpSpatialIndex(
        const std::vector<Ionp>& ionpInstances,
        double minimum_cell_size = 0.0);

    void rebuild(
        const std::vector<Ionp>& ionpInstances,
        double minimum_cell_size = 0.0);
    void translate(double delta_x, double delta_y, double delta_z);
    bool empty() const;
    double max_outer_radius() const;
    double max_core_radius() const;
    bool has_coatings() const;
    const std::vector<unsigned int>& all_indices() const;
    std::vector<unsigned int> point_candidates(const Point3D& point, double search_radius) const;
    std::vector<unsigned int> segment_candidates(const Point3D& start, const Point3D& end, double search_radius) const;
    std::vector<unsigned int> variable_step_candidates(const Point3D& point) const;
    void point_candidates(const Point3D& point,
        double search_radius,
        std::vector<unsigned int>& candidates,
        double translation_x = 0.0,
        double translation_y = 0.0,
        double translation_z = 0.0) const;
    void segment_candidates(const Point3D& start,
        const Point3D& end,
        double search_radius,
        std::vector<unsigned int>& candidates,
        double translation_x = 0.0,
        double translation_y = 0.0,
        double translation_z = 0.0,
        bool sort_candidates = true) const;
    void variable_step_candidates(const Point3D& point,
        std::vector<unsigned int>& candidates,
        double translation_x = 0.0,
        double translation_y = 0.0,
        double translation_z = 0.0,
        bool sort_candidates = true) const;

private:
    struct CellKey {
        int x;
        int y;
        int z;

        bool operator==(const CellKey& other) const {
            return x == other.x && y == other.y && z == other.z;
        }
    };

    struct CellKeyHash {
        std::size_t operator()(const CellKey& key) const;
    };

    CellKey key_for_point(const Point3D& point) const;
    void append_candidates_in_box(double min_x,
        double max_x,
        double min_y,
        double max_y,
        double min_z,
        double max_z,
        std::vector<unsigned int>& candidates,
        double translation_x,
        double translation_y,
        double translation_z,
        bool sort_candidates) const;

    double m_cell_size = 1.0;
    double m_max_outer_radius = 0.0;
    double m_max_core_radius = 0.0;
    bool m_has_coatings = false;
    double m_variable_step_cell_size = 1.0;
    double m_translation_x = 0.0;
    double m_translation_y = 0.0;
    double m_translation_z = 0.0;
    std::vector<unsigned int> m_all_indices;
    std::vector<CellKey> m_cell_keys;
    std::unordered_map<CellKey, std::vector<unsigned int>, CellKeyHash> m_cells;
    std::unordered_map<CellKey, std::vector<unsigned int>, CellKeyHash> m_variable_step_cells;
};

// Cell lookup includes periodic image proxies; proxies for one physical cell
// retain the same cell index so membrane crossings preserve its identity.
class CellSpatialIndex {
public:
    struct Proxy {
        Point3D center;
        std::size_t cell_index;
        bool primary_image;
    };

    CellSpatialIndex() = default;
    CellSpatialIndex(const std::vector<Point3D>& cell_origins,
        const Point3D& fallback_origin,
        double cell_radius,
        double voxel_x_length,
        double voxel_y_length,
        double voxel_z_length,
        bool periodic);

    bool empty() const;
    bool has_disjoint_interiors() const;
    double cell_radius() const;
    const Proxy& proxy(std::size_t proxy_index) const;
    const std::vector<std::size_t>& all_proxy_indices() const;
    const std::vector<std::size_t>& proxies_for_cell(std::size_t cell_index) const;
    void segment_candidates(const Point3D& start,
        const Point3D& end,
        double expansion,
        std::vector<std::size_t>& candidates) const;
    void point_candidates(const Point3D& point,
        double expansion,
        std::vector<std::size_t>& candidates) const;

private:
    struct CellKey {
        int x;
        int y;
        int z;

        bool operator==(const CellKey& other) const {
            return x == other.x && y == other.y && z == other.z;
        }
    };

    struct CellKeyHash {
        std::size_t operator()(const CellKey& key) const;
    };

    void append_candidates_in_box(double min_x,
        double max_x,
        double min_y,
        double max_y,
        double min_z,
        double max_z,
        std::vector<std::size_t>& candidates) const;

    double m_cell_radius = 0.0;
    double m_bin_size = 1.0;
    bool m_disjoint_interiors = true;
    std::vector<Proxy> m_proxies;
    std::vector<std::size_t> m_all_proxy_indices;
    std::vector<std::vector<std::size_t>> m_proxy_indices_by_cell;
    std::unordered_map<CellKey, std::vector<std::size_t>, CellKeyHash> m_bins;
};

class MRsequence
{
public:
    MRsequence(
        SimMatrix& simM,
        int sequence_n,
        int prot_n,
        int bound_prot_n,
        double r_i,
        double e_time,
        double time_s,
        double e_spacing,
        int numThread,
        int m_nxnodes,
        double m_xlength,
        Point3D& dirIONP,
        double velIONP,
        double thickness_IONPcoat,
        double coating_diffusion_coeff,
        double cell_diffusion_coeff,
        bool check_coating_perm,
        double coating_perm,
        double back_T2,
        double cell_T2,
        double magM,
        Point3D& Cellorig,
        double Cellrad,
        double cell_perm,
        Point3D& Cubeorig,
        Point3D& Cubesize,
        bool BGrid,
        bool var_Step,
        bool record_positions_flag,
        ImportanceSamplingConfig importance_sampling_cfg,
        IonpCheckMode ionp_check_mode_cfg);

    void replaceProton(
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
        const Point3D* ionp_translation = nullptr);
    void replaceIONP( std::vector<Ionp>& ionpInstances, IonpSpatialIndex* ionp_index = nullptr );
    double computeBFieldAt(
        const Point3D& proton_pos,
        const std::vector<Ionp>& ionpInstances,
        const std::vector<double>& BfieldGrid,
        bool useGrid,
        const Point3D* ionp_translation = nullptr,
        const std::vector<double>* field_amplitudes = nullptr) const;
    void writeData( const std::vector<std::complex<double>>& signalShellProtons, const std::vector<double>& timeV, const std::vector<Ionp>& ionpInstances_orig ) const;
    void set_output_context( std::filesystem::path output_dir );
    void set_bfield_solver( BFieldSolver solver );
    void set_surrounding_voxels( bool enabled );
    void set_surrounding_voxel_superparticles( std::size_t count );
    // Run the configured sequence and fill signal/time outputs. Optional null
    // diagnostic pointers disable those outputs. Proton and IONP vectors are
    // mutable simulation state; geometry used by indices must stay consistent.
    void Sequence(
        std::vector<Proton>& protonInstances,
        std::vector<Ionp>& ionpInstances_orig,
        std::vector<std::complex<double>>& signalVector,
        std::vector<double>& timeV,
        std::vector<std::complex<double>>& signalVectorMSE,
        std::vector<double>& timeVMSE,
        std::vector<Point3D>& cell_origins,
        ImportanceRegionSignals* regionSignalVectors = nullptr,
        ImportanceRegionSignals* regionSignalVectorsMSE = nullptr,
        ImportanceProposalSignalDiagnostics* proposalSignalDiagnostics = nullptr,
        ImportanceProposalSignalDiagnostics* proposalSignalDiagnosticsMSE = nullptr);

private:
    SimMatrix simMatrix;
    int sequence_num, protons_n, bound_protons_n;
    double r_ionp, echo_time, time_step, echo_spacing;
    int numThreads, matrix_nxnodes;
    double matrix_xlength;
    std::filesystem::path output_tmp_dir;
    Point3D directionIONP;
    double vIONP, thickness_IONPcoating, coating_diffusion_coefficient, cell_diffusion_coefficient;
    double outside_diffusion_step, cell_diffusion_step, coating_diffusion_step;
    double outside_diffusion_substep, cell_diffusion_substep, coating_diffusion_substep;
    double coating_permeability, background_T2, incell_T2, magnetic_moment;
    Point3D cell_o, cube_o, cube_s;
    double cell_r, cell_permeability;
    bool Grid, var_Steps;
    bool record_positions;
    bool check_coating_permeability;
    IonpCheckMode ionp_check_mode;
    BFieldSolver bfield_solver;
    bool surrounding_voxels;
    std::size_t surrounding_voxel_superparticles;
    ImportanceSamplingConfig importance_sampling;
    IonpPopulationMap ionp_populations;

    void writeProtonPositions(const std::vector<std::vector<TrajectoryEvent>>& events, std::size_t ionp_count) const;

};

inline MRsequence::MRsequence(
    SimMatrix& simM,
    int sequence_n,
    int prot_n,
    int bound_prot_n,
    double r_i,
    double e_time,
    double time_s,
    double e_spacing,
    int numThread,
    int m_nxnodes,
    double m_xlength,
    Point3D& dirIONP,
    double velIONP,
    double thickness_IONPcoat,
    double coating_diffusion_coeff,
    double cell_diffusion_coeff,
    bool check_coating_perm,
    double coating_perm,
    double back_T2,
    double cell_T2,
    double magM,
    Point3D& Cellorig,
    double Cellrad,
    double cell_perm,
    Point3D& Cubeorig,
    Point3D& Cubesize,
    bool BGrid,
    bool var_Step,
    bool record_positions_flag,
    ImportanceSamplingConfig importance_sampling_cfg,
    IonpCheckMode ionp_check_mode_cfg):
  simMatrix(simM),
  sequence_num(sequence_n),
  protons_n(prot_n),
  bound_protons_n(bound_prot_n),
  r_ionp(r_i),
  echo_time(e_time),
  time_step(time_s),
  echo_spacing(e_spacing),
  numThreads(numThread),
  matrix_nxnodes(m_nxnodes),
  matrix_xlength(m_xlength),
  directionIONP(dirIONP),
  vIONP(velIONP),
  thickness_IONPcoating(thickness_IONPcoat),
  coating_diffusion_coefficient(coating_diffusion_coeff),
  cell_diffusion_coefficient(cell_diffusion_coeff),
  outside_diffusion_step(std::sqrt(6 * D_D_CONST * time_s)),
  cell_diffusion_step(std::sqrt(6 * cell_diffusion_coeff * time_s)),
  coating_diffusion_step(std::sqrt(6 * coating_diffusion_coeff * time_s)),
  outside_diffusion_substep(std::sqrt(6 * D_D_CONST * (time_s / 64.0))),
  cell_diffusion_substep(std::sqrt(6 * cell_diffusion_coeff * (time_s / 64.0))),
  coating_diffusion_substep(std::sqrt(6 * coating_diffusion_coeff * (time_s / 64.0))),
  coating_permeability(coating_perm),
  background_T2(back_T2),
  incell_T2(cell_T2),
  magnetic_moment(magM),
  cell_o(Cellorig),
  cube_o(Cubeorig),
  cube_s(Cubesize),
  cell_r(Cellrad),
  cell_permeability(cell_perm),
  Grid(BGrid),
  var_Steps(var_Step),
  record_positions(record_positions_flag),
  check_coating_permeability(check_coating_perm),
  ionp_check_mode(ionp_check_mode_cfg),
  bfield_solver(BFieldSolver::Current),
  surrounding_voxels(false),
  surrounding_voxel_superparticles(0),
  importance_sampling(importance_sampling_cfg){
  importance_sampling.validate_mode_exclusivity();
  if(!importance_sampling.overlapping_proposals_enabled()) {
    importance_sampling.materialize_legacy_regions();
  }
}

#endif
