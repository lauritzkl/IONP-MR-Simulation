/**
 * @file SimSpace.hpp
 *
 * Build physical IONP, aggregate and cell geometry from loaded or synthetic
 * inputs. Vectors/references are caller-owned inputs and outputs; distances
 * use metres. Surrounding-voxel helpers create magnetic sources separately
 * from the central voxel's physical collision geometry.
 */

/*****************************************************************************
*Class to place IONP in SimVolume in different spacial distributions

*17-05-2022 (Lauritz Klünder)
******************************************************************************/

#ifndef SIMSPACE_HPP_
#define SIMSPACE_HPP_

#include <cstddef>
#include <filesystem>
#include <optional>
#include <vector>

class Aggregate;
class Ionp;
class Point3D;
class SimMatrix;

namespace SimSpace {

    // Historical placement is the default, matching the executable. Callers
    // can explicitly pass true to use the overlap-grid implementation.
    void simulationSpace_inCell(
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
        bool use_overlap_grid = false,
        std::optional<std::filesystem::path> cell_position_path = std::nullopt,
        std::optional<std::filesystem::path> cell_load_path = std::nullopt);
    void simulationSpace_outCell(
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
        bool use_overlap_grid = false);
    void simulationSpace_outCellAggregates(
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
        const std::vector<Point3D>& c_cell_vec);
    std::vector<Ionp> surroundingSpace( const std::vector<Ionp>& ionpInstances, const SimMatrix& simMatrix );
    std::vector<Ionp> surroundingSpace( const std::vector<Ionp>& ionpInstances, const SimMatrix& simMatrix, std::size_t superparticle_cap, const std::vector<Point3D>& cell_origins, const Point3D& fallback_cell_origin, double cell_radius );
    void simulationSpace(
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
        std::optional<std::filesystem::path> cell_position_path = std::nullopt);
    void simulationSpace(
        std::vector<Aggregate>& aggregateInstances,
        std::vector<Ionp>& ionpInstances,
        SimMatrix& simMatrix,
        int agg_n,
        int ionp_n,
        double r_ionp);

}

#endif
