/**
 * @file SimMatrix.cpp
 *
 * Dipole B_z evaluation, full-source grid construction and trilinear
 * interpolation. A grid contains every magnetic source; spatial collision
 * indices do not truncate the magnetic field. See SimMatrix.hpp for layout.
 */

/*****************************************************************************
*Class to define simulation matrix (cpp)

*3-9-2021
******************************************************************************/

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <vector>

#include <omp.h>

#include "mptmacros.h"

#include "Ionp.hpp"
#include "Point3D.hpp"
#include "SimMatrix.hpp"

namespace {

inline std::size_t gridIndex(unsigned int x, unsigned int y, unsigned int z, unsigned int nx, unsigned int ny) {
    return static_cast<std::size_t>(z) * nx * ny + static_cast<std::size_t>(x) * ny + y;
}

} // namespace

// B_z = amplitude * (3*cos(theta)^2 - 1) / distance^3, with the moment
// aligned along z. Core volume sets amplitude in the radius-based overload.
double SimMatrix::calculateBfield(const Point3D& p1, const Point3D& p2, const double r) const {
    const double dx = p2.x() - p1.x();
    const double dy = p2.y() - p1.y();
    const double dz = p2.z() - p1.z();

    const double distance_sq = dx * dx + dy * dy + dz * dz;
    const double inv_distance = 1.0 / std::sqrt(distance_sq);
    const double inv_distance_cubed = inv_distance * inv_distance * inv_distance;

    const double cos_theta = dz * inv_distance;
    const double angular_component = (3.0 * square(cos_theta) - 1.0);
    const double radius_cubed = r * r * r;

    return B_fixed * radius_cubed * inv_distance_cubed * angular_component;
}

// Explicit magnetic moment (A*m^2) replaces the radius-derived amplitude;
// the r parameter is retained for the existing interface.
double SimMatrix::calculateBfield(const Point3D& p1, const Point3D& p2, const double r, const double magM) const {
    (void) r;

    const double dx = p2.x() - p1.x();
    const double dy = p2.y() - p1.y();
    const double dz = p2.z() - p1.z();

    const double distance_sq = dx * dx + dy * dy + dz * dz;
    const double inv_distance = 1.0 / std::sqrt(distance_sq);
    const double inv_distance_cubed = inv_distance * inv_distance * inv_distance;

    const double cos_theta = dz * inv_distance;
    const double angular_component = (3.0 * square(cos_theta) - 1.0);
    const double B_fixed_new = D_MO / (4.0 * D_PI) * magM;

    return B_fixed_new * inv_distance_cubed * angular_component;
}

double SimMatrix::calculateBfieldWithAmplitude(
    const Point3D& p1,
    const Point3D& p2,
    double amplitude) const {
    const double dx = p2.x() - p1.x();
    const double dy = p2.y() - p1.y();
    const double dz = p2.z() - p1.z();

    const double distance_sq = dx * dx + dy * dy + dz * dz;
    const double inv_distance = 1.0 / std::sqrt(distance_sq);
    const double inv_distance_cubed =
        inv_distance * inv_distance * inv_distance;
    const double cos_theta = dz * inv_distance;
    const double angular_component = (3.0 * square(cos_theta) - 1.0);

    return amplitude * inv_distance_cubed * angular_component;
}

void SimMatrix::calculateBfield(std::vector<double>& BfieldGrid, const std::vector<Ionp>& ionpInstances, double magnetic_moment) const {

    std::cout << "Starting to calculate the B-Field Grid...." << std::endl;

    if (m_nx_nodes <= 0 || m_ny_nodes <= 0 || m_nz_nodes <= 0) {
        BfieldGrid.clear();
        return;
    }

    const std::size_t grid_size = static_cast<std::size_t>(m_nx_nodes) * static_cast<std::size_t>(m_ny_nodes) * static_cast<std::size_t>(m_nz_nodes);
    BfieldGrid.assign(grid_size, 0.0);

    const double start_x = -m_x_length / 2.0;
    const double start_y = -m_y_length / 2.0;
    const double start_z = -m_z_length / 2.0;

    const double stepsize_x = (m_nx_nodes > 1) ? m_x_length / static_cast<double>(m_nx_nodes - 1) : 0.0;
    const double stepsize_y = (m_ny_nodes > 1) ? m_y_length / static_cast<double>(m_ny_nodes - 1) : 0.0;
    const double stepsize_z = (m_nz_nodes > 1) ? m_z_length / static_cast<double>(m_nz_nodes - 1) : 0.0;

    const bool use_magnetic_moment = (magnetic_moment != 0.0);
    if (ionpInstances.empty()) {
        std::fill(BfieldGrid.begin(), BfieldGrid.end(), 0.0);
        return;
    }

    // The full all-IONP grid is the reference grid calculation. Nearby
    // interpolation error is corrected later by the optional hybrid evaluator.
    #pragma omp parallel for
    for (unsigned int z = 0; z < SC_UI(m_nz_nodes); ++z) {
        for (unsigned int x = 0; x < SC_UI(m_nx_nodes); ++x) {
            for (unsigned int y = 0; y < SC_UI(m_ny_nodes); ++y) {
                Point3D gridPoint;
                gridPoint.set(start_x + stepsize_x * x,
                              start_y + stepsize_y * y,
                              start_z + stepsize_z * z);

                double bfield_tot = 0.0;
                if (use_magnetic_moment) {
                    for (unsigned int index = 0;
                         index < SC_UI(ionpInstances.size());
                         ++index) {
                        bfield_tot += calculateBfield(
                            gridPoint,
                            ionpInstances[index].xt(),
                            ionpInstances[index].radius(),
                            magnetic_moment);
                    }
                }
                else {
                    for (unsigned int index = 0;
                         index < SC_UI(ionpInstances.size());
                         ++index) {
                        bfield_tot += calculateBfield(
                            gridPoint,
                            ionpInstances[index].xt(),
                            ionpInstances[index].radius());
                    }
                }

                BfieldGrid[gridIndex(
                    x, y, z, SC_UI(m_nx_nodes), SC_UI(m_ny_nodes))] =
                    bfield_tot;
            }
        }
    }

    if (m_save_grid && !m_save_grid_path.empty()) {
        std::error_code ec;
        std::filesystem::create_directories(m_save_grid_path.parent_path(), ec);
        std::ofstream out(m_save_grid_path, std::ios::binary);
        if (!out.is_open()) {
            std::cerr << "Warning: unable to write B-field grid to " << m_save_grid_path << std::endl;
        } else {
            // Native binary layout: three uint32 dimensions, then contiguous
            // doubles in z/x/y order. No coordinates, units or endian marker.
            const std::uint32_t nx = static_cast<std::uint32_t>(m_nx_nodes);
            const std::uint32_t ny = static_cast<std::uint32_t>(m_ny_nodes);
            const std::uint32_t nz = static_cast<std::uint32_t>(m_nz_nodes);
            out.write(reinterpret_cast<const char*>(&nx), sizeof(nx));
            out.write(reinterpret_cast<const char*>(&ny), sizeof(ny));
            out.write(reinterpret_cast<const char*>(&nz), sizeof(nz));
            out.write(reinterpret_cast<const char*>(BfieldGrid.data()), static_cast<std::streamsize>(BfieldGrid.size() * sizeof(double)));
        }
    }

    std::cout << "B-Field Grid calculation finished...." << std::endl;
}

// Cache the eight corners and tensor-product linear weights. At an upper
// grid edge the clamped corner indices may coincide; their weights handle it.
GridInterpolationStencil SimMatrix::interpolationStencil(const Point3D& p1) const {

    const double step_x = m_grid_step_x;
    const double step_y = m_grid_step_y;
    const double step_z = m_grid_step_z;

    const double origin_x = m_grid_origin_x;
    const double origin_y = m_grid_origin_y;
    const double origin_z = m_grid_origin_z;

    const auto clamp_index = [](double value, unsigned int max_index) {
        const double clamped = std::max(0.0, std::min(value, static_cast<double>(max_index)));
        return static_cast<unsigned int>(clamped);
    };

    const unsigned int max_x_index = m_grid_max_x_index;
    const unsigned int max_y_index = m_grid_max_y_index;
    const unsigned int max_z_index = m_grid_max_z_index;

    const double scaled_x = (step_x != 0.0) ? (p1.x() - origin_x) / step_x : 0.0;
    const double scaled_y = (step_y != 0.0) ? (p1.y() - origin_y) / step_y : 0.0;
    const double scaled_z = (step_z != 0.0) ? (p1.z() - origin_z) / step_z : 0.0;

    const unsigned int x_indexDOWN = clamp_index(std::floor(scaled_x), max_x_index);
    const unsigned int y_indexDOWN = clamp_index(std::floor(scaled_y), max_y_index);
    const unsigned int z_indexDOWN = clamp_index(std::floor(scaled_z), max_z_index);

    const unsigned int x_indexUP = std::min(x_indexDOWN + 1, max_x_index);
    const unsigned int y_indexUP = std::min(y_indexDOWN + 1, max_y_index);
    const unsigned int z_indexUP = std::min(z_indexDOWN + 1, max_z_index);

    const double xd = (x_indexUP == x_indexDOWN || step_x == 0.0) ? 0.0 : (p1.x() - (origin_x + step_x * x_indexDOWN)) / step_x;
    const double yd = (y_indexUP == y_indexDOWN || step_y == 0.0) ? 0.0 : (p1.y() - (origin_y + step_y * y_indexDOWN)) / step_y;
    const double zd = (z_indexUP == z_indexDOWN || step_z == 0.0) ? 0.0 : (p1.z() - (origin_z + step_z * z_indexDOWN)) / step_z;

    GridInterpolationStencil stencil;
    const std::array<unsigned int, 8> xs = {
        x_indexDOWN, x_indexUP, x_indexDOWN, x_indexUP,
        x_indexDOWN, x_indexUP, x_indexDOWN, x_indexUP};
    const std::array<unsigned int, 8> ys = {
        y_indexDOWN, y_indexDOWN, y_indexUP, y_indexUP,
        y_indexDOWN, y_indexDOWN, y_indexUP, y_indexUP};
    const std::array<unsigned int, 8> zs = {
        z_indexDOWN, z_indexDOWN, z_indexDOWN, z_indexDOWN,
        z_indexUP, z_indexUP, z_indexUP, z_indexUP};
    stencil.weights = {
        (1. - xd) * (1. - yd) * (1. - zd),
        xd * (1. - yd) * (1. - zd),
        (1. - xd) * yd * (1. - zd),
        xd * yd * (1. - zd),
        (1. - xd) * (1. - yd) * zd,
        xd * (1. - yd) * zd,
        (1. - xd) * yd * zd,
        xd * yd * zd};
    for (std::size_t corner = 0; corner < stencil.indices.size(); ++corner) {
        stencil.indices[corner] = gridIndex(
            xs[corner], ys[corner], zs[corner],
            SC_UI(m_nx_nodes), SC_UI(m_ny_nodes));
        stencil.points[corner].set(
            origin_x + step_x * xs[corner],
            origin_y + step_y * ys[corner],
            origin_z + step_z * zs[corner]);
    }
    return stencil;
}

double SimMatrix::interpolate(
    const Point3D& p1,
    const std::vector<double>& BfieldGrid) const {
    return interpolate(interpolationStencil(p1), BfieldGrid);
}

double SimMatrix::interpolate(
    const GridInterpolationStencil& stencil,
    const std::vector<double>& BfieldGrid) const {
    double result = 0.0;
    for (std::size_t corner = 0; corner < stencil.indices.size(); ++corner) {
        result += stencil.weights[corner] * BfieldGrid[stencil.indices[corner]];
    }
    return result;
}

double SimMatrix::interpolateBfieldContribution(
    const Point3D& p1,
    const Point3D& ionp_position,
    double amplitude) const {
    return interpolateBfieldContribution(
        interpolationStencil(p1), ionp_position, amplitude);
}

double SimMatrix::interpolateBfieldContribution(
    const GridInterpolationStencil& stencil,
    const Point3D& ionp_position,
    double amplitude) const {
    double result = 0.0;
    for (std::size_t corner = 0; corner < stencil.points.size(); ++corner) {
        if (stencil.weights[corner] == 0.0) {
            continue;
        }
        result += stencil.weights[corner] * calculateBfieldWithAmplitude(
            stencil.points[corner], ionp_position, amplitude);
    }
    return result;
}

double SimMatrix::maxGridSpacing() const {
    return std::max({m_grid_step_x, m_grid_step_y, m_grid_step_z});
}

void SimMatrix::set_save_grid_path(const std::optional<std::filesystem::path>& path) {
    m_save_grid = path.has_value() && !path->empty();
    if (m_save_grid) {
        m_save_grid_path = *path;
    }
}
