/**
 * @file SimMatrix.hpp
 *
 * Voxel dimensions and scalar B_z field-grid operations. Positions/spacing
 * use metres and field values use tesla. Grid storage is supplied by callers.
 * Grid coordinates are centred on zero even when x0 stores another origin.
 */

/*****************************************************************************
*Class to define the simulation matrix (header .hpp)

*3-9-2021
******************************************************************************/

#ifndef SIMMATRIX_HPP_
#define SIMMATRIX_HPP_

#include <array>
#include <cstddef>
#include <vector>
#include <filesystem>
#include <optional>
#include "Point3D.hpp"

class Ionp;

// Eight flattened indices, positions and trilinear weights, reusable for the
// full grid and individual-source contributions at the same query point.
struct GridInterpolationStencil {
    std::array<std::size_t, 8> indices{};
    std::array<Point3D, 8> points{};
    std::array<double, 8> weights{};
};

class SimMatrix {
public:

    SimMatrix(const Point3D& p, double x_l, double y_l, double z_l, int nx, int ny, int nz);
    SimMatrix(double x0, double y0, double z0, double x_l, double y_l, double z_l, int nx, int ny, int nz);
    SimMatrix(double x_l, double y_l, double z_l, int nx, int ny, int nz);
    SimMatrix(double l, int n);

    SimMatrix(const SimMatrix& m); //copy
    ~SimMatrix();

    inline const Point3D& x0() const;
    inline double xLength() const;
    inline double yLength() const;
    inline double zLength() const;
    inline int nxNodes() const;
    inline int nyNodes() const;
    inline int nzNodes() const;

    inline void addBValue(double b);
    inline void addPosValue(const Point3D& p);

    // p1 is the query, p2 the source centre. Callers must avoid coincidence
    // with the dipole centre. r is core radius; magM overrides the moment.
    double calculateBfield(const Point3D& p1, const Point3D& p2, const double r) const; //calculate B between 2 points
    double calculateBfield(const Point3D& p1, const Point3D& p2, const double r, const double magM) const; //calculate B between 2 points
    double calculateBfieldWithAmplitude(const Point3D& p1, const Point3D& p2, double amplitude) const;

    // Fill the caller's grid in z*nx*ny + x*ny + y order (y varies fastest).
    void calculateBfield(std::vector<double>& BfieldGrid, const std::vector<Ionp>& ionpInstances, double magnetic_moment) const;

    // Clamp grid indices at the edges. Periodic wrapping is the caller's job.
    double interpolate(const Point3D& p1, const std::vector<double>& BfieldGrid) const;
    double interpolate(const GridInterpolationStencil& stencil, const std::vector<double>& BfieldGrid) const;
    GridInterpolationStencil interpolationStencil(const Point3D& p1) const;
    double interpolateBfieldContribution(const Point3D& p1, const Point3D& ionp_position, double amplitude) const;
    double interpolateBfieldContribution(const GridInterpolationStencil& stencil, const Point3D& ionp_position, double amplitude) const;
    double maxGridSpacing() const;

    void set_save_grid_path(const std::optional<std::filesystem::path>& path);

private:
    Point3D m_x0; //origin
    const double m_x_length; //fixed
    const double m_y_length;
    const double m_z_length;
    const int m_nx_nodes;
    const int m_ny_nodes;
    const int m_nz_nodes;
    const double m_grid_step_x;
    const double m_grid_step_y;
    const double m_grid_step_z;
    const double m_grid_origin_x;
    const double m_grid_origin_y;
    const double m_grid_origin_z;
    const unsigned int m_grid_max_x_index;
    const unsigned int m_grid_max_y_index;
    const unsigned int m_grid_max_z_index;

public:
    // Legacy auxiliary containers; the active grid is passed explicitly.
    // The copy constructor copies matrix metadata, not these containers.
    std::vector<double> b_field;
    std::vector<int>::iterator b_field_iter;
    std::vector<Point3D> pos_nodes;
    std::vector<int>::iterator Iter;

private:
    bool m_save_grid = false;
    std::filesystem::path m_save_grid_path;
};

//inline definitions
inline SimMatrix::SimMatrix(const Point3D& p, double x_l, double y_l, double z_l, int nx, int ny, int nz) // most general(point)
    : m_x0(p), m_x_length(x_l), m_y_length(y_l), m_z_length(z_l), m_nx_nodes(nx), m_ny_nodes(ny), m_nz_nodes(nz),
      m_grid_step_x(nx > 1 ? x_l / static_cast<double>(nx - 1) : 0.0),
      m_grid_step_y(ny > 1 ? y_l / static_cast<double>(ny - 1) : 0.0),
      m_grid_step_z(nz > 1 ? z_l / static_cast<double>(nz - 1) : 0.0),
      m_grid_origin_x(-x_l / 2.0), m_grid_origin_y(-y_l / 2.0), m_grid_origin_z(-z_l / 2.0),
      m_grid_max_x_index(static_cast<unsigned int>(nx > 0 ? nx - 1 : 0)),
      m_grid_max_y_index(static_cast<unsigned int>(ny > 0 ? ny - 1 : 0)),
      m_grid_max_z_index(static_cast<unsigned int>(nz > 0 ? nz - 1 : 0)) {}
inline SimMatrix::SimMatrix(
    double x0,
    double y0,
    double z0,
    double x_l,
    double y_l,
    double z_l,
    int nx,
    int ny,
    int nz) // most general(exactly the numbers)
    : SimMatrix(Point3D(x0, y0, z0), x_l, y_l, z_l, nx, ny, nz) {}
inline SimMatrix::SimMatrix(double x_l, double y_l, double z_l, int nx, int ny, int nz) //origin(0,0,0), different lengths, nNodes
    : SimMatrix(Point3D(0., 0., 0.), x_l, y_l, z_l, nx, ny, nz) {}
inline SimMatrix::SimMatrix(double l, int n) // origin(0,0,0), same length, nNodes
    : SimMatrix(Point3D(0., 0., 0.), l, l, l, n, n, n) {}

inline SimMatrix::SimMatrix(const SimMatrix& m)
    : m_x0(m.m_x0), m_x_length(m.m_x_length), m_y_length(m.m_y_length), m_z_length(m.m_z_length), m_nx_nodes(m.m_nx_nodes), m_ny_nodes(m.m_ny_nodes), m_nz_nodes(m.m_nz_nodes),
      m_grid_step_x(m.m_grid_step_x), m_grid_step_y(m.m_grid_step_y), m_grid_step_z(m.m_grid_step_z),
      m_grid_origin_x(m.m_grid_origin_x), m_grid_origin_y(m.m_grid_origin_y), m_grid_origin_z(m.m_grid_origin_z),
      m_grid_max_x_index(m.m_grid_max_x_index), m_grid_max_y_index(m.m_grid_max_y_index), m_grid_max_z_index(m.m_grid_max_z_index),
      m_save_grid(m.m_save_grid), m_save_grid_path(m.m_save_grid_path) {}
inline SimMatrix::~SimMatrix() {}

inline const Point3D& SimMatrix::x0() const {
    return m_x0;
}
inline double SimMatrix::xLength() const {
    return m_x_length;
}
inline double SimMatrix::yLength() const {
    return m_y_length;
}
inline double SimMatrix::zLength() const {
    return m_z_length;
}
inline int SimMatrix::nxNodes() const {
    return m_nx_nodes;
}
inline int SimMatrix::nyNodes() const {
    return m_ny_nodes;
}
inline int SimMatrix::nzNodes() const {
    return m_nz_nodes;
}

inline void SimMatrix::addBValue(double b) {
   b_field.push_back(b);
}
inline void SimMatrix::addPosValue(const Point3D& p) {
   pos_nodes.push_back(p);
}

#endif
