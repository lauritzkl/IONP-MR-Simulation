/**
 * @file CorrectedHybridField.cpp
 *
 * Build a spatial lookup for local dipole corrections, then reuse one
 * eight-corner grid stencil per query. Subtracting each source's interpolated
 * contribution before adding its exact field avoids counting that source twice.
 */

#include "CorrectedHybridField.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "Ionp.hpp"
#include "SimMatrix.hpp"

namespace {

constexpr double kFieldRangeFactor = 8.0;

// Cubic smoothstep complement: full correction inside start, zero outside
// end, with zero slope at both ends of the transition band.
double correction_weight(double distance, double start, double end) {
    if (distance <= start) {
        return 1.0;
    }
    if (distance >= end || end <= start) {
        return 0.0;
    }
    const double t = (distance - start) / (end - start);
    return 1.0 - 3.0 * t * t + 2.0 * t * t * t;
}

inline double field_from_displacement(
    double dx,
    double dy,
    double dz,
    double amplitude) {
    const double distance_squared = dx * dx + dy * dy + dz * dz;
    const double inverse_distance = 1.0 / std::sqrt(distance_squared);
    const double inverse_distance_cubed =
        inverse_distance * inverse_distance * inverse_distance;
    const double cos_theta = dz * inverse_distance;
    return amplitude * inverse_distance_cubed *
        (3.0 * cos_theta * cos_theta - 1.0);
}

} // namespace

// Sort occupied bins lexicographically so each bin owns one contiguous range.
bool CorrectedHybridField::cellKeyLess(
    const CellKey& lhs,
    const CellKey& rhs) {
    if (lhs.x != rhs.x) {
        return lhs.x < rhs.x;
    }
    if (lhs.y != rhs.y) {
        return lhs.y < rhs.y;
    }
    return lhs.z < rhs.z;
}

std::size_t CorrectedHybridField::cellKeyHash(const CellKey& key) {
    std::size_t hash = static_cast<std::size_t>(key.x) * 73856093u;
    hash ^= static_cast<std::size_t>(key.y) * 19349663u;
    hash ^= static_cast<std::size_t>(key.z) * 83492791u;
    return hash;
}

CorrectedHybridField::CorrectedHybridField(
    const SimMatrix& matrix,
    const std::vector<Ionp>& ionps,
    const std::vector<double>& amplitudes,
    const std::vector<double>& grid)
    : m_matrix(matrix), m_grid(grid) {
    if (amplitudes.size() != ionps.size()) {
        throw std::invalid_argument(
            "Corrected hybrid field requires one amplitude per IONP");
    }
    if (grid.empty()) {
        throw std::invalid_argument(
            "Corrected hybrid field requires a non-empty B-field grid");
    }
    if (matrix.nxNodes() <= 1 ||
        matrix.nyNodes() <= 1 ||
        matrix.nzNodes() <= 1) {
        throw std::invalid_argument(
            "Corrected hybrid field requires at least two grid nodes per axis");
    }
    const std::size_t expected_grid_size =
        static_cast<std::size_t>(matrix.nxNodes()) *
        static_cast<std::size_t>(matrix.nyNodes()) *
        static_cast<std::size_t>(matrix.nzNodes());
    if (grid.size() != expected_grid_size) {
        throw std::invalid_argument(
            "Corrected hybrid field grid size does not match its matrix");
    }
    const double grid_spacing = matrix.maxGridSpacing();
    if (grid_spacing <= 0.0) {
        throw std::invalid_argument(
            "Corrected hybrid field requires at least two grid nodes per axis");
    }

    m_x.reserve(ionps.size());
    m_y.reserve(ionps.size());
    m_z.reserve(ionps.size());
    m_amplitude = amplitudes;
    m_correction_start.reserve(ionps.size());
    m_correction_end.reserve(ionps.size());

    for (const Ionp& ionp : ionps) {
        m_x.push_back(ionp.xt().x());
        m_y.push_back(ionp.xt().y());
        m_z.push_back(ionp.xt().z());
        const double start = std::max(
            kFieldRangeFactor * ionp.outer_radius(),
            3.0 * grid_spacing);
        const double end = start + grid_spacing;
        m_correction_start.push_back(start);
        m_correction_end.push_back(end);
        m_maximum_correction_end =
            std::max(m_maximum_correction_end, end);
    }

    if (m_maximum_correction_end > 0.0) {
        m_cell_size = m_maximum_correction_end;
    }
    std::vector<CellEntry> entries;
    entries.reserve(ionps.size());
    for (unsigned int index = 0;
         index < static_cast<unsigned int>(ionps.size());
         ++index) {
        entries.push_back({
            keyForPosition(m_x[index], m_y[index], m_z[index]),
            index});
    }
    std::sort(
        entries.begin(),
        entries.end(),
        [](const CellEntry& lhs, const CellEntry& rhs) {
            return cellKeyLess(lhs.key, rhs.key);
        });
    // Keep particles from each occupied cell contiguous. The open-addressed
    // table below indexes only occupied cells and avoids per-query allocation.
    m_cell_particles.reserve(entries.size());
    for (const CellEntry& entry : entries) {
        if (m_cell_ranges.empty() ||
            !(m_cell_ranges.back().key == entry.key)) {
            m_cell_ranges.push_back({
                entry.key,
                m_cell_particles.size(),
                m_cell_particles.size()});
        }
        m_cell_particles.push_back(entry.particle_index);
        m_cell_ranges.back().end = m_cell_particles.size();
    }
    if (!m_cell_ranges.empty()) {
        std::size_t slot_count = 1;
        while (slot_count < 2 * m_cell_ranges.size()) {
            slot_count *= 2;
        }
        m_cell_hash_slots.assign(
            slot_count, std::numeric_limits<std::size_t>::max());
        const std::size_t mask = slot_count - 1;
        for (std::size_t range_index = 0;
             range_index < m_cell_ranges.size();
             ++range_index) {
            std::size_t slot =
                cellKeyHash(m_cell_ranges[range_index].key) & mask;
            while (m_cell_hash_slots[slot] !=
                   std::numeric_limits<std::size_t>::max()) {
                slot = (slot + 1) & mask;
            }
            m_cell_hash_slots[slot] = range_index;
        }
    }
}

CorrectedHybridField::CellKey CorrectedHybridField::keyForPosition(
    double x,
    double y,
    double z) const {
    return {
        static_cast<int>(std::floor(x / m_cell_size)),
        static_cast<int>(std::floor(y / m_cell_size)),
        static_cast<int>(std::floor(z / m_cell_size))};
}

const CorrectedHybridField::CellRange* CorrectedHybridField::findCell(
    const CellKey& key) const {
    if (m_cell_hash_slots.empty()) {
        return nullptr;
    }
    const std::size_t empty = std::numeric_limits<std::size_t>::max();
    const std::size_t mask = m_cell_hash_slots.size() - 1;
    std::size_t slot = cellKeyHash(key) & mask;
    while (m_cell_hash_slots[slot] != empty) {
        const CellRange& range = m_cell_ranges[m_cell_hash_slots[slot]];
        if (range.key == key) {
            return &range;
        }
        slot = (slot + 1) & mask;
    }
    return nullptr;
}

// Use the full grid as the baseline and correct only nearby sources. The
// correction is weight * (exact source field - its grid interpolation).
double CorrectedHybridField::evaluate(const Point3D& position) const {
    const GridInterpolationStencil stencil =
        m_matrix.interpolationStencil(position);
    double field = m_matrix.interpolate(stencil, m_grid);
    if (m_cell_ranges.empty() || m_maximum_correction_end <= 0.0) {
        return field;
    }
    const CellKey minimum = keyForPosition(
        position.x() - m_maximum_correction_end,
        position.y() - m_maximum_correction_end,
        position.z() - m_maximum_correction_end);
    const CellKey maximum = keyForPosition(
        position.x() + m_maximum_correction_end,
        position.y() + m_maximum_correction_end,
        position.z() + m_maximum_correction_end);
    const auto distance_to_cell_axis = [this](double coordinate, int cell) {
        const double lower = static_cast<double>(cell) * m_cell_size;
        const double upper = lower + m_cell_size;
        if (coordinate < lower) {
            return lower - coordinate;
        }
        if (coordinate > upper) {
            return coordinate - upper;
        }
        return 0.0;
    };
    const double maximum_distance_squared =
        m_maximum_correction_end * m_maximum_correction_end;
    for (int z = minimum.z; z <= maximum.z; ++z) {
        const double cell_dz = distance_to_cell_axis(position.z(), z);
        for (int x = minimum.x; x <= maximum.x; ++x) {
            const double cell_dx = distance_to_cell_axis(position.x(), x);
            for (int y = minimum.y; y <= maximum.y; ++y) {
                const double cell_dy = distance_to_cell_axis(position.y(), y);
                if (cell_dx * cell_dx + cell_dy * cell_dy +
                    cell_dz * cell_dz >= maximum_distance_squared) {
                    continue;
                }
                const CellRange* range = findCell(CellKey{x, y, z});
                if (range == nullptr) {
                    continue;
                }
                for (std::size_t particle = range->begin;
                     particle < range->end;
                     ++particle) {
                    const unsigned int index = m_cell_particles[particle];
                    const double dx = m_x[index] - position.x();
                    const double dy = m_y[index] - position.y();
                    const double dz = m_z[index] - position.z();
                    const double distance_squared =
                        dx * dx + dy * dy + dz * dz;
                    const double end = m_correction_end[index];
                    if (distance_squared >= end * end) {
                        continue;
                    }
                    const double distance = std::sqrt(distance_squared);
                    const double weight = correction_weight(
                        distance, m_correction_start[index], end);
                    if (weight == 0.0) {
                        continue;
                    }
                    const double exact = field_from_displacement(
                        dx, dy, dz, m_amplitude[index]);
                    double interpolated = 0.0;
                    for (std::size_t corner = 0;
                         corner < stencil.points.size();
                         ++corner) {
                        if (stencil.weights[corner] == 0.0) {
                            continue;
                        }
                        interpolated += stencil.weights[corner] *
                            field_from_displacement(
                                m_x[index] - stencil.points[corner].x(),
                                m_y[index] - stencil.points[corner].y(),
                                m_z[index] - stencil.points[corner].z(),
                                m_amplitude[index]);
                    }
                    field += weight * (exact - interpolated);
                }
            }
        }
    }
    return field;
}

double CorrectedHybridField::correctionStart(std::size_t ionp_index) const {
    return m_correction_start.at(ionp_index);
}

double CorrectedHybridField::correctionEnd(std::size_t ionp_index) const {
    return m_correction_end.at(ionp_index);
}

double CorrectedHybridField::maximumCorrectionEnd() const {
    return m_maximum_correction_end;
}
