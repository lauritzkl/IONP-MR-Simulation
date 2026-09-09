/**
 * @file CorrectedHybridField.hpp
 *
 * Stationary-source field evaluator: full-grid interpolation plus a smooth
 * near-source correction. Matrix and grid references must outlive this object
 * and remain unchanged during queries. Source positions/amplitudes are copied;
 * rebuild the evaluator when the underlying geometry changes.
 */

#ifndef CORRECTED_HYBRID_FIELD_HPP_
#define CORRECTED_HYBRID_FIELD_HPP_

#include <cstddef>
#include <vector>

#include "Point3D.hpp"

class Ionp;
class SimMatrix;

class CorrectedHybridField {
public:
    CorrectedHybridField(const SimMatrix& matrix,
                         const std::vector<Ionp>& ionps,
                         const std::vector<double>& amplitudes,
                         const std::vector<double>& grid);

    // Read-only query; concurrent calls may share this immutable evaluator.
    // Returns scalar B_z in tesla; a dipole centre is a singular query point.
    double evaluate(const Point3D& position) const;
    double correctionStart(std::size_t ionp_index) const;
    double correctionEnd(std::size_t ionp_index) const;
    double maximumCorrectionEnd() const;

private:
    struct CellKey {
        int x;
        int y;
        int z;

        bool operator==(const CellKey& other) const {
            return x == other.x && y == other.y && z == other.z;
        }
    };

    struct CellEntry {
        CellKey key;
        unsigned int particle_index;
    };

    struct CellRange {
        CellKey key;
        std::size_t begin;
        std::size_t end;
    };

    static bool cellKeyLess(const CellKey& lhs, const CellKey& rhs);
    static std::size_t cellKeyHash(const CellKey& key);
    CellKey keyForPosition(double x, double y, double z) const;
    const CellRange* findCell(const CellKey& key) const;

    // Borrowed grid geometry/data; the remaining vectors own source/index data.
    const SimMatrix& m_matrix;
    const std::vector<double>& m_grid;
    std::vector<double> m_x;
    std::vector<double> m_y;
    std::vector<double> m_z;
    std::vector<double> m_amplitude;
    std::vector<double> m_correction_start;
    std::vector<double> m_correction_end;
    double m_maximum_correction_end = 0.0;
    double m_cell_size = 1.0;
    std::vector<CellRange> m_cell_ranges;
    std::vector<unsigned int> m_cell_particles;
    std::vector<std::size_t> m_cell_hash_slots;
};

#endif
