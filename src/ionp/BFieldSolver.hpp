/**
 * @file BFieldSolver.hpp
 *
 * Field-solver selection shared by configuration parsing and MRsequence.
 * The enum and its names describe existing runtime choices; they do not
 * change the dipole model or collision geometry.
 */

#ifndef BFIELD_SOLVER_HPP_
#define BFIELD_SOLVER_HPP_

enum class BFieldSolver {
    Current,
    CorrectedHybrid
};

inline const char* bfield_solver_name(BFieldSolver solver) {
    return solver == BFieldSolver::CorrectedHybrid
        ? "corrected_hybrid"
        : "current";
}

#endif
