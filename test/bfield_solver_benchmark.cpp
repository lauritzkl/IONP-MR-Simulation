#include <chrono>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

#include "CorrectedHybridField.hpp"
#include "Ionp.hpp"
#include "SimMatrix.hpp"
#include "mptmacros.h"

namespace {

using Clock = std::chrono::steady_clock;

constexpr std::size_t kParticleCount = 75000;
constexpr std::size_t kQueryCount = 256;
constexpr double kMinimumSpeedup = 100.0;

double fractional_part(double value) {
	return value - std::floor(value);
}

} // namespace

int main() {
	const double length = 171e-6;
	const int nodes = 250;
	const double core_radius = 12.5e-9;
	const double coating_thickness = 37.5e-9;
	SimMatrix matrix(length, nodes);

	std::vector<Ionp> ionps;
	std::vector<double> amplitudes;
	ionps.reserve(kParticleCount);
	amplitudes.reserve(kParticleCount);
	const double amplitude =
		B_fixed * core_radius * core_radius * core_radius;
	for (std::size_t index = 0; index < kParticleCount; ++index) {
		const double x = length *
			(fractional_part(0.5 + 0.7548776662466927 * index) - 0.5);
		const double y = length *
			(fractional_part(0.25 + 0.5698402909980532 * index) - 0.5);
		const double z = length *
			(fractional_part(0.75 + 0.4385790213902124 * index) - 0.5);
		ionps.emplace_back(
			x, y, z, 0.0, core_radius, coating_thickness);
		amplitudes.push_back(amplitude);
	}

	// Grid construction is intentionally excluded from this benchmark. Its
	// size and node spacing are realistic, while values do not affect query cost.
	const std::size_t grid_size =
		static_cast<std::size_t>(nodes) * nodes * nodes;
	const std::vector<double> grid(grid_size, 0.0);
	const CorrectedHybridField hybrid(matrix, ionps, amplitudes, grid);

	std::vector<Point3D> queries;
	queries.reserve(kQueryCount);
	for (const Ionp& ionp : ionps) {
		if (queries.size() == kQueryCount) {
			break;
		}
		if (std::abs(ionp.xt().x()) >= 0.4 * length ||
			std::abs(ionp.xt().y()) >= 0.4 * length ||
			std::abs(ionp.xt().z()) >= 0.4 * length) {
			continue;
		}
		// Force every timed query into the near field of a particle without
		// evaluating inside its magnetic core.
		queries.emplace_back(
			ionp.xt().x() + 4.0 * ionp.outer_radius(),
			ionp.xt().y(),
			ionp.xt().z());
	}
	if (queries.size() != kQueryCount) {
		std::cerr << "Unable to construct the requested near-field queries"
			<< std::endl;
		return 1;
	}
	std::size_t corrected_source_count = 0;
	const double correction_end = hybrid.maximumCorrectionEnd();
	for (const Point3D& query : queries) {
		for (const Ionp& ionp : ionps) {
			const double dx = ionp.xt().x() - query.x();
			const double dy = ionp.xt().y() - query.y();
			const double dz = ionp.xt().z() - query.z();
			if (dx * dx + dy * dy + dz * dz <
				correction_end * correction_end) {
				++corrected_source_count;
			}
		}
	}

	volatile double sink = hybrid.evaluate(queries.front());
	const auto direct_start = Clock::now();
	for (const Point3D& query : queries) {
		double field = 0.0;
		for (std::size_t index = 0; index < ionps.size(); ++index) {
			field += matrix.calculateBfieldWithAmplitude(
				query, ionps[index].xt(), amplitudes[index]);
		}
		sink = sink + field;
	}
	const auto direct_end = Clock::now();

	const auto hybrid_start = Clock::now();
	for (const Point3D& query : queries) {
		sink = sink + hybrid.evaluate(query);
	}
	const auto hybrid_end = Clock::now();

	const double direct_seconds =
		std::chrono::duration<double>(direct_end - direct_start).count();
	const double hybrid_seconds =
		std::chrono::duration<double>(hybrid_end - hybrid_start).count();
	const double speedup = direct_seconds / hybrid_seconds;
	std::cout << "Particles: " << kParticleCount << "\n";
	std::cout << "Grid nodes per axis: " << nodes << "\n";
	std::cout << "Maximum grid spacing: " << matrix.maxGridSpacing() << " m\n";
	std::cout << "Direct query time: " << direct_seconds << " s\n";
	std::cout << "Corrected-hybrid query time: " << hybrid_seconds << " s\n";
	std::cout << "Average corrected sources per query: "
		<< static_cast<double>(corrected_source_count) / kQueryCount << "\n";
	std::cout << "Corrected-hybrid speedup: " << speedup << "x\n";
	if (!std::isfinite(sink) || speedup < kMinimumSpeedup) {
		std::cerr << "Required corrected-hybrid speedup is "
			<< kMinimumSpeedup << "x" << std::endl;
		return 1;
	}
	return 0;
}
