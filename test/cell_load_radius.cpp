#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iostream>
#include <vector>

#include "Aggregate.hpp"
#include "Ionp.hpp"
#include "Point3D.hpp"
#include "RandNumberBetween.hpp"
#include "SimMatrix.hpp"
#include "SimSpace.hpp"

namespace {

bool approximately_equal(double actual, double expected, double relative_tolerance) {
	return std::abs(actual - expected) <=
		relative_tolerance * std::max(std::abs(expected), 1e-30);
}

} // namespace

int main() {
	rng_config::set_seed(20260812u);

	const int intracellular_particle_count = 10;
	const int extracellular_particle_count = 4;
	const double cell_load_fg = 5428.134;
	const double expected_equivalent_radius = std::pow(
		(6.299e-20 * cell_load_fg) /
			(1000. * static_cast<double>(intracellular_particle_count)),
		1. / 3.);
	const double primary_radius = 4e-9;
	const double secondary_radius = 12.5e-9;
	const double configured_outer_radius = 50e-9;

	SimMatrix matrix(
		0., 0., 0.,
		188.207205776206e-6,
		188.207205776206e-6,
		188.207205776206e-6,
		0, 0, 0);
	Point3D cell_origin;
	Point3D voxel_position;
	std::vector<Point3D> cell_origins;
	std::vector<Ionp> ionps;

	SimSpace::simulationSpace_inCell(
		ionps,
		matrix,
		intracellular_particle_count,
		0.,
		primary_radius,
		secondary_radius,
		0.,
		1e-12,
		1e-12,
		configured_outer_radius,
		cell_origin,
		10e-6,
		1,
		cell_load_fg,
		0.,
		1,
		cell_origins,
		voxel_position,
		false,
		std::filesystem::current_path());

	if (ionps.size() != static_cast<std::size_t>(intracellular_particle_count)) {
		std::cerr << "Unexpected intracellular particle count\n";
		return 1;
	}
	for (const Ionp& ionp : ionps) {
		if (!approximately_equal(ionp.radius(), expected_equivalent_radius, 1e-12)) {
			std::cerr << "Cell-load-derived radius did not override the configured distribution\n";
			return 1;
		}
		if (ionp.coating() != 0.) {
			std::cerr << "A configured outer radius below the equivalent core radius must not create a shell\n";
			return 1;
		}
	}

	SimSpace::simulationSpace_outCell(
		ionps,
		matrix,
		extracellular_particle_count,
		0.,
		primary_radius,
		secondary_radius,
		0.,
		1e-12,
		1e-12,
		configured_outer_radius,
		cell_origin,
		10e-6,
		1,
		cell_origins,
		voxel_position,
		false,
		std::filesystem::current_path());

	const std::size_t expected_total =
		static_cast<std::size_t>(intracellular_particle_count +
			extracellular_particle_count);
	if (ionps.size() != expected_total) {
		std::cerr << "Unexpected total particle count\n";
		return 1;
	}
	for (std::size_t index = intracellular_particle_count;
		 index < ionps.size();
		 ++index) {
		if (!approximately_equal(ionps[index].radius(), secondary_radius, 1e-9)) {
			std::cerr << "Extracellular placement did not preserve the configured radius distribution\n";
			return 1;
		}
		if (!approximately_equal(
				ionps[index].outer_radius(), configured_outer_radius, 1e-12)) {
			std::cerr << "Extracellular coating outer radius changed unexpectedly\n";
			return 1;
		}
	}

	return 0;
}
