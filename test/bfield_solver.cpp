#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "CorrectedHybridField.hpp"
#include "Ionp.hpp"
#include "SimMatrix.hpp"
#include "SimSpace.hpp"
#include "mptmacros.h"

namespace {

void expect(bool condition, const char* message) {
	if (!condition) {
		throw std::runtime_error(message);
	}
}

void expect_close(double actual,
				  double expected,
				  const char* message,
				  double relative_tolerance = 1e-11,
				  double absolute_tolerance = 1e-20) {
	const double scale = std::max(std::abs(actual), std::abs(expected));
	if (std::abs(actual - expected) >
		absolute_tolerance + relative_tolerance * scale) {
		throw std::runtime_error(message);
	}
}

double amplitude_for(const Ionp& ionp) {
	return B_fixed * ionp.radius() * ionp.radius() * ionp.radius();
}

std::vector<double> amplitudes_for(const std::vector<Ionp>& ionps) {
	std::vector<double> amplitudes;
	amplitudes.reserve(ionps.size());
	for (const Ionp& ionp : ionps) {
		amplitudes.push_back(amplitude_for(ionp));
	}
	return amplitudes;
}

void expect_particle_copy_with_translation(
	const Ionp& actual,
	const Ionp& source,
	const Point3D& translation,
	double coordinate_tolerance = 0.0) {
	expect_close(actual.x0().x(), source.x0().x() + translation.x(),
				 "Translated IONP has the wrong initial x coordinate", 0.0,
				 coordinate_tolerance);
	expect_close(actual.x0().y(), source.x0().y() + translation.y(),
				 "Translated IONP has the wrong initial y coordinate", 0.0,
				 coordinate_tolerance);
	expect_close(actual.x0().z(), source.x0().z() + translation.z(),
				 "Translated IONP has the wrong initial z coordinate", 0.0,
				 coordinate_tolerance);
	expect_close(actual.xt().x(), source.xt().x() + translation.x(),
				 "Translated IONP has the wrong current x coordinate", 0.0,
				 coordinate_tolerance);
	expect_close(actual.xt().y(), source.xt().y() + translation.y(),
				 "Translated IONP has the wrong current y coordinate", 0.0,
				 coordinate_tolerance);
	expect_close(actual.xt().z(), source.xt().z() + translation.z(),
				 "Translated IONP has the wrong current z coordinate", 0.0,
				 coordinate_tolerance);
	expect_close(actual.Particle::diffTime(), source.Particle::diffTime(),
				 "Translated IONP lost its particle diffusion time", 0.0, 0.0);
	expect_close(actual.diffTime(), source.diffTime(),
				 "Translated IONP lost its IONP diffusion time", 0.0, 0.0);
	expect_close(actual.radius(), source.radius(),
				 "Translated IONP has the wrong core radius", 0.0, 0.0);
	expect_close(actual.coating(), source.coating(),
				 "Translated IONP has the wrong coating", 0.0, 0.0);
}

double explicit_surrounding_field(
	const SimMatrix& matrix,
	const std::vector<Ionp>& physical_ionps,
	const Point3D& position) {
	double result = 0.0;
	for (int x_offset = -1; x_offset <= 1; ++x_offset) {
		for (int y_offset = -1; y_offset <= 1; ++y_offset) {
			for (int z_offset = -1; z_offset <= 1; ++z_offset) {
				const Point3D translation(
					x_offset * matrix.xLength(),
					y_offset * matrix.yLength(),
					z_offset * matrix.zLength());
				for (const Ionp& ionp : physical_ionps) {
					result += matrix.calculateBfield(
						position, ionp.xt() + translation, ionp.radius());
				}
			}
		}
	}
	return result;
}

void test_surrounding_space_copies_complete_realized_particles() {
	const SimMatrix matrix(12.0, 18.0, 24.0, 5, 7, 9);
	Ionp coated(-4.2, 2.3, -7.1, 0.0, 0.23, 0.65);
	coated.setX0(Point3D(-4.7, 2.1, -7.4));
	coated.Particle::setDiffTime(1.25);
	Ionp uncoated(3.4, -5.6, 8.2, 0.0, 0.31, 0.0);
	uncoated.setX0(Point3D(3.0, -5.1, 8.0));
	uncoated.Particle::setDiffTime(2.5);
	const std::vector<Ionp> physical_ionps = {coated, uncoated};
	const std::vector<Ionp> original_ionps = physical_ionps;

	const std::vector<Ionp> field_sources =
		SimSpace::surroundingSpace(physical_ionps, matrix);
	expect(field_sources.size() == 27 * physical_ionps.size(),
		   "Surrounding space did not create exactly 27 source sets");

	for (std::size_t particle = 0; particle < physical_ionps.size(); ++particle) {
		expect_particle_copy_with_translation(
			field_sources[particle], physical_ionps[particle], Point3D());
		expect_particle_copy_with_translation(
			physical_ionps[particle], original_ionps[particle], Point3D());
	}

	std::size_t source_index = physical_ionps.size();
	for (int x_offset = -1; x_offset <= 1; ++x_offset) {
		for (int y_offset = -1; y_offset <= 1; ++y_offset) {
			for (int z_offset = -1; z_offset <= 1; ++z_offset) {
				if (x_offset == 0 && y_offset == 0 && z_offset == 0) {
					continue;
				}
				const Point3D translation(
					x_offset * matrix.xLength(),
					y_offset * matrix.yLength(),
					z_offset * matrix.zLength());
				for (const Ionp& ionp : physical_ionps) {
					expect_particle_copy_with_translation(
						field_sources.at(source_index++), ionp, translation);
				}
			}
		}
	}
	expect(source_index == field_sources.size(),
		   "Surrounding space produced a missing or duplicate translation");

	const std::vector<Ionp> repeated_sources =
		SimSpace::surroundingSpace(physical_ionps, matrix);
	for (std::size_t particle = 0; particle < field_sources.size(); ++particle) {
		expect_particle_copy_with_translation(
			repeated_sources[particle], field_sources[particle], Point3D());
	}
}

void test_surrounding_space_rejects_non_positive_dimensions() {
	const std::vector<Ionp> ionps = {
		Ionp(-0.35, 0.45, -0.2, 0.0, 0.23, 0.0)};
	const std::vector<SimMatrix> invalid_matrices = {
		SimMatrix(0.0, 18.0, 24.0, 2, 2, 2),
		SimMatrix(12.0, 0.0, 24.0, 2, 2, 2),
		SimMatrix(12.0, 18.0, -1.0, 2, 2, 2)};
	for (const SimMatrix& matrix : invalid_matrices) {
		bool threw = false;
		try {
			(void) SimSpace::surroundingSpace(ionps, matrix);
		} catch (const std::invalid_argument&) {
			threw = true;
		}
		expect(threw,
			   "Surrounding space accepted a non-positive voxel dimension");
	}
}

double core_radius_cubed_sum(
	const std::vector<Ionp>& ionps,
	const std::vector<std::size_t>& indices) {
	double sum = 0.0;
	for (std::size_t index : indices) {
		const double radius = ionps.at(index).radius();
		sum += radius * radius * radius;
	}
	return sum;
}

void test_surrounding_superparticles_conserve_each_group_and_repeat_template() {
	const SimMatrix matrix(30.0, 30.0, 30.0, 4, 4, 4);
	const std::vector<Point3D> cells = {
		Point3D(-8.0, 0.0, 0.0),
		Point3D(8.0, 0.0, 0.0)};
	constexpr double cell_radius = 5.0;
	const std::vector<Ionp> physical_ionps = {
		Ionp(-9.0, 0.0, 0.0, 0.0, 0.40, 0.10),
		Ionp(-7.5, 1.0, 0.0, 0.0, 0.55, 0.00),
		Ionp(-8.0, -1.0, 1.0, 0.0, 0.65, 0.00),
		Ionp(7.0, 0.0, 0.0, 0.0, 0.70, 0.00),
		Ionp(8.0, 1.0, 0.0, 0.0, 0.80, 0.00),
		Ionp(9.0, -1.0, 0.0, 0.0, 0.90, 0.00),
		Ionp(0.0, 8.0, 0.0, 0.0, 0.30, 0.20),
		Ionp(0.0, -8.0, 0.0, 0.0, 0.40, 0.00),
		Ionp(0.0, 0.0, 8.0, 0.0, 0.50, 0.00)};
	constexpr std::size_t cap = 2;
	constexpr std::size_t expected_template_count = 6;

	rng_config::set_seed(0x5a17u);
	const std::vector<Ionp> field_sources = SimSpace::surroundingSpace(
		physical_ionps,
		matrix,
		cap,
		cells,
		Point3D(),
		cell_radius);
	expect(
		field_sources.size() ==
			physical_ionps.size() + 26 * expected_template_count,
		"Per-group compression produced the wrong magnetic source count");
	for (std::size_t index = 0; index < physical_ionps.size(); ++index) {
		expect_particle_copy_with_translation(
			field_sources[index], physical_ionps[index], Point3D());
	}

	const Point3D first_translation(
		-matrix.xLength(), -matrix.yLength(), -matrix.zLength());
	std::vector<std::size_t> template_group_counts(3, 0);
	std::vector<double> template_group_iron(3, 0.0);
	for (std::size_t template_index = 0;
		 template_index < expected_template_count; ++template_index) {
		const Ionp& proxy =
			field_sources[physical_ionps.size() + template_index];
		const Point3D position = proxy.xt() - first_translation;
		expect_close(proxy.coating(), 0.0,
			"Compressed surrounding proxy retained a coating", 0.0, 0.0);
		expect(std::abs(position.x()) + proxy.radius() <=
			   matrix.xLength() / 2.0 + 1e-12 &&
			   std::abs(position.y()) + proxy.radius() <=
			   matrix.yLength() / 2.0 + 1e-12 &&
			   std::abs(position.z()) + proxy.radius() <=
			   matrix.zLength() / 2.0 + 1e-12,
			"A surrounding proxy core extends outside the centered voxel");

		std::size_t group = 2;
		std::size_t containing_cell_count = 0;
		for (std::size_t cell_index = 0; cell_index < cells.size(); ++cell_index) {
			if (distance(position, cells[cell_index]) + proxy.radius() <=
				cell_radius + 1e-12) {
				group = cell_index;
				++containing_cell_count;
			}
		}
		expect(containing_cell_count <= 1,
			"A generated proxy is ambiguously contained by multiple cells");
		if (group == 2) {
			for (const Point3D& cell : cells) {
				expect(distance(position, cell) + 1e-12 >=
					   cell_radius + proxy.radius(),
					"An extracellular proxy core intersects a cell");
			}
		}
		++template_group_counts[group];
		template_group_iron[group] +=
			proxy.radius() * proxy.radius() * proxy.radius();
	}
	for (std::size_t group = 0; group < template_group_counts.size(); ++group) {
		expect(template_group_counts[group] == cap,
			"A compressed group did not use the configured per-group cap");
	}
	expect_close(
		template_group_iron[0],
		core_radius_cubed_sum(physical_ionps, {0, 1, 2}),
		"First cell's proxy template changed its iron amount");
	expect_close(
		template_group_iron[1],
		core_radius_cubed_sum(physical_ionps, {3, 4, 5}),
		"Second cell's proxy template changed its iron amount");
	expect_close(
		template_group_iron[2],
		core_radius_cubed_sum(physical_ionps, {6, 7, 8}),
		"Extracellular proxy template changed its iron amount");

	for (std::size_t left = 0; left < expected_template_count; ++left) {
		const Ionp& left_proxy =
			field_sources[physical_ionps.size() + left];
		for (std::size_t right = left + 1;
			 right < expected_template_count; ++right) {
			const Ionp& right_proxy =
				field_sources[physical_ionps.size() + right];
			expect(distance(left_proxy.xt(), right_proxy.xt()) + 1e-12 >=
				   left_proxy.radius() + right_proxy.radius(),
				"Generated surrounding proxies overlap");
		}
	}

	std::size_t source_index = physical_ionps.size();
	for (int x_offset = -1; x_offset <= 1; ++x_offset) {
		for (int y_offset = -1; y_offset <= 1; ++y_offset) {
			for (int z_offset = -1; z_offset <= 1; ++z_offset) {
				if (x_offset == 0 && y_offset == 0 && z_offset == 0) {
					continue;
				}
				const Point3D translation_difference(
					x_offset * matrix.xLength() - first_translation.x(),
					y_offset * matrix.yLength() - first_translation.y(),
					z_offset * matrix.zLength() - first_translation.z());
				for (std::size_t template_index = 0;
					 template_index < expected_template_count; ++template_index) {
					expect_particle_copy_with_translation(
						field_sources.at(source_index++),
						field_sources[physical_ionps.size() + template_index],
						translation_difference,
						1e-12);
				}
			}
		}
	}
	expect(source_index == field_sources.size(),
		"Compressed surrounding voxels did not repeat one complete template");

	rng_config::set_seed(0x5a17u);
	const std::vector<Ionp> repeated_sources = SimSpace::surroundingSpace(
		physical_ionps,
		matrix,
		cap,
		cells,
		Point3D(),
		cell_radius);
	expect(repeated_sources.size() == field_sources.size(),
		"Seeded superparticle generation changed its source count");
	for (std::size_t index = 0; index < field_sources.size(); ++index) {
		expect_particle_copy_with_translation(
			repeated_sources[index], field_sources[index], Point3D());
	}
}

void test_surrounding_superparticles_keep_small_groups_exact() {
	const SimMatrix matrix(20.0, 20.0, 20.0, 3, 3, 3);
	Ionp intracellular(0.5, 0.0, 0.0, 0.0, 0.4, 0.2);
	intracellular.setX0(Point3D(0.3, 0.1, 0.0));
	Ionp extracellular(7.0, 0.0, 0.0, 0.0, 0.6, 0.1);
	extracellular.setX0(Point3D(6.8, 0.0, 0.1));
	const std::vector<Ionp> physical_ionps = {
		intracellular, extracellular};
	const std::vector<Point3D> cells = {Point3D()};

	const std::vector<Ionp> exact_sources =
		SimSpace::surroundingSpace(physical_ionps, matrix);
	const std::vector<Ionp> capped_sources = SimSpace::surroundingSpace(
		physical_ionps, matrix, 2, cells, Point3D(), 4.0);
	expect(capped_sources.size() == exact_sources.size(),
		"Groups at or below the cap were unexpectedly compressed");
	for (std::size_t index = 0; index < exact_sources.size(); ++index) {
		expect_particle_copy_with_translation(
			capped_sources[index], exact_sources[index], Point3D());
	}

	const std::vector<Ionp> straddling = {
		Ionp(4.0, 0.0, 0.0, 0.0, 0.5, 0.0)};
	const std::vector<Ionp> zero_cap = SimSpace::surroundingSpace(
		straddling, matrix, 0, cells, Point3D(), 4.0);
	const std::vector<Ionp> zero_cap_reference =
		SimSpace::surroundingSpace(straddling, matrix);
	expect(zero_cap.size() == zero_cap_reference.size(),
		"A zero superparticle cap did not retain exact replication");
	for (std::size_t index = 0; index < zero_cap.size(); ++index) {
		expect_particle_copy_with_translation(
			zero_cap[index], zero_cap_reference[index], Point3D());
	}
}

void test_surrounding_superparticles_support_no_cell_and_zero_iron_groups() {
	const SimMatrix matrix(20.0, 20.0, 20.0, 3, 3, 3);
	const std::vector<Ionp> physical_ionps = {
		Ionp(-6.0, -5.0, -4.0, 0.0, 0.35, 0.0),
		Ionp(6.0, 5.0, 4.0, 0.0, 0.45, 0.0),
		Ionp(-5.0, 5.0, 4.0, 0.0, 0.55, 0.0),
		Ionp(5.0, -5.0, -4.0, 0.0, 0.65, 0.0)};
	rng_config::set_seed(0x31u);
	const std::vector<Ionp> field_sources = SimSpace::surroundingSpace(
		physical_ionps, matrix, 2, {}, Point3D(100.0, 0.0, 0.0), 0.0);
	expect(field_sources.size() == physical_ionps.size() + 26 * 2,
		"No-cell compression did not create one extracellular template group");
	double proxy_iron = 0.0;
	for (std::size_t index = 0; index < 2; ++index) {
		const Ionp& proxy = field_sources[physical_ionps.size() + index];
		proxy_iron += proxy.radius() * proxy.radius() * proxy.radius();
		const Point3D position = proxy.xt() + Point3D(
			matrix.xLength(), matrix.yLength(), matrix.zLength());
		expect(std::abs(position.x()) + proxy.radius() <= 10.0 + 1e-12 &&
			   std::abs(position.y()) + proxy.radius() <= 10.0 + 1e-12 &&
			   std::abs(position.z()) + proxy.radius() <= 10.0 + 1e-12,
			"A no-cell extracellular proxy extends outside the voxel");
	}
	expect_close(
		proxy_iron,
		core_radius_cubed_sum(physical_ionps, {0, 1, 2, 3}),
		"No-cell extracellular compression changed the iron amount");

	const std::vector<Ionp> zero_iron = {
		Ionp(-1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
		Ionp(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
		Ionp(1.0, 0.0, 0.0, 0.0, 0.0, 0.0)};
	const std::vector<Ionp> zero_iron_sources = SimSpace::surroundingSpace(
		zero_iron, matrix, 2, {}, Point3D(), 0.0);
	expect(zero_iron_sources.size() == zero_iron.size(),
		"A compressed zero-iron group created meaningless zero-radius proxies");
}

void test_surrounding_superparticles_reject_invalid_geometry() {
	const SimMatrix matrix(20.0, 20.0, 20.0, 3, 3, 3);
	bool threw = false;
	try {
		(void) SimSpace::surroundingSpace(
			{Ionp(5.0, 0.0, 0.0, 0.0, 0.5, 0.0)},
			matrix,
			1,
			{Point3D()},
			Point3D(),
			5.0);
	}
	catch (const std::runtime_error& error) {
		threw = std::string(error.what()).find("straddles") !=
			std::string::npos;
	}
	expect(threw,
		"Surrounding compression accepted a membrane-straddling IONP");

	threw = false;
	try {
		(void) SimSpace::surroundingSpace(
			{Ionp(0.5, 0.0, 0.0, 0.0, 0.5, 0.0)},
			matrix,
			1,
			{Point3D(), Point3D(1.0, 0.0, 0.0)},
			Point3D(),
			5.0);
	}
	catch (const std::runtime_error& error) {
		threw = std::string(error.what()).find("more than one cell") !=
			std::string::npos;
	}
	expect(threw,
		"Surrounding compression accepted an ambiguously contained IONP");

	threw = false;
	try {
		(void) SimSpace::surroundingSpace(
			{Ionp(-0.1, 0.0, 0.0, 0.0, 0.8, 0.0),
			 Ionp(0.1, 0.0, 0.0, 0.0, 0.8, 0.0)},
			matrix,
			1,
			{Point3D()},
			Point3D(),
			1.0);
	}
	catch (const std::runtime_error& error) {
		threw = std::string(error.what()).find(
			"Increase Surrounding_voxel_superparticles") != std::string::npos;
	}
	expect(threw,
		"An impossible proxy geometry did not recommend increasing the cap");
}

void test_surrounding_sources_match_explicit_field_reference() {
	const SimMatrix matrix(12.0, 18.0, 24.0, 4, 5, 6);
	const std::vector<Ionp> physical_ionps = {
		Ionp(-5.65, 0.35, -1.15, 0.0, 0.23, 0.4),
		Ionp(1.45, -4.75, 7.35, 0.0, 0.31, 0.0)};
	const std::vector<Ionp> field_sources =
		SimSpace::surroundingSpace(physical_ionps, matrix);
	std::vector<double> grid;
	matrix.calculateBfield(grid, field_sources, 0.0);

	const double x_step = matrix.xLength() / (matrix.nxNodes() - 1.0);
	const double y_step = matrix.yLength() / (matrix.nyNodes() - 1.0);
	const double z_step = matrix.zLength() / (matrix.nzNodes() - 1.0);
	std::size_t grid_index = 0;
	for (int z = 0; z < matrix.nzNodes(); ++z) {
		for (int x = 0; x < matrix.nxNodes(); ++x) {
			for (int y = 0; y < matrix.nyNodes(); ++y) {
				const Point3D node(
					-matrix.xLength() / 2.0 + x_step * x,
					-matrix.yLength() / 2.0 + y_step * y,
					-matrix.zLength() / 2.0 + z_step * z);
				expect_close(
					grid.at(grid_index++),
					explicit_surrounding_field(matrix, physical_ionps, node),
					"Surrounding-source grid differs from the explicit 27-voxel sum",
					2e-13);
			}
		}
	}

	const Point3D edge_query(5.8, 0.2, -1.0);
	double direct_field = 0.0;
	for (const Ionp& ionp : field_sources) {
		direct_field += matrix.calculateBfield(
			edge_query, ionp.xt(), ionp.radius());
	}
	const double explicit_field =
		explicit_surrounding_field(matrix, physical_ionps, edge_query);
	expect_close(
		direct_field,
		explicit_field,
		"Direct edge field differs from the explicit 27-voxel sum",
		2e-13);
	double central_field = 0.0;
	for (const Ionp& ionp : physical_ionps) {
		central_field += matrix.calculateBfield(
			edge_query, ionp.xt(), ionp.radius());
	}
	expect(std::abs(direct_field - central_field) >
		   1e-6 * std::max(std::abs(direct_field), std::abs(central_field)),
		   "Edge reference did not exercise a surrounding-voxel contribution");
}

void test_corrected_hybrid_includes_near_boundary_replica() {
	const SimMatrix matrix(12.0, 18.0, 24.0, 13, 19, 25);
	const std::vector<Ionp> physical_ionps = {
		Ionp(-5.7, 0.2, -0.1, 0.0, 0.23, 0.0)};
	const std::vector<Ionp> field_sources =
		SimSpace::surroundingSpace(physical_ionps, matrix);
	const std::vector<double> amplitudes = amplitudes_for(field_sources);
	std::vector<double> grid;
	matrix.calculateBfield(grid, field_sources, 0.0);
	const CorrectedHybridField field(
		matrix, field_sources, amplitudes, grid);
	const Point3D query(5.7, 0.1, -0.2);
	const Point3D near_replica_position(
		physical_ionps.front().xt().x() + matrix.xLength(),
		physical_ionps.front().xt().y(),
		physical_ionps.front().xt().z());

	std::size_t near_replica = field_sources.size();
	for (std::size_t index = 0; index < field_sources.size(); ++index) {
		if (field_sources[index].xt() == near_replica_position) {
			near_replica = index;
			break;
		}
	}
	expect(near_replica < field_sources.size(),
		   "Could not find the neighboring image near the voxel boundary");
	expect(distance(query, near_replica_position) <=
		   field.correctionStart(near_replica),
		   "Boundary image is not in the exact hybrid correction region");

	double expected = matrix.calculateBfieldWithAmplitude(
		query, near_replica_position, amplitudes[near_replica]);
	for (std::size_t index = 0; index < field_sources.size(); ++index) {
		if (index == near_replica) {
			continue;
		}
		expect(distance(query, field_sources[index].xt()) >=
			   field.correctionEnd(index),
			   "Hybrid boundary reference has an unintended second near source");
		expected += matrix.interpolateBfieldContribution(
			query, field_sources[index].xt(), amplitudes[index]);
	}
	expect_close(
		field.evaluate(query),
		expected,
		"Corrected hybrid field did not correct the neighboring image",
		2e-13);
}

void test_full_grid_is_all_particle_reference() {
	SimMatrix matrix(24.0, 13);
	const std::vector<Ionp> ionps = {
		Ionp(-0.35, 0.45, -0.2, 0.0, 0.23, 0.6),
		Ionp(7.1, -6.3, 5.4, 0.0, 0.31, 0.0),
		Ionp(-8.2, 7.4, 6.7, 0.0, 0.19, 0.2)};
	std::vector<double> grid;
	matrix.calculateBfield(grid, ionps, 0.0);

	const double step = matrix.xLength() /
		static_cast<double>(matrix.nxNodes() - 1);
	std::size_t index = 0;
	for (int z = 0; z < matrix.nzNodes(); ++z) {
		for (int x = 0; x < matrix.nxNodes(); ++x) {
			for (int y = 0; y < matrix.nyNodes(); ++y) {
				const Point3D node(
					-matrix.xLength() / 2.0 + step * x,
					-matrix.yLength() / 2.0 + step * y,
					-matrix.zLength() / 2.0 + step * z);
				double reference = 0.0;
				for (const Ionp& ionp : ionps) {
					reference += matrix.calculateBfield(
						node, ionp.xt(), ionp.radius());
				}
				expect_close(
					grid.at(index++),
					reference,
					"Full grid differs from the all-IONP reference at a node",
					0.0,
					0.0);
			}
		}
	}
}

void test_coating_does_not_change_grid_amplitude() {
	SimMatrix matrix(24.0, 13);
	const std::vector<Ionp> uncoated = {
		Ionp(-0.35, 0.45, -0.2, 0.0, 0.23, 0.0)};
	const std::vector<Ionp> coated = {
		Ionp(-0.35, 0.45, -0.2, 0.0, 0.23, 1.7)};
	std::vector<double> uncoated_grid;
	std::vector<double> coated_grid;
	matrix.calculateBfield(uncoated_grid, uncoated, 0.0);
	matrix.calculateBfield(coated_grid, coated, 0.0);
	expect(uncoated_grid == coated_grid,
		   "Coating radius changed the magnetic grid amplitude");
}

void test_corrected_hybrid_single_particle() {
	SimMatrix matrix(24.0, 13);
	const std::vector<Ionp> ionps = {
		Ionp(-0.35, 0.45, -0.2, 0.0, 0.23, 0.6)};
	const std::vector<double> amplitudes = amplitudes_for(ionps);
	std::vector<double> grid;
	matrix.calculateBfield(grid, ionps, 0.0);
	const CorrectedHybridField field(matrix, ionps, amplitudes, grid);

	const std::vector<Point3D> exact_queries = {
		// Inside the coating.
		Point3D(0.15, 0.65, -0.1),
		// Between grid nodes.
		Point3D(2.85, 1.55, -0.9),
		// Exactly at a grid node.
		Point3D(2.0, 2.0, 2.0)};
	for (const Point3D& query : exact_queries) {
		expect_close(
			field.evaluate(query),
			matrix.calculateBfieldWithAmplitude(
				query, ionps.front().xt(), amplitudes.front()),
			"Corrected hybrid did not reproduce the exact near field");
	}

	const double r0 = field.correctionStart(0);
	const double r1 = field.correctionEnd(0);
	const Point3D at_r0(
		ionps.front().xt().x() + r0,
		ionps.front().xt().y(),
		ionps.front().xt().z());
	const Point3D at_r1(
		ionps.front().xt().x() + r1,
		ionps.front().xt().y(),
		ionps.front().xt().z());
	const Point3D at_midpoint(
		ionps.front().xt().x() + (r0 + r1) / 2.0,
		ionps.front().xt().y(),
		ionps.front().xt().z());

	expect_close(
		field.evaluate(at_r0),
		matrix.calculateBfieldWithAmplitude(
			at_r0, ionps.front().xt(), amplitudes.front()),
		"Correction was not fully exact at R0");
	expect_close(
		field.evaluate(at_r1),
		matrix.interpolate(at_r1, grid),
		"Correction did not finish at R1");
	const double midpoint_exact = matrix.calculateBfieldWithAmplitude(
		at_midpoint, ionps.front().xt(), amplitudes.front());
	const double midpoint_grid = matrix.interpolate(at_midpoint, grid);
	expect_close(
		field.evaluate(at_midpoint),
		midpoint_grid + 0.5 * (midpoint_exact - midpoint_grid),
		"Correction did not use cubic smoothstep in the R0/R1 transition");
}

void test_corrected_hybrid_preserves_far_sources() {
	SimMatrix matrix(24.0, 13);
	const std::vector<Ionp> ionps = {
		Ionp(-0.35, 0.45, -0.2, 0.0, 0.23, 0.25),
		Ionp(9.1, 8.8, 9.3, 0.0, 0.31, 0.0)};
	const std::vector<double> amplitudes = amplitudes_for(ionps);
	std::vector<double> grid;
	matrix.calculateBfield(grid, ionps, 0.0);
	const CorrectedHybridField field(matrix, ionps, amplitudes, grid);
	const Point3D query(0.75, 0.85, 0.1);

	const double expected = matrix.calculateBfieldWithAmplitude(
		query, ionps[0].xt(), amplitudes[0]) +
		matrix.interpolateBfieldContribution(
			query, ionps[1].xt(), amplitudes[1]);
	expect_close(
		field.evaluate(query),
		expected,
		"Corrected hybrid lost or double-counted a distant source");
}

void test_corrected_hybrid_at_volume_boundary() {
	SimMatrix matrix(24.0, 13);
	const std::vector<Ionp> ionps = {
		Ionp(-10.4, -10.1, -9.8, 0.0, 0.23, 0.4)};
	const std::vector<double> amplitudes = amplitudes_for(ionps);
	std::vector<double> grid;
	matrix.calculateBfield(grid, ionps, 0.0);
	const CorrectedHybridField field(matrix, ionps, amplitudes, grid);
	const Point3D query(-12.0, -11.4, -11.2);
	expect_close(
		field.evaluate(query),
		matrix.calculateBfieldWithAmplitude(
			query, ionps.front().xt(), amplitudes.front()),
		"Corrected hybrid failed at a volume boundary");
}

void test_corrected_hybrid_rejects_invalid_grids() {
	const std::vector<Ionp> ionps = {
		Ionp(-0.35, 0.45, -0.2, 0.0, 0.23, 0.0)};
	const std::vector<double> amplitudes = amplitudes_for(ionps);
	bool threw = false;
	try {
		const SimMatrix invalid_matrix(24.0, 1);
		const std::vector<double> grid(1, 0.0);
		const CorrectedHybridField field(
			invalid_matrix, ionps, amplitudes, grid);
		(void) field;
	} catch (const std::invalid_argument&) {
		threw = true;
	}
	expect(threw, "Corrected hybrid accepted a one-node grid");

	threw = false;
	try {
		const SimMatrix matrix(24.0, 13);
		const std::vector<double> wrong_size_grid(8, 0.0);
		const CorrectedHybridField field(
			matrix, ionps, amplitudes, wrong_size_grid);
		(void) field;
	} catch (const std::invalid_argument&) {
		threw = true;
	}
	expect(threw, "Corrected hybrid accepted a mismatched grid size");
}

} // namespace

int main() {
	try {
		test_surrounding_space_copies_complete_realized_particles();
		test_surrounding_space_rejects_non_positive_dimensions();
		test_surrounding_superparticles_conserve_each_group_and_repeat_template();
		test_surrounding_superparticles_keep_small_groups_exact();
		test_surrounding_superparticles_support_no_cell_and_zero_iron_groups();
		test_surrounding_superparticles_reject_invalid_geometry();
		test_surrounding_sources_match_explicit_field_reference();
		test_corrected_hybrid_includes_near_boundary_replica();
		test_full_grid_is_all_particle_reference();
		test_coating_does_not_change_grid_amplitude();
		test_corrected_hybrid_single_particle();
		test_corrected_hybrid_preserves_far_sources();
		test_corrected_hybrid_at_volume_boundary();
		test_corrected_hybrid_rejects_invalid_grids();
	} catch (const std::exception& error) {
		std::cerr << error.what() << std::endl;
		return 1;
	}
	return 0;
}
