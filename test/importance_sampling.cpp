#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "ImportanceSampling.hpp"
#include "Ionp.hpp"
#include "MRsequence.hpp"
#include "mptmacros.h"
#include "Point3D.hpp"
#include "Proton.hpp"
#include "RandomList.hpp"
#include "RandNumberBetween.hpp"
#include "SimMatrix.hpp"
#include "SimSpace.hpp"

bool should_use_variable_time_steps(bool var_steps,
									bool importance_sampling_enabled,
									ProtonSampleRegion sample_region);
bool should_use_field_grid_for_step(bool normal_grid_selection,
									bool variable_step_active);
double diffusion_coefficient_for_position(const Point3D& position,
										  const std::vector<Ionp>& ionpInstances,
										  const IonpSpatialIndex& ionp_index,
										  IonpCheckMode ionp_check_mode,
										  double coating_diffusion_coefficient,
										  bool in_cell,
										  double cell_diffusion_coefficient);
int coating_index_for_position(const Point3D& position,
							   const std::vector<Ionp>& ionpInstances,
							   const IonpSpatialIndex& ionp_index,
							   IonpCheckMode ionp_check_mode);
bool has_near_ionp_for_variable_step(const Point3D& position,
									 const std::vector<Ionp>& ionpInstances,
									 const IonpSpatialIndex& ionp_index,
									 IonpCheckMode ionp_check_mode);
bool enforce_ionp_boundaries(Proton& proton,
							 Point3D& lastPosition,
							 const std::vector<Ionp>& ionpInstances,
							 const IonpSpatialIndex& ionp_index,
							 IonpCheckMode ionp_check_mode,
							 bool force_all_ionp_candidates,
							 int& coating_ionp_index,
							 bool& coating_state_known,
							 double current_time,
							 bool check_coating_permeability,
							 double coating_permeability,
							 RandomList& random_unit,
							 std::vector<TrajectoryEvent>* event_log);
bool enforce_cell_boundaries(Proton& proton,
							 Point3D& lastPosition,
							 const std::vector<Point3D>& cell_origins,
							 const Point3D& fallback_origin,
							 double cell_radius,
							 double permeability,
							 RandomList& random_unit,
							 double current_time,
							 std::vector<TrajectoryEvent>* event_log);
bool enforce_cell_boundaries_periodic(Proton& proton,
								  Point3D& lastPosition,
								  const std::vector<Point3D>& cell_origins,
								  const Point3D& fallback_origin,
								  double cell_radius,
								  double permeability,
								  double voxel_x_length,
								  double voxel_y_length,
								  double voxel_z_length,
								  RandomList& random_unit,
								  double current_time,
								  std::vector<TrajectoryEvent>* event_log);

namespace {

constexpr double kPi = 3.14159265358979323846;

bool approx_equal(double lhs, double rhs, double tolerance = 1e-10) {
	return std::abs(lhs - rhs) <= tolerance;
}

void expect(bool condition, const char* message) {
	if (!condition) {
		throw std::runtime_error(message);
	}
}

void expect_indices_equal(const std::vector<unsigned int>& actual,
						  const std::vector<unsigned int>& expected,
						  const char* message) {
	if (actual != expected) {
		throw std::runtime_error(message);
	}
}

void expect_complex_close(const std::complex<double>& actual,
						  const std::complex<double>& expected,
						  double tolerance = 1e-10) {
	expect(std::abs(actual.real() - expected.real()) <= tolerance, "Complex real part mismatch");
	expect(std::abs(actual.imag() - expected.imag()) <= tolerance, "Complex imaginary part mismatch");
}

MRsequence make_sequence(const SimMatrix& sim,
						 int free_protons,
						 int bound_protons,
						 double echo_time,
						 double time_step,
						 double ionp_radius,
						 double coating_thickness,
						 const Point3D& cell_origin,
						 double cell_radius,
						 ImportanceSamplingConfig importance_sampling,
						 bool variable_steps = false,
						 IonpCheckMode ionp_check_mode = IonpCheckMode::Spatial,
						 double magnetic_moment = 0.0) {
	Point3D direction(1.0, 0.0, 0.0);
	Point3D cube_origin(0.0, 0.0, 0.0);
	Point3D cube_size(0.0, 0.0, 0.0);
	SimMatrix sim_copy = sim;
	Point3D cell_origin_copy = cell_origin;
	return MRsequence(sim_copy,
					  0,
					  free_protons,
					  bound_protons,
					  ionp_radius,
					  echo_time,
					  time_step,
					  time_step,
					  1,
					  0,
					  sim.xLength(),
					  direction,
					  0.0,
					  coating_thickness,
					  1e-10,
					  D_D_CONST,
					  false,
					  1.0,
					  0.0,
					  0.0,
					  magnetic_moment,
					  cell_origin_copy,
					  cell_radius,
					  0.0,
					  cube_origin,
					  cube_size,
					  false,
					  variable_steps,
					  false,
					  importance_sampling,
					  ionp_check_mode);
}

double point_distance(const Point3D& lhs, const Point3D& rhs) {
	const double dx = lhs.x() - rhs.x();
	const double dy = lhs.y() - rhs.y();
	const double dz = lhs.z() - rhs.z();
	return std::sqrt(dx * dx + dy * dy + dz * dz);
}

bool proposal_contains_membership(
	FreeSamplingProposal proposal,
	const JointSampleRegion& membership) {
	switch (proposal) {
		case FreeSamplingProposal::Uniform:
			return true;
		case FreeSamplingProposal::CellInside:
			return membership.cell == CellSampleRegion::Inside;
		case FreeSamplingProposal::CellNear:
			return membership.cell == CellSampleRegion::Near;
		case FreeSamplingProposal::CellFar:
			return membership.cell == CellSampleRegion::Far;
		case FreeSamplingProposal::IntracellularIonpNear:
			return membership.intracellular_ionp == IonpDistanceRegion::Near;
		case FreeSamplingProposal::IntracellularIonpIntermediate:
			return membership.intracellular_ionp ==
				IonpDistanceRegion::Intermediate;
		case FreeSamplingProposal::IntracellularIonpFar:
			return membership.intracellular_ionp == IonpDistanceRegion::Far;
		case FreeSamplingProposal::ExtracellularIonpNear:
			return membership.extracellular_ionp == IonpDistanceRegion::Near;
		case FreeSamplingProposal::ExtracellularIonpIntermediate:
			return membership.extracellular_ionp ==
				IonpDistanceRegion::Intermediate;
		case FreeSamplingProposal::ExtracellularIonpFar:
			return membership.extracellular_ionp == IonpDistanceRegion::Far;
	}
	return false;
}

struct BoundaryCheckResult {
	bool repeat;
	Point3D proton_position;
	Point3D last_position;
	int coating_ionp_index;
	bool coating_state_known;
};

struct CellBoundaryCheckResult {
	bool reflected;
	Point3D proton_position;
	Point3D last_position;
	bool in_cell;
	int cell_number;
	std::vector<TrajectoryEvent> events;
};

CellBoundaryCheckResult run_cell_boundary_check(const Point3D& start,
												 const Point3D& end,
												 const std::vector<Point3D>& cell_origins,
												 const Point3D& fallback_origin,
												 double cell_radius,
												 double permeability,
												 double random_draw,
												 bool initially_in_cell = false,
												 int initial_cell_number = 0) {
	RandomList random_unit(random_draw, random_draw, false);
	Point3D last_position = start;
	Proton proton(0.0, 0.0);
	proton.setXt(end);
	proton.inCell(initially_in_cell);
	proton.Cellnum(initial_cell_number);
	std::vector<TrajectoryEvent> events;

	const bool reflected = enforce_cell_boundaries(proton,
											   last_position,
											   cell_origins,
											   fallback_origin,
											   cell_radius,
											   permeability,
											   random_unit,
											   1.25,
											   &events);

	return {reflected, proton.xt(), last_position, proton.inCell(), proton.Cellnum(), events};
}

CellBoundaryCheckResult run_periodic_cell_boundary_check(const Point3D& start,
													  const Point3D& end,
													  const std::vector<Point3D>& cell_origins,
													  const Point3D& fallback_origin,
													  double cell_radius,
													  double permeability,
													  double random_draw) {
	RandomList random_unit(random_draw, random_draw, false);
	Point3D last_position = start;
	Proton proton(0.0, 0.0);
	proton.setXt(end);
	std::vector<TrajectoryEvent> events;

	const bool reflected = enforce_cell_boundaries_periodic(proton,
													last_position,
													cell_origins,
													fallback_origin,
													cell_radius,
													permeability,
													10.0,
													10.0,
													10.0,
													random_unit,
													2.5,
													&events);

	return {reflected, proton.xt(), last_position, proton.inCell(), proton.Cellnum(), events};
}

BoundaryCheckResult run_boundary_check(IonpCheckMode ionp_check_mode,
									   const std::vector<Ionp>& ionps,
									   const Point3D& start,
									   const Point3D& end,
									   bool check_coating_permeability,
									   double coating_permeability) {
	const IonpSpatialIndex ionp_index(ionps);
	RandomList random_unit(0.0, 1.0, false);
	Point3D last_position = start;
	Proton proton(0.0, 0.0);
	proton.setXt(end);
	int coating_ionp_index = -1;
	bool coating_state_known = false;

	const bool repeat = enforce_ionp_boundaries(proton,
												last_position,
												ionps,
												ionp_index,
												ionp_check_mode,
												false,
												coating_ionp_index,
												coating_state_known,
												0.0,
												check_coating_permeability,
												coating_permeability,
												random_unit,
												nullptr);

	return {repeat, proton.xt(), last_position, coating_ionp_index, coating_state_known};
}

void expect_boundary_results_equal(const BoundaryCheckResult& lhs, const BoundaryCheckResult& rhs) {
	expect(lhs.repeat == rhs.repeat, "IONP check modes disagreed on repeat state");
	expect(point_distance(lhs.proton_position, rhs.proton_position) <= 1e-12, "IONP check modes produced different proton positions");
	expect(point_distance(lhs.last_position, rhs.last_position) <= 1e-12, "IONP check modes produced different last positions");
	expect(lhs.coating_ionp_index == rhs.coating_ionp_index, "IONP check modes produced different coating indices");
	expect(lhs.coating_state_known == rhs.coating_state_known, "IONP check modes produced different coating state flags");
}

Point3D simulate_single_free_proton_final_position(ImportanceSamplingConfig config,
													 bool variable_steps,
													 unsigned int seed) {
	rng_config::set_seed(seed);

	SimMatrix sim(0.0, 0.0, 0.0, 20.0, 20.0, 20.0, 0, 0, 0);
	const Point3D cell_origin(0.0, 0.0, 0.0);
	MRsequence sequence = make_sequence(sim,
										 config.total_free_protons(),
										 0,
										 0.01,
										 0.01,
										 0.5,
										 0.0,
										 cell_origin,
										 1.0,
										 config,
										 variable_steps);

	std::vector<Proton> protons;
	std::vector<Ionp> ionps = {Ionp(0.0, 0.0, 0.0, 0.0, 0.5, 0.0)};
	std::vector<std::complex<double>> signal;
	std::vector<double> timeV;
	std::vector<std::complex<double>> signal_mse;
	std::vector<double> timeV_mse;
	std::vector<Point3D> cell_origins;

	sequence.Sequence(protons, ionps, signal, timeV, signal_mse, timeV_mse, cell_origins);
	expect(protons.size() == static_cast<std::size_t>(config.total_free_protons()), "Unexpected proton count in variable-step test");
	return protons.front().xt();
}

void test_variable_step_eligibility_helper() {
	for (std::size_t region_index = 0;
		 region_index < kFreeSampleRegionCount;
		 ++region_index) {
		const ProtonSampleRegion region =
			static_cast<ProtonSampleRegion>(region_index);
		expect(should_use_variable_time_steps(true, false, region), "Non-IS runs should keep variable steps enabled for every region");
		expect(!should_use_variable_time_steps(false, false, region), "var_Steps=false should disable variable steps without IS");
		expect(!should_use_variable_time_steps(false, true, region), "var_Steps=false should disable variable steps with IS");
	}

	expect(should_use_variable_time_steps(true, true, ProtonSampleRegion::Inside), "IS inside region should keep variable steps");
	expect(should_use_variable_time_steps(true, true, ProtonSampleRegion::Near), "IS near region should use dynamic variable steps");
	expect(should_use_variable_time_steps(true, true, ProtonSampleRegion::Far), "IS far region should use dynamic variable steps");
	expect(!should_use_variable_time_steps(true, true, ProtonSampleRegion::BoundIntracellular), "IS intracellular bound region should not use variable steps");
	expect(!should_use_variable_time_steps(true, true, ProtonSampleRegion::BoundExtracellular), "IS extracellular bound region should not use variable steps");
	expect(!should_use_variable_time_steps(true, false, ProtonSampleRegion::BoundIntracellular), "Intracellular bound samples should remain fixed-step without IS");
	expect(!should_use_variable_time_steps(true, false, ProtonSampleRegion::BoundExtracellular), "Extracellular bound samples should remain fixed-step without IS");

	expect(!should_use_field_grid_for_step(true, true), "Variable substeps must override an enabled field grid");
	expect(!should_use_field_grid_for_step(false, true), "Variable substeps must remain direct when the field grid is disabled");
	expect(should_use_field_grid_for_step(true, false), "Normal steps should resume using an enabled field grid");
	expect(!should_use_field_grid_for_step(false, false), "Normal steps must remain direct when the field grid is disabled");
}

void test_surrounding_sources_drive_variable_step_proximity() {
	const SimMatrix sim(
		0.0, 0.0, 0.0, 20.0, 22.0, 24.0, 0, 0, 0);
	const std::vector<Ionp> physical_ionps = {
		Ionp(-9.5, 0.0, 0.0, 0.0, 0.1, 0.0)};
	const Point3D edge_position(9.8, 0.0, 0.0);
	const IonpSpatialIndex physical_index(physical_ionps);
	expect(!has_near_ionp_for_variable_step(
			   edge_position,
			   physical_ionps,
			   physical_index,
			   IonpCheckMode::Spatial),
		   "Central-only proximity unexpectedly saw an opposite-face IONP");

	const std::vector<Ionp> field_sources =
		SimSpace::surroundingSpace(physical_ionps, sim);
	const IonpSpatialIndex field_index(field_sources);
	expect(has_near_ionp_for_variable_step(
			   edge_position,
			   field_sources,
			   field_index,
			   IonpCheckMode::Spatial),
		   "Variable-step proximity missed the neighboring voxel image");
}

void test_surrounding_sources_move_together_in_direct_field() {
	ImportanceSamplingConfig config;
	const SimMatrix sim(
		0.0, 0.0, 0.0, 20.0, 22.0, 24.0, 0, 0, 0);
	const Point3D cell_origin;
	MRsequence sequence = make_sequence(
		sim, 0, 0, 0.0, 0.01, 0.1, 0.0, cell_origin, 0.0, config);
	const std::vector<Ionp> physical_ionps = {
		Ionp(-9.5, 0.7, -1.1, 0.0, 0.1, 0.0)};
	const std::vector<Ionp> field_sources =
		SimSpace::surroundingSpace(physical_ionps, sim);
	const Point3D translation(1.25, -0.75, 0.5);
	const Point3D query(9.7, -0.1, -0.4);
	const std::vector<double> no_grid;

	const double actual = sequence.computeBFieldAt(
		query, field_sources, no_grid, false, &translation);
	double expected = 0.0;
	for (const Ionp& source : field_sources) {
		expected += sim.calculateBfield(
			query, source.xt() + translation, source.radius());
	}
	expect(approx_equal(actual, expected, 1e-20),
		   "Moving direct-field evaluation did not translate every magnetic image");
}

void test_superparticle_configuration_invariants() {
	ImportanceSamplingConfig config;
	const SimMatrix sim(
		0.0, 0.0, 0.0, 20.0, 22.0, 24.0, 0, 0, 0);
	const Point3D cell_origin;
	std::vector<Proton> protons;
	std::vector<Ionp> ionps;
	std::vector<std::complex<double>> signal;
	std::vector<double> timeV;
	std::vector<std::complex<double>> signal_mse;
	std::vector<double> timeV_mse;
	std::vector<Point3D> cell_origins;

	MRsequence disabled_surrounding = make_sequence(
		sim, 0, 0, 0.0, 0.01, 0.1, 0.0, cell_origin, 0.0, config);
	disabled_surrounding.set_surrounding_voxel_superparticles(1);
	bool rejected_without_surrounding = false;
	try {
		disabled_surrounding.Sequence(
			protons, ionps, signal, timeV, signal_mse, timeV_mse, cell_origins);
	}
	catch (const std::invalid_argument& error) {
		rejected_without_surrounding =
			std::string(error.what()).find("Surrounding_voxels = Yes") !=
			std::string::npos;
	}
	expect(rejected_without_surrounding,
		   "A superparticle cap without surrounding voxels was not rejected");

	MRsequence fixed_moment = make_sequence(
		sim,
		0,
		0,
		0.0,
		0.01,
		0.1,
		0.0,
		cell_origin,
		0.0,
		config,
		false,
		IonpCheckMode::Spatial,
		1.0);
	fixed_moment.set_surrounding_voxels(true);
	fixed_moment.set_surrounding_voxel_superparticles(1);
	bool rejected_fixed_moment = false;
	try {
		fixed_moment.Sequence(
			protons, ionps, signal, timeV, signal_mse, timeV_mse, cell_origins);
	}
	catch (const std::invalid_argument& error) {
		rejected_fixed_moment =
			std::string(error.what()).find("Magnetic Moment = 0") !=
			std::string::npos;
	}
	expect(rejected_fixed_moment,
		   "Superparticles with a fixed magnetic moment were not rejected");
}

void test_signal_aggregation_helper() {
	ImportanceSamplingConfig config;
	config.enabled = true;
	config.protons_inside = 1;
	config.protons_near = 1;
	config.protons_far = 1;
	config.weight_inside = 1.0;
	config.weight_near = 2.0;
	config.weight_far = 3.0;
	config.weight_bound = 4.0;

	const auto normalized_weights = normalized_sample_weights(config, true);
	std::array<std::size_t, kSampleRegionCount> sample_counts{};
	sample_counts[sample_region_index(ProtonSampleRegion::Inside)] = 1;
	sample_counts[sample_region_index(ProtonSampleRegion::Near)] = 1;
	sample_counts[sample_region_index(ProtonSampleRegion::Far)] = 1;
	sample_counts[sample_region_index(ProtonSampleRegion::Bound)] = 1;
	const std::vector<ProtonSampleRegion> sample_regions = {
		ProtonSampleRegion::Inside,
		ProtonSampleRegion::Near,
		ProtonSampleRegion::Far,
		ProtonSampleRegion::Bound
	};
	const std::vector<std::vector<double>> phases = {
		{0.0, 0.0},
		{kPi / 2.0, kPi / 2.0},
		{kPi, kPi},
		{-kPi / 2.0, -kPi / 2.0}
	};

	const WeightedSignalAggregation aggregation =
		aggregate_weighted_signals(phases, sample_regions, sample_counts, normalized_weights);

	const std::size_t inside_index =
		sample_region_index(ProtonSampleRegion::Inside);
	const std::size_t near_index =
		sample_region_index(ProtonSampleRegion::Near);
	const std::size_t far_index =
		sample_region_index(ProtonSampleRegion::Far);
	const std::size_t bound_index =
		sample_region_index(ProtonSampleRegion::Bound);
	const std::complex<double> expected =
		normalized_weights[inside_index] * std::complex<double>(1.0, 0.0) +
		normalized_weights[near_index] * std::complex<double>(0.0, 1.0) +
		normalized_weights[far_index] * std::complex<double>(-1.0, 0.0) +
		normalized_weights[bound_index] * std::complex<double>(0.0, -1.0);

	expect(aggregation.combined_signal.size() == 2, "Unexpected aggregation output size");
	expect_complex_close(aggregation.combined_signal[0], expected);
	expect_complex_close(aggregation.combined_signal[1], expected);
}

void test_overlapping_proposal_names_and_counts() {
	const std::array<const char*, kFreeSamplingProposalCount> expected_names = {
		"uniform",
		"cell_inside",
		"cell_near",
		"cell_far",
		"intracellular_ionp_near",
		"intracellular_ionp_intermediate",
		"intracellular_ionp_far",
		"extracellular_ionp_near",
		"extracellular_ionp_intermediate",
		"extracellular_ionp_far"
	};

	ImportanceSamplingConfig config;
	config.enabled = true;
	config.mode = ImportanceSamplingMode::OverlappingProposal;
	int expected_total = 0;
	for (std::size_t index = 0; index < kFreeSamplingProposalCount; ++index) {
		const FreeSamplingProposal proposal =
			static_cast<FreeSamplingProposal>(index);
		expect(free_sampling_proposal_index(proposal) == index,
			   "Overlapping proposal index is not stable");
		expect(free_sampling_proposal_name(proposal) == expected_names[index],
			   "Overlapping proposal name is unclear or unstable");
		config.proposal_protons[index] = static_cast<int>(index + 1);
		expected_total += static_cast<int>(index + 1);
	}

	expect(config.overlapping_proposals_enabled(),
		   "Configured overlapping proposals were not enabled");
	expect(config.total_proposal_protons() == expected_total,
		   "Overlapping proposal counts did not sum correctly");
	expect(config.total_free_protons() == expected_total,
		   "Total free count did not select overlapping proposal counts");

	ImportanceSamplingConfig overflowing;
	overflowing.mode = ImportanceSamplingMode::OverlappingProposal;
	overflowing.proposal_protons[0] = std::numeric_limits<int>::max();
	overflowing.proposal_protons[1] = 1;
	bool overflow_threw = false;
	try {
		(void) overflowing.total_proposal_protons();
	}
	catch (const std::runtime_error&) {
		overflow_threw = true;
	}
	expect(overflow_threw,
		   "An overflowing total proposal count was accepted");

	ImportanceSamplingConfig contradictory = config;
	contradictory.joint_regions_configured = true;
	bool contradictory_mode_threw = false;
	try {
		(void) contradictory.total_free_protons();
	}
	catch (const std::runtime_error&) {
		contradictory_mode_threw = true;
	}
	expect(contradictory_mode_threw,
		   "Contradictory legacy and overlapping mode state was accepted");

	ImportanceSamplingConfig missing_mode;
	missing_mode.proposal_protons[free_sampling_proposal_index(
		FreeSamplingProposal::Uniform)] = 1;
	bool missing_mode_threw = false;
	try {
		(void) missing_mode.total_free_protons();
	}
	catch (const std::runtime_error&) {
		missing_mode_threw = true;
	}
	expect(missing_mode_threw,
		   "Proposal counts without OverlappingProposal mode were accepted");
}

void test_overlapping_effective_count_redistribution() {
	std::array<int, kFreeSamplingProposalCount> requested = {
		5, 3, 4, 6, 2, 7, 8, 9, 10, 11
	};
	std::array<double, kFreeSamplingProposalCount> fractions = {
		1.0, 0.0, 0.1, 0.9, 0.02, 0.0, 0.98, 0.0, 0.03, 0.97
	};

	const EffectiveProposalCounts effective =
		resolve_effective_proposal_counts(requested, fractions);
	const std::size_t requested_total = std::accumulate(
		requested.begin(), requested.end(), std::size_t{0},
		[](std::size_t total, int count) {
			return total + static_cast<std::size_t>(count);
		});
	const std::size_t effective_total = std::accumulate(
		effective.counts.begin(), effective.counts.end(), std::size_t{0});

	expect(effective_total == requested_total,
		   "Zero-support redistribution changed the total trajectory count");
	expect(effective.redistributed_to_uniform == 19,
		   "Unexpected number of zero-support samples was redistributed");
	expect(effective.counts[free_sampling_proposal_index(
		   FreeSamplingProposal::Uniform)] == 24,
		   "Zero-support samples were not moved to the uniform proposal");
	expect(effective.counts[free_sampling_proposal_index(
		   FreeSamplingProposal::CellInside)] == 0 &&
		   effective.counts[free_sampling_proposal_index(
		   FreeSamplingProposal::IntracellularIonpIntermediate)] == 0 &&
		   effective.counts[free_sampling_proposal_index(
		   FreeSamplingProposal::ExtracellularIonpNear)] == 0,
		   "A zero-support proposal retained samples");
	expect(effective.redistributed[free_sampling_proposal_index(
		   FreeSamplingProposal::CellInside)] &&
		   effective.redistributed[free_sampling_proposal_index(
		   FreeSamplingProposal::IntracellularIonpIntermediate)] &&
		   effective.redistributed[free_sampling_proposal_index(
		   FreeSamplingProposal::ExtracellularIonpNear)],
		   "Redistributed proposals were not identified for diagnostics");
	expect(effective.counts[free_sampling_proposal_index(
		   FreeSamplingProposal::IntracellularIonpFar)] == 8 &&
		   effective.counts[free_sampling_proposal_index(
		   FreeSamplingProposal::ExtracellularIonpFar)] == 11,
		   "Supported far proposals were incorrectly redistributed");

	bool missing_uniform_threw = false;
	try {
		auto no_uniform = requested;
		no_uniform[free_sampling_proposal_index(
			FreeSamplingProposal::Uniform)] = 0;
		(void) resolve_effective_proposal_counts(no_uniform, fractions);
	}
	catch (const std::runtime_error&) {
		missing_uniform_threw = true;
	}
	expect(missing_uniform_threw,
		   "Overlapping proposals accepted a missing uniform safety proposal");

	bool negative_count_threw = false;
	try {
		auto negative = requested;
		negative[free_sampling_proposal_index(
			FreeSamplingProposal::CellNear)] = -1;
		(void) resolve_effective_proposal_counts(negative, fractions);
	}
	catch (const std::runtime_error&) {
		negative_count_threw = true;
	}
	expect(negative_count_threw,
		   "Overlapping proposals accepted a negative proton count");
}

void test_overlapping_mis_weight_formula() {
	std::array<std::size_t, kFreeSamplingProposalCount> counts{};
	counts[free_sampling_proposal_index(FreeSamplingProposal::Uniform)] = 10;
	counts[free_sampling_proposal_index(FreeSamplingProposal::CellInside)] = 20;
	counts[free_sampling_proposal_index(
		FreeSamplingProposal::IntracellularIonpNear)] = 30;
	counts[free_sampling_proposal_index(
		FreeSamplingProposal::ExtracellularIonpFar)] = 40;

	std::array<double, kFreeSamplingProposalCount> fractions{};
	fractions[free_sampling_proposal_index(FreeSamplingProposal::Uniform)] = 1.0;
	fractions[free_sampling_proposal_index(FreeSamplingProposal::CellInside)] = 0.2;
	fractions[free_sampling_proposal_index(
		FreeSamplingProposal::IntracellularIonpNear)] = 0.05;
	fractions[free_sampling_proposal_index(
		FreeSamplingProposal::ExtracellularIonpFar)] = 0.8;

	const JointSampleRegion overlap{
		CellSampleRegion::Inside,
		IonpDistanceRegion::Near,
		IonpDistanceRegion::Far
	};
	const JointSampleRegion background{
		CellSampleRegion::Far,
		IonpDistanceRegion::Far,
		IonpDistanceRegion::Near
	};
	const JointSampleRegion cell_only{
		CellSampleRegion::Inside,
		IonpDistanceRegion::Far,
		IonpDistanceRegion::Near
	};

	const double overlap_expected =
		1.0 / (10.0 / 1.0 + 20.0 / 0.2 + 30.0 / 0.05 + 40.0 / 0.8);
	const double background_expected = 1.0 / (10.0 / 1.0);
	const double cell_only_expected =
		1.0 / (10.0 / 1.0 + 20.0 / 0.2);

	expect(approx_equal(
		   overlapping_mis_raw_weight(overlap, counts, fractions),
		   overlap_expected),
		   "MIS weight did not include every proposal covering an overlap point");
	expect(approx_equal(
		   overlapping_mis_raw_weight(background, counts, fractions),
		   background_expected),
		   "Uniform-only point received an incorrect MIS weight");
	expect(approx_equal(
		   overlapping_mis_raw_weight(cell_only, counts, fractions),
		   cell_only_expected),
		   "Single-marginal point received an incorrect MIS weight");

	const std::vector<JointSampleRegion> memberships = {
		overlap, background, cell_only
	};
	const std::vector<double> raw = overlapping_mis_raw_weights(
		memberships, counts, fractions);
	const std::vector<double> normalized =
		overlapping_mis_self_normalized_weights(
			memberships, counts, fractions);
	expect(raw.size() == memberships.size() &&
		   normalized.size() == memberships.size(),
		   "Vector MIS weighting changed the sample count");
	const double raw_total = overlap_expected + background_expected +
		cell_only_expected;
	expect(approx_equal(normalized[0], overlap_expected / raw_total) &&
		   approx_equal(normalized[1], background_expected / raw_total) &&
		   approx_equal(normalized[2], cell_only_expected / raw_total),
		   "Self-normalized MIS weights are incorrect");
	expect(approx_equal(
		   std::accumulate(normalized.begin(), normalized.end(), 0.0), 1.0),
		   "Self-normalized MIS weights do not sum to one");
}

void test_overlapping_signal_aggregation_helper() {
	const std::vector<std::vector<double>> phases = {
		{0.0, 0.0},
		{kPi / 2.0, kPi / 2.0},
		{kPi, kPi},
		{-kPi / 2.0, -kPi / 2.0}
	};
	const std::vector<FreeSamplingProposal> proposals = {
		FreeSamplingProposal::Uniform,
		FreeSamplingProposal::CellInside,
		FreeSamplingProposal::IntracellularIonpNear,
		FreeSamplingProposal::IntracellularIonpNear
	};
	const std::vector<double> normalized_weights = {0.1, 0.2, 0.3, 0.4};

	const OverlappingWeightedSignalAggregation aggregation =
		aggregate_overlapping_weighted_signals(
			phases, proposals, normalized_weights);
	const std::complex<double> expected(-0.2, -0.2);
	expect(aggregation.combined_signal.size() == 2,
		   "Overlapping aggregation returned the wrong time dimension");
	expect_complex_close(aggregation.combined_signal[0], expected);
	expect_complex_close(aggregation.combined_signal[1], expected);

	const std::size_t uniform = free_sampling_proposal_index(
		FreeSamplingProposal::Uniform);
	const std::size_t cell_inside = free_sampling_proposal_index(
		FreeSamplingProposal::CellInside);
	const std::size_t intracellular_near = free_sampling_proposal_index(
		FreeSamplingProposal::IntracellularIonpNear);
	expect_complex_close(
		aggregation.proposal_means[uniform][0],
		std::complex<double>(1.0, 0.0));
	expect_complex_close(
		aggregation.proposal_means[cell_inside][0],
		std::complex<double>(0.0, 1.0));
	expect_complex_close(
		aggregation.proposal_means[intracellular_near][0],
		std::complex<double>(-0.5, -0.5));
	expect_complex_close(
		aggregation.proposal_contributions[uniform][0],
		std::complex<double>(0.1, 0.0));
	expect_complex_close(
		aggregation.proposal_contributions[cell_inside][0],
		std::complex<double>(0.0, 0.2));
	expect_complex_close(
		aggregation.proposal_contributions[intracellular_near][0],
		std::complex<double>(-0.3, -0.4));

	const std::vector<std::vector<double>> zero_phases(4,
		std::vector<double>{0.0});
	const OverlappingWeightedSignalAggregation initial =
		aggregate_overlapping_weighted_signals(
			zero_phases, proposals, normalized_weights);
	expect(initial.combined_signal.size() == 1 &&
		   approx_equal(initial.combined_signal.front().real(), 1.0) &&
		   approx_equal(initial.combined_signal.front().imag(), 0.0),
		   "Self-normalized overlapping signal is not one at initial time");
}

void test_legacy_automatic_weights_survive_rematerialization() {
	ImportanceSamplingConfig config;
	config.enabled = true;
	config.protons_inside = 10000;
	config.protons_near = 6000;
	config.protons_far = 4000;
	config.materialize_legacy_regions();

	const std::size_t inside_index =
		sample_region_index(ProtonSampleRegion::Inside);
	const std::size_t near_index =
		sample_region_index(ProtonSampleRegion::Near);
	const std::size_t far_index =
		sample_region_index(ProtonSampleRegion::Far);
	config.joint_weights[inside_index] = 0.001;
	config.joint_weights[near_index] = 0.039;
	config.joint_weights[far_index] = 0.96;

	// MRsequence materializes its copy again after main has estimated weights.
	config.materialize_legacy_regions();
	const auto weights = normalized_sample_weights(config, 0, 0);

	expect(approx_equal(weights[inside_index], 0.001),
		   "Rematerialization discarded the automatic inside weight");
	expect(approx_equal(weights[near_index], 0.039),
		   "Rematerialization discarded the automatic near weight");
	expect(approx_equal(weights[far_index], 0.96),
		   "Rematerialization discarded the automatic far weight");
}

void test_all_joint_and_bound_buckets() {
	ImportanceSamplingConfig config;
	config.enabled = true;
	config.joint_regions_configured = true;
	config.split_bound_configured = true;
	std::vector<ProtonSampleRegion> sample_regions;
	std::vector<std::vector<double>> phases;
	for (std::size_t index = 0; index < kFreeSampleRegionCount; ++index) {
		config.joint_protons[index] = 1;
		config.joint_weights[index] = static_cast<double>(index + 1);
		sample_regions.push_back(static_cast<ProtonSampleRegion>(index));
		phases.push_back({0.0});
		const JointSampleRegion decoded =
			decode_joint_sample_region(
				static_cast<ProtonSampleRegion>(index));
		expect(joint_sample_region_index(
				   decoded.cell,
				   decoded.intracellular_ionp,
				   decoded.extracellular_ionp) == index,
			   "Joint sample-region encoding is not reversible");
	}
	config.bound_protons_intracellular = 1;
	config.bound_protons_extracellular = 1;
	config.weight_bound_intracellular = 28.0;
	config.weight_bound_extracellular = 29.0;
	sample_regions.push_back(ProtonSampleRegion::BoundIntracellular);
	sample_regions.push_back(ProtonSampleRegion::BoundExtracellular);
	phases.push_back({0.0});
	phases.push_back({0.0});

	const auto counts = sample_region_counts(config, 1, 1);
	const auto weights = normalized_sample_weights(config, 1, 1);
	const WeightedSignalAggregation aggregation =
		aggregate_weighted_signals(
			phases, sample_regions, counts, weights);

	expect(config.total_free_protons() == 27,
		   "All joint buckets did not contribute to the effective free count");
	expect(sample_regions.size() == kSampleRegionCount,
		   "The full joint/bound label set has the wrong size");
	expect(aggregation.combined_signal.size() == 1,
		   "Full joint aggregation returned the wrong time dimension");
	expect_complex_close(
		aggregation.combined_signal.front(),
		std::complex<double>(1.0, 0.0));
}

void test_clear_region_names() {
	const ProtonSampleRegion free_region = make_joint_sample_region(
		CellSampleRegion::Inside,
		IonpDistanceRegion::Intermediate,
		IonpDistanceRegion::Far);
	expect(
		sample_region_name(free_region) ==
			"cell_inside__intracellular_ionp_intermediate__extracellular_ionp_far",
		"Free-region name does not clearly separate the three sampling axes");
	expect(
		sample_region_name(ProtonSampleRegion::BoundIntracellular) ==
			"bound_intracellular_ionp",
		"Intracellular bound-region name is unclear");
	expect(
		sample_region_name(ProtonSampleRegion::BoundExtracellular) ==
			"bound_extracellular_ionp",
		"Extracellular bound-region name is unclear");
}

void test_realized_outer_radius_semantics() {
	expect(approx_equal(ionp_shell_thickness(1.0, 0.0), 0.0),
		   "Coating=0 should produce no stored shell");
	expect(approx_equal(ionp_shell_thickness(1.0, 1.5), 0.5),
		   "Configured coating radius did not produce the expected shell thickness");
	expect(approx_equal(ionp_shell_thickness(1.5, 1.5), 0.0),
		   "Core equal to coating radius should have no shell");
	expect(approx_equal(ionp_shell_thickness(2.0, 1.5), 0.0),
		   "Core larger than coating radius should clamp the shell to zero");

	const Ionp coated(
		0.0, 0.0, 0.0, 0.0, 1.0,
		ionp_shell_thickness(1.0, 1.5));
	const Ionp oversized_core(
		0.0, 0.0, 0.0, 0.0, 2.0,
		ionp_shell_thickness(2.0, 1.5));
	expect(approx_equal(coated.outer_radius(), 1.5),
		   "Coated particle outer radius is not core plus stored shell");
	expect(approx_equal(oversized_core.outer_radius(), 2.0),
		   "Oversized core should remain the realized outer radius");

	const std::array<double, 4> distributed_cores = {0.5, 1.0, 1.5, 2.0};
	const std::array<double, 4> expected_outer = {1.5, 1.5, 1.5, 2.0};
	for (std::size_t index = 0; index < distributed_cores.size(); ++index) {
		const Ionp particle(
			0.0,
			0.0,
			0.0,
			0.0,
			distributed_cores[index],
			ionp_shell_thickness(distributed_cores[index], 1.5));
		expect(approx_equal(particle.outer_radius(), expected_outer[index]),
			   "Distributed core realization produced the wrong outer radius");
	}
}

void test_ionp_population_and_joint_region_classification() {
	const Point3D cell_origin(0.0, 0.0, 0.0);
	const std::vector<Point3D> cells = {cell_origin};
	const std::vector<Ionp> ionps = {
		Ionp(0.0, 0.0, 0.0, 0.0, 0.5, 0.5),
		Ionp(8.0, 0.0, 0.0, 0.0, 0.5, 0.5)
	};
	const IonpPopulationMap populations =
		classify_ionp_populations(ionps, cells, cell_origin, 5.0);
	expect(populations.intracellular_indices.size() == 1 &&
		   populations.intracellular_indices.front() == 0,
		   "Intracellular IONP population classification failed");
	expect(populations.extracellular_indices.size() == 1 &&
		   populations.extracellular_indices.front() == 1,
		   "Extracellular IONP population classification failed");
	expect(populations.coated_intracellular_count == 1 &&
		   populations.coated_extracellular_count == 1,
		   "Coated population counts are incorrect");

	ImportanceSamplingConfig config;
	config.enabled = true;
	config.cell_near_cutoff = 7.0;
	config.intracellular_ionp_near_cutoff = 3.0;
	config.intracellular_ionp_far_cutoff = 5.0;
	config.extracellular_ionp_near_cutoff = 3.0;
	config.extracellular_ionp_far_cutoff = 5.0;

	ProtonSampleRegion region;
	expect(classify_joint_sample_region(
			   Point3D(4.0, 0.0, 0.0),
			   cells,
			   cell_origin,
			   5.0,
			   ionps,
			   populations,
			   config,
			   region),
		   "Valid overlap point was rejected");
	expect(region == make_joint_sample_region(
					 CellSampleRegion::Inside,
					 IonpDistanceRegion::Near,
					 IonpDistanceRegion::Near),
		   "Simultaneous cell/internal/external membership was not represented by one joint tuple");
	const IonpSpatialIndex particle_index(ionps);
	const std::vector<unsigned int> candidates =
		particle_index.point_candidates(Point3D(4.0, 0.0, 0.0), 6.0);
	ProtonSampleRegion candidate_region;
	expect(classify_joint_sample_region_from_ionp_candidates(
			   Point3D(4.0, 0.0, 0.0),
			   cells,
			   cell_origin,
			   5.0,
			   ionps,
			   populations,
			   candidates,
			   config,
			   candidate_region) &&
		   candidate_region == region,
		   "Spatial-candidate and full joint classifiers disagreed");

	expect(classify_joint_sample_region(
			   Point3D(6.0, 0.0, 0.0),
			   cells,
			   cell_origin,
			   5.0,
			   ionps,
			   populations,
			   config,
			   region),
		   "Valid near-cell overlap point was rejected");
	expect(region == make_joint_sample_region(
					 CellSampleRegion::Near,
					 IonpDistanceRegion::Intermediate,
					 IonpDistanceRegion::Near),
		   "Joint boundary classification produced the wrong tuple");

	expect(classify_joint_sample_region(
			   Point3D(6.1, 0.0, 0.0),
			   cells,
			   cell_origin,
			   5.0,
			   ionps,
			   populations,
			   config,
			   region) &&
		   region == make_joint_sample_region(
						 CellSampleRegion::Near,
						 IonpDistanceRegion::Far,
						 IonpDistanceRegion::Near),
		   "Distance immediately beyond Y should classify as far");

	expect(!classify_joint_sample_region(
			   Point3D(0.5, 0.0, 0.0),
			   cells,
			   cell_origin,
			   5.0,
			   ionps,
			   populations,
			   config,
			   region),
		   "Point inside an IONP outer radius should be excluded");

	bool straddling_threw = false;
	try {
		const std::vector<Ionp> straddling = {
			Ionp(5.0, 0.0, 0.0, 0.0, 0.5, 0.5)
		};
		(void) classify_ionp_populations(
			straddling, cells, cell_origin, 5.0);
	}
	catch (const std::runtime_error&) {
		straddling_threw = true;
	}
	expect(straddling_threw,
		   "Membrane-straddling IONP should be rejected");

	const IonpPopulationMap no_cell_populations =
		classify_ionp_populations(ionps, {}, Point3D(), 0.0);
	expect(no_cell_populations.intracellular_indices.empty() &&
		   no_cell_populations.extracellular_indices.size() == ionps.size(),
		   "No-cell geometry should classify every IONP as extracellular");

	bool overlapping_cell_straddle_threw = false;
	try {
		const std::vector<Point3D> overlapping_cells = {
			Point3D(0.0, 0.0, 0.0),
			Point3D(4.5, 0.0, 0.0)
		};
		const std::vector<Ionp> ambiguous = {
			Ionp(0.0, 0.0, 0.0, 0.0, 0.5, 0.5)
		};
		(void) classify_ionp_populations(
			ambiguous, overlapping_cells, Point3D(), 5.0);
	}
	catch (const std::runtime_error&) {
		overlapping_cell_straddle_threw = true;
	}
	expect(overlapping_cell_straddle_threw,
		   "An IONP contained by one cell but straddling another should be rejected");
}

void test_joint_particle_auto_volume_weights() {
	rng_config::set_seed(8765);

	ImportanceSamplingConfig config;
	config.enabled = true;
	config.joint_regions_configured = true;
	config.extracellular_ionp_near_cutoff = 1.0;
	config.extracellular_ionp_far_cutoff = 2.0;
	const std::size_t near_index = joint_sample_region_index(
		CellSampleRegion::Far,
		IonpDistanceRegion::Far,
		IonpDistanceRegion::Near);
	const std::size_t intermediate_index = joint_sample_region_index(
		CellSampleRegion::Far,
		IonpDistanceRegion::Far,
		IonpDistanceRegion::Intermediate);
	const std::size_t far_index = joint_sample_region_index(
		CellSampleRegion::Far,
		IonpDistanceRegion::Far,
		IonpDistanceRegion::Far);
	config.joint_protons[near_index] = 1;
	config.joint_protons[intermediate_index] = 1;
	config.joint_protons[far_index] = 1;

	const std::vector<Ionp> ionps = {
		Ionp(0.0, 0.0, 0.0, 0.0, 1.0, 0.0)
	};
	const IonpPopulationMap populations =
		classify_ionp_populations(ionps, {}, Point3D(), 0.0);
	const auto weights = estimate_active_joint_region_volume_weights(
		config,
		{},
		Point3D(),
		0.0,
		ionps,
		populations,
		20.0,
		20.0,
		20.0,
		300000);

	const double box_volume = 8000.0;
	const double near_exact =
		(4.0 / 3.0) * kPi * (std::pow(2.0, 3) - std::pow(1.0, 3)) /
		box_volume;
	const double intermediate_exact =
		(4.0 / 3.0) * kPi * (std::pow(3.0, 3) - std::pow(2.0, 3)) /
		box_volume;
	const double far_exact =
		1.0 - (4.0 / 3.0) * kPi * std::pow(3.0, 3) / box_volume;
	expect(std::abs(weights[near_index] - near_exact) < 0.0004,
		   "Automatic near-IONP joint weight is inaccurate");
	expect(std::abs(weights[intermediate_index] - intermediate_exact) < 0.0007,
		   "Automatic intermediate-IONP joint weight is inaccurate");
	expect(std::abs(weights[far_index] - far_exact) < 0.008,
		   "Automatic far-water joint weight is inaccurate");

	const std::vector<Ionp> overlapping_ionps = {
		Ionp(-1.5, 0.0, 0.0, 0.0, 1.0, 0.0),
		Ionp(1.5, 0.0, 0.0, 0.0, 1.0, 0.0)
	};
	const IonpPopulationMap overlapping_populations =
		classify_ionp_populations(
			overlapping_ionps, {}, Point3D(), 0.0);
	ImportanceSamplingConfig overlap_config;
	overlap_config.enabled = true;
	overlap_config.joint_regions_configured = true;
	overlap_config.extracellular_ionp_near_cutoff = 1.0;
	overlap_config.extracellular_ionp_far_cutoff = 2.0;
	overlap_config.joint_protons[near_index] = 1;
	const auto overlap_weights =
		estimate_active_joint_region_volume_weights(
			overlap_config,
			{},
			Point3D(),
			0.0,
			overlapping_ionps,
			overlapping_populations,
			20.0,
			20.0,
			20.0,
			300000);
	const double overlap_lens = 11.0 * kPi / 12.0;
	const double overlap_union_shell =
		2.0 * (4.0 / 3.0) * kPi *
			(std::pow(2.0, 3) - std::pow(1.0, 3)) -
		overlap_lens;
	expect(std::abs(
			   overlap_weights[near_index] -
			   overlap_union_shell / box_volume) < 0.0005,
		   "Overlapping-shell automatic weight did not correct proposal coverage");
}

void test_overlapping_proposal_volume_estimation() {
	ImportanceSamplingConfig config;
	config.enabled = true;
	config.mode = ImportanceSamplingMode::OverlappingProposal;
	config.extracellular_ionp_near_cutoff = 1.0;
	config.extracellular_ionp_far_cutoff = 2.0;
	config.proposal_protons[free_sampling_proposal_index(
		FreeSamplingProposal::Uniform)] = 1;
	config.proposal_protons[free_sampling_proposal_index(
		FreeSamplingProposal::ExtracellularIonpNear)] = 1;
	config.proposal_protons[free_sampling_proposal_index(
		FreeSamplingProposal::ExtracellularIonpIntermediate)] = 1;
	config.proposal_protons[free_sampling_proposal_index(
		FreeSamplingProposal::ExtracellularIonpFar)] = 1;

	const std::vector<Ionp> ionps = {
		Ionp(0.0, 0.0, 0.0, 0.0, 1.0, 0.0)
	};
	const IonpPopulationMap populations =
		classify_ionp_populations(ionps, {}, Point3D(), 0.0);
	const OverlappingProposalVolumeEstimate estimate =
		estimate_overlapping_proposal_volumes(
			config,
			{},
			Point3D(),
			0.0,
			ionps,
			populations,
			20.0,
			20.0,
			20.0,
			300000);

	const double box_volume = 8000.0;
	const double excluded_volume = (4.0 / 3.0) * kPi;
	const double accessible_volume = box_volume - excluded_volume;
	const double near_fraction =
		(4.0 / 3.0) * kPi * (std::pow(2.0, 3) - std::pow(1.0, 3)) /
		accessible_volume;
	const double intermediate_fraction =
		(4.0 / 3.0) * kPi * (std::pow(3.0, 3) - std::pow(2.0, 3)) /
		accessible_volume;
	const double far_fraction = 1.0 - near_fraction - intermediate_fraction;

	const auto fraction = [&](FreeSamplingProposal proposal) {
		return estimate.region_fractions[free_sampling_proposal_index(proposal)];
	};
	expect(estimate.pilot_samples == 300000,
		   "Overlapping volume estimator did not report its pilot size");
	expect(std::abs(estimate.accessible_volume - accessible_volume) < 35.0,
		   "Accessible free-water volume estimate is inaccurate");
	expect(approx_equal(fraction(FreeSamplingProposal::Uniform), 1.0),
		   "Uniform proposal does not cover all accessible water");
	expect(approx_equal(fraction(FreeSamplingProposal::CellInside), 0.0) &&
		   approx_equal(fraction(FreeSamplingProposal::CellNear), 0.0) &&
		   approx_equal(fraction(FreeSamplingProposal::CellFar), 1.0),
		   "No-cell geometry produced incorrect cell proposal volumes");
	expect(approx_equal(
		   fraction(FreeSamplingProposal::IntracellularIonpNear), 0.0) &&
		   approx_equal(
		   fraction(FreeSamplingProposal::IntracellularIonpIntermediate), 0.0) &&
		   approx_equal(
		   fraction(FreeSamplingProposal::IntracellularIonpFar), 1.0),
		   "Absent intracellular population did not collapse to the far proposal");
	expect(std::abs(
		   fraction(FreeSamplingProposal::ExtracellularIonpNear) -
		   near_fraction) < 0.0005,
		   "Extracellular near marginal volume is inaccurate");
	expect(std::abs(
		   fraction(FreeSamplingProposal::ExtracellularIonpIntermediate) -
		   intermediate_fraction) < 0.0008,
		   "Extracellular intermediate marginal volume is inaccurate");
	expect(std::abs(
		   fraction(FreeSamplingProposal::ExtracellularIonpFar) -
		   far_fraction) < 0.003,
		   "Extracellular far marginal volume is inaccurate");
}

void test_joint_region_initialization() {
	rng_config::set_seed(7712);

	ImportanceSamplingConfig config;
	config.enabled = true;
	config.joint_regions_configured = true;
	config.cell_near_cutoff = 7.0;
	config.intracellular_ionp_near_cutoff = 5.0;
	config.intracellular_ionp_far_cutoff = 6.0;
	config.extracellular_ionp_near_cutoff = 5.0;
	config.extracellular_ionp_far_cutoff = 6.0;
	const ProtonSampleRegion target = make_joint_sample_region(
		CellSampleRegion::Inside,
		IonpDistanceRegion::Near,
		IonpDistanceRegion::Near);
	const std::size_t target_index = sample_region_index(target);
	config.joint_protons[target_index] = 8;
	config.joint_weights[target_index] = 1.0;

	SimMatrix sim(
		0.0, 0.0, 0.0, 20.0, 20.0, 20.0, 0, 0, 0);
	const Point3D cell_origin(0.0, 0.0, 0.0);
	MRsequence sequence = make_sequence(
		sim,
		config.total_free_protons(),
		0,
		0.0,
		0.01,
		0.5,
		0.0,
		cell_origin,
		5.0,
		config);

	std::vector<Proton> protons;
	std::vector<Ionp> ionps = {
		Ionp(0.0, 0.0, 0.0, 0.0, 0.5, 0.5),
		Ionp(7.0, 0.0, 0.0, 0.0, 0.5, 0.5)
	};
	std::vector<std::complex<double>> signal;
	std::vector<double> timeV;
	std::vector<std::complex<double>> signal_mse;
	std::vector<double> timeV_mse;
	std::vector<Point3D> cells = {cell_origin};
	ImportanceRegionSignals region_signals;
	sequence.Sequence(
		protons,
		ionps,
		signal,
		timeV,
		signal_mse,
		timeV_mse,
		cells,
		&region_signals);

	const IonpPopulationMap populations =
		classify_ionp_populations(ionps, cells, cell_origin, 5.0);
	expect(protons.size() == 8,
		   "Unexpected joint-region initialized proton count");
	for (const Proton& proton : protons) {
		ProtonSampleRegion classified;
		expect(classify_joint_sample_region(
				   proton.xt(),
				   cells,
				   cell_origin,
				   5.0,
				   ionps,
				   populations,
				   config,
				   classified) &&
			   classified == target,
			   "Proton was not initialized in its requested joint region");
	}
	expect(signal.size() == 1 && approx_equal(signal.front().real(), 1.0),
		   "Joint-region zero-time signal should be normalized");
	expect(region_signals[target_index].size() == 1 &&
		   approx_equal(region_signals[target_index].front().real(), 1.0),
		   "Requested joint-region signal was not returned");
}

void test_overlapping_proposal_initialization() {
	rng_config::set_seed(1984);

	ImportanceSamplingConfig config;
	config.enabled = true;
	config.mode = ImportanceSamplingMode::OverlappingProposal;
	config.cell_near_cutoff = 7.0;
	config.intracellular_ionp_near_cutoff = 1.0;
	config.intracellular_ionp_far_cutoff = 2.0;
	config.extracellular_ionp_near_cutoff = 1.0;
	config.extracellular_ionp_far_cutoff = 2.0;
	for (std::size_t index = 0; index < kFreeSamplingProposalCount; ++index) {
		config.proposal_protons[index] = 3;
	}

	SimMatrix sim(
		0.0, 0.0, 0.0, 20.0, 20.0, 20.0, 0, 0, 0);
	const Point3D cell_origin(0.0, 0.0, 0.0);
	MRsequence sequence = make_sequence(
		sim,
		config.total_free_protons(),
		0,
		0.0,
		0.01,
		0.5,
		0.0,
		cell_origin,
		5.0,
		config);

	std::vector<Proton> protons;
	std::vector<Ionp> ionps = {
		// Complete outer sphere is inside the cell.
		Ionp(4.0, 0.0, 0.0, 0.0, 0.5, 0.0),
		// Complete outer sphere is outside the cell; near shells overlap at
		// the membrane without either particle straddling it.
		Ionp(6.0, 0.0, 0.0, 0.0, 0.5, 0.0)
	};
	std::vector<std::complex<double>> signal;
	std::vector<double> timeV;
	std::vector<std::complex<double>> signal_mse;
	std::vector<double> timeV_mse;
	std::vector<Point3D> cells = {cell_origin};
	ImportanceProposalSignalDiagnostics diagnostics;

	sequence.Sequence(
		protons,
		ionps,
		signal,
		timeV,
		signal_mse,
		timeV_mse,
		cells,
		nullptr,
		nullptr,
		&diagnostics,
		nullptr);

	expect(protons.size() ==
		   static_cast<std::size_t>(config.total_free_protons()),
		   "Overlapping initialization changed the requested total count");
	const IonpPopulationMap populations =
		classify_ionp_populations(ionps, cells, cell_origin, 5.0);
	std::size_t proton_index = 0;
	for (std::size_t proposal_index = 0;
		 proposal_index < kFreeSamplingProposalCount;
		 ++proposal_index) {
		const FreeSamplingProposal proposal =
			static_cast<FreeSamplingProposal>(proposal_index);
		for (int sample = 0;
			 sample < config.proposal_protons[proposal_index];
			 ++sample, ++proton_index) {
			ProtonSampleRegion classified;
			expect(classify_joint_sample_region(
					   protons[proton_index].xt(),
					   cells,
					   cell_origin,
					   5.0,
					   ionps,
					   populations,
					   config,
					   classified),
				   "Overlapping proposal initialized inside an IONP outer radius");
			const JointSampleRegion membership =
				decode_joint_sample_region(classified);
			expect(proposal_contains_membership(proposal, membership),
				   "A proposal-generated proton does not satisfy its marginal region");
		}
	}
	expect(signal.size() == 1 &&
		   approx_equal(signal.front().real(), 1.0) &&
		   approx_equal(signal.front().imag(), 0.0),
		   "Overlapping proposal signal is not normalized at initial time");
	double contribution_sum = 0.;
	for (std::size_t proposal_index = 0;
		 proposal_index < kFreeSamplingProposalCount;
		 ++proposal_index) {
		expect(diagnostics.means[proposal_index].size() == 1 &&
			   approx_equal(diagnostics.means[proposal_index].front().real(), 1.0),
			   "Active proposal mean was not returned at initial time");
		expect(diagnostics.contributions[proposal_index].size() == 1,
			   "Active proposal MIS contribution was not returned");
		contribution_sum += diagnostics.contributions[proposal_index].front().real();
	}
	expect(approx_equal(contribution_sum, 1.0),
		   "Proposal MIS contributions do not reconstruct the free signal");
}

void test_overlapping_zero_support_falls_back_to_uniform() {
	rng_config::set_seed(1985);

	ImportanceSamplingConfig config;
	config.enabled = true;
	config.mode = ImportanceSamplingMode::OverlappingProposal;
	config.extracellular_ionp_near_cutoff = 1.0;
	config.extracellular_ionp_far_cutoff = 2.0;
	config.proposal_protons[free_sampling_proposal_index(
		FreeSamplingProposal::Uniform)] = 3;
	config.proposal_protons[free_sampling_proposal_index(
		FreeSamplingProposal::ExtracellularIonpNear)] = 4;

	SimMatrix sim(
		0.0, 0.0, 0.0, 20.0, 20.0, 20.0, 0, 0, 0);
	const Point3D cell_origin;
	MRsequence sequence = make_sequence(
		sim,
		config.total_free_protons(),
		0,
		0.0,
		0.01,
		0.0,
		0.0,
		cell_origin,
		0.0,
		config);

	std::vector<Proton> protons;
	std::vector<Ionp> ionps;
	std::vector<std::complex<double>> signal;
	std::vector<double> timeV;
	std::vector<std::complex<double>> signal_mse;
	std::vector<double> timeV_mse;
	std::vector<Point3D> cells;
	ImportanceProposalSignalDiagnostics diagnostics;
	sequence.Sequence(
		protons,
		ionps,
		signal,
		timeV,
		signal_mse,
		timeV_mse,
		cells,
		nullptr,
		nullptr,
		&diagnostics,
		nullptr);

	expect(protons.size() == 7,
		   "Zero-support fallback did not preserve the total proton budget");
	expect(signal.size() == 1 && approx_equal(signal.front().real(), 1.0),
		   "Zero-support fallback did not produce a normalized signal");
	const std::size_t uniform = free_sampling_proposal_index(
		FreeSamplingProposal::Uniform);
	const std::size_t extracellular_near = free_sampling_proposal_index(
		FreeSamplingProposal::ExtracellularIonpNear);
	expect(diagnostics.means[uniform].size() == 1 &&
		   diagnostics.means[extracellular_near].empty(),
		   "Zero-support samples were not reported under their effective uniform proposal");
}

void test_ionp_spatial_index_helper() {
	{
		const std::vector<Ionp> ionps;
		const IonpSpatialIndex index(ionps);
		expect(index.empty(), "Empty IONP index should report empty");
		expect(index.point_candidates(Point3D(0.0, 0.0, 0.0), 1.0).empty(), "Empty IONP index returned point candidates");
		expect(index.segment_candidates(Point3D(0.0, 0.0, 0.0), Point3D(1.0, 0.0, 0.0), 1.0).empty(), "Empty IONP index returned segment candidates");
	}

	const std::vector<Ionp> ionps = {
		Ionp(4.0, 0.0, 0.0, 0.0, 1.0, 0.5),
		Ionp(0.0, 0.0, 0.0, 0.0, 1.0, 0.5),
		Ionp(2.0, 0.0, 0.0, 0.0, 0.25, 0.25)
	};
	const IonpSpatialIndex index(ionps);

	expect(!index.empty(), "Non-empty IONP index should not report empty");
	expect(approx_equal(index.max_outer_radius(), 1.5), "Unexpected max outer radius");
	expect_indices_equal(index.point_candidates(Point3D(1.25, 0.0, 0.0), index.max_outer_radius()),
						 {1, 2},
						 "Point candidate order should follow original IONP order");
	expect(coating_index_for_position(Point3D(1.25, 0.0, 0.0), ionps, index, IonpCheckMode::Spatial) == 1, "Point inside coating was not identified");
	expect(coating_index_for_position(Point3D(0.5, 0.0, 0.0), ionps, index, IonpCheckMode::Spatial) == -1, "Core point should not be classified as coating");
	expect(coating_index_for_position(Point3D(10.0, 0.0, 0.0), ionps, index, IonpCheckMode::Spatial) == -1, "Outside point should not be classified as coating");

	expect_indices_equal(index.segment_candidates(Point3D(-2.0, 0.0, 0.0), Point3D(1.25, 0.0, 0.0), index.max_outer_radius()),
						 {1, 2},
						 "Segment crossing coating/core returned unexpected candidates");
	expect_indices_equal(index.segment_candidates(Point3D(3.25, 0.0, 0.0), Point3D(4.75, 0.0, 0.0), index.max_outer_radius()),
						 {0, 2},
						 "Segment crossing second coating returned unexpected candidates");

	expect(has_near_ionp_for_variable_step(Point3D(10.0, 0.0, 0.0), ionps, index, IonpCheckMode::Spatial), "Variable-step broad phase missed an exact near IONP");
	expect(!has_near_ionp_for_variable_step(Point3D(100.0, 0.0, 0.0), ionps, index, IonpCheckMode::Spatial), "Variable-step broad phase reported a far-away IONP");

	std::vector<Ionp> moved = {Ionp(100.0, 0.0, 0.0, 0.0, 1.0, 0.5)};
	IonpSpatialIndex moved_index(moved);
	expect(coating_index_for_position(Point3D(1.25, 0.0, 0.0), moved, moved_index, IonpCheckMode::Spatial) == -1, "Moved index should not contain stale old position");
	moved[0].setXt(Point3D(0.0, 0.0, 0.0));
	moved_index.rebuild(moved);
	expect(coating_index_for_position(Point3D(1.25, 0.0, 0.0), moved, moved_index, IonpCheckMode::Spatial) == 0, "Rebuilt index did not reflect moved IONP position");
}

void test_translated_ionp_spatial_index_helper() {
	const std::vector<Ionp> base_ionps = {
		Ionp(-2.0, 0.0, 0.0, 0.0, 0.5, 0.25),
		Ionp(1.0, 1.0, 0.0, 0.0, 0.25, 0.0),
		Ionp(4.0, -1.0, 0.5, 0.0, 0.75, 0.1)
	};
	const Point3D translation(7.25, -3.5, 2.0);
	const IonpSpatialIndex static_index(base_ionps);
	IonpSpatialIndex internally_translated_index(base_ionps);
	internally_translated_index.translate(
		translation.x(),
		translation.y(),
		translation.z());

	std::vector<Ionp> moved_ionps = base_ionps;
	for (std::size_t index = 0; index < moved_ionps.size(); ++index) {
		const Point3D center = base_ionps[index].xt();
		moved_ionps[index].setXt(
			center.x() + translation.x(),
			center.y() + translation.y(),
			center.z() + translation.z());
	}
	const IonpSpatialIndex rebuilt_index(moved_ionps);

	const Point3D point(translation.x() - 1.25,
						translation.y(),
						translation.z());
	std::vector<unsigned int> external_candidates;
	static_index.point_candidates(
		point,
		static_index.max_outer_radius(),
		external_candidates,
		translation.x(),
		translation.y(),
		translation.z());
	expect_indices_equal(
		external_candidates,
		rebuilt_index.point_candidates(point, rebuilt_index.max_outer_radius()),
		"Externally translated point query disagreed with a rebuilt index");
	expect_indices_equal(
		external_candidates,
		internally_translated_index.point_candidates(
			point,
			internally_translated_index.max_outer_radius()),
		"External and internal IONP translations disagreed for a point query");

	const Point3D segment_start(translation.x() - 4.0,
								translation.y(),
								translation.z());
	const Point3D segment_end(translation.x() + 5.0,
							  translation.y(),
							  translation.z());
	static_index.segment_candidates(
		segment_start,
		segment_end,
		static_index.max_outer_radius(),
		external_candidates,
		translation.x(),
		translation.y(),
		translation.z());
	expect_indices_equal(
		external_candidates,
		internally_translated_index.segment_candidates(
			segment_start,
			segment_end,
			internally_translated_index.max_outer_radius()),
		"External and internal IONP translations disagreed for a segment query");

	static_index.variable_step_candidates(
		point,
		external_candidates,
		translation.x(),
		translation.y(),
		translation.z());
	expect_indices_equal(
		external_candidates,
		internally_translated_index.variable_step_candidates(point),
		"External and internal IONP translations disagreed for a variable-step query");
}

void test_explicit_random_list_seed() {
	RandomList first_uniform(0.0, 1.0, false, 123456U);
	RandomList second_uniform(0.0, 1.0, false, 123456U);
	for (unsigned int draw = 0; draw < 10050; ++draw) {
		expect(first_uniform.getNumber() == second_uniform.getNumber(),
			   "Explicitly seeded uniform RandomList streams diverged");
	}

	RandomList first_gaussian(0.0, 1.0, true, 987654U);
	RandomList second_gaussian(0.0, 1.0, true, 987654U);
	for (unsigned int draw = 0; draw < 10050; ++draw) {
		expect(first_gaussian.getNumber() == second_gaussian.getNumber(),
			   "Explicitly seeded Gaussian RandomList streams diverged");
	}
}

void test_coating_diffusion_coefficient_helper() {
	const std::vector<Ionp> ionps = {Ionp(0.0, 0.0, 0.0, 0.0, 1.0, 0.5)};
	const IonpSpatialIndex ionp_index(ionps);
	const double coating_diffusion = 1e-10;
	const double cell_diffusion = 1e-9;

	expect(approx_equal(diffusion_coefficient_for_position(Point3D(1.25, 0.0, 0.0), ionps, ionp_index, IonpCheckMode::Spatial, coating_diffusion, false, cell_diffusion), coating_diffusion), "Coating proton did not use the configured diffusion coefficient");
	expect(approx_equal(diffusion_coefficient_for_position(Point3D(2.0, 0.0, 0.0), ionps, ionp_index, IonpCheckMode::Spatial, coating_diffusion, false, cell_diffusion), D_D_CONST), "Outside proton should use the default diffusion coefficient");
	expect(approx_equal(diffusion_coefficient_for_position(Point3D(2.0, 0.0, 0.0), ionps, ionp_index, IonpCheckMode::Spatial, coating_diffusion, true, cell_diffusion), cell_diffusion), "Inside-cell proton should use the configured cell diffusion coefficient");
	expect(approx_equal(diffusion_coefficient_for_position(Point3D(0.5, 0.0, 0.0), ionps, ionp_index, IonpCheckMode::Spatial, coating_diffusion, false, cell_diffusion), D_D_CONST), "Core positions should not be treated as coating positions");
	expect(approx_equal(diffusion_coefficient_for_position(Point3D(1.25, 0.0, 0.0), ionps, ionp_index, IonpCheckMode::Spatial, coating_diffusion, true, cell_diffusion), coating_diffusion), "Coating diffusion should take precedence over cell diffusion");
	expect(approx_equal(diffusion_coefficient_for_position(Point3D(1.25, 0.0, 0.0), ionps, ionp_index, IonpCheckMode::Spatial, 0.0, true, cell_diffusion), 0.0), "Explicit zero coating diffusion coefficient was not preserved");
	expect(approx_equal(diffusion_coefficient_for_position(Point3D(2.0, 0.0, 0.0), ionps, ionp_index, IonpCheckMode::Spatial, coating_diffusion, true, 0.0), 0.0), "Explicit zero cell diffusion coefficient was not preserved");
}

void test_ionp_check_modes_equivalent() {
	const std::vector<Ionp> ionps = {
		Ionp(100.0, 0.0, 0.0, 0.0, 1.0, 0.5),
		Ionp(0.0, 0.0, 0.0, 0.0, 1.0, 0.5)
	};
	const IonpSpatialIndex ionp_index(ionps);
	const Point3D coating_position(1.25, 0.0, 0.0);
	const Point3D far_position(50.0, 0.0, 0.0);
	const double coating_diffusion = 1e-10;

	expect(coating_index_for_position(coating_position, ionps, ionp_index, IonpCheckMode::Spatial) ==
		   coating_index_for_position(coating_position, ionps, ionp_index, IonpCheckMode::All),
		   "Spatial and all-IONP coating lookup disagreed");
	expect(approx_equal(diffusion_coefficient_for_position(coating_position, ionps, ionp_index, IonpCheckMode::Spatial, coating_diffusion, false, D_D_CONST),
						diffusion_coefficient_for_position(coating_position, ionps, ionp_index, IonpCheckMode::All, coating_diffusion, false, D_D_CONST)),
		   "Spatial and all-IONP coating diffusion disagreed");
	expect(has_near_ionp_for_variable_step(coating_position, ionps, ionp_index, IonpCheckMode::Spatial) ==
		   has_near_ionp_for_variable_step(coating_position, ionps, ionp_index, IonpCheckMode::All),
		   "Spatial and all-IONP near-IONP checks disagreed near target");
	expect(has_near_ionp_for_variable_step(far_position, ionps, ionp_index, IonpCheckMode::Spatial) ==
		   has_near_ionp_for_variable_step(far_position, ionps, ionp_index, IonpCheckMode::All),
		   "Spatial and all-IONP near-IONP checks disagreed far from target");

	expect_boundary_results_equal(
		run_boundary_check(IonpCheckMode::Spatial, ionps, Point3D(1.25, 0.0, 0.0), Point3D(0.75, 0.0, 0.0), false, 1.0),
		run_boundary_check(IonpCheckMode::All, ionps, Point3D(1.25, 0.0, 0.0), Point3D(0.75, 0.0, 0.0), false, 1.0));
	expect_boundary_results_equal(
		run_boundary_check(IonpCheckMode::Spatial, ionps, Point3D(2.0, 0.0, 0.0), Point3D(0.5, 0.0, 0.0), true, 0.0),
		run_boundary_check(IonpCheckMode::All, ionps, Point3D(2.0, 0.0, 0.0), Point3D(0.5, 0.0, 0.0), true, 0.0));
}

void test_coating_permeability_boundary_helper() {
	const std::vector<Ionp> ionps = {Ionp(0.0, 0.0, 0.0, 0.0, 1.0, 0.5)};
	const IonpSpatialIndex ionp_index(ionps);
	RandomList random_unit(0.0, 1.0, false);

	{
		Point3D last_position(2.0, 0.0, 0.0);
		Proton proton(0.0, 0.0);
		proton.setXt(Point3D(1.25, 0.0, 0.0));
		int coating_ionp_index = -1;
		bool coating_state_known = false;
		const bool repeat = enforce_ionp_boundaries(proton, last_position, ionps, ionp_index, IonpCheckMode::Spatial, false, coating_ionp_index, coating_state_known, 0.0, false, 1.0, random_unit, nullptr);
		expect(!repeat, "Fully permeable coating should not request a repeated step");
		expect(approx_equal(proton.xt().x(), 1.25), "Fully permeable coating should allow outside-to-coating crossing");
		expect(coating_state_known && coating_ionp_index == 0, "Accepted outside-to-coating crossing did not cache coating membership");
	}

	{
		Point3D last_position(2.0, 0.0, 0.0);
		Proton proton(0.0, 0.0);
		proton.setXt(Point3D(1.25, 0.0, 0.0));
		int coating_ionp_index = -1;
		bool coating_state_known = false;
		const bool repeat = enforce_ionp_boundaries(proton, last_position, ionps, ionp_index, IonpCheckMode::Spatial, false, coating_ionp_index, coating_state_known, 0.0, true, 1.0, random_unit, nullptr);
		expect(!repeat, "Permeability 1.0 coating should not request a repeated step");
		expect(approx_equal(proton.xt().x(), 1.25), "Permeability 1.0 coating should allow outside-to-coating crossing");
		expect(coating_state_known && coating_ionp_index == 0, "Permeability 1.0 crossing did not cache coating membership");
	}

	{
		Point3D last_position(2.0, 0.0, 0.0);
		Proton proton(0.0, 0.0);
		proton.setXt(Point3D(1.25, 0.0, 0.0));
		int coating_ionp_index = -1;
		bool coating_state_known = false;
		const bool repeat = enforce_ionp_boundaries(proton, last_position, ionps, ionp_index, IonpCheckMode::Spatial, false, coating_ionp_index, coating_state_known, 0.0, true, 0.0, random_unit, nullptr);
		expect(!repeat, "Denied outside-to-coating crossing should reflect, not repeat");
		expect(proton.xt().x() > 1.5, "Denied outside-to-coating crossing did not reflect at the coating boundary");
		expect(coating_state_known && coating_ionp_index == -1, "Denied outside-to-coating crossing should cache outside membership");
	}

	{
		Point3D last_position(1.25, 0.0, 0.0);
		Proton proton(0.0, 0.0);
		proton.setXt(Point3D(2.0, 0.0, 0.0));
		int coating_ionp_index = -1;
		bool coating_state_known = false;
		const bool repeat = enforce_ionp_boundaries(proton, last_position, ionps, ionp_index, IonpCheckMode::Spatial, false, coating_ionp_index, coating_state_known, 0.0, true, 0.0, random_unit, nullptr);
		expect(!repeat, "Denied coating-to-outside crossing should reflect, not repeat");
		expect(proton.xt().x() < 1.5 && proton.xt().x() >= 1.0 - 1e-9, "Denied coating-to-outside crossing did not reflect back into the coating");
		expect(coating_state_known && coating_ionp_index == 0, "Denied coating-to-outside crossing should cache coating membership");
	}

	{
		Point3D last_position(1.25, 0.0, 0.0);
		Proton proton(0.0, 0.0);
		proton.setXt(Point3D(0.75, 0.0, 0.0));
		int coating_ionp_index = -1;
		bool coating_state_known = false;
		const bool repeat = enforce_ionp_boundaries(proton, last_position, ionps, ionp_index, IonpCheckMode::Spatial, false, coating_ionp_index, coating_state_known, 0.0, false, 1.0, random_unit, nullptr);
		expect(!repeat, "Core crossing should reflect, not repeat");
		expect(proton.xt().x() > 1.0, "Core crossing did not reflect at the impermeable IONP core");
		expect(coating_state_known && coating_ionp_index == 0, "Core reflection should cache final coating membership");
	}

	{
		Point3D last_position(2.0, 0.0, 0.0);
		Proton proton(0.0, 0.0);
		proton.setXt(Point3D(0.5, 0.0, 0.0));
		int coating_ionp_index = -1;
		bool coating_state_known = false;
		const bool repeat = enforce_ionp_boundaries(proton, last_position, ionps, ionp_index, IonpCheckMode::Spatial, false, coating_ionp_index, coating_state_known, 0.0, true, 0.0, random_unit, nullptr);
		expect(!repeat, "Denied outside-to-core path should reflect at the coating boundary, not repeat");
		expect(proton.xt().x() > 1.5, "Outside-to-core path with zero coating permeability did not stop at the coating boundary first");
		expect(coating_state_known && coating_ionp_index == -1, "Denied outside-to-core path should cache outside membership");
	}
}

void test_ionp_path_ordering_and_tunneling() {
	{
		const std::vector<Ionp> ionps = {Ionp(0.0, 0.0, 0.0, 0.0, 1.0, 0.0)};
		const BoundaryCheckResult result =
			run_boundary_check(IonpCheckMode::Spatial,
							   ionps,
							   Point3D(-2.0, 0.0, 0.0),
							   Point3D(2.0, 0.0, 0.0),
							   false,
							   1.0);
		expect(!result.repeat, "Outside-to-outside core passage should reflect, not retry");
		expect(approx_equal(result.proton_position.x(), -4.0),
			   "Outside-to-outside core passage tunneled through the IONP");
	}

	{
		// Deliberately store the farther IONP first. Collision handling must follow
		// traveled-path order, not the order of IONPs in the input vector.
		const std::vector<Ionp> ionps = {
			Ionp(2.0, 0.0, 0.0, 0.0, 0.5, 0.0),
			Ionp(-2.0, 0.0, 0.0, 0.0, 0.5, 0.0)
		};
		const BoundaryCheckResult result =
			run_boundary_check(IonpCheckMode::Spatial,
							   ionps,
							   Point3D(-5.0, 0.0, 0.0),
							   Point3D(5.0, 0.0, 0.0),
							   false,
							   1.0);
		expect(!result.repeat, "Nearest-core collision should reflect, not retry");
		expect(approx_equal(result.proton_position.x(), -10.0),
			   "A farther IONP was processed before the nearest collision");
	}

	{
		// This chord intersects only the coating, with both endpoints outside.
		const std::vector<Ionp> ionps = {Ionp(0.0, 0.0, 0.0, 0.0, 0.25, 0.75)};
		const BoundaryCheckResult result =
			run_boundary_check(IonpCheckMode::Spatial,
							   ionps,
							   Point3D(-2.0, 0.5, 0.0),
							   Point3D(2.0, 0.5, 0.0),
							   true,
							   0.0);
		expect(!result.repeat, "Denied outside-to-outside coating passage should reflect, not retry");
		expect(result.proton_position.x() < -2.2 && result.proton_position.y() > 2.9,
			   "Outside-to-outside coating passage skipped the first membrane crossing");
	}
}

void test_cell_permeability_boundary_helper() {
	const Point3D fallback_origin(0.0, 0.0, 0.0);
	const std::vector<Point3D> fallback_cell;

	{
		const CellBoundaryCheckResult result = run_cell_boundary_check(
			Point3D(-2.0, 0.0, 0.0), Point3D(0.0, 0.0, 0.0), fallback_cell, fallback_origin, 1.0, 1.0, 0.75);
		expect(!result.reflected, "Permeability 1.0 cell entry should not reflect");
		expect(approx_equal(result.proton_position.x(), 0.0), "Accepted cell entry changed the endpoint");
		expect(result.in_cell && result.cell_number == 1, "Accepted fallback-cell entry did not update membership");
		expect(result.events.empty(), "Accepted cell entry should not record a reflection event");
	}

	{
		const CellBoundaryCheckResult result = run_cell_boundary_check(
			Point3D(-2.0, 0.0, 0.0), Point3D(2.0, 0.0, 0.0), fallback_cell, fallback_origin, 1.0, 0.0, 0.25);
		expect(result.reflected, "Permeability 0.0 cell entry should reflect");
		expect(approx_equal(result.proton_position.x(), -4.0), "Cell reflection did not preserve unused path length");
		expect(!result.in_cell && result.cell_number == 0, "Denied cell entry changed membership");
		expect(result.events.size() == 1, "Denied cell entry should record one reflection event");
		expect(approx_equal(result.events.front().position.x(), -1.0, 2e-6), "Cell entry reflection recorded the wrong hit point");
		expect(result.events.front().kind == 1 && approx_equal(result.events.front().time, 1.25), "Cell reflection event metadata is incorrect");
	}

	{
		const CellBoundaryCheckResult result = run_cell_boundary_check(
			Point3D(0.0, 0.0, 0.0), Point3D(2.0, 0.0, 0.0), fallback_cell, fallback_origin, 1.0, 1.0, 0.75, true, 1);
		expect(!result.reflected, "Permeability 1.0 cell exit should not reflect");
		expect(approx_equal(result.proton_position.x(), 2.0), "Accepted cell exit changed the endpoint");
		expect(!result.in_cell && result.cell_number == 0, "Accepted cell exit did not clear membership");
	}

	{
		const CellBoundaryCheckResult result = run_cell_boundary_check(
			Point3D(0.0, 0.0, 0.0), Point3D(2.0, 0.0, 0.0), fallback_cell, fallback_origin, 1.0, 0.0, 0.25, true, 1);
		expect(result.reflected, "Permeability 0.0 cell exit should reflect");
		expect(approx_equal(result.proton_position.x(), 0.0), "Denied cell exit did not continue its remaining length inside");
		expect(result.in_cell && result.cell_number == 1, "Denied cell exit did not retain membership");
	}

	{
		const CellBoundaryCheckResult accepted = run_cell_boundary_check(
			Point3D(-2.0, 0.0, 0.0), Point3D(0.0, 0.0, 0.0), fallback_cell, fallback_origin, 1.0, 0.5, 0.25);
		const CellBoundaryCheckResult denied = run_cell_boundary_check(
			Point3D(-2.0, 0.0, 0.0), Point3D(0.0, 0.0, 0.0), fallback_cell, fallback_origin, 1.0, 0.5, 0.5);
		expect(!accepted.reflected && accepted.in_cell, "Fractional permeability should accept when draw < probability");
		expect(denied.reflected && !denied.in_cell, "Fractional permeability should deny when draw >= probability");
	}

	{
		const CellBoundaryCheckResult result = run_cell_boundary_check(
			Point3D(-2.0, 0.0, 0.0), Point3D(2.0, 0.0, 0.0), fallback_cell, fallback_origin, 1.0, 1.0, 0.75);
		expect(!result.reflected, "Fully permeable outside-to-outside passage should not reflect");
		expect(approx_equal(result.proton_position.x(), 2.0), "Full cell passage changed the endpoint");
		expect(!result.in_cell && result.cell_number == 0, "Full cell passage left stale membership");
	}

	{
		const std::vector<Point3D> cells = {Point3D(-2.0, 0.0, 0.0), Point3D(2.0, 0.0, 0.0)};
		const CellBoundaryCheckResult result = run_cell_boundary_check(
			Point3D(-4.0, 0.0, 0.0), Point3D(2.0, 0.0, 0.0), cells, fallback_origin, 0.5, 1.0, 0.75);
		expect(!result.reflected, "Fully permeable multicell passage should not reflect");
		expect(result.in_cell && result.cell_number == 2, "Multicell passage did not retain the one-based destination Cellnum");
	}

	{
		const std::vector<Point3D> touching_cells = {Point3D(-1.0, 0.0, 0.0), Point3D(1.0, 0.0, 0.0)};
		const CellBoundaryCheckResult result = run_cell_boundary_check(
			Point3D(-1.0, 0.0, 0.0), Point3D(1.0, 0.0, 0.0), touching_cells, fallback_origin, 1.0, 1.0, 0.75, true, 1);
		expect(!result.reflected, "Fully permeable touching-cell transfer should not reflect");
		expect(result.in_cell && result.cell_number == 2, "Touching-cell transfer skipped the destination membrane");
	}

	{
		const std::vector<Point3D> touching_cells = {Point3D(-1.0, 0.0, 0.0), Point3D(1.0, 0.0, 0.0)};
		rng_config::set_seed(24680);
		RandomList reference_random(0.0, 1.0, false);
		const double first_draw = reference_random.getNumber();
		const double second_draw = reference_random.getNumber();
		const double expected_next_draw = reference_random.getNumber();
		const double permeability = (std::max(first_draw, second_draw) + 1.0) / 2.0;

		rng_config::set_seed(24680);
		RandomList crossing_random(0.0, 1.0, false);
		Point3D last_position(1.0, 0.0, 0.0);
		Proton proton(0.0, 0.0);
		proton.setXt(Point3D(-1.0, 0.0, 0.0));
		proton.inCell(true);
		proton.Cellnum(2);
		enforce_cell_boundaries(proton,
			last_position,
			touching_cells,
			fallback_origin,
			1.0,
			permeability,
			crossing_random,
			0.0,
			nullptr);
		expect(proton.inCell() && proton.Cellnum() == 1, "Reverse touching-cell transfer assigned the wrong destination");
		expect(approx_equal(crossing_random.getNumber(), expected_next_draw),
			"Touching-cell transfer did not consume independent exit and entry decisions");
	}

	{
		const std::vector<Point3D> cells = {Point3D(-2.0, 0.0, 0.0), Point3D(2.0, 0.0, 0.0)};
		const CellBoundaryCheckResult result = run_cell_boundary_check(
			Point3D(0.0, 0.0, 0.0), Point3D(4.0, 0.0, 0.0), cells, fallback_origin, 1.0, 0.0, 0.25);
		expect(result.reflected, "Impermeable cells should report reflected crossings");
		expect(approx_equal(result.proton_position.x(), 0.0), "Two cell reflections did not preserve total path length");
		expect(!result.in_cell && result.cell_number == 0, "Two cell reflections left stale membership");
		expect(result.events.size() == 2, "Path walker did not process both cell reflections in path order");
		expect(approx_equal(result.events[0].position.x(), 1.0, 2e-6) &&
			   approx_equal(result.events[1].position.x(), -1.0, 2e-6),
			   "Multiple cell reflection events were not ordered along the traveled path");
	}

	{
		const std::vector<Point3D> loaded_cell = {Point3D(3.0, 0.0, 0.0)};
		const CellBoundaryCheckResult result = run_cell_boundary_check(
			Point3D(1.5, 0.0, 0.0), Point3D(3.0, 0.0, 0.0), loaded_cell, fallback_origin, 1.0, 1.0, 0.75);
		expect(result.in_cell && result.cell_number == 1, "A single loaded cell center was not authoritative over the fallback origin");
	}

	{
		const CellBoundaryCheckResult result = run_cell_boundary_check(
			Point3D(-2.0, 1.0, 0.0), Point3D(2.0, 1.0, 0.0), fallback_cell, fallback_origin, 1.0, 0.0, 0.25);
		expect(!result.reflected, "Tangential cell contact should not count as a crossing");
		expect(point_distance(result.proton_position, Point3D(2.0, 1.0, 0.0)) <= 1e-12, "Tangential cell contact changed the endpoint");
		expect(!result.in_cell && result.cell_number == 0 && result.events.empty(), "Tangential contact changed state or logged a reflection");
	}

	{
		const std::vector<Point3D> cells = {fallback_origin};
		const CellBoundaryCheckResult result = run_cell_boundary_check(
			Point3D(-2.0, 0.0, 0.0), Point3D(2.0, 0.0, 0.0), cells, fallback_origin, 0.0, 0.0, 0.25);
		expect(!result.reflected, "Zero-radius cell geometry should not reflect");
		expect(approx_equal(result.proton_position.x(), 2.0), "Zero-radius cell geometry changed the endpoint");
		expect(!result.in_cell && result.cell_number == 0, "Zero-radius cell geometry changed membership");
	}

	// Exact endpoint probabilities are deterministic and must not advance the
	// diffusion RNG stream.
	auto next_random_after_crossing = [&](double permeability) {
		rng_config::set_seed(8675309);
		RandomList random_unit(0.0, 1.0, false);
		Point3D last_position(-2.0, 0.0, 0.0);
		Proton proton(0.0, 0.0);
		proton.setXt(Point3D(0.0, 0.0, 0.0));
		enforce_cell_boundaries(proton, last_position, fallback_cell, fallback_origin, 1.0, permeability, random_unit, 0.0, nullptr);
		return random_unit.getNumber();
	};
	auto first_random = []() {
		rng_config::set_seed(8675309);
		RandomList random_unit(0.0, 1.0, false);
		return random_unit.getNumber();
	};
	expect(approx_equal(next_random_after_crossing(0.0), first_random()), "Permeability 0.0 consumed a random draw");
	expect(approx_equal(next_random_after_crossing(1.0), first_random()), "Permeability 1.0 consumed a random draw");
}

void test_periodic_cell_boundary_helper() {
	const Point3D fallback_origin(0.0, 0.0, 0.0);

	{
		const std::vector<Point3D> central_cell = {fallback_origin};
		const CellBoundaryCheckResult result = run_periodic_cell_boundary_check(
			Point3D(4.0, 0.0, 0.0), Point3D(6.0, 0.0, 0.0), central_cell, fallback_origin, 1.0, 0.0, 0.25);
		expect(!result.reflected, "Periodic motion was checked as a false chord through the voxel");
		expect(approx_equal(result.proton_position.x(), -4.0), "Periodic path did not preserve and wrap its full length");
		expect(!result.in_cell && result.cell_number == 0, "Periodic false-chord test changed cell membership");
	}

	{
		const std::vector<Point3D> wrapped_cell = {Point3D(-4.0, 0.0, 0.0)};
		const CellBoundaryCheckResult accepted = run_periodic_cell_boundary_check(
			Point3D(4.0, 0.0, 0.0), Point3D(6.0, 0.0, 0.0), wrapped_cell, fallback_origin, 0.5, 1.0, 0.75);
		expect(!accepted.reflected, "Accepted post-wrap cell entry should not reflect");
		expect(approx_equal(accepted.proton_position.x(), -4.0), "Accepted post-wrap crossing changed the endpoint");
		expect(accepted.in_cell && accepted.cell_number == 1, "Post-wrap cell entry did not update membership");

		const CellBoundaryCheckResult denied = run_periodic_cell_boundary_check(
			Point3D(4.0, 0.0, 0.0), Point3D(6.0, 0.0, 0.0), wrapped_cell, fallback_origin, 0.5, 0.0, 0.25);
		expect(denied.reflected, "Denied post-wrap cell entry should reflect");
		expect(approx_equal(std::abs(denied.proton_position.x()), 5.0, 2e-6), "Post-wrap reflection did not preserve unused path length");
		expect(!denied.in_cell && denied.cell_number == 0, "Denied post-wrap cell entry changed membership");
		expect(denied.events.size() == 1 && approx_equal(denied.events.front().position.x(), -4.5, 2e-6),
			   "Post-wrap reflection recorded the wrong membrane hit");
	}

	{
		const std::vector<Point3D> face_touching_cell = {Point3D(-4.0, 0.0, 0.0)};
		const CellBoundaryCheckResult accepted = run_periodic_cell_boundary_check(
			Point3D(4.0, 0.0, 0.0), Point3D(6.0, 0.0, 0.0), face_touching_cell, fallback_origin, 1.0, 1.0, 0.75);
		expect(!accepted.reflected, "Accepted membrane crossing at a periodic face should not reflect");
		expect(approx_equal(accepted.proton_position.x(), -4.0), "Face-coincident crossing did not wrap its remaining path");
		expect(accepted.in_cell && accepted.cell_number == 1, "Face-coincident entry did not update membership");

		const CellBoundaryCheckResult denied = run_periodic_cell_boundary_check(
			Point3D(4.0, 0.0, 0.0), Point3D(6.0, 0.0, 0.0), face_touching_cell, fallback_origin, 1.0, 0.0, 0.25);
		expect(denied.reflected, "Denied membrane crossing at a periodic face should reflect");
		expect(approx_equal(denied.proton_position.x(), 4.0), "Face-coincident reflection did not preserve unused path length");
		expect(!denied.in_cell && denied.cell_number == 0, "Denied face-coincident entry changed membership");
	}

	{
		const std::vector<Point3D> seam_cell = {Point3D(-4.8, 0.0, 0.0)};
		const CellBoundaryCheckResult result = run_periodic_cell_boundary_check(
			Point3D(4.5, 0.0, 0.0), Point3D(4.7, 0.0, 0.0), seam_cell, fallback_origin, 1.0, 0.0, 0.25);
		expect(!result.reflected, "Motion inside a periodic cell image should not create a crossing");
		expect(approx_equal(result.proton_position.x(), 4.7), "Motion inside a periodic cell image changed the endpoint");
		expect(result.in_cell && result.cell_number == 1, "Periodic-image membership was not synchronized");
	}
}

void test_region_initialization() {
	ImportanceSamplingConfig config;
	config.enabled = true;
	config.cell_near_cutoff = 3.0;
	config.protons_inside = 3;
	config.protons_near = 4;
	config.protons_far = 5;
	config.weight_inside = 1.0;
	config.weight_near = 1.0;
	config.weight_far = 1.0;

	SimMatrix sim(0.0, 0.0, 0.0, 20.0, 20.0, 20.0, 0, 0, 0);
	const Point3D cell_origin(0.0, 0.0, 0.0);
	MRsequence sequence = make_sequence(sim, config.total_free_protons(), 0, 0.0, 0.01, 0.0, 0.0, cell_origin, 1.0, config);

	std::vector<Proton> protons;
	std::vector<Ionp> ionps;
	std::vector<std::complex<double>> signal;
	std::vector<double> timeV;
	std::vector<std::complex<double>> signal_mse;
	std::vector<double> timeV_mse;
	std::vector<Point3D> cell_origins;
	ImportanceRegionSignals region_signals;

	sequence.Sequence(
		protons,
		ionps,
		signal,
		timeV,
		signal_mse,
		timeV_mse,
		cell_origins,
		&region_signals);

	expect(protons.size() == static_cast<std::size_t>(config.total_free_protons()), "Unexpected proton count after initialization");
	for (int proton_index = 0; proton_index < config.protons_inside; ++proton_index) {
		expect(distance(protons[static_cast<std::size_t>(proton_index)].xt(), cell_origin) <= 1.0 + 1e-12, "Inside proton was initialized outside the cell");
	}
	for (int proton_index = config.protons_inside; proton_index < config.protons_inside + config.protons_near; ++proton_index) {
		const double proton_distance = distance(protons[static_cast<std::size_t>(proton_index)].xt(), cell_origin);
		expect(proton_distance > 1.0 - 1e-12, "Near proton was initialized inside the cell");
		expect(proton_distance <= config.cell_near_cutoff + 1e-12, "Near proton was initialized outside the near region");
	}
	for (int proton_index = config.protons_inside + config.protons_near; proton_index < config.total_free_protons(); ++proton_index) {
		const double proton_distance = distance(protons[static_cast<std::size_t>(proton_index)].xt(), cell_origin);
		expect(proton_distance > config.cell_near_cutoff - 1e-12, "Far proton was initialized inside the near region");
	}
	expect(signal.size() == 1, "Unexpected zero-time signal length");
	expect(approx_equal(signal.front().real(), 1.0), "Zero-field signal magnitude should stay at one");
	expect(approx_equal(signal.front().imag(), 0.0), "Zero-field signal should be purely real");
	const std::array<ProtonSampleRegion, 3> legacy_regions = {
		ProtonSampleRegion::Inside,
		ProtonSampleRegion::Near,
		ProtonSampleRegion::Far
	};
	for (ProtonSampleRegion region : legacy_regions) {
		const auto& regional_signal =
			region_signals[sample_region_index(region)];
		expect(regional_signal.size() == 1 &&
			   approx_equal(regional_signal.front().real(), 1.0),
			   "Legacy cell importance-sampling region signal was not returned");
	}
}

void test_loaded_cell_region_initialization() {
	rng_config::set_seed(4123);

	ImportanceSamplingConfig config;
	config.enabled = true;
	config.protons_inside = 6;
	config.weight_inside = 1.0;

	SimMatrix sim(0.0, 0.0, 0.0, 20.0, 20.0, 20.0, 0, 0, 0);
	const Point3D fallback_origin(-3.0, 0.0, 0.0);
	const Point3D loaded_origin(3.0, 0.0, 0.0);
	MRsequence sequence = make_sequence(sim, config.total_free_protons(), 0, 0.0, 0.01, 0.0, 0.0, fallback_origin, 1.0, config);

	std::vector<Proton> protons;
	std::vector<Ionp> ionps;
	std::vector<std::complex<double>> signal;
	std::vector<double> timeV;
	std::vector<std::complex<double>> signal_mse;
	std::vector<double> timeV_mse;
	std::vector<Point3D> cell_origins = {loaded_origin};

	sequence.Sequence(protons, ionps, signal, timeV, signal_mse, timeV_mse, cell_origins);

	expect(protons.size() == static_cast<std::size_t>(config.protons_inside), "Unexpected loaded-cell proton count");
	for (const Proton& proton : protons) {
		expect(distance(proton.xt(), loaded_origin) <= 1.0 + 1e-12, "Inside proton ignored the authoritative loaded cell center");
		expect(proton.inCell() && proton.Cellnum() == 1, "Loaded-cell initialization assigned incorrect live membership");
	}
}

void test_auto_volume_weights() {
	rng_config::set_seed(1234);

	ImportanceSamplingConfig config;
	config.enabled = true;
	config.cell_near_cutoff = 2.0;
	config.protons_inside = 1;
	config.protons_near = 1;
	config.protons_far = 1;

	const std::array<double, 3> auto_weights =
		estimate_active_free_region_volume_weights(config, {}, Point3D(0.0, 0.0, 0.0), 1.0, 10.0, 10.0, 10.0, 300000);

	const double inside_exact = (4.0 / 3.0) * kPi * std::pow(1.0, 3) / 1000.0;
	const double near_exact = (4.0 / 3.0) * kPi * (std::pow(2.0, 3) - std::pow(1.0, 3)) / 1000.0;
	const double far_exact = 1.0 - inside_exact - near_exact;

	expect(std::abs(auto_weights[0] - inside_exact) < 0.0015, "Auto inside weight is not close to the expected single-cell volume fraction");
	expect(std::abs(auto_weights[1] - near_exact) < 0.0025, "Auto near weight is not close to the expected single-cell shell fraction");
	expect(std::abs(auto_weights[2] - far_exact) < 0.003, "Auto far weight is not close to the expected remaining fraction");
}

void test_bound_weight_at_initial_time() {
	ImportanceSamplingConfig config;
	config.enabled = true;
	config.cell_near_cutoff = 3.0;
	config.protons_inside = 2;
	config.protons_near = 2;
	config.protons_far = 2;
	config.weight_inside = 2.0;
	config.weight_near = 3.0;
	config.weight_far = 4.0;
	config.weight_bound = 5.0;

	SimMatrix sim(0.0, 0.0, 0.0, 20.0, 20.0, 20.0, 0, 0, 0);
	const Point3D cell_origin(0.0, 0.0, 0.0);
	MRsequence sequence = make_sequence(sim, config.total_free_protons(), 2, 0.0, 0.01, 0.5, 0.1, cell_origin, 1.0, config);

	std::vector<Proton> protons;
	std::vector<Ionp> ionps = {
		Ionp(4.0, 0.0, 0.0, 0.0, 0.5, 0.25)
	};
	std::vector<std::complex<double>> signal;
	std::vector<double> timeV;
	std::vector<std::complex<double>> signal_mse;
	std::vector<double> timeV_mse;
	std::vector<Point3D> cell_origins;

	sequence.Sequence(protons, ionps, signal, timeV, signal_mse, timeV_mse, cell_origins);

	expect(signal.size() == 1, "Unexpected initial-time signal length with bound protons");
	expect(approx_equal(signal.front().real(), 1.0), "Weighted initial-time signal magnitude should stay at one");
	expect(approx_equal(signal.front().imag(), 0.0), "Weighted initial-time signal should be purely real");
}

void test_surrounding_voxels_keep_bound_jobs_physical_only() {
	rng_config::set_seed(7321);

	ImportanceSamplingConfig config;
	const SimMatrix sim(
		0.0, 0.0, 0.0, 20.0, 22.0, 24.0, 0, 0, 0);
	const Point3D cell_origin;
	MRsequence sequence = make_sequence(
		sim,
		0,
		2,
		0.0,
		0.01,
		0.5,
		0.25,
		cell_origin,
		0.0,
		config);
	sequence.set_surrounding_voxels(true);
	sequence.set_surrounding_voxel_superparticles(1);

	std::vector<Proton> protons;
	std::vector<Ionp> ionps = {
		Ionp(-3.0, 0.0, 0.0, 0.0, 0.5, 0.25),
		Ionp(3.0, 0.0, 0.0, 0.0, 0.7, 0.4)};
	const std::vector<Ionp> original_ionps = ionps;
	std::vector<std::complex<double>> signal;
	std::vector<double> timeV;
	std::vector<std::complex<double>> signal_mse;
	std::vector<double> timeV_mse;
	std::vector<Point3D> cell_origins;

	sequence.Sequence(
		protons,
		ionps,
		signal,
		timeV,
		signal_mse,
		timeV_mse,
		cell_origins);

	expect(ionps.size() == original_ionps.size(),
		   "Surrounding magnetic replicas changed the physical IONP count");
	for (std::size_t index = 0; index < ionps.size(); ++index) {
		expect(point_distance(ionps[index].x0(), original_ionps[index].x0()) == 0.0 &&
			   point_distance(ionps[index].xt(), original_ionps[index].xt()) == 0.0 &&
			   ionps[index].radius() == original_ionps[index].radius() &&
			   ionps[index].coating() == original_ionps[index].coating(),
			   "Surrounding magnetic replicas mutated a physical IONP");
	}
	expect(protons.size() == 4,
		   "Surrounding magnetic replicas incorrectly created bound-proton jobs");
	expect(signal.size() == 1 && approx_equal(signal.front().real(), 1.0),
		   "Physical-only bound-proton run produced an invalid initial signal");
}

void test_overlapping_bound_weight_at_initial_time() {
	ImportanceSamplingConfig config;
	config.enabled = true;
	config.mode = ImportanceSamplingMode::OverlappingProposal;
	config.proposal_protons[free_sampling_proposal_index(
		FreeSamplingProposal::Uniform)] = 4;
	config.weight_bound = 0.25;

	SimMatrix sim(
		0.0, 0.0, 0.0, 20.0, 20.0, 20.0, 0, 0, 0);
	const Point3D cell_origin;
	MRsequence sequence = make_sequence(
		sim,
		config.total_free_protons(),
		2,
		0.0,
		0.01,
		0.5,
		0.25,
		cell_origin,
		0.0,
		config);

	std::vector<Proton> protons;
	std::vector<Ionp> ionps = {
		Ionp(0.0, 0.0, 0.0, 0.0, 0.5, 0.25)
	};
	std::vector<std::complex<double>> signal;
	std::vector<double> timeV;
	std::vector<std::complex<double>> signal_mse;
	std::vector<double> timeV_mse;
	std::vector<Point3D> cells;
	sequence.Sequence(
		protons,
		ionps,
		signal,
		timeV,
		signal_mse,
		timeV_mse,
		cells);

	expect(protons.size() == 6,
		   "Overlapping free/bound run produced the wrong sample count");
	expect(signal.size() == 1 &&
		   approx_equal(signal.front().real(), 1.0) &&
		   approx_equal(signal.front().imag(), 0.0),
		   "Overlapping free/bound signal is not normalized at initial time");
}

void test_variable_step_available_to_all_free_samples() {
	ImportanceSamplingConfig inside_config;
	inside_config.enabled = true;
	inside_config.cell_near_cutoff = 3.0;
	inside_config.protons_inside = 1;
	inside_config.weight_inside = 1.0;

	const Point3D inside_fixed_step = simulate_single_free_proton_final_position(inside_config, false, 2026);
	const Point3D inside_variable_step = simulate_single_free_proton_final_position(inside_config, true, 2026);
	expect(point_distance(inside_fixed_step, inside_variable_step) > 1e-6, "IS inside samples should still use the subdivided variable-step path");

	ImportanceSamplingConfig near_config;
	near_config.enabled = true;
	near_config.cell_near_cutoff = 3.0;
	near_config.protons_near = 1;
	near_config.weight_near = 1.0;

	const Point3D near_fixed_step = simulate_single_free_proton_final_position(near_config, false, 2026);
	const Point3D near_variable_step = simulate_single_free_proton_final_position(near_config, true, 2026);
	expect(point_distance(near_fixed_step, near_variable_step) > 1e-6, "IS near samples should use dynamic variable steps when enabled");
}

void test_zero_time_step_fails_fast() {
	ImportanceSamplingConfig config;
	SimMatrix sim(0.0, 0.0, 0.0, 20.0, 20.0, 20.0, 0, 0, 0);
	const Point3D cell_origin(0.0, 0.0, 0.0);
	MRsequence sequence = make_sequence(sim,
										 1,
										 0,
										 0.01,
										 0.0,
										 0.5,
										 0.0,
										 cell_origin,
										 0.0,
										 config,
										 true);

	std::vector<Proton> protons;
	std::vector<Ionp> ionps = {Ionp(0.0, 0.0, 0.0, 0.0, 0.5, 0.0)};
	std::vector<std::complex<double>> signal;
	std::vector<double> timeV;
	std::vector<std::complex<double>> signal_mse;
	std::vector<double> timeV_mse;
	std::vector<Point3D> cell_origins;

	bool threw = false;
	try {
		sequence.Sequence(protons, ionps, signal, timeV, signal_mse, timeV_mse, cell_origins);
	} catch (const std::runtime_error&) {
		threw = true;
	}

	expect(threw, "Zero time_step should fail fast instead of entering the propagation loop");
}

void test_replace_proton_records_current_time() {
	ImportanceSamplingConfig config;
	SimMatrix sim(0.0, 0.0, 0.0, 20.0, 20.0, 20.0, 0, 0, 0);
	const Point3D cell_origin(0.0, 0.0, 0.0);
	MRsequence sequence = make_sequence(sim,
										1,
										0,
										1.0,
										0.25,
										0.5,
										0.0,
										cell_origin,
										0.0,
										config);

	std::vector<Proton> protons = {Proton(0.0, 0.0)};
	const std::vector<Ionp> ionps;
	const IonpSpatialIndex ionp_index(ionps);
	const std::vector<Point3D> cell_origins;
	const CellSpatialIndex cell_index(cell_origins,
									 cell_origin,
									 0.0,
									 sim.xLength(),
									 sim.yLength(),
									 sim.zLength(),
									 true);
	int coating_ionp_index = -1;
	bool coating_state_known = false;
	RandomList random_x(-1.0, 1.0, false);
	RandomList random_y(-1.0, 1.0, false);
	RandomList random_z(-1.0, 1.0, false);
	RandomList random_diff(0.25, 0.25, false);

	sequence.replaceProton(protons,
						   ionps,
						   ionp_index,
						   cell_index,
						   coating_ionp_index,
						   coating_state_known,
						   cell_origins,
						   0.25,
						   0,
						   2,
						   0.75,
						   ProtonSampleRegion::Far,
						   random_x,
						   random_y,
						   random_z,
						   random_diff,
						   nullptr);
	expect(approx_equal(protons.front().diffTime(), 0.75),
		   "replaceProton did not retain the physical current time");
}

void test_combined_cell_ionp_path_ordering() {
	ImportanceSamplingConfig config;
	SimMatrix sim(0.0, 0.0, 0.0, 20.0, 20.0, 20.0, 0, 0, 0);
	const Point3D cell_origin(-1.0, 0.0, 0.0);
	const double diffusion_time = 2e9; // sqrt(6 * D * dt) = 6
	MRsequence sequence = make_sequence(sim,
										1,
										0,
										diffusion_time,
										diffusion_time,
										0.25,
										0.0,
										cell_origin,
										0.5,
										config);

	std::vector<Proton> protons = {
		Proton(Point3D(3.0, 0.0, 0.0), 0.0, 0.0)
	};
	const std::vector<Ionp> ionps = {
		Ionp(1.5, 0.0, 0.0, 0.0, 0.25, 0.0)
	};
	const IonpSpatialIndex ionp_index(ionps);
	const std::vector<Point3D> cell_origins = {cell_origin};
	const CellSpatialIndex cell_index(cell_origins,
									 cell_origin,
									 0.5,
									 sim.xLength(),
									 sim.yLength(),
									 sim.zLength(),
									 true);
	int coating_ionp_index = -1;
	bool coating_state_known = false;
	RandomList random_x(-1.0, 1.0, false);
	RandomList random_y(-1.0, 1.0, false);
	RandomList random_z(-1.0, 1.0, false);
	RandomList random_diff(0.5, 0.5, false); // direction along -x

	sequence.replaceProton(protons,
						   ionps,
						   ionp_index,
						   cell_index,
						   coating_ionp_index,
						   coating_state_known,
						   cell_origins,
						   diffusion_time,
						   0,
						   1,
						   diffusion_time,
						   ProtonSampleRegion::Far,
						   random_x,
						   random_y,
						   random_z,
						   random_diff,
						   nullptr);

	expect(approx_equal(protons.front().xt().x(), 6.5),
		   "Cell and IONP collisions were not processed in traveled-path order");
	expect(!protons.front().inCell(),
		   "An unreachable later cell collision changed proton membership");
}

} // namespace

int main() {
	try {
		test_variable_step_eligibility_helper();
		test_surrounding_sources_drive_variable_step_proximity();
		test_surrounding_sources_move_together_in_direct_field();
		test_superparticle_configuration_invariants();
		test_signal_aggregation_helper();
		test_overlapping_proposal_names_and_counts();
		test_overlapping_effective_count_redistribution();
		test_overlapping_mis_weight_formula();
		test_overlapping_signal_aggregation_helper();
		test_legacy_automatic_weights_survive_rematerialization();
		test_all_joint_and_bound_buckets();
		test_clear_region_names();
		test_realized_outer_radius_semantics();
		test_ionp_population_and_joint_region_classification();
		test_joint_particle_auto_volume_weights();
		test_overlapping_proposal_volume_estimation();
		test_joint_region_initialization();
		test_overlapping_proposal_initialization();
		test_overlapping_zero_support_falls_back_to_uniform();
		test_ionp_spatial_index_helper();
		test_translated_ionp_spatial_index_helper();
		test_explicit_random_list_seed();
		test_coating_diffusion_coefficient_helper();
		test_ionp_check_modes_equivalent();
		test_coating_permeability_boundary_helper();
		test_ionp_path_ordering_and_tunneling();
		test_cell_permeability_boundary_helper();
		test_periodic_cell_boundary_helper();
		test_region_initialization();
		test_loaded_cell_region_initialization();
		test_auto_volume_weights();
		test_bound_weight_at_initial_time();
		test_surrounding_voxels_keep_bound_jobs_physical_only();
		test_overlapping_bound_weight_at_initial_time();
		test_variable_step_available_to_all_free_samples();
		test_zero_time_step_fails_fast();
		test_replace_proton_records_current_time();
		test_combined_cell_ionp_path_ordering();
	} catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		return 1;
	}
	return 0;
}
