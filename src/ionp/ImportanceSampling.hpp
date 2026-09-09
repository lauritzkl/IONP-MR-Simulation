/**
 * @file ImportanceSampling.hpp
 *
 * Sampling labels, configuration, geometric volume estimates and signal
 * aggregation. Exclusive joint regions and overlapping proposals use different
 * weighting rules. Labels record the initial sampling assignment; diffusion
 * across a boundary does not reassign a proton's estimator weight.
 */

#ifndef IMPORTANCE_SAMPLING_HPP_
#define IMPORTANCE_SAMPLING_HPP_

#include <array>
#include <complex>
#include <cstddef>
#include <string>
#include <vector>

#include "Ionp.hpp"
#include "Point3D.hpp"

enum class CellSampleRegion : unsigned char {
    Inside = 0,
    Near = 1,
    Far = 2
};

enum class IonpDistanceRegion : unsigned char {
    Near = 0,
    Intermediate = 1,
    Far = 2
};

enum class IonpPopulation : unsigned char {
    Intracellular = 0,
    Extracellular = 1
};

enum class ImportanceSamplingMode : unsigned char {
    LegacyJoint = 0,
    OverlappingProposal = 1
};

// Overlapping importance-sampling proposals. Unlike ProtonSampleRegion, these
// labels identify the distribution used to draw a proton, not an exclusive
// three-axis spatial stratum.
enum class FreeSamplingProposal : unsigned char {
    Uniform = 0,
    CellInside = 1,
    CellNear = 2,
    CellFar = 3,
    IntracellularIonpNear = 4,
    IntracellularIonpIntermediate = 5,
    IntracellularIonpFar = 6,
    ExtracellularIonpNear = 7,
    ExtracellularIonpIntermediate = 8,
    ExtracellularIonpFar = 9
};

constexpr std::size_t kCellSampleRegionCount = 3;
constexpr std::size_t kIonpDistanceRegionCount = 3;
constexpr std::size_t kFreeSamplingProposalCount = 10;
constexpr std::size_t kFreeSampleRegionCount =
    kCellSampleRegionCount * kIonpDistanceRegionCount * kIonpDistanceRegionCount;
constexpr std::size_t kBoundIntracellularSampleRegionIndex = kFreeSampleRegionCount;
constexpr std::size_t kBoundExtracellularSampleRegionIndex = kFreeSampleRegionCount + 1;
constexpr std::size_t kSampleRegionCount = kFreeSampleRegionCount + 2;

constexpr std::size_t joint_sample_region_index(CellSampleRegion cell_region,
                                             IonpDistanceRegion intracellular_ionp_region,
                                             IonpDistanceRegion extracellular_ionp_region) {
    return static_cast<std::size_t>(cell_region) * 9u +
        static_cast<std::size_t>(intracellular_ionp_region) * 3u +
        static_cast<std::size_t>(extracellular_ionp_region);
}

enum class ProtonSampleRegion : unsigned char {
    BoundIntracellular = static_cast<unsigned char>(kBoundIntracellularSampleRegionIndex),
    BoundExtracellular = static_cast<unsigned char>(kBoundExtracellularSampleRegionIndex),

    // Legacy aliases map cell-only importance sampling onto both IONP-far axes.
    Inside = static_cast<unsigned char>(
        joint_sample_region_index(CellSampleRegion::Inside,
                                  IonpDistanceRegion::Far,
                                  IonpDistanceRegion::Far)),
    Near = static_cast<unsigned char>(
        joint_sample_region_index(CellSampleRegion::Near,
                                  IonpDistanceRegion::Far,
                                  IonpDistanceRegion::Far)),
    Far = static_cast<unsigned char>(
        joint_sample_region_index(CellSampleRegion::Far,
                                  IonpDistanceRegion::Far,
                                  IonpDistanceRegion::Far)),
    Bound = BoundExtracellular
};

struct JointSampleRegion {
    CellSampleRegion cell = CellSampleRegion::Far;
    IonpDistanceRegion intracellular_ionp = IonpDistanceRegion::Far;
    IonpDistanceRegion extracellular_ionp = IonpDistanceRegion::Far;
};

// Classification of physical IONPs against cells; these population labels
// are independent of the evolving cell membership of a diffusing proton.
struct IonpPopulationMap {
    std::vector<IonpPopulation> populations;
    std::vector<std::size_t> intracellular_indices;
    std::vector<std::size_t> extracellular_indices;
    std::size_t coated_intracellular_count = 0;
    std::size_t coated_extracellular_count = 0;
};

struct ImportanceSamplingConfig {
    bool enabled = false;

    // Existing cell-only interface.
    // Cell cutoff is distance from the cell CENTRE (metres), unlike IONP
    // cutoffs below, which measure distance from the outer particle surface.
    double cell_near_cutoff = 0.;
    int protons_inside = 0;
    int protons_near = 0;
    int protons_far = 0;
    double weight_inside = 0.;
    double weight_near = 0.;
    double weight_far = 0.;
    double weight_bound = 0.;

    // Population-specific IONP surface-distance cutoffs.
    double intracellular_ionp_near_cutoff = 0.;
    double intracellular_ionp_far_cutoff = 0.;
    double extracellular_ionp_near_cutoff = 0.;
    double extracellular_ionp_far_cutoff = 0.;

    std::array<int, kFreeSampleRegionCount> joint_protons{};
    std::array<double, kFreeSampleRegionCount> joint_weights{};
    bool joint_regions_configured = false;

    // Overlapping-proposal interface. Proposal counts control how many starts
    // are drawn from each marginal region; signal weights are derived
    // automatically with deterministic-mixture MIS.
    std::array<int, kFreeSamplingProposalCount> proposal_protons{};
    // Accessible proposal volumes relative to accessible free water, used to
    // convert proposal counts into mixture densities for MIS weights.
    std::array<double, kFreeSamplingProposalCount> proposal_volume_fractions{};
    double free_water_volume_fraction = 0.;
    bool proposal_volumes_estimated = false;
    ImportanceSamplingMode mode = ImportanceSamplingMode::LegacyJoint;

    int bound_protons_intracellular = 0;
    int bound_protons_extracellular = 0;
    double weight_bound_intracellular = 0.;
    double weight_bound_extracellular = 0.;
    bool split_bound_configured = false;

    int total_free_protons() const;
    int total_proposal_protons() const;
    bool overlapping_proposals_enabled() const;
    void validate_mode_exclusivity() const;
    bool intracellular_ionp_sampling_enabled() const;
    bool extracellular_ionp_sampling_enabled() const;
    void materialize_legacy_regions();
};

constexpr std::size_t free_sampling_proposal_index(
    FreeSamplingProposal proposal) {
    return static_cast<std::size_t>(proposal);
}

std::string free_sampling_proposal_name(FreeSamplingProposal proposal);
FreeSamplingProposal cell_sampling_proposal(CellSampleRegion region);
FreeSamplingProposal ionp_sampling_proposal(
    IonpPopulation population,
    IonpDistanceRegion region);
std::array<FreeSamplingProposal, 3> marginal_sampling_proposals(
    const JointSampleRegion& membership);

std::size_t sample_region_index(ProtonSampleRegion region);
ProtonSampleRegion make_joint_sample_region(CellSampleRegion cell_region,
                                            IonpDistanceRegion intracellular_ionp_region,
                                            IonpDistanceRegion extracellular_ionp_region);
JointSampleRegion decode_joint_sample_region(ProtonSampleRegion region);
std::string sample_region_name(ProtonSampleRegion region);

bool has_explicit_free_region_weights(const ImportanceSamplingConfig& config);

std::array<std::size_t, kSampleRegionCount> sample_region_counts(
    const ImportanceSamplingConfig& config,
    std::size_t bound_intracellular_sample_count,
    std::size_t bound_extracellular_sample_count);

// Legacy convenience overload: common bound samples are mapped to the extracellular bucket.
std::array<std::size_t, kSampleRegionCount> sample_region_counts(
    const ImportanceSamplingConfig& config,
    std::size_t bound_sample_count);

IonpPopulationMap classify_ionp_populations(const std::vector<Ionp>& ionps,
                                             const std::vector<Point3D>& cell_origins,
                                             const Point3D& cell_origin,
                                             double cell_radius);

ProtonSampleRegion classify_sample_region(const Point3D& point,
                                          const std::vector<Point3D>& cell_origins,
                                          const Point3D& cell_origin,
                                          double cell_r,
                                          const ImportanceSamplingConfig& importance_sampling);

bool classify_joint_sample_region(const Point3D& point,
                                  const std::vector<Point3D>& cell_origins,
                                  const Point3D& cell_origin,
                                  double cell_r,
                                  const std::vector<Ionp>& ionps,
                                  const IonpPopulationMap& populations,
                                  const ImportanceSamplingConfig& importance_sampling,
                                  ProtonSampleRegion& region);

bool classify_joint_sample_region_from_ionp_candidates(
    const Point3D& point,
    const std::vector<Point3D>& cell_origins,
    const Point3D& cell_origin,
    double cell_r,
    const std::vector<Ionp>& ionps,
    const IonpPopulationMap& populations,
    const std::vector<unsigned int>& ionp_candidates,
    const ImportanceSamplingConfig& importance_sampling,
    ProtonSampleRegion& region);

std::array<double, kFreeSampleRegionCount> estimate_active_joint_region_volume_weights(
    const ImportanceSamplingConfig& config,
    const std::vector<Point3D>& cell_origins,
    const Point3D& cell_origin,
    double cell_r,
    const std::vector<Ionp>& ionps,
    const IonpPopulationMap& populations,
    double sim_x_length,
    double sim_y_length,
    double sim_z_length,
    std::size_t sample_count = 200000);

// Existing cell-only estimator retained for compatibility and regression tests.
std::array<double, 3> estimate_active_free_region_volume_weights(
    const ImportanceSamplingConfig& config,
    const std::vector<Point3D>& cell_origins,
    const Point3D& cell_origin,
    double cell_r,
    double sim_x_length,
    double sim_y_length,
    double sim_z_length,
    std::size_t sample_count = 200000);

std::array<double, kSampleRegionCount> normalized_sample_weights(
    const ImportanceSamplingConfig& config,
    std::size_t bound_intracellular_sample_count,
    std::size_t bound_extracellular_sample_count);

// Legacy convenience overload.
std::array<double, kSampleRegionCount> normalized_sample_weights(
    const ImportanceSamplingConfig& config,
    bool include_bound);

using ImportanceRegionSignals =
    std::array<std::vector<std::complex<double>>, kSampleRegionCount>;

struct WeightedSignalAggregation {
    ImportanceRegionSignals region_signals;
    std::vector<std::complex<double>> combined_signal;
};

WeightedSignalAggregation aggregate_weighted_signals(
    const std::vector<std::vector<double>>& phase_vector,
    const std::vector<ProtonSampleRegion>& sample_regions,
    const std::array<std::size_t, kSampleRegionCount>& sample_counts,
    const std::array<double, kSampleRegionCount>& normalized_weights);

struct OverlappingProposalVolumeEstimate {
    double accessible_volume = 0.;
    double box_volume = 0.;
    std::array<double, kFreeSamplingProposalCount> region_volumes{};
    std::array<double, kFreeSamplingProposalCount> region_fractions{};
    std::size_t pilot_samples = 0;
};

// Estimate all marginal proposal volumes in one deterministic pilot pass. The
// returned fractions are relative to accessible free-proton volume, so the
// uniform proposal fraction is one whenever accessible volume is positive.
OverlappingProposalVolumeEstimate estimate_overlapping_proposal_volumes(
    const ImportanceSamplingConfig& config,
    const std::vector<Point3D>& cell_origins,
    const Point3D& cell_origin,
    double cell_r,
    const std::vector<Ionp>& ionps,
    const IonpPopulationMap& populations,
    double sim_x_length,
    double sim_y_length,
    double sim_z_length,
    std::size_t sample_count = 200000);

struct EffectiveProposalCounts {
    std::array<std::size_t, kFreeSamplingProposalCount> counts{};
    std::array<bool, kFreeSamplingProposalCount> redistributed{};
    std::size_t redistributed_to_uniform = 0;
};

// Resolve geometry-dependent proposal counts. A requested proposal with zero
// estimated volume is disabled and its complete budget is moved to the
// mandatory uniform proposal, keeping the total trajectory count unchanged.
EffectiveProposalCounts resolve_effective_proposal_counts(
    const std::array<int, kFreeSamplingProposalCount>& requested_counts,
    const std::array<double, kFreeSamplingProposalCount>& region_fractions);

// Deterministic-mixture balance-heuristic coefficient for one accessible
// starting point. `membership` describes the point, independent of which
// proposal generated it.
double overlapping_mis_raw_weight(
    const JointSampleRegion& membership,
    const std::array<std::size_t, kFreeSamplingProposalCount>& effective_counts,
    const std::array<double, kFreeSamplingProposalCount>& region_fractions);

std::vector<double> overlapping_mis_raw_weights(
    const std::vector<JointSampleRegion>& memberships,
    const std::array<std::size_t, kFreeSamplingProposalCount>& effective_counts,
    const std::array<double, kFreeSamplingProposalCount>& region_fractions);

std::vector<double> self_normalize_importance_weights(
    const std::vector<double>& raw_weights);

std::vector<double> overlapping_mis_self_normalized_weights(
    const std::vector<JointSampleRegion>& memberships,
    const std::array<std::size_t, kFreeSamplingProposalCount>& effective_counts,
    const std::array<double, kFreeSamplingProposalCount>& region_fractions);

using ImportanceProposalSignals =
    std::array<std::vector<std::complex<double>>, kFreeSamplingProposalCount>;

struct OverlappingWeightedSignalAggregation {
    ImportanceProposalSignals proposal_means;
    ImportanceProposalSignals proposal_contributions;
    std::vector<std::complex<double>> combined_signal;
};

// Aggregate the normalized free signal with per-proton MIS weights. Proposal
// signals are unweighted conditional means and are diagnostics only; because
// proposals overlap, they must not be summed as exclusive compartments.
OverlappingWeightedSignalAggregation aggregate_overlapping_weighted_signals(
    const std::vector<std::vector<double>>& phase_vector,
    const std::vector<FreeSamplingProposal>& sample_proposals,
    const std::vector<double>& normalized_sample_weights);

#endif
