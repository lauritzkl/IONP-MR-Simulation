/**
 * @file Data.cpp
 *
 * Configuration parsing proceeds through key normalization, typed conversion,
 * validation and default resolution. Auxiliary geometry/histogram readers and
 * output serialization follow. Parsing deliberately retains the existing
 * permissive numeric-prefix and stream-extraction behavior.
 */

#include <algorithm>
#include <cctype>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "mptmacros.h"
#include "Point3D.hpp"
#include "Ionp.hpp"
#include "Aggregate.hpp"
#include "SimMatrix.hpp"
#include "Data.hpp"
#include "RandNumberBetween.hpp"
#include "OutputUtils.hpp"

namespace {

enum class ConfigKey {
	IonpRadius,
	AggRadius,
	MatrixX0,
	MatrixY0,
	MatrixZ0,
	Origin,
	MatrixXLength,
	MatrixYLength,
	MatrixZLength,
	Length,
	MatrixXNodes,
	MatrixYNodes,
	MatrixZNodes,
	Nodes,
	IonpNumber,
	AggNumber,
	IonpPositionFile,
	AggPositionFile,
	CellPositionFile,
	CellLoadFile,
	ProtonNumber,
	EchoTime,
	TimeStep,
	VariableTime,
	SurroundingVoxels,
	SurroundingVoxelSuperparticles,
	EchoSpacing,
	Threads,
	VelocityIonp,
	Direction,
	Sequence,
	CellHist,
	RadiusHist,
	Coating,
	CoatingDiffusionCoefficient,
	CellDiffusionCoefficient,
	CoatingPermeability,
	CellPermeability,
	BoundProtonNumber,
	BackgroundT2,
	CellT2,
	IonpRadius1,
	IonpPer,
	IonpRadiusStd,
	IonpRadiusStd1,
	AggRadiusStd,
	Concentration,
	ConcentrationExternal,
	MagneticMoment,
	CellOrigin,
	VoxelPosition,
	CubeOrigin,
	CubeSize,
	CellRadius,
	CellConcentration,
	CellConcentrationStd,
	CellConcentrationDist,
	CellNumber,
	BFieldSolver,
	ImportanceSamplingEnabled,
	ImportanceSamplingCellNearCutoff,
	ImportanceSamplingIntracellularIonpNearCutoff,
	ImportanceSamplingIntracellularIonpFarCutoff,
	ImportanceSamplingExtracellularIonpNearCutoff,
	ImportanceSamplingExtracellularIonpFarCutoff,
	ImportanceSamplingProposalProtonsUniform,
	ImportanceSamplingProposalProtonsCellInside,
	ImportanceSamplingProposalProtonsCellNear,
	ImportanceSamplingProposalProtonsCellFar,
	ImportanceSamplingProposalProtonsIntracellularIonpNear,
	ImportanceSamplingProposalProtonsIntracellularIonpIntermediate,
	ImportanceSamplingProposalProtonsIntracellularIonpFar,
	ImportanceSamplingProposalProtonsExtracellularIonpNear,
	ImportanceSamplingProposalProtonsExtracellularIonpIntermediate,
	ImportanceSamplingProposalProtonsExtracellularIonpFar,
	ImportanceSamplingProtonsInside,
	ImportanceSamplingProtonsNear,
	ImportanceSamplingProtonsFar,
	ImportanceSamplingWeightInside,
	ImportanceSamplingWeightNear,
	ImportanceSamplingWeightFar,
	ImportanceSamplingWeightBound,
	ImportanceSamplingBoundProtonsIntracellular,
	ImportanceSamplingBoundProtonsExtracellular,
	ImportanceSamplingWeightBoundIntracellular,
	ImportanceSamplingWeightBoundExtracellular
};

std::string trim(std::string_view input) {
	const auto first = input.find_first_not_of(" \t\r\n");
	if (first == std::string_view::npos) {
		return {};
	}
	const auto last = input.find_last_not_of(" \t\r\n");
	return std::string(input.substr(first, last - first + 1));
}

// Aliases ignore case, whitespace and underscores. Keep this mapping in sync
// with the configuration reference when a future functional change adds keys.
std::string normalize_key(std::string_view key) {
	std::string normalized;
	normalized.reserve(key.size());
	for (char c : key) {
		if (std::isspace(static_cast<unsigned char>(c)) || c == '_') {
			continue;
		}
		normalized.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
	}
	return normalized;
}

const std::unordered_map<std::string, ConfigKey> kKeyLookup = {
	{"ionpr", ConfigKey::IonpRadius}, {"ionpradius", ConfigKey::IonpRadius}, {"rionp", ConfigKey::IonpRadius}, {"radiusionp", ConfigKey::IonpRadius},
	{"aggr", ConfigKey::AggRadius}, {"aggradius", ConfigKey::AggRadius}, {"aggregateradius", ConfigKey::AggRadius}, {"aggregater", ConfigKey::AggRadius}, {"ragg", ConfigKey::AggRadius}, {"radiusagg", ConfigKey::AggRadius},
	{"matrixx0", ConfigKey::MatrixX0}, {"matrixy0", ConfigKey::MatrixY0}, {"matrixz0", ConfigKey::MatrixZ0},
	{"origin", ConfigKey::Origin},
	{"matrixxlength", ConfigKey::MatrixXLength}, {"matrixylength", ConfigKey::MatrixYLength}, {"matrixzlength", ConfigKey::MatrixZLength},
	{"xlength", ConfigKey::MatrixXLength}, {"ylength", ConfigKey::MatrixYLength}, {"zlength", ConfigKey::MatrixZLength},
	{"length", ConfigKey::Length},
	{"matrixxnodes", ConfigKey::MatrixXNodes}, {"matrixynodes", ConfigKey::MatrixYNodes}, {"matrixznodes", ConfigKey::MatrixZNodes},
	{"xnodes", ConfigKey::MatrixXNodes}, {"ynodes", ConfigKey::MatrixYNodes}, {"znodes", ConfigKey::MatrixZNodes},
	{"nodes", ConfigKey::Nodes},
	{"ionpn", ConfigKey::IonpNumber}, {"ionpnum", ConfigKey::IonpNumber}, {"ionpnumber", ConfigKey::IonpNumber},
	{"aggn", ConfigKey::AggNumber}, {"aggnum", ConfigKey::AggNumber}, {"aggnumber", ConfigKey::AggNumber},
	{"ionppositionfile", ConfigKey::IonpPositionFile},
	{"aggpositionfile", ConfigKey::AggPositionFile},
	{"cellpositionfile", ConfigKey::CellPositionFile},
	{"cellloadfile", ConfigKey::CellLoadFile},
	{"protn", ConfigKey::ProtonNumber}, {"protnum", ConfigKey::ProtonNumber}, {"protnumber", ConfigKey::ProtonNumber},
	{"protonn", ConfigKey::ProtonNumber}, {"protonnum", ConfigKey::ProtonNumber}, {"protonnumber", ConfigKey::ProtonNumber},
	{"echotime", ConfigKey::EchoTime}, {"te", ConfigKey::EchoTime}, {"techo", ConfigKey::EchoTime},
	{"timestep", ConfigKey::TimeStep}, {"dt", ConfigKey::TimeStep}, {"deltat", ConfigKey::TimeStep},
	{"variabletime", ConfigKey::VariableTime}, {"variabletimestep", ConfigKey::VariableTime},
	{"surroundingvoxels", ConfigKey::SurroundingVoxels},
	{"surroundingvoxelsuperparticles", ConfigKey::SurroundingVoxelSuperparticles},
	{"echospacing", ConfigKey::EchoSpacing}, {"setime", ConfigKey::EchoSpacing},
	{"threads", ConfigKey::Threads}, {"numthreads", ConfigKey::Threads}, {"numberthreads", ConfigKey::Threads},
	{"vionp", ConfigKey::VelocityIonp}, {"velocityionp", ConfigKey::VelocityIonp},
	{"direction", ConfigKey::Direction}, {"movedirection", ConfigKey::Direction}, {"movementdirection", ConfigKey::Direction},
	{"sequence", ConfigKey::Sequence}, {"mrsequence", ConfigKey::Sequence},
	{"cellhist", ConfigKey::CellHist}, {"cellhistogram", ConfigKey::CellHist},
	{"radiushist", ConfigKey::RadiusHist}, {"radiushistogram", ConfigKey::RadiusHist},
	{"coating", ConfigKey::Coating}, {"shell", ConfigKey::Coating}, {"coatingthickness", ConfigKey::Coating},
	{"coatingd", ConfigKey::CoatingDiffusionCoefficient}, {"dcoating", ConfigKey::CoatingDiffusionCoefficient},
	{"coatingdiffusion", ConfigKey::CoatingDiffusionCoefficient}, {"coatingdiffusioncoefficient", ConfigKey::CoatingDiffusionCoefficient},
	{"celld", ConfigKey::CellDiffusionCoefficient}, {"dcell", ConfigKey::CellDiffusionCoefficient},
	{"celldiffusion", ConfigKey::CellDiffusionCoefficient}, {"celldiffusioncoefficient", ConfigKey::CellDiffusionCoefficient},
	{"coatingpermeability", ConfigKey::CoatingPermeability}, {"coatingp", ConfigKey::CoatingPermeability},
	{"coatingprobability", ConfigKey::CoatingPermeability}, {"coatingpermeationprobability", ConfigKey::CoatingPermeability},
	{"cellpermeability", ConfigKey::CellPermeability}, {"cellp", ConfigKey::CellPermeability},
	{"cellprobability", ConfigKey::CellPermeability}, {"cellpermeationprobability", ConfigKey::CellPermeability},
	{"cellularpermeability", ConfigKey::CellPermeability},
	{"boundprotn", ConfigKey::BoundProtonNumber}, {"boundprotnum", ConfigKey::BoundProtonNumber}, {"boundprotnumber", ConfigKey::BoundProtonNumber},
	{"boundprotonn", ConfigKey::BoundProtonNumber}, {"boundprotonnum", ConfigKey::BoundProtonNumber}, {"boundprotonnumber", ConfigKey::BoundProtonNumber},
	{"t2out", ConfigKey::BackgroundT2}, {"backgroundt2", ConfigKey::BackgroundT2},
	{"t2in", ConfigKey::CellT2}, {"cellt2", ConfigKey::CellT2},
	{"ionpr1", ConfigKey::IonpRadius1}, {"ionpradius1", ConfigKey::IonpRadius1}, {"rionp1", ConfigKey::IonpRadius1}, {"radiusionp1", ConfigKey::IonpRadius1},
	{"ionpper", ConfigKey::IonpPer}, {"ionpprobability", ConfigKey::IonpPer}, {"ionpdistributionprobability", ConfigKey::IonpPer},
	{"ionprstd", ConfigKey::IonpRadiusStd}, {"ionpradiusstd", ConfigKey::IonpRadiusStd}, {"rionpstd", ConfigKey::IonpRadiusStd}, {"radiusionpstd", ConfigKey::IonpRadiusStd},
	{"ionprstd1", ConfigKey::IonpRadiusStd1}, {"ionpradiusstd1", ConfigKey::IonpRadiusStd1}, {"rionpstd1", ConfigKey::IonpRadiusStd1}, {"radiusionpstd1", ConfigKey::IonpRadiusStd1},
	{"aggrstd", ConfigKey::AggRadiusStd}, {"aggradiusstd", ConfigKey::AggRadiusStd}, {"aggregateradiusstd", ConfigKey::AggRadiusStd}, {"aggregaterstd", ConfigKey::AggRadiusStd}, {"raggstd", ConfigKey::AggRadiusStd}, {"radiusaggstd", ConfigKey::AggRadiusStd},
	{"concentration", ConfigKey::Concentration}, {"conc", ConfigKey::Concentration},
	{"concentrationextern", ConfigKey::ConcentrationExternal}, {"concextern", ConfigKey::ConcentrationExternal}, {"concentrationext", ConfigKey::ConcentrationExternal}, {"concext", ConfigKey::ConcentrationExternal},
	{"magneticmoment", ConfigKey::MagneticMoment},
	{"cellorigin", ConfigKey::CellOrigin}, {"cellorig", ConfigKey::CellOrigin},
	{"voxelposition", ConfigKey::VoxelPosition}, {"voxelpos", ConfigKey::VoxelPosition},
	{"cubeorigin", ConfigKey::CubeOrigin}, {"cubeorig", ConfigKey::CubeOrigin},
	{"cubesize", ConfigKey::CubeSize},
	{"cellradius", ConfigKey::CellRadius}, {"cellrad", ConfigKey::CellRadius},
	{"cellconcentration", ConfigKey::CellConcentration}, {"cellconc", ConfigKey::CellConcentration},
	{"cellconcentrationstd", ConfigKey::CellConcentrationStd}, {"cellconcstd", ConfigKey::CellConcentrationStd},
	{"cellconcentrationdistribution", ConfigKey::CellConcentrationDist}, {"cellconcdist", ConfigKey::CellConcentrationDist},
	{"cellnumber", ConfigKey::CellNumber}, {"cellnum", ConfigKey::CellNumber},
	{"bfieldsolver", ConfigKey::BFieldSolver},
	{"importancesampling", ConfigKey::ImportanceSamplingEnabled}, {"isenabled", ConfigKey::ImportanceSamplingEnabled},
	{"isnearcutoff", ConfigKey::ImportanceSamplingCellNearCutoff},
	{"iscellnearcutoff", ConfigKey::ImportanceSamplingCellNearCutoff},
	{"isintracellularionpnearcutoff", ConfigKey::ImportanceSamplingIntracellularIonpNearCutoff},
	{"isintracellularionpfarcutoff", ConfigKey::ImportanceSamplingIntracellularIonpFarCutoff},
	{"isextracellularionpnearcutoff", ConfigKey::ImportanceSamplingExtracellularIonpNearCutoff},
	{"isextracellularionpfarcutoff", ConfigKey::ImportanceSamplingExtracellularIonpFarCutoff},
	{"isprotonsuniform", ConfigKey::ImportanceSamplingProposalProtonsUniform},
	{"isprotonscellinside", ConfigKey::ImportanceSamplingProposalProtonsCellInside},
	{"isprotonscellnear", ConfigKey::ImportanceSamplingProposalProtonsCellNear},
	{"isprotonscellfar", ConfigKey::ImportanceSamplingProposalProtonsCellFar},
	{"isprotonsintracellularionpnear", ConfigKey::ImportanceSamplingProposalProtonsIntracellularIonpNear},
	{"isprotonsintracellularionpintermediate", ConfigKey::ImportanceSamplingProposalProtonsIntracellularIonpIntermediate},
	{"isprotonsintracellularionpfar", ConfigKey::ImportanceSamplingProposalProtonsIntracellularIonpFar},
	{"isprotonsextracellularionpnear", ConfigKey::ImportanceSamplingProposalProtonsExtracellularIonpNear},
	{"isprotonsextracellularionpintermediate", ConfigKey::ImportanceSamplingProposalProtonsExtracellularIonpIntermediate},
	{"isprotonsextracellularionpfar", ConfigKey::ImportanceSamplingProposalProtonsExtracellularIonpFar},
	// Legacy population names retained for existing INI files.
	{"isionpinsidenearcutoff", ConfigKey::ImportanceSamplingIntracellularIonpNearCutoff},
	{"isionpinsidefarcutoff", ConfigKey::ImportanceSamplingIntracellularIonpFarCutoff},
	{"isionpoutsidenearcutoff", ConfigKey::ImportanceSamplingExtracellularIonpNearCutoff},
	{"isionpoutsidefarcutoff", ConfigKey::ImportanceSamplingExtracellularIonpFarCutoff},
	{"isprotonsinside", ConfigKey::ImportanceSamplingProtonsInside},
	{"isprotonsnear", ConfigKey::ImportanceSamplingProtonsNear},
	{"isprotonsfar", ConfigKey::ImportanceSamplingProtonsFar},
	{"isweightinside", ConfigKey::ImportanceSamplingWeightInside},
	{"isweightnear", ConfigKey::ImportanceSamplingWeightNear},
	{"isweightfar", ConfigKey::ImportanceSamplingWeightFar},
	{"isweightbound", ConfigKey::ImportanceSamplingWeightBound},
	{"isboundprotonsintracellularionp", ConfigKey::ImportanceSamplingBoundProtonsIntracellular},
	{"isboundprotonsextracellularionp", ConfigKey::ImportanceSamplingBoundProtonsExtracellular},
	{"isweightboundintracellularionp", ConfigKey::ImportanceSamplingWeightBoundIntracellular},
	{"isweightboundextracellularionp", ConfigKey::ImportanceSamplingWeightBoundExtracellular},
	// Legacy population names retained for existing INI files.
	{"isboundprotonsionpinside", ConfigKey::ImportanceSamplingBoundProtonsIntracellular},
	{"isboundprotonsionpoutside", ConfigKey::ImportanceSamplingBoundProtonsExtracellular},
	{"isweightboundionpinside", ConfigKey::ImportanceSamplingWeightBoundIntracellular},
	{"isweightboundionpoutside", ConfigKey::ImportanceSamplingWeightBoundExtracellular}
};

// std::stod accepts a valid numeric prefix; no consumed-length check is made.
// Preserve that behavior here, including its treatment of trailing text.
double parse_double(std::string_view value, const std::string& key, std::size_t line_no) {
	try {
		return std::stod(std::string(value));
	} catch (const std::exception& e) {
		throw std::runtime_error("Unable to parse numeric value for '" + key + "' at line " + std::to_string(line_no) + ": " + std::string(value));
	}
}

double parse_unit_interval_probability(std::string_view value, const std::string& key, std::size_t line_no) {
	const std::string text(value);
	std::size_t consumed = 0;
	double parsed = 0.;
	try {
		parsed = std::stod(text, &consumed);
	} catch (const std::exception&) {
		throw std::runtime_error("Unable to parse probability for '" + key + "' at line " +
			std::to_string(line_no) + ": " + text);
	}
	if (consumed != text.size() || !std::isfinite(parsed) || parsed < 0. || parsed > 1.) {
		throw std::runtime_error("'" + key + "' must be a probability in [0, 1] without percent notation at line " +
			std::to_string(line_no));
	}
	return parsed;
}

double parse_probability_or_percent(std::string_view value, const std::string& key, std::size_t line_no) {
	const double parsed = parse_double(value, key, line_no);
	if (parsed >= 0. && parsed <= 1.) {
		return parsed;
	}
	if (parsed > 1. && parsed <= 100.) {
		return parsed / 100.;
	}
	throw std::runtime_error("'" + key + "' must be a probability in [0, 1] or a percent in [0, 100] at line " + std::to_string(line_no));
}

int parse_int(std::string_view value, const std::string& key, std::size_t line_no) {
	try {
		return std::stoi(std::string(value));
	} catch (const std::exception& e) {
		throw std::runtime_error("Unable to parse integer value for '" + key + "' at line " + std::to_string(line_no) + ": " + std::string(value));
	}
}

std::size_t parse_non_negative_size(std::string_view value,
									const std::string& key,
									std::size_t line_no) {
	const std::string text(value);
	std::size_t consumed = 0;
	long long parsed = 0;
	try {
		parsed = std::stoll(text, &consumed);
	} catch (const std::exception&) {
		throw std::runtime_error("'" + key + "' must be a non-negative integer at line " +
			std::to_string(line_no) + ": " + text);
	}
	if (consumed != text.size() || parsed < 0 ||
		static_cast<unsigned long long>(parsed) >
			std::numeric_limits<std::size_t>::max()) {
		throw std::runtime_error("'" + key + "' must be a non-negative integer at line " +
			std::to_string(line_no) + ": " + text);
	}
	return static_cast<std::size_t>(parsed);
}

bool parse_yes_no(std::string_view value, const std::string& key, std::size_t line_no) {
	const auto normalized = normalize_key(value);
	if (normalized == "yes" || normalized == "ja" || normalized == "true") {
		return true;
	}
	if (normalized == "no" || normalized == "nein" || normalized == "false") {
		return false;
	}
	throw std::runtime_error("Expected Yes/No value for '" + key + "' at line " + std::to_string(line_no) + ": " + std::string(value));
}

// Decode the three-axis region suffix before mapping the value to its
// allocation/weight slot: cell, intracellular IONP, extracellular IONP.
bool parse_joint_importance_key(const std::string& normalized_key,
									bool& is_weight,
									std::size_t& region_index) {
	static const std::array<std::string, 3> cell_names = {
		"inside", "near", "far"
	};
	static const std::array<std::string, 3> proximity_names = {
		"near", "intermediate", "far"
	};
	static const std::array<std::string, 3> legacy_proximity_names = {
		"near", "middle", "far"
	};

	for (std::size_t cell = 0; cell < cell_names.size(); ++cell) {
		for (std::size_t intracellular = 0;
			 intracellular < proximity_names.size();
			 ++intracellular) {
			for (std::size_t extracellular = 0;
				 extracellular < proximity_names.size();
				 ++extracellular) {
				const std::string clear_suffix =
					"cell" + cell_names[cell] +
					"intracellularionp" + proximity_names[intracellular] +
					"extracellularionp" + proximity_names[extracellular];
				const std::string legacy_suffix =
					"cell" + cell_names[cell] +
					"ionpinside" + legacy_proximity_names[intracellular] +
					"ionpoutside" + legacy_proximity_names[extracellular];
				const std::size_t parsed_region_index =
					cell * 9u + intracellular * 3u + extracellular;
				if (normalized_key == "isprotons" + clear_suffix ||
					normalized_key == "isprotons" + legacy_suffix) {
					is_weight = false;
					region_index = parsed_region_index;
					return true;
				}
				if (normalized_key == "isweight" + clear_suffix ||
					normalized_key == "isweight" + legacy_suffix) {
					is_weight = true;
					region_index = parsed_region_index;
					return true;
				}
			}
		}
	}
	return false;
}

std::vector<double> parse_vector(std::string_view value, std::size_t expected, const std::string& key, std::size_t line_no) {
	std::istringstream iss{std::string(value)};
	std::vector<double> result;
	double element;
	while (iss >> element) {
		result.push_back(element);
	}
	if (result.size() != expected) {
		throw std::runtime_error("Expected " + std::to_string(expected) + " values for '" + key + "' at line " + std::to_string(line_no) + " but got " + std::to_string(result.size()));
	}
	return result;
}

ConfigKey to_config_key(const std::string& raw_key, std::size_t line_no) {
	const auto normalized = normalize_key(raw_key);
	const auto it = kKeyLookup.find(normalized);
	if (it == kKeyLookup.end()) {
		throw std::runtime_error("Unknown configuration key '" + raw_key + "' at line " + std::to_string(line_no));
	}
	return it->second;
}

// Validate counts, cutoffs, active weights and mutually exclusive sampling
// modes together; individual scalar conversion cannot check these relationships.
void validate_importance_sampling_config(const ImportanceSamplingConfig& config,
										 double cell_radius,
										 int bound_proton_number) {
	if (!config.enabled) {
		return;
	}

	const bool proposal_mode = config.overlapping_proposals_enabled();
	auto proposal_count = [&](FreeSamplingProposal proposal) {
		return config.proposal_protons[free_sampling_proposal_index(proposal)];
	};
	if (proposal_mode) {
		const bool has_legacy_free_values =
			config.protons_inside != 0 || config.protons_near != 0 ||
			config.protons_far != 0 || config.weight_inside != 0. ||
			config.weight_near != 0. || config.weight_far != 0. ||
			std::any_of(
				config.joint_protons.begin(),
				config.joint_protons.end(),
				[](int count) { return count != 0; }) ||
			std::any_of(
				config.joint_weights.begin(),
				config.joint_weights.end(),
				[](double weight) { return weight != 0.; });
		if (config.joint_regions_configured || has_legacy_free_values) {
			throw std::runtime_error(
				"Overlapping free-proton proposals cannot be mixed with legacy "
				"cell-only counts, joint free-region keys, or manual free-region weights");
		}
		for (std::size_t index = 0;
			 index < kFreeSamplingProposalCount;
			 ++index) {
			if (config.proposal_protons[index] < 0) {
				throw std::runtime_error(
					"Overlapping importance-sampling proposal counts must be non-negative");
			}
		}
		if (proposal_count(FreeSamplingProposal::Uniform) <= 0) {
			throw std::runtime_error(
				"Overlapping importance sampling requires a positive IS_protons_uniform count");
		}
	}

	bool uses_non_far_cell_region = proposal_mode &&
		(proposal_count(FreeSamplingProposal::CellInside) > 0 ||
		 proposal_count(FreeSamplingProposal::CellNear) > 0);
	bool uses_near_cell_region = proposal_mode &&
		proposal_count(FreeSamplingProposal::CellNear) > 0;
	bool uses_intracellular_ionp_region = proposal_mode &&
		(proposal_count(FreeSamplingProposal::IntracellularIonpNear) > 0 ||
		 proposal_count(FreeSamplingProposal::IntracellularIonpIntermediate) > 0 ||
		 proposal_count(FreeSamplingProposal::IntracellularIonpFar) > 0);
	bool uses_extracellular_ionp_region = proposal_mode &&
		(proposal_count(FreeSamplingProposal::ExtracellularIonpNear) > 0 ||
		 proposal_count(FreeSamplingProposal::ExtracellularIonpIntermediate) > 0 ||
		 proposal_count(FreeSamplingProposal::ExtracellularIonpFar) > 0);
	for (std::size_t index = 0; index < kFreeSampleRegionCount; ++index) {
		if (config.joint_protons[index] < 0) {
			throw std::runtime_error(
				"Importance-sampling proton counts must be non-negative");
		}
		if (config.joint_weights[index] < 0.) {
			throw std::runtime_error(
				"Importance-sampling weights must be non-negative");
		}
		if (config.joint_weights[index] > 0. &&
			config.joint_protons[index] == 0) {
			throw std::runtime_error(
				"Positive importance-sampling weight has no samples for region " +
				sample_region_name(static_cast<ProtonSampleRegion>(index)));
		}
		const JointSampleRegion region =
			decode_joint_sample_region(
				static_cast<ProtonSampleRegion>(index));
		if (config.joint_protons[index] > 0 ||
			config.joint_weights[index] > 0.) {
			uses_non_far_cell_region =
				uses_non_far_cell_region ||
				region.cell != CellSampleRegion::Far;
			uses_near_cell_region =
				uses_near_cell_region ||
				region.cell == CellSampleRegion::Near;
			uses_intracellular_ionp_region =
				uses_intracellular_ionp_region ||
				region.intracellular_ionp != IonpDistanceRegion::Far;
			uses_extracellular_ionp_region =
				uses_extracellular_ionp_region ||
				region.extracellular_ionp != IonpDistanceRegion::Far;
		}
	}

	if (uses_non_far_cell_region && cell_radius <= 0.) {
		throw std::runtime_error(
			"Inside/near-cell importance strata require a positive Cell Radius");
	}
	const bool cell_near_cutoff_is_active =
		uses_near_cell_region || config.cell_near_cutoff > 0.;
	if (cell_radius > 0. &&
		cell_near_cutoff_is_active &&
		config.cell_near_cutoff <= cell_radius) {
		throw std::runtime_error(
			"IS_cell_near_cutoff must be larger than Cell Radius");
	}
	auto validate_ionp_cutoffs = [](bool used,
									bool enabled,
									double near_cutoff,
									double far_cutoff,
									const char* population) {
		if (used && !enabled) {
			throw std::runtime_error(
				std::string("Requested ") + population +
				" IONP regions require their cutoff pair");
		}
		if (enabled &&
			(near_cutoff <= 0. || far_cutoff <= near_cutoff)) {
			throw std::runtime_error(
				std::string("Importance-sampling ") + population +
				" IONP cutoffs must satisfy 0 < near < far");
		}
	};
	validate_ionp_cutoffs(
		uses_intracellular_ionp_region,
		config.intracellular_ionp_sampling_enabled(),
		config.intracellular_ionp_near_cutoff,
		config.intracellular_ionp_far_cutoff,
		"intracellular");
	validate_ionp_cutoffs(
		uses_extracellular_ionp_region,
		config.extracellular_ionp_sampling_enabled(),
		config.extracellular_ionp_near_cutoff,
		config.extracellular_ionp_far_cutoff,
		"extracellular");

	if (config.total_free_protons() <= 0) {
		throw std::runtime_error("Importance-sampling mode requires at least one free proton sample");
	}
	if (config.weight_bound < 0. ||
		config.weight_bound_intracellular < 0. ||
		config.weight_bound_extracellular < 0. ||
		config.bound_protons_intracellular < 0 ||
		config.bound_protons_extracellular < 0) {
		throw std::runtime_error(
			"Importance-sampling bound counts and weights must be non-negative");
	}
	if (config.split_bound_configured) {
		if (bound_proton_number != 0 || config.weight_bound != 0.) {
			throw std::runtime_error(
				"Bound_Prot_n/IS_weight_bound cannot be mixed with population-specific bound importance-sampling keys");
		}
		if ((config.bound_protons_intracellular > 0 &&
			 config.weight_bound_intracellular <= 0.) ||
			(config.bound_protons_extracellular > 0 &&
			 config.weight_bound_extracellular <= 0.)) {
			throw std::runtime_error(
				"Each active population-specific bound sample group requires a positive weight");
		}
		if ((config.weight_bound_intracellular > 0. &&
			 config.bound_protons_intracellular == 0) ||
			(config.weight_bound_extracellular > 0. &&
			 config.bound_protons_extracellular == 0)) {
			throw std::runtime_error(
				"A positive population-specific bound weight requires bound samples in that population");
		}
	}
	else if (bound_proton_number > 0 && config.weight_bound <= 0.) {
		throw std::runtime_error(
			"IS_weight_bound must be positive when Bound_Prot_n > 0 in importance-sampling mode");
	}
	else if (config.weight_bound > 0. && bound_proton_number == 0) {
		throw std::runtime_error(
			"A positive IS_weight_bound requires Bound_Prot_n > 0");
	}

	const double bound_weight_total = config.split_bound_configured
		? config.weight_bound_intracellular + config.weight_bound_extracellular
		: (bound_proton_number > 0 ? config.weight_bound : 0.);
	if (!has_explicit_free_region_weights(config) &&
		bound_weight_total == 0.) {
		return;
	}
	const double active_weight_total = std::accumulate(
		config.joint_weights.begin(),
		config.joint_weights.end(),
		0.0) + bound_weight_total;
	if (active_weight_total <= 0.) {
		throw std::runtime_error("Importance-sampling weights must sum to a positive value");
	}
}

} // namespace

void Data::set_output_context(std::filesystem::path tmp_dir, std::filesystem::path log_dir) {
	output_tmp_dir = std::move(tmp_dir);
	output_log_dir = std::move(log_dir);
	first_tmp_output_index = output::next_file_index(output_tmp_dir);
}

// Read key/value pairs first. Dependent defaults and cross-option validation
// follow the parsing loop so they use the complete configuration.
void Data::readData (  )
{
  std::cout << "Entering readFunc...." << std::endl;
	std::ifstream file_in(filename, std::ios::in);

	if (!file_in.is_open()) {
		throw std::runtime_error("Unable to open input file: " + filename);
	}

	std::string line_in;
	std::size_t line_no = 0;
	bool legacy_free_count_key_configured = false;
	bool manual_free_weight_key_configured = false;
	while (std::getline(file_in, line_in)) {
		++line_no;
		const std::string trimmed_line = trim(line_in);
		if (trimmed_line.empty() || trimmed_line[0] == '#') {
			continue;
		}

		const std::size_t pos = trimmed_line.find("=");
		if (pos == std::string::npos) {
			throw std::runtime_error("Missing '=' in config file at line " + std::to_string(line_no) + ": " + trimmed_line);
		}

		const std::string raw_key = trim(trimmed_line.substr(0, pos));
		const std::string value_str = trim(trimmed_line.substr(pos + 1));

		if (value_str.empty()) {
			throw std::runtime_error("Missing value for key '" + raw_key + "' at line " + std::to_string(line_no));
		}

		bool joint_weight = false;
		std::size_t joint_region_index = 0;
		if (parse_joint_importance_key(
				normalize_key(raw_key),
				joint_weight,
				joint_region_index)) {
			importance_sampling.joint_regions_configured = true;
			if (joint_weight) {
				importance_sampling.joint_weights[joint_region_index] =
					parse_double(value_str, raw_key, line_no);
			}
			else {
				importance_sampling.joint_protons[joint_region_index] =
					parse_int(value_str, raw_key, line_no);
			}
			continue;
		}

		const ConfigKey key = to_config_key(raw_key, line_no);

		switch (key) {
			case ConfigKey::IonpRadius:
				radius_ionp = parse_double(value_str, raw_key, line_no);
				break;
			case ConfigKey::AggRadius:
				radius_agg = parse_double(value_str, raw_key, line_no);
				break;
			case ConfigKey::MatrixX0:
				matrix_x0 = parse_double(value_str, raw_key, line_no);
				break;
			case ConfigKey::MatrixY0:
				matrix_y0 = parse_double(value_str, raw_key, line_no);
				break;
			case ConfigKey::MatrixZ0:
				matrix_z0 = parse_double(value_str, raw_key, line_no);
				break;
			case ConfigKey::Origin: {
				const double origin = parse_double(value_str, raw_key, line_no);
				matrix_x0 = origin;
				matrix_y0 = origin;
				matrix_z0 = origin;
				break;
			}
			case ConfigKey::MatrixXLength:
				matrix_xlength = parse_double(value_str, raw_key, line_no);
				break;
			case ConfigKey::MatrixYLength:
				matrix_ylength = parse_double(value_str, raw_key, line_no);
				break;
			case ConfigKey::MatrixZLength:
				matrix_zlength = parse_double(value_str, raw_key, line_no);
				break;
			case ConfigKey::Length: {
				const double length = parse_double(value_str, raw_key, line_no);
				matrix_xlength = length;
				matrix_ylength = length;
				matrix_zlength = length;
				break;
			}
			case ConfigKey::MatrixXNodes:
				matrix_nxnodes = parse_int(value_str, raw_key, line_no);
				break;
			case ConfigKey::MatrixYNodes:
				matrix_nynodes = parse_int(value_str, raw_key, line_no);
				break;
			case ConfigKey::MatrixZNodes:
				matrix_nznodes = parse_int(value_str, raw_key, line_no);
				break;
			case ConfigKey::Nodes: {
				const int nodes = parse_int(value_str, raw_key, line_no);
				matrix_nxnodes = nodes;
				matrix_nynodes = nodes;
				matrix_nznodes = nodes;
				break;
			}
			case ConfigKey::IonpNumber:
				ionp_number = parse_int(value_str, raw_key, line_no);
				break;
			case ConfigKey::AggNumber:
				agg_number = parse_int(value_str, raw_key, line_no);
				break;
			case ConfigKey::IonpPositionFile:
				ionp_position_file = value_str;
				break;
			case ConfigKey::AggPositionFile:
				agg_position_file = value_str;
				break;
			case ConfigKey::CellPositionFile:
				cell_position_file = value_str;
				break;
			case ConfigKey::CellLoadFile:
				cell_load_file = value_str;
				break;
			case ConfigKey::ProtonNumber:
				proton_number = parse_int(value_str, raw_key, line_no);
				break;
			case ConfigKey::EchoTime:
				echo_time = parse_double(value_str, raw_key, line_no);
				break;
			case ConfigKey::TimeStep:
				time_step = parse_double(value_str, raw_key, line_no);
				break;
			case ConfigKey::VariableTime:
				variable_time = parse_yes_no(value_str, raw_key, line_no);
				variable_time_explicit = true;
				break;
			case ConfigKey::SurroundingVoxels:
				surrounding_voxels = parse_yes_no(value_str, raw_key, line_no);
				break;
			case ConfigKey::SurroundingVoxelSuperparticles:
				surrounding_voxel_superparticles =
					parse_non_negative_size(value_str, raw_key, line_no);
				break;
			case ConfigKey::EchoSpacing:
				echo_spacing = parse_double(value_str, raw_key, line_no);
				break;
			case ConfigKey::Threads:
				numThreads = parse_int(value_str, raw_key, line_no);
				break;
			case ConfigKey::VelocityIonp:
				velocityIONP = parse_double(value_str, raw_key, line_no);
				break;
			case ConfigKey::Direction: {
				const auto values = parse_vector(value_str, 3, raw_key, line_no);
				directionIONP.set(values[0], values[1], values[2]);
				break;
			}
			case ConfigKey::Sequence: {
				const auto seq_value = normalize_key(value_str);
				if (seq_value == "fid") {
					sequence_num = 0;
					std::cout << "Sequence FID is used" << std::endl;
				}
				else if (seq_value == "se" || seq_value == "spinecho") {
					sequence_num = 1;
					std::cout << "Sequence Spin Echo is used" << std::endl;
				}
				else if (seq_value == "mse" || seq_value == "multispinecho") {
					sequence_num = 2;
					std::cout << "Sequence Multi Spin Echo is used" << std::endl;
				}
				else{
					throw std::runtime_error("Unknown sequence: " + value_str + " at line " + std::to_string(line_no));
				}
				break;
			}
			case ConfigKey::CellHist: {
				const auto hist_value = normalize_key(value_str);
				if (hist_value == "yes" || hist_value == "ja") {
					Cell_hist = true;
				}
				else if (hist_value == "no" || hist_value == "nein") {
					Cell_hist = false;
				}
				else{
					throw std::runtime_error("Neither Yes nor No provided for Cell histogram at line " + std::to_string(line_no) + ": " + value_str);
				}
				break;
			}
			case ConfigKey::RadiusHist: {
				const auto hist_value = normalize_key(value_str);
				if (hist_value == "yes" || hist_value == "ja") {
					Radius_hist = true;
				}
				else if (hist_value == "no" || hist_value == "nein") {
					Radius_hist = false;
				}
				else{
					throw std::runtime_error("Neither Yes nor No provided for Radius histogram at line " + std::to_string(line_no) + ": " + value_str);
				}
				break;
			}
			case ConfigKey::Coating:
				thickness_IONPcoating = parse_double(value_str, raw_key, line_no);
				if (thickness_IONPcoating < 0.) {
					throw std::runtime_error(
						"Configured IONP coating outer radius must be non-negative at line " +
						std::to_string(line_no));
				}
				break;
			case ConfigKey::CoatingDiffusionCoefficient:
				coating_diffusion_coefficient = parse_double(value_str, raw_key, line_no);
				if (coating_diffusion_coefficient < 0.) {
					throw std::runtime_error("Coating diffusion coefficient must be non-negative at line " + std::to_string(line_no));
				}
				coating_diffusion_coefficient_explicit = true;
				break;
			case ConfigKey::CellDiffusionCoefficient:
				cell_diffusion_coefficient = parse_double(value_str, raw_key, line_no);
				if (cell_diffusion_coefficient < 0.) {
					throw std::runtime_error("Cell diffusion coefficient must be non-negative at line " + std::to_string(line_no));
				}
				cell_diffusion_coefficient_explicit = true;
				break;
			case ConfigKey::CoatingPermeability:
				coating_permeability = parse_double(value_str, raw_key, line_no);
				if (coating_permeability < 0. || coating_permeability > 1.) {
					throw std::runtime_error("Coating permeability must be a probability in [0, 1] at line " + std::to_string(line_no));
				}
				coating_permeability_explicit = true;
				break;
			case ConfigKey::CellPermeability:
				cell_permeability = parse_unit_interval_probability(value_str, raw_key, line_no);
				cell_permeability_explicit = true;
				break;
			case ConfigKey::BoundProtonNumber:
				bound_proton_number = parse_int(value_str, raw_key, line_no);
				if (bound_proton_number < 0) {
					throw std::runtime_error(
						"Bound_Prot_n must be non-negative at line " +
						std::to_string(line_no));
				}
				break;
			case ConfigKey::BackgroundT2:
				background_T2 = parse_double(value_str, raw_key, line_no);
				break;
			case ConfigKey::CellT2:
				cell_T2 = parse_double(value_str, raw_key, line_no);
				break;
			case ConfigKey::IonpRadius1:
				radius_ionp_1 = parse_double(value_str, raw_key, line_no);
				break;
			case ConfigKey::IonpPer:
				ionp_per = parse_probability_or_percent(value_str, raw_key, line_no);
				break;
			case ConfigKey::IonpRadiusStd:
				radius_ionp_std = parse_double(value_str, raw_key, line_no);
				break;
			case ConfigKey::IonpRadiusStd1:
				radius_ionp_std_1 = parse_double(value_str, raw_key, line_no);
				break;
			case ConfigKey::AggRadiusStd:
				radius_agg_std = parse_double(value_str, raw_key, line_no);
				break;
			case ConfigKey::Concentration:
				concentration = parse_double(value_str, raw_key, line_no);
				break;
			case ConfigKey::ConcentrationExternal:
				concentration_ext = parse_double(value_str, raw_key, line_no);
				break;
			case ConfigKey::MagneticMoment:
				mag_moment = parse_double(value_str, raw_key, line_no);
				break;
			case ConfigKey::CellOrigin: {
				const auto values = parse_vector(value_str, 3, raw_key, line_no);
				Cell_Origin.set(values[0], values[1], values[2]);
				break;
			}
			case ConfigKey::VoxelPosition: {
				const auto values = parse_vector(value_str, 3, raw_key, line_no);
				Voxel_Pos.set(values[0], values[1], values[2]);
				break;
			}
			case ConfigKey::CubeOrigin: {
				const auto values = parse_vector(value_str, 3, raw_key, line_no);
				Cube_Origin.set(values[0], values[1], values[2]);
				break;
			}
			case ConfigKey::CubeSize: {
				const auto values = parse_vector(value_str, 3, raw_key, line_no);
				Cube_size.set(values[0], values[1], values[2]);
				break;
			}
			case ConfigKey::CellRadius:
				Cell_radius = parse_double(value_str, raw_key, line_no);
				break;
			case ConfigKey::CellConcentration:
				Cell_concentration = parse_double(value_str, raw_key, line_no);
				break;
			case ConfigKey::CellConcentrationStd:
				Cell_concentration_std = parse_double(value_str, raw_key, line_no);
				break;
			case ConfigKey::CellConcentrationDist: {
				const auto dist_value = normalize_key(value_str);
				if (dist_value == "lognormal" || dist_value == "lognorm") {
					Cell_concentration_dist = 0;
					std::cout << "Lognormal distribution for cell concentration used!" << std::endl;
				}
				else if (dist_value == "normal" || dist_value == "norm") {
					Cell_concentration_dist = 1;
					std::cout << "Normal distribution for cell concentration used!" << std::endl;
				}
				else{
					throw std::runtime_error("Unknown distribution for cell concentration at line " + std::to_string(line_no) + ": " + value_str);
				}
				break;
			}
			case ConfigKey::CellNumber:
				Cell_number = parse_int(value_str, raw_key, line_no);
				break;
			case ConfigKey::BFieldSolver: {
				const std::string solver = normalize_key(value_str);
				if (solver == "current") {
					bfield_solver = BFieldSolver::Current;
				}
				else if (solver == "correctedhybrid" || solver == "hybrid") {
					bfield_solver = BFieldSolver::CorrectedHybrid;
				}
				else {
					throw std::runtime_error(
						"Unknown BField_solver '" + value_str +
						"' at line " + std::to_string(line_no) +
						"; expected current or corrected_hybrid");
				}
				break;
			}
			case ConfigKey::ImportanceSamplingEnabled:
				importance_sampling.enabled = parse_yes_no(value_str, raw_key, line_no);
				break;
			case ConfigKey::ImportanceSamplingCellNearCutoff:
				importance_sampling.cell_near_cutoff = parse_double(value_str, raw_key, line_no);
				break;
			case ConfigKey::ImportanceSamplingIntracellularIonpNearCutoff:
				importance_sampling.intracellular_ionp_near_cutoff =
					parse_double(value_str, raw_key, line_no);
				break;
			case ConfigKey::ImportanceSamplingIntracellularIonpFarCutoff:
				importance_sampling.intracellular_ionp_far_cutoff =
					parse_double(value_str, raw_key, line_no);
				break;
			case ConfigKey::ImportanceSamplingExtracellularIonpNearCutoff:
				importance_sampling.extracellular_ionp_near_cutoff =
					parse_double(value_str, raw_key, line_no);
				break;
			case ConfigKey::ImportanceSamplingExtracellularIonpFarCutoff:
				importance_sampling.extracellular_ionp_far_cutoff =
					parse_double(value_str, raw_key, line_no);
				break;
			case ConfigKey::ImportanceSamplingProposalProtonsUniform:
				importance_sampling.proposal_protons[
					free_sampling_proposal_index(FreeSamplingProposal::Uniform)] =
					parse_int(value_str, raw_key, line_no);
				importance_sampling.mode =
					ImportanceSamplingMode::OverlappingProposal;
				break;
			case ConfigKey::ImportanceSamplingProposalProtonsCellInside:
				importance_sampling.proposal_protons[
					free_sampling_proposal_index(FreeSamplingProposal::CellInside)] =
					parse_int(value_str, raw_key, line_no);
				importance_sampling.mode =
					ImportanceSamplingMode::OverlappingProposal;
				break;
			case ConfigKey::ImportanceSamplingProposalProtonsCellNear:
				importance_sampling.proposal_protons[
					free_sampling_proposal_index(FreeSamplingProposal::CellNear)] =
					parse_int(value_str, raw_key, line_no);
				importance_sampling.mode =
					ImportanceSamplingMode::OverlappingProposal;
				break;
			case ConfigKey::ImportanceSamplingProposalProtonsCellFar:
				importance_sampling.proposal_protons[
					free_sampling_proposal_index(FreeSamplingProposal::CellFar)] =
					parse_int(value_str, raw_key, line_no);
				importance_sampling.mode =
					ImportanceSamplingMode::OverlappingProposal;
				break;
			case ConfigKey::ImportanceSamplingProposalProtonsIntracellularIonpNear:
				importance_sampling.proposal_protons[
					free_sampling_proposal_index(
						FreeSamplingProposal::IntracellularIonpNear)] =
					parse_int(value_str, raw_key, line_no);
				importance_sampling.mode =
					ImportanceSamplingMode::OverlappingProposal;
				break;
			case ConfigKey::ImportanceSamplingProposalProtonsIntracellularIonpIntermediate:
				importance_sampling.proposal_protons[
					free_sampling_proposal_index(
						FreeSamplingProposal::IntracellularIonpIntermediate)] =
					parse_int(value_str, raw_key, line_no);
				importance_sampling.mode =
					ImportanceSamplingMode::OverlappingProposal;
				break;
			case ConfigKey::ImportanceSamplingProposalProtonsIntracellularIonpFar:
				importance_sampling.proposal_protons[
					free_sampling_proposal_index(
						FreeSamplingProposal::IntracellularIonpFar)] =
					parse_int(value_str, raw_key, line_no);
				importance_sampling.mode =
					ImportanceSamplingMode::OverlappingProposal;
				break;
			case ConfigKey::ImportanceSamplingProposalProtonsExtracellularIonpNear:
				importance_sampling.proposal_protons[
					free_sampling_proposal_index(
						FreeSamplingProposal::ExtracellularIonpNear)] =
					parse_int(value_str, raw_key, line_no);
				importance_sampling.mode =
					ImportanceSamplingMode::OverlappingProposal;
				break;
			case ConfigKey::ImportanceSamplingProposalProtonsExtracellularIonpIntermediate:
				importance_sampling.proposal_protons[
					free_sampling_proposal_index(
						FreeSamplingProposal::ExtracellularIonpIntermediate)] =
					parse_int(value_str, raw_key, line_no);
				importance_sampling.mode =
					ImportanceSamplingMode::OverlappingProposal;
				break;
			case ConfigKey::ImportanceSamplingProposalProtonsExtracellularIonpFar:
				importance_sampling.proposal_protons[
					free_sampling_proposal_index(
						FreeSamplingProposal::ExtracellularIonpFar)] =
					parse_int(value_str, raw_key, line_no);
				importance_sampling.mode =
					ImportanceSamplingMode::OverlappingProposal;
				break;
			case ConfigKey::ImportanceSamplingProtonsInside:
				importance_sampling.protons_inside = parse_int(value_str, raw_key, line_no);
				legacy_free_count_key_configured = true;
				break;
			case ConfigKey::ImportanceSamplingProtonsNear:
				importance_sampling.protons_near = parse_int(value_str, raw_key, line_no);
				legacy_free_count_key_configured = true;
				break;
			case ConfigKey::ImportanceSamplingProtonsFar:
				importance_sampling.protons_far = parse_int(value_str, raw_key, line_no);
				legacy_free_count_key_configured = true;
				break;
			case ConfigKey::ImportanceSamplingWeightInside:
				importance_sampling.weight_inside = parse_double(value_str, raw_key, line_no);
				manual_free_weight_key_configured = true;
				break;
			case ConfigKey::ImportanceSamplingWeightNear:
				importance_sampling.weight_near = parse_double(value_str, raw_key, line_no);
				manual_free_weight_key_configured = true;
				break;
			case ConfigKey::ImportanceSamplingWeightFar:
				importance_sampling.weight_far = parse_double(value_str, raw_key, line_no);
				manual_free_weight_key_configured = true;
				break;
			case ConfigKey::ImportanceSamplingWeightBound:
				importance_sampling.weight_bound = parse_double(value_str, raw_key, line_no);
				break;
			case ConfigKey::ImportanceSamplingBoundProtonsIntracellular:
				importance_sampling.bound_protons_intracellular =
					parse_int(value_str, raw_key, line_no);
				importance_sampling.split_bound_configured = true;
				break;
			case ConfigKey::ImportanceSamplingBoundProtonsExtracellular:
				importance_sampling.bound_protons_extracellular =
					parse_int(value_str, raw_key, line_no);
				importance_sampling.split_bound_configured = true;
				break;
			case ConfigKey::ImportanceSamplingWeightBoundIntracellular:
				importance_sampling.weight_bound_intracellular =
					parse_double(value_str, raw_key, line_no);
				importance_sampling.split_bound_configured = true;
				break;
			case ConfigKey::ImportanceSamplingWeightBoundExtracellular:
				importance_sampling.weight_bound_extracellular =
					parse_double(value_str, raw_key, line_no);
				importance_sampling.split_bound_configured = true;
				break;
		}
	}
	file_in.close();

	if (!variable_time_explicit) {
		variable_time = (time_step == 0.);
	}

	if (surrounding_voxel_superparticles > 0 && !surrounding_voxels) {
		throw std::runtime_error(
			"Surrounding_voxel_superparticles > 0 requires Surrounding_voxels = Yes");
	}
	if (surrounding_voxel_superparticles > 0 && mag_moment != 0.) {
		throw std::runtime_error(
			"Surrounding_voxel_superparticles > 0 requires Magnetic Moment = 0");
	}

	if (importance_sampling.mode == ImportanceSamplingMode::OverlappingProposal &&
		(legacy_free_count_key_configured ||
		 manual_free_weight_key_configured ||
		 importance_sampling.joint_regions_configured)) {
		throw std::runtime_error(
			"Overlapping IS_protons_* proposal keys cannot be mixed with legacy "
			"cell-only counts, joint free-region keys, or manual free-region weights");
	}

	importance_sampling.materialize_legacy_regions();
	validate_importance_sampling_config(importance_sampling, Cell_radius, bound_proton_number);
	effective_proton_number = proton_number;

  std::cout<<"\n"<<"IONP radius: "<< radius_ionp <<std::endl;
  std::cout<<"IONP number: "<< ionp_number <<std::endl;
  if (ionp_position_file) {
    std::cout<<"IONP_position_file: "<< *ionp_position_file <<std::endl;
  }
  if (agg_position_file) {
    std::cout<<"Agg_position_file: "<< *agg_position_file <<std::endl;
  }
  if (cell_position_file) {
    std::cout<<"Cell_position_file: "<< *cell_position_file <<std::endl;
  }
  if (cell_load_file) {
    std::cout<<"Cell_load_file: "<< *cell_load_file <<std::endl;
  }
  std::cout<<"Configured IONP coating outer radius: "<< thickness_IONPcoating <<std::endl;
  if (coating_diffusion_coefficient_explicit) {
    std::cout<<"IONP coating diffusion coefficient: "<< coating_diffusion_coefficient <<std::endl;
  }
  else {
    std::cout<<"IONP coating diffusion coefficient: same as outside ("<< D_D_CONST <<")"<<std::endl;
  }
  if (coating_permeability_explicit) {
    std::cout<<"IONP coating permeability: "<< coating_permeability <<std::endl;
  }
  else {
    std::cout<<"IONP coating permeability: fully permeable"<<std::endl;
  }
  if (cell_diffusion_coefficient_explicit) {
    std::cout<<"Cell diffusion coefficient: "<< cell_diffusion_coefficient <<std::endl;
  }
  else {
    std::cout<<"Cell diffusion coefficient: same as outside ("<< D_D_CONST <<")"<<std::endl;
  }
  std::cout<<"IONP radius std (var. radius): "<< radius_ionp_std <<std::endl;
  std::cout<<"Aggregate radius: "<< radius_agg <<std::endl;
  std::cout<<"Aggregate number: "<< agg_number <<std::endl;
  std::cout<<"Aggregate radius std (var. radius): "<< radius_agg_std <<std::endl;
  std::cout<<"Proton number: "<< proton_number <<std::endl;
  std::cout<<"Bound Proton number: "<< bound_proton_number <<"\n" << std::endl;
  if (importance_sampling.enabled) {
    std::cout << "Importance sampling enabled" << std::endl;
	    std::cout << "IS cell near cutoff: " << importance_sampling.cell_near_cutoff << std::endl;
	    if (importance_sampling.intracellular_ionp_sampling_enabled()) {
	      std::cout << "IS intracellular IONP surface cutoffs (near / far): "
	                << importance_sampling.intracellular_ionp_near_cutoff << " / "
	                << importance_sampling.intracellular_ionp_far_cutoff << std::endl;
	    }
	    if (importance_sampling.extracellular_ionp_sampling_enabled()) {
	      std::cout << "IS extracellular IONP surface cutoffs (near / far): "
	                << importance_sampling.extracellular_ionp_near_cutoff << " / "
	                << importance_sampling.extracellular_ionp_far_cutoff << std::endl;
    }
	    std::cout << "IS effective free proton count: "
	              << importance_sampling.total_free_protons() << std::endl;
	    if (importance_sampling.overlapping_proposals_enabled()) {
	      std::cout << "IS requested overlapping proposal counts:" << std::endl;
	      for (std::size_t proposal_index = 0;
	           proposal_index < kFreeSamplingProposalCount;
	           ++proposal_index) {
	        const int count = importance_sampling.proposal_protons[proposal_index];
	        if (count == 0) {
	          continue;
	        }
	        std::cout << "  "
	                  << free_sampling_proposal_name(
	                       static_cast<FreeSamplingProposal>(proposal_index))
	                  << ": requested samples=" << count << std::endl;
	      }
	    }
	    else {
	      for (std::size_t region_index = 0;
	           region_index < kFreeSampleRegionCount;
	           ++region_index) {
	        if (importance_sampling.joint_protons[region_index] == 0 &&
	            importance_sampling.joint_weights[region_index] == 0.) {
	          continue;
	        }
	        std::cout << "  "
	                  << sample_region_name(static_cast<ProtonSampleRegion>(region_index))
	                  << ": samples=" << importance_sampling.joint_protons[region_index]
	                  << ", weight=" << importance_sampling.joint_weights[region_index]
	                  << std::endl;
	      }
	    }
    if (importance_sampling.split_bound_configured) {
	      std::cout << "IS bound samples per coated IONP (intracellular / extracellular): "
	                << importance_sampling.bound_protons_intracellular << " / "
	                << importance_sampling.bound_protons_extracellular << std::endl;
	      std::cout << "IS bound weights (intracellular / extracellular): "
	                << importance_sampling.weight_bound_intracellular << " / "
	                << importance_sampling.weight_bound_extracellular << std::endl;
    }
    else {
      std::cout << "IS common bound weight: "
                << importance_sampling.weight_bound << std::endl;
    }
    std::cout << std::endl;
  }

  if ( mag_moment != 0. ){
    std::cout<<"Magnetic Moment: "<< mag_moment << std::endl;
  }

  if ( radius_ionp_1 != 0. && radius_ionp != 0. ){
    std::cout<<"Bimodal radii distribution -> here for Resovist" << std::endl;
    std::cout<<"Primary radius probability: " << ionp_per << ", secondary radius probability: " << (1. - ionp_per) << std::endl;
    std::cout<<"IONP radius 1: " << radius_ionp_1 << std::endl;
    std::cout<<"IONP radius 1 std: " << radius_ionp_std_1 << "\n" << std::endl;
  }

  std::cout<<"Concentration: "<< concentration <<" µmol/l\n" << std::endl;

  if ( concentration_ext != 0. ){
    std::cout<<"Concentration of IONPs not in cells: "<< concentration_ext <<" µmol/l\n" << std::endl;
  }

  std::cout<<"Set T2 outside of cell: "<< background_T2 << std::endl;
  std::cout<<"Set T2 inside of cell:  "<< cell_T2 <<"\n" << std::endl;

  std::cout<<"Origin: "<<matrix_x0<<" "<<matrix_y0<<" "<<matrix_z0<<std::endl;
  std::cout<<"Nodes: "<<matrix_nxnodes<<" "<<matrix_nynodes<<" "<<matrix_nznodes<<std::endl;
  std::cout<<"Dimensions given: "<<matrix_xlength<<" "<<matrix_ylength<<" "<<matrix_zlength<<std::endl;

  if( matrix_xlength == 0. || matrix_ylength == 0. || matrix_zlength == 0.){
    std::cout<<"No Simulation Volume size given ... using constant volume fraction "<<std::endl;
    std::cout<<"Volume fraction: "<<D_VOL_FRAC<<std::endl;

    double mtx_dimensions;

    if (agg_number == 0. || radius_agg == 0.){
      mtx_dimensions = pow(ionp_number*(4./3.)*(D_PI)*(1/D_VOL_FRAC), 1./3.)*radius_ionp;
    }
    else {
      mtx_dimensions = pow(ionp_number*agg_number*(4./3.)*(D_PI)*(1/D_VOL_FRAC), 1./3.)*radius_ionp;
    }

    std::cout<<"Dimensions given: "<<mtx_dimensions<<"\n" <<std::endl;

    matrix_xlength = mtx_dimensions;
    matrix_ylength = mtx_dimensions;
    matrix_zlength = mtx_dimensions;
  }

  std::cout << "\nCell Origin: " << Cell_Origin.x() << " " << Cell_Origin.y() << " " << Cell_Origin.z() << std::endl;
  std::cout << "Cell Radius: " << Cell_radius << "\n";
  std::cout << "Cell Number: " << Cell_number << "\n";
  if (cell_permeability_explicit) {
    std::cout << "Cell permeability: " << cell_permeability << "\n\n";
  }
  else {
    std::cout << "Cell permeability: impermeable (default 0)\n\n";
  }
  if ( Cell_concentration != 0. ){
    std::cout << "median Cell Concentration: " << Cell_concentration << " fg/cell" << std::endl;
    std::cout << "Cell Concentration Std: " << Cell_concentration_std << std::endl;
  }

  std::cout << "Center of the voxel at position: " << Voxel_Pos.x() << " " << Voxel_Pos.y() << " " << Voxel_Pos.z() << std::endl;

  if ( Cube_size.x() != 0. ){
    std::cout << "Cube Origin: " << Cube_Origin.x() << " " << Cube_Origin.y() << " " << Cube_Origin.z() << std::endl;
    std::cout << "Cube Size: " << Cube_size.x() << " " << Cube_size.y() << " " << Cube_size.z() << "\n\n";
  }

  if (Cell_hist == true){
    std::cout << "An iron-load per cell histogram is used from the configured INI root Histogram directory" << "\n" << std::endl;
  }
  if (Radius_hist == true){
	std::cout << "A radius histogram for the IONPs is used from the configured INI root Radius_Histogram directory" << "\n" << std::endl;
  }

  /*if (time_step == 0.){

    //Vuong method
    time_step = square(radius_ionp)/(6 * D_D_CONST);
  }*/
  /*if (numThreads == 0){
    numThreads = 10;
  }*/
  if (velocityIONP != 0){
    std::cout << "The IONPs are moving during the Simulation!" << std::endl;
    std::cout << "Velocity: " << velocityIONP << "m/s" << std::endl;
    std::cout << "Direction: " << directionIONP.x() << " " << directionIONP.y() << " " << directionIONP.z() << std::endl;
  }

  std::cout<<"Echo Time: "<< echo_time <<std::endl;
  std::cout<<"Time Step: "<< time_step <<std::endl;
  std::cout<<"Variable Time Step: "<< (variable_time ? "Yes" : "No") <<std::endl;
  std::cout<<"Surrounding voxels: "<< (surrounding_voxels ? "Yes" : "No") <<std::endl;
  std::cout<<"Surrounding voxel superparticles per group: "
           << surrounding_voxel_superparticles
           << (surrounding_voxel_superparticles == 0 ? " (exact replicas)" : "")
           <<std::endl;
  std::cout<<"Echo Spacing (SE/MSE): "<< echo_spacing <<std::endl;

  std::cout<<"Threads: "<< numThreads <<std::endl;
  std::cout<<"B-field solver: "<< bfield_solver_name(bfield_solver) <<std::endl;
}

// Read value/frequency pairs and replace randomHist with 2000 sampled bins.
// The histogram sampler consumes its own seed from the shared run counter.
void Data::readHist ( const std::string& path_to_ionp_data, std::vector<double>& randomHist)
{
  std::cout << "Beginning to read Histogram File" << std::endl;

  std::ifstream file_in(path_to_ionp_data, std::ios::in);

	if (!file_in.is_open()) {
		throw std::runtime_error("Unable to open input file: " + path_to_ionp_data);
	}

  std::vector<double> bins;
  std::vector<int> frequencies;

  double bin;
  double frequency;

  while (file_in >> bin >> frequency) {
    bins.push_back(bin);
    frequencies.push_back(SC_I(frequency));
  }

  RandNumberHist rand(bins, frequencies);

  int lenVector = 2000;

  randomHist = rand.randomVector(lenVector);

  std::cout << "Random Vector according to Histogram: " << randomHist.size() << std::endl;

}

// Append coordinate triples as IONPs, converting configured outer radius to
// shell thickness. Stream extraction terminates at malformed/incomplete input.
void Data::readIONP (
    const std::string& path_to_ionp_data,
    std::vector<Ionp>& ionpInstances,
    SimMatrix& simMatrix,
    double r_ionp,
    double thickness_ionp_coating)
{
  std::cout << "Beginning to read IONP position file" << std::endl;

  (void) ionpInstances;
  (void) simMatrix;

  std::ifstream file_in(path_to_ionp_data, std::ios::in);

	if (!file_in.is_open()) {
		throw std::runtime_error("Unable to open input file: " + path_to_ionp_data);
	}

  double x, y, z;
  const double coating = ionp_shell_thickness(r_ionp, thickness_ionp_coating);
  while (file_in >> x >> y >> z) {
    ionpInstances.push_back(Ionp(x,y,z, 0., r_ionp, coating));
  }

  for (unsigned int i = 0; i < SC_UI(ionpInstances.size()); ++i) {
		for (unsigned int k = 0; k < SC_UI(ionpInstances.size()); ++k) {
      if (i != k){
        if (distance(ionpInstances[i].xt(), ionpInstances[k].xt()) <
            ionpInstances[i].outer_radius() +
              ionpInstances[k].outer_radius()) {
					throw std::runtime_error("IONPs are placed too close to each other: " + std::to_string(i) + " and " + std::to_string(k));
				}
			}
    }
    if (std::abs(ionpInstances[i].xt().x()) +
          ionpInstances[i].outer_radius() > simMatrix.xLength() / 2.) {
      throw std::runtime_error("IONP " + std::to_string(i) + " is outside of Simulation Volume in x direction: " + std::to_string(ionpInstances[i].xt().x()) + " vs. " + std::to_string(simMatrix.xLength() / 2.));
		}
		if (std::abs(ionpInstances[i].xt().y()) +
          ionpInstances[i].outer_radius() > simMatrix.yLength() / 2.) {
			throw std::runtime_error("IONP " + std::to_string(i) + " is outside of Simulation Volume in y direction: " + std::to_string(ionpInstances[i].xt().y()) + " vs. " + std::to_string(simMatrix.yLength() / 2.));
		}
		if (std::abs(ionpInstances[i].xt().z()) +
          ionpInstances[i].outer_radius() > simMatrix.zLength() / 2.) {
			throw std::runtime_error("IONP " + std::to_string(i) + " is outside of Simulation Volume in z direction: " + std::to_string(ionpInstances[i].xt().z()) + " vs. " + std::to_string(simMatrix.zLength() / 2.));
		}

    std::cout << ionpInstances[i].xt().x() << " " << ionpInstances[i].xt().y() << " " << ionpInstances[i].xt().z() << std::endl;
  }

  std::cout << "IONPs created: " << ionpInstances.size() << std::endl;

}

// Append cell centres in metres; the configured common radius is separate.
void Data::readCell ( const std::string& path_to_cell_data, std::vector<Point3D>& cell_origins, double r_cell )
{
  std::cout << "Beginning to read Cell position file" << std::endl;

  (void) cell_origins;

  std::ifstream file_in(path_to_cell_data, std::ios::in);

	if (!file_in.is_open()) {
		throw std::runtime_error("Unable to open input file: " + path_to_cell_data);
	}

  double x, y, z;
  while (file_in >> x >> y >> z) {
    cell_origins.push_back(Point3D(x,y,z));
  }

  for (unsigned int i = 0; i < SC_UI(cell_origins.size()); ++i) {
		for (unsigned int k = 0; k < SC_UI(cell_origins.size()); ++k) {
      if (i != k){
        if (distance(cell_origins[i], cell_origins[k]) < (2. * r_cell)) {
					throw std::runtime_error("Cells are placed too close to each other: " + std::to_string(i) + " and " + std::to_string(k));
				}
			}
    }

  }

  std::cout << "Cells created: " << cell_origins.size() << std::endl;

}

// Append scalar cell loads in the input order (fg per cell in configuration use).
void Data::readCell ( const std::string& path_to_cell_data, std::vector<double>& cell_loads)
{
  std::cout << "Beginning to read Cell iron load file" << std::endl;

  (void) cell_loads;

  std::ifstream file_in(path_to_cell_data, std::ios::in);

	if (!file_in.is_open()) {
		throw std::runtime_error("Unable to open input file: " + path_to_cell_data);
	}

  double x;
  while (file_in >> x) {
    cell_loads.push_back(x);
  }

  std::cout << "Iron load for " << cell_loads.size() << " Cells" << std::endl;

}

// Append aggregate centres; their configured radius controls subsequent placement.
void Data::readAgg ( const std::string& path_to_ionp_data, std::vector<Aggregate>& ionpInstances, SimMatrix& simMatrix, double r_ionp )
{
  std::cout << "Beginning to read Agg position file" << std::endl;

  (void) ionpInstances;
  (void) simMatrix;

  std::ifstream file_in(path_to_ionp_data, std::ios::in);

	if (!file_in.is_open()) {
		throw std::runtime_error("Unable to open input file: " + path_to_ionp_data);
	}

  double x, y, z;
  while (file_in >> x >> y >> z) {
    ionpInstances.push_back(Aggregate(x,y,z, 0., r_ionp, 0.));
  }

  for (unsigned int i = 0; i < SC_UI(ionpInstances.size()); ++i) {
		for (unsigned int k = 0; k < SC_UI(ionpInstances.size()); ++k) {
      if (i != k){
        if (distance(ionpInstances[i].xt(), ionpInstances[k].xt()) < (2. * r_ionp)) {
					throw std::runtime_error("Aggregates are placed too close to each other: " + std::to_string(i) + " and " + std::to_string(k));
				}
			}
    }
    if (std::abs(ionpInstances[i].xt().x()) > simMatrix.xLength() / 2.) {
      throw std::runtime_error("Aggregate " + std::to_string(i) + " is outside of Simulation Volume in x direction: " + std::to_string(ionpInstances[i].xt().x()) + " vs. " + std::to_string(simMatrix.xLength() / 2.));
		}
		if (std::abs(ionpInstances[i].xt().y()) > simMatrix.yLength() / 2.) {
			throw std::runtime_error("Aggregate " + std::to_string(i) + " is outside of Simulation Volume in y direction: " + std::to_string(ionpInstances[i].xt().y()) + " vs. " + std::to_string(simMatrix.yLength() / 2.));
		}
		if (std::abs(ionpInstances[i].xt().z()) > simMatrix.zLength() / 2.) {
			throw std::runtime_error("Aggregate " + std::to_string(i) + " is outside of Simulation Volume in z direction: " + std::to_string(ionpInstances[i].xt().z()) + " vs. " + std::to_string(simMatrix.zLength() / 2.));
		}
  }

}

void Data::writeData ( const std::vector<std::complex<double>>& Signal_vec, const std::vector<double>& timeV ) const
{
	writeSignalData(Signal_vec, timeV, "", false);
}

void Data::writeImportanceRegionData(
	const std::vector<std::complex<double>>& Signal_vec,
	const std::vector<double>& timeV,
	ProtonSampleRegion region) const {
	writeSignalData(
		Signal_vec,
		timeV,
		"_region_" + sample_region_name(region),
		true);
}

void Data::writeImportanceProposalData(
	const std::vector<std::complex<double>>& Signal_vec,
	const std::vector<double>& timeV,
	FreeSamplingProposal proposal,
	bool weighted_contribution) const {
	writeSignalData(
		Signal_vec,
		timeV,
		"_proposal_" + free_sampling_proposal_name(proposal) +
			(weighted_contribution
				? "_mis_contribution"
				: "_mean"),
		true);
}

// Common signal serializer. Main signals use magnitude; requested complex
// diagnostics also preserve real/imaginary parts for contribution reconstruction.
void Data::writeSignalData(
	const std::vector<std::complex<double>>& Signal_vec,
	const std::vector<double>& timeV,
	const std::string& region_suffix,
	bool complex_output) const
{

  if (output_tmp_dir.empty()) {
		throw std::runtime_error("Output directory is not configured for signal output");
	}
	std::filesystem::create_directories(output_tmp_dir);

	const std::size_t file_index = output::next_file_index(output_tmp_dir);
	const std::string prefix = std::to_string(file_index);

	  std::string magnetization_out_name;
	const std::string signal_type =
		complex_output ? "_signal_complex_" : "_signal_";

	if ( sequence_num == 0 ){
		magnetization_out_name = prefix + signal_type + "FID" + region_suffix + "_" + std::to_string(int (radius_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionp_number) + "IONPs_" + std::to_string(effective_proton_number) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te.dat";
	}
	else if ( sequence_num == 1 ){
		magnetization_out_name = prefix + signal_type + "SE" + region_suffix + "_" + std::to_string(int (radius_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionp_number) + "IONPs_" + std::to_string(effective_proton_number) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te.dat";
	}
	else if ( sequence_num == 2 ){
		magnetization_out_name = prefix + signal_type + "MSE" + region_suffix + "_" + std::to_string(int (radius_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionp_number) + "IONPs_" + std::to_string(effective_proton_number) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te_" + std::to_string(int (echo_spacing*pow(10., 3.))) + "ms_spacing.dat";
	  }

	const auto out_path = output_tmp_dir / magnetization_out_name;

  std::cout << out_path.string() << std::endl;

	std::ofstream outfile(out_path, std::ios::out | std::ios::binary);
  outfile.precision(15);
	if (!outfile.is_open()) {
		throw std::runtime_error("Unable to open output file: " + out_path.string());
	}
	for (unsigned int i = 0; i < Signal_vec.size(); ++i) {
		if (complex_output) {
			outfile << Signal_vec[i].real() << "\t"
				<< Signal_vec[i].imag() << "\t"
				<< timeV[i] << std::endl;
		}
		else {
			outfile << Signal_vec[i].real() << "\t"
				<< timeV[i] << std::endl;
		}
	}
	outfile.close();

}

void Data::writeData ( const std::vector<double>& Phase_vec ) const
{

  if (output_tmp_dir.empty()) {
		throw std::runtime_error("Output directory is not configured for phase output");
	}
	std::filesystem::create_directories(output_tmp_dir);

	const std::size_t file_index = output::next_file_index(output_tmp_dir);
	const std::string prefix = std::to_string(file_index);

  std::string magnetization_out_name;

  if ( sequence_num == 0 ){
    magnetization_out_name = prefix + "_Phase_FID_" + std::to_string(int (radius_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionp_number) + "IONPs_" + std::to_string(effective_proton_number) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te.dat";
  }
  else if ( sequence_num == 1 ){
    magnetization_out_name = prefix + "_Phase_SE_" + std::to_string(int (radius_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionp_number) + "IONPs_" + std::to_string(effective_proton_number) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te.dat";
  }
  else if ( sequence_num == 2 ){
    magnetization_out_name = prefix + "_Phase_MSE_" + std::to_string(int (radius_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionp_number) + "IONPs_" + std::to_string(effective_proton_number) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te_" + std::to_string(int (echo_spacing*pow(10., 3.))) + "ms_spacing.dat";
  }

	const auto out_path = output_tmp_dir / magnetization_out_name;

  std::cout << out_path.string() << std::endl;

	std::ofstream outfile(out_path, std::ios::out | std::ios::binary);
  outfile.precision(15);
	if (!outfile.is_open()) {
		throw std::runtime_error("Unable to open output file: " + out_path.string());
	}
	for (unsigned int i = 0; i < Phase_vec.size(); ++i) {
		outfile << Phase_vec[i] << std::endl;
	}
	outfile.close();

}

void Data::writeData ( const std::vector<Ionp>& IONP_vec ) const
{

  if (output_tmp_dir.empty()) {
		throw std::runtime_error("Output directory is not configured for IONP output");
	}
	std::filesystem::create_directories(output_tmp_dir);

	const std::size_t file_index = output::next_file_index(output_tmp_dir);
	const std::string prefix = std::to_string(file_index);

  std::string ionps_out_name;

  if ( sequence_num == 0 ){
    ionps_out_name = prefix + "_IONPs_FID_" + std::to_string(int (radius_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionp_number) + "IONPs_" + std::to_string(effective_proton_number) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te.dat";
  }
  else if ( sequence_num == 1 ){
    ionps_out_name = prefix + "_IONPs_SE_" + std::to_string(int (radius_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionp_number) + "IONPs_" + std::to_string(effective_proton_number) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te.dat";
  }
  else if ( sequence_num == 2 ){
    ionps_out_name = prefix + "_IONPs_MSE_" + std::to_string(int (radius_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionp_number) + "IONPs_" + std::to_string(effective_proton_number) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te_" + std::to_string(int (echo_spacing*pow(10., 3.))) + "ms_spacing.dat";
  }

	const auto out_path = output_tmp_dir / ionps_out_name;

  std::cout << out_path.string() << std::endl;

	std::ofstream outfile(out_path, std::ios::out | std::ios::binary);
  outfile.precision(15);
	if (!outfile.is_open()) {
		throw std::runtime_error("Unable to open output file: " + out_path.string());
	}
  for (unsigned int i = 0; i < IONP_vec.size(); ++i) {
    outfile << IONP_vec[i].xt().x()<<"\t"<<IONP_vec[i].xt().y()<<"\t"<<IONP_vec[i].xt().z()<< std::endl;
  }
	outfile.close();

}

void Data::writeData ( const std::vector<Aggregate>& Agg_vec ) const
{

  if (output_tmp_dir.empty()) {
		throw std::runtime_error("Output directory is not configured for aggregate output");
	}
	std::filesystem::create_directories(output_tmp_dir);

	const std::size_t file_index = output::next_file_index(output_tmp_dir);
	const std::string prefix = std::to_string(file_index);

  std::string agg_out_name;

  if ( sequence_num == 0 ){
    agg_out_name = prefix +  "_Aggs_FID_" + std::to_string(int (radius_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionp_number) + "IONPs_" + std::to_string(effective_proton_number) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te.dat";
  }
  else if ( sequence_num == 1 ){
    agg_out_name = prefix +  "_Aggs_SE_" + std::to_string(int (radius_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionp_number) + "IONPs_" + std::to_string(effective_proton_number) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te.dat";
  }
  else if ( sequence_num == 2 ){
    agg_out_name = prefix +  "_Aggs_MSE_" + std::to_string(int (radius_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionp_number) + "IONPs_" + std::to_string(effective_proton_number) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te_" + std::to_string(int (echo_spacing*pow(10., 3.))) + "ms_spacing.dat";
  }

	const auto out_path = output_tmp_dir / agg_out_name;

  std::cout << out_path.string() << std::endl;

	std::ofstream outfile(out_path, std::ios::out | std::ios::binary);
  outfile.precision(15);
	if (!outfile.is_open()) {
		throw std::runtime_error("Unable to open output file: " + out_path.string());
	}
  for (unsigned int i = 0; i < Agg_vec.size(); ++i) {
    outfile << Agg_vec[i].xt().x()<<"\t"<<Agg_vec[i].xt().y()<<"\t"<<Agg_vec[i].xt().z()<< std::endl;
  }
	outfile.close();

}

void Data::writeData ( const std::vector<Proton>& Proton_vec ) const
{

  if (output_tmp_dir.empty()) {
		throw std::runtime_error("Output directory is not configured for proton output");
	}
	std::filesystem::create_directories(output_tmp_dir);

	const std::size_t file_index = output::next_file_index(output_tmp_dir);
	const std::string prefix = std::to_string(file_index);

  std::string agg_out_name;

  if ( sequence_num == 0 ){
    agg_out_name = prefix +  "_Protons_FID_" + std::to_string(int (radius_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionp_number) + "IONPs_" + std::to_string(effective_proton_number) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te.dat";
  }
  else if ( sequence_num == 1 ){
    agg_out_name = prefix +  "_Protons_SE_" + std::to_string(int (radius_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionp_number) + "IONPs_" + std::to_string(effective_proton_number) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te.dat";
  }
  else if ( sequence_num == 2 ){
    agg_out_name = prefix +  "_Protons_MSE_" + std::to_string(int (radius_ionp*pow(10., 9.))) + "nm_" + std::to_string(ionp_number) + "IONPs_" + std::to_string(effective_proton_number) + "Protons_" + std::to_string(int (time_step*pow(10., 9.))) + "ns_dt_" + std::to_string(int (echo_time*pow(10., 3.))) + "ms_te_" + std::to_string(int (echo_spacing*pow(10., 3.))) + "ms_spacing.dat";
  }

	const auto out_path = output_tmp_dir / agg_out_name;

  std::cout << out_path.string() << std::endl;

	std::ofstream outfile(out_path, std::ios::out | std::ios::binary);
  outfile.precision(15);
	if (!outfile.is_open()) {
		throw std::runtime_error("Unable to open output file: " + out_path.string());
	}
  for (unsigned int i = 0; i < Proton_vec.size(); ++i) {
    outfile << Proton_vec[i].xt().x()<<"\t"<<Proton_vec[i].xt().y()<<"\t"<<Proton_vec[i].xt().z()<< std::endl;
  }
	outfile.close();

}

// Redirect stdout to the numbered log, or restore its original stream buffer.
void Data::log ( bool reset )
{

  static std::streambuf* previous_buf = nullptr;
  static std::ofstream out;

  if (!reset){
		if (output_log_dir.empty()) {
			throw std::runtime_error("Log directory is not configured");
		}
		if (output_tmp_dir.empty()) {
			throw std::runtime_error("Tmp directory is not configured");
		}
		std::filesystem::create_directories(output_log_dir);

		const std::size_t log_index = output::next_file_index(output_log_dir);
		const std::size_t tmp_index = (first_tmp_output_index != 0) ? first_tmp_output_index : output::next_file_index(output_tmp_dir);
		const auto log_out_path = output_log_dir / (std::to_string(log_index) + "_" + std::to_string(tmp_index) + "_log.txt");

		out.open(log_out_path);
		if (!out.is_open()) {
			throw std::runtime_error("Unable to open log file: " + log_out_path.string());
		}

		previous_buf = std::cout.rdbuf(out.rdbuf()); //redirect std::cout to log.txt!
  }
  else if (previous_buf){
    std::cout.rdbuf(previous_buf); //reset to standard output again
		previous_buf = nullptr;
		out.close();
  }

}

// Configuration accessors: these expose stored/resolved values without unit conversion.
double Data::get_matrix_x0() const{
  return matrix_x0;
}
double Data::get_matrix_y0() const{
  return matrix_y0;
}
double Data::get_matrix_z0() const{
  return matrix_z0;
}
double Data::get_matrix_xlength() const{
  return matrix_xlength;
}
double Data::get_matrix_ylength() const{
  return matrix_ylength;
}
double Data::get_matrix_zlength() const{
  return matrix_zlength;
}
int Data::get_matrix_nxnodes() const{
  return matrix_nxnodes;
}
int Data::get_matrix_nynodes() const{
  return matrix_nynodes;
}
int Data::get_matrix_nznodes() const{
  return matrix_nznodes;
}
double Data::get_radius_ionp() const{
  return radius_ionp;
}
double Data::get_radius_agg() const{
  return radius_agg;
}
int Data::get_ionp_number() const{
  return ionp_number;
}
int Data::get_agg_number() const{
  return agg_number;
}
std::optional<std::string> Data::get_ionp_position_file() const{
  return ionp_position_file;
}
std::optional<std::string> Data::get_agg_position_file() const{
  return agg_position_file;
}
std::optional<std::string> Data::get_cell_position_file() const{
  return cell_position_file;
}
std::optional<std::string> Data::get_cell_load_file() const{
  return cell_load_file;
}
int Data::get_proton_number() const{
  return proton_number;
}
double Data::get_echo_time() const{
  return echo_time;
}
double Data::get_time_step() const{
  return time_step;
}
bool Data::get_variable_time() const{
  return variable_time;
}
bool Data::get_surrounding_voxels() const{
  return surrounding_voxels;
}
std::size_t Data::get_surrounding_voxel_superparticles() const{
  return surrounding_voxel_superparticles;
}
double Data::get_echo_spacing() const{
  return echo_spacing;
}
int Data::get_numThreads() const{
  return numThreads;
}
int Data::get_sequence_num() const{
  return sequence_num;
}
double Data::get_velocity_IONP() const{
  return velocityIONP;
}
Point3D Data::get_direction_IONP() const{
  return directionIONP;
}
double Data::get_thickness_IONPcoating() const{
  return thickness_IONPcoating;
}
bool Data::has_coating_diffusion_coefficient() const{
  return coating_diffusion_coefficient_explicit;
}
double Data::get_coating_diffusion_coefficient() const{
  return coating_diffusion_coefficient;
}
bool Data::has_cell_diffusion_coefficient() const{
  return cell_diffusion_coefficient_explicit;
}
double Data::get_cell_diffusion_coefficient() const{
  return cell_diffusion_coefficient;
}
bool Data::has_coating_permeability() const{
  return coating_permeability_explicit;
}
double Data::get_coating_permeability() const{
  return coating_permeability;
}
bool Data::has_cell_permeability() const{
  return cell_permeability_explicit;
}
double Data::get_cell_permeability() const{
  return cell_permeability;
}
int Data::get_bound_proton_number() const{
  return bound_proton_number;
}
double Data::get_background_T2() const{
  return background_T2;
}
double Data::get_cell_T2() const{
  return cell_T2;
}
double Data::get_radius_ionp_1() const{
  return radius_ionp_1;
}
double Data::get_ionp_per() const{
  return ionp_per;
}
double Data::get_radius_ionp_std() const{
  return radius_ionp_std;
}
double Data::get_radius_ionp_std_1() const{
  return radius_ionp_std_1;
}
double Data::get_radius_agg_std() const{
  return radius_agg_std;
}
double Data::get_concentration() const{
  return concentration;
}
double Data::get_concentration_ext() const{
  return concentration_ext;
}
double Data::get_mag_moment() const{
  return mag_moment;
}
Point3D Data::get_Cell_Origin() const{
  return Cell_Origin;
}
Point3D Data::get_Cube_Origin() const{
  return Cube_Origin;
}
Point3D Data::get_Cube_size() const{
  return Cube_size;
}
Point3D Data::get_Voxel_position() const{
  return Voxel_Pos;
}
double Data::get_Cell_radius() const{
  return Cell_radius;
}
double Data::get_Cell_concentration() const{
  return Cell_concentration;
}
double Data::get_Cell_concentration_std() const{
  return Cell_concentration_std;
}
int Data::get_Cell_concentration_dist() const{
  return Cell_concentration_dist;
}
int Data::get_Cell_number() const{
  return Cell_number;
}
bool Data::get_Cell_hist() const{
  return Cell_hist;
}
bool Data::get_Radius_hist() const{
  return Radius_hist;
}
BFieldSolver Data::get_bfield_solver() const{
  return bfield_solver;
}
ImportanceSamplingConfig Data::get_importance_sampling() const{
  return importance_sampling;
}

void Data::set_time_step(double dt){
  time_step = dt;
}

void Data::set_effective_proton_number(int proton_count){
  effective_proton_number = proton_count;
}
