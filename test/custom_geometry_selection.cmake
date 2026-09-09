cmake_minimum_required(VERSION 3.16)

if(NOT DEFINED IONP_EXECUTABLE OR IONP_EXECUTABLE STREQUAL "")
  message(FATAL_ERROR "IONP_EXECUTABLE must name the simulator executable")
endif()
if(NOT EXISTS "${IONP_EXECUTABLE}")
  message(FATAL_ERROR "Simulator executable does not exist: ${IONP_EXECUTABLE}")
endif()
if(NOT DEFINED TEST_ROOT OR TEST_ROOT STREQUAL "")
  message(FATAL_ERROR "TEST_ROOT must name a build-local scratch directory")
endif()
get_filename_component(test_root_name "${TEST_ROOT}" NAME)
if(NOT test_root_name MATCHES "^custom_geometry_selection_")
  message(FATAL_ERROR
    "Refusing to clear a scratch directory with an unexpected name: ${TEST_ROOT}")
endif()

file(REMOVE_RECURSE "${TEST_ROOT}")
file(MAKE_DIRECTORY "${TEST_ROOT}")

function(write_test_config config_path ionp_number)
  set(config_text [=[
IONP_radius = 1e-7
IONP_radius_std = 0
Coating = 0
Coating_D = 0
Coating_permeability = 1
Agg_radius = 0
Agg_num = 0
Cell_radius = 0
Cell_num = 0
Cell_hist = No
Conc = 0
Conc_ext = 0
IONP_num = @IONP_NUMBER@
Proton_num = 1
Bound_Prot_n = 0
T_E = 1e-6
T2_out = 0.05
T2_in = 0
Length = 2e-5
time_step = 1e-6
variable_Time = No
Nodes = 0
Threads = 1
Sequence = FID
]=])
  string(REPLACE "@IONP_NUMBER@" "${ionp_number}" config_text "${config_text}")
  file(WRITE "${config_path}" "${config_text}")
  foreach(config_line IN LISTS ARGN)
    file(APPEND "${config_path}" "${config_line}\n")
  endforeach()
endfunction()

function(prepare_case case_name ionp_number)
  set(case_dir "${TEST_ROOT}/${case_name}")
  file(MAKE_DIRECTORY
    "${case_dir}/ini/IONP_pos"
    "${case_dir}/ini/Agg_pos"
    "${case_dir}/output"
    "${case_dir}/logs"
  )
  write_test_config("${case_dir}/ini/config.txt" "${ionp_number}" ${ARGN})
  set(CASE_DIR "${case_dir}" PARENT_SCOPE)
endfunction()

function(run_simulator case_dir output_variable result_variable)
  execute_process(
    COMMAND
      "${IONP_EXECUTABLE}"
      --ini-dir "${case_dir}/ini"
      --ini-root "${case_dir}/ini"
      --output-dir "${case_dir}/output"
      --log-dir "${case_dir}/logs"
      --no-log
      --seed 314159
    RESULT_VARIABLE run_result
    OUTPUT_VARIABLE run_stdout
    ERROR_VARIABLE run_stderr
    TIMEOUT 30
  )
  set(${output_variable} "${run_stdout}\n${run_stderr}" PARENT_SCOPE)
  set(${result_variable} "${run_result}" PARENT_SCOPE)
endfunction()

function(assert_contains variable_name expected context)
  string(FIND "${${variable_name}}" "${expected}" found_at)
  if(found_at EQUAL -1)
    message(FATAL_ERROR
      "${context}: expected output to contain:\n${expected}\n"
      "Complete simulator output:\n${${variable_name}}")
  endif()
endfunction()

function(assert_success result_value output_variable context)
  if(NOT "${result_value}" STREQUAL "0")
    message(FATAL_ERROR
      "${context}: simulator process returned ${result_value}\n"
      "Complete simulator output:\n${${output_variable}}")
  endif()
endfunction()

function(find_geometry_outputs output_dir output_variable)
  file(GLOB geometry_files "${output_dir}/*_IONPs_*.dat")
  list(SORT geometry_files)
  set(${output_variable} "${geometry_files}" PARENT_SCOPE)
endfunction()

function(assert_geometry_count output_dir expected_count context)
  find_geometry_outputs("${output_dir}" geometry_files)
  list(LENGTH geometry_files actual_count)
  if(NOT actual_count EQUAL expected_count)
    message(FATAL_ERROR
      "${context}: expected ${expected_count} IONP geometry output(s), "
      "found ${actual_count}: ${geometry_files}")
  endif()
endfunction()

function(assert_record_count geometry_file expected_count context)
  file(READ "${geometry_file}" geometry_text)
  string(REGEX MATCHALL "[^\r\n]+" geometry_records "${geometry_text}")
  list(LENGTH geometry_records actual_count)
  if(NOT actual_count EQUAL expected_count)
    message(FATAL_ERROR
      "${context}: expected ${expected_count} geometry record(s) in "
      "${geometry_file}, found ${actual_count}:\n${geometry_text}")
  endif()
endfunction()

function(assert_invalid_selector case_name selector_line expected_diagnostic)
  prepare_case("${case_name}" 0 "${selector_line}")
  run_simulator("${CASE_DIR}" invalid_output invalid_result)

  # main intentionally catches errors for each config and continues the batch.
  assert_success("${invalid_result}" invalid_output "${case_name}")
  assert_contains(invalid_output "Simulation failed for" "${case_name}")
  assert_contains(invalid_output "${expected_diagnostic}" "${case_name}")
  assert_geometry_count("${CASE_DIR}/output" 0 "${case_name}")
endfunction()

# Two configs in one --ini-dir must independently select their geometry. The
# first selector also verifies that spaces in a leaf filename are preserved.
set(batch_dir "${TEST_ROOT}/batch")
file(MAKE_DIRECTORY
  "${batch_dir}/ini/IONP_pos"
  "${batch_dir}/output"
  "${batch_dir}/logs"
)
file(WRITE "${batch_dir}/ini/IONP_pos/first positions.txt" "-2e-6 0 0\n")
file(WRITE "${batch_dir}/ini/IONP_pos/second.dat" "2e-6 0 0\n")
write_test_config(
  "${batch_dir}/ini/01_first.txt" 0
  "IONP_position_file = first positions.txt"
)
write_test_config(
  "${batch_dir}/ini/02_second.txt" 0
  "IONP_position_file = second.dat"
)
run_simulator("${batch_dir}" batch_output batch_result)
assert_success("${batch_result}" batch_output "per-config selection batch")
assert_contains(batch_output "Using explicit IONP position file:" "per-config selection batch")
assert_contains(batch_output "first positions.txt" "per-config selection batch")
assert_contains(batch_output "second.dat" "per-config selection batch")
find_geometry_outputs("${batch_dir}/output" batch_geometry_files)
list(LENGTH batch_geometry_files batch_geometry_count)
if(NOT batch_geometry_count EQUAL 2)
  message(FATAL_ERROR
    "per-config selection batch: expected two geometry outputs, found "
    "${batch_geometry_count}: ${batch_geometry_files}\n${batch_output}")
endif()
list(GET batch_geometry_files 0 first_geometry_file)
list(GET batch_geometry_files 1 second_geometry_file)
assert_record_count("${first_geometry_file}" 1 "first selected geometry")
assert_record_count("${second_geometry_file}" 1 "second selected geometry")
file(READ "${first_geometry_file}" first_geometry)
file(READ "${second_geometry_file}" second_geometry)
string(SUBSTRING "${first_geometry}" 0 1 first_x_prefix)
string(SUBSTRING "${second_geometry}" 0 1 second_x_prefix)
if(NOT first_x_prefix STREQUAL "-")
  message(FATAL_ERROR
    "first config did not emit its negative-x geometry:\n${first_geometry}")
endif()
if(NOT second_x_prefix STREQUAL "2")
  message(FATAL_ERROR
    "second config did not emit its positive-x geometry:\n${second_geometry}")
endif()
if(first_geometry STREQUAL second_geometry)
  message(FATAL_ERROR "the two selected geometry outputs unexpectedly match")
endif()

# Surrounding voxels must expand only the magnetic source collection. The
# physical geometry output remains the single loaded central IONP.
prepare_case(
  "surrounding_voxels"
  0
  "IONP_position_file = central.txt"
  "Surrounding_voxels = Yes"
)
file(WRITE "${CASE_DIR}/ini/IONP_pos/central.txt" "-2e-6 0 0\n")
run_simulator("${CASE_DIR}" surrounding_output surrounding_result)
assert_success("${surrounding_result}" surrounding_output "surrounding voxels")
assert_contains(surrounding_output "Surrounding voxels: Yes" "surrounding voxels")
assert_contains(surrounding_output "Physical IONP count: 1" "surrounding voxels")
assert_contains(
  surrounding_output
  "Surrounding voxel source mode: exact"
  "surrounding voxels"
)
assert_contains(
  surrounding_output
  "Surrounding voxel template source count: 1"
  "surrounding voxels"
)
assert_contains(
  surrounding_output
  "Surrounding voxel replica count: 26"
  "surrounding voxels"
)
assert_contains(
  surrounding_output
  "Surrounding magnetic replica count: 26"
  "surrounding voxels"
)
assert_contains(
  surrounding_output
  "Total magnetic field-source count: 27"
  "surrounding voxels"
)
find_geometry_outputs("${CASE_DIR}/output" surrounding_geometry_files)
list(LENGTH surrounding_geometry_files surrounding_geometry_count)
if(NOT surrounding_geometry_count EQUAL 1)
  message(FATAL_ERROR
    "surrounding voxels: expected one physical geometry output, found "
    "${surrounding_geometry_count}: ${surrounding_geometry_files}\n"
    "${surrounding_output}")
endif()
list(GET surrounding_geometry_files 0 surrounding_geometry_file)
assert_record_count(
  "${surrounding_geometry_file}"
  1
  "surrounding voxels physical output"
)

# A positive per-group cap replaces only each neighboring voxel's magnetic
# template. The two central particles remain physical and appear in output.
prepare_case(
  "surrounding_voxel_superparticles"
  0
  "IONP_position_file = central.txt"
  "Surrounding_voxels = Yes"
  "Surrounding_voxel_superparticles = 1"
)
file(WRITE
  "${CASE_DIR}/ini/IONP_pos/central.txt"
  "-2e-6 0 0\n2e-6 0 0\n"
)
run_simulator(
  "${CASE_DIR}"
  superparticle_output
  superparticle_result
)
assert_success(
  "${superparticle_result}"
  superparticle_output
  "surrounding voxel superparticles"
)
assert_contains(
  superparticle_output
  "Physical IONP count: 2"
  "surrounding voxel superparticles"
)
assert_contains(
  superparticle_output
  "Surrounding voxel source mode: superparticles"
  "surrounding voxel superparticles"
)
assert_contains(
  superparticle_output
  "Surrounding voxel superparticle cap: 1"
  "surrounding voxel superparticles"
)
assert_contains(
  superparticle_output
  "Surrounding voxel template source count: 1"
  "surrounding voxel superparticles"
)
assert_contains(
  superparticle_output
  "Physical IONP sum(r^3):"
  "surrounding voxel superparticles"
)
assert_contains(
  superparticle_output
  "Surrounding template sum(r^3):"
  "surrounding voxel superparticles"
)
assert_contains(
  superparticle_output
  "Surrounding template iron difference sum(r^3):"
  "surrounding voxel superparticles"
)
assert_contains(
  superparticle_output
  "Surrounding voxel replica count: 26"
  "surrounding voxel superparticles"
)
assert_contains(
  superparticle_output
  "Surrounding magnetic replica count: 26"
  "surrounding voxel superparticles"
)
assert_contains(
  superparticle_output
  "Total magnetic field-source count: 28"
  "surrounding voxel superparticles"
)
find_geometry_outputs(
  "${CASE_DIR}/output"
  superparticle_geometry_files
)
list(LENGTH superparticle_geometry_files superparticle_geometry_count)
if(NOT superparticle_geometry_count EQUAL 1)
  message(FATAL_ERROR
    "surrounding voxel superparticles: expected one physical geometry output, "
    "found ${superparticle_geometry_count}: ${superparticle_geometry_files}\n"
    "${superparticle_output}")
endif()
list(GET superparticle_geometry_files 0 superparticle_geometry_file)
assert_record_count(
  "${superparticle_geometry_file}"
  2
  "surrounding voxel superparticles physical output"
)

# Omitting the selector retains the legacy scan-and-concatenate behavior.
prepare_case("legacy_scan" 0)
file(WRITE "${CASE_DIR}/ini/IONP_pos/left.txt" "-2e-6 0 0\n")
file(WRITE "${CASE_DIR}/ini/IONP_pos/right.dat" "2e-6 0 0\n")
run_simulator("${CASE_DIR}" legacy_output legacy_result)
assert_success("${legacy_result}" legacy_output "legacy selector omission")
assert_contains(
  legacy_output
  "Using legacy IONP position directory scan:"
  "legacy selector omission"
)
assert_contains(legacy_output "left.txt" "legacy selector omission")
assert_contains(legacy_output "right.dat" "legacy selector omission")
find_geometry_outputs("${CASE_DIR}/output" legacy_geometry_files)
list(LENGTH legacy_geometry_files legacy_geometry_count)
if(NOT legacy_geometry_count EQUAL 1)
  message(FATAL_ERROR
    "legacy selector omission: expected one geometry output, found "
    "${legacy_geometry_count}: ${legacy_geometry_files}\n${legacy_output}")
endif()
list(GET legacy_geometry_files 0 legacy_geometry_file)
assert_record_count("${legacy_geometry_file}" 2 "legacy selector omission")

# Aggregate selection uses the same leaf-file contract and must not consume a
# second supported file sitting beside the selected one.
prepare_case(
  "aggregate_selector" 1
  "Agg_radius = 1e-6"
  "Agg_position_file = chosen aggregate.dat"
)
file(WRITE "${CASE_DIR}/ini/Agg_pos/chosen aggregate.dat" "2e-6 0 0\n")
file(WRITE "${CASE_DIR}/ini/Agg_pos/unselected.txt" "-2e-6 0 0\n")
run_simulator("${CASE_DIR}" aggregate_output aggregate_result)
assert_success("${aggregate_result}" aggregate_output "aggregate selector")
assert_contains(
  aggregate_output
  "Using explicit aggregate position file:"
  "aggregate selector"
)
assert_contains(aggregate_output "chosen aggregate.dat" "aggregate selector")
string(FIND "${aggregate_output}" "unselected.txt" unselected_log_position)
if(NOT unselected_log_position EQUAL -1)
  message(FATAL_ERROR
    "aggregate selector unexpectedly consumed unselected.txt:\n${aggregate_output}")
endif()
file(GLOB aggregate_geometry_files "${CASE_DIR}/output/*_Aggs_*.dat")
list(LENGTH aggregate_geometry_files aggregate_geometry_count)
if(NOT aggregate_geometry_count EQUAL 1)
  message(FATAL_ERROR
    "aggregate selector: expected one aggregate geometry output, found "
    "${aggregate_geometry_count}: ${aggregate_geometry_files}\n${aggregate_output}")
endif()
list(GET aggregate_geometry_files 0 aggregate_geometry_file)
assert_record_count("${aggregate_geometry_file}" 1 "aggregate selector")
file(READ "${aggregate_geometry_file}" aggregate_geometry)
string(SUBSTRING "${aggregate_geometry}" 0 1 aggregate_x_prefix)
if(NOT aggregate_x_prefix STREQUAL "2")
  message(FATAL_ERROR
    "aggregate selector did not emit the chosen positive-x geometry:\n"
    "${aggregate_geometry}")
endif()

# An inactive selector is warned about but never validated. The deliberately
# invalid path must therefore allow the normal generated-position branch to run.
prepare_case(
  "inactive_selector" 1
  "IONP_position_file = ../not-a-valid-selection.csv"
)
run_simulator("${CASE_DIR}" inactive_output inactive_result)
assert_success("${inactive_result}" inactive_output "inactive selector")
assert_contains(
  inactive_output
  "Warning: ignoring IONP_position_file="
  "inactive selector"
)
assert_geometry_count("${CASE_DIR}/output" 1 "inactive selector")

# Active selectors validate the target before simulation. Although the process
# exits zero after catching each per-config error, none may emit geometry.
assert_invalid_selector(
  "missing_target"
  "IONP_position_file = missing.txt"
  "does not resolve to a regular file:"
)

prepare_case(
  "directory_target" 0
  "IONP_position_file = geometry.txt"
)
file(MAKE_DIRECTORY "${CASE_DIR}/ini/IONP_pos/geometry.txt")
run_simulator("${CASE_DIR}" directory_output directory_result)
assert_success("${directory_result}" directory_output "directory target")
assert_contains(directory_output "Simulation failed for" "directory target")
assert_contains(
  directory_output
  "does not resolve to a regular file:"
  "directory target"
)
assert_geometry_count("${CASE_DIR}/output" 0 "directory target")

prepare_case(
  "unsupported_extension" 0
  "IONP_position_file = geometry.csv"
)
file(WRITE "${CASE_DIR}/ini/IONP_pos/geometry.csv" "0 0 0\n")
run_simulator("${CASE_DIR}" extension_output extension_result)
assert_success("${extension_result}" extension_output "unsupported extension")
assert_contains(extension_output "Simulation failed for" "unsupported extension")
assert_contains(
  extension_output
  "must name a .txt or .dat file:"
  "unsupported extension"
)
assert_geometry_count("${CASE_DIR}/output" 0 "unsupported extension")

assert_invalid_selector(
  "path_component"
  "IONP_position_file = nested/geometry.txt"
  "must be a leaf filename without directory components:"
)

prepare_case("absolute_path" 0)
file(WRITE "${CASE_DIR}/outside.txt" "0 0 0\n")
file(APPEND
  "${CASE_DIR}/ini/config.txt"
  "IONP_position_file = ${CASE_DIR}/outside.txt\n"
)
run_simulator("${CASE_DIR}" absolute_output absolute_result)
assert_success("${absolute_result}" absolute_output "absolute path")
assert_contains(absolute_output "Simulation failed for" "absolute path")
assert_contains(
  absolute_output
  "must be a leaf filename without directory components:"
  "absolute path"
)
assert_geometry_count("${CASE_DIR}/output" 0 "absolute path")

assert_invalid_selector(
  "empty_selector"
  "IONP_position_file ="
  "Missing value for key 'IONP_position_file'"
)

prepare_case(
  "empty_geometry_file" 0
  "IONP_position_file = empty.dat"
)
file(WRITE "${CASE_DIR}/ini/IONP_pos/empty.dat" "")
run_simulator("${CASE_DIR}" empty_file_output empty_file_result)
assert_success("${empty_file_result}" empty_file_output "empty geometry file")
assert_contains(empty_file_output "Simulation failed for" "empty geometry file")
assert_contains(
  empty_file_output
  "Selected IONP position file produced no records:"
  "empty geometry file"
)
assert_geometry_count("${CASE_DIR}/output" 0 "empty geometry file")
