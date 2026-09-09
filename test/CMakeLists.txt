# Keep assertions active in test executables in every build configuration.
# These directory options do not affect ionp_core or the simulator executable.
if(MSVC)
  add_compile_options(/UNDEBUG)
else()
  add_compile_options(-UNDEBUG)
endif()

add_executable(basic_parsing basic_parsing.cpp)

target_link_libraries(
  basic_parsing
  PRIVATE ionp_core
)

add_test(NAME basic_parsing COMMAND basic_parsing)
set_tests_properties(basic_parsing PROPERTIES ENVIRONMENT "LSAN_OPTIONS=detect_leaks=0")

add_executable(cell_load_radius cell_load_radius.cpp)
target_link_libraries(
  cell_load_radius
  PRIVATE ionp_core
)
add_test(NAME cell_load_radius COMMAND cell_load_radius)
set_tests_properties(cell_load_radius PROPERTIES ENVIRONMENT "LSAN_OPTIONS=detect_leaks=0")

add_executable(elastic_reflection elastic_reflection.cpp)
target_link_libraries(
  elastic_reflection
  PRIVATE ionp_core
)
add_test(NAME elastic_reflection COMMAND elastic_reflection)
set_tests_properties(elastic_reflection PROPERTIES ENVIRONMENT "LSAN_OPTIONS=detect_leaks=0")

add_executable(importance_sampling importance_sampling.cpp)
target_link_libraries(
  importance_sampling
  PRIVATE ionp_core
)
add_test(NAME importance_sampling COMMAND importance_sampling)
set_tests_properties(importance_sampling PROPERTIES ENVIRONMENT "LSAN_OPTIONS=detect_leaks=0")

add_executable(bfield_solver bfield_solver.cpp)
target_link_libraries(
  bfield_solver
  PRIVATE ionp_core
)
add_test(NAME bfield_solver COMMAND bfield_solver)
set_tests_properties(bfield_solver PROPERTIES ENVIRONMENT "LSAN_OPTIONS=detect_leaks=0")

add_executable(bfield_solver_benchmark bfield_solver_benchmark.cpp)
target_link_libraries(
  bfield_solver_benchmark
  PRIVATE ionp_core
)
add_test(NAME bfield_solver_benchmark COMMAND bfield_solver_benchmark)
set_tests_properties(
  bfield_solver_benchmark
  PROPERTIES
    ENVIRONMENT "LSAN_OPTIONS=detect_leaks=0"
    LABELS "benchmark"
)

add_test(
  NAME custom_geometry_selection
  COMMAND
    "${CMAKE_COMMAND}"
    "-DIONP_EXECUTABLE=$<TARGET_FILE:IONP>"
    "-DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}/custom_geometry_selection_tmp"
    -P "${CMAKE_CURRENT_SOURCE_DIR}/custom_geometry_selection.cmake"
)
set_tests_properties(
  custom_geometry_selection
  PROPERTIES
    ENVIRONMENT "LSAN_OPTIONS=detect_leaks=0"
    TIMEOUT 60
)
