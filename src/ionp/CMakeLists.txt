set( PROJNAME IONP )

set( IONP_SOURCES
	  Aggregate.cpp
	  Ionp.cpp
	  Particle.cpp
	  Point3D.cpp
	  Proton.cpp
	  SimMatrix.cpp
	  Vector3D.cpp
	  Data.cpp
	  CorrectedHybridField.cpp
	  ImportanceSampling.cpp
    RandomList.cpp
    MRsequence.cpp
	  SimSpace.cpp
)

find_package(OpenMP REQUIRED COMPONENTS CXX)

add_library( ionp_core ${IONP_SOURCES} )

target_link_libraries(
  ionp_core
  PUBLIC project_options
         project_warnings
         OpenMP::OpenMP_CXX
)

target_include_directories(
  ionp_core
  PUBLIC "${CMAKE_BINARY_DIR}/configured_files/include"
          "${CMAKE_SOURCE_DIR}/include"
          "${CMAKE_SOURCE_DIR}/src/ionp"
)

add_executable( ${PROJNAME} main.cpp )

target_link_libraries(
  ${PROJNAME}
  PRIVATE ionp_core
)
