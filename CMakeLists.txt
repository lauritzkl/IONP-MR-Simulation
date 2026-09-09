cmake_minimum_required(VERSION 3.16...3.23)

# Set the project name and language
project(
  Simulation
  VERSION 0.0.1
  DESCRIPTION "Monte Carlo simulation of proton diffusion and MR signal around iron oxide nanoparticles"
  LANGUAGES CXX C)

message(STATUS "Started CMake for ${PROJECT_NAME} v${PROJECT_VERSION}...\n")

string( TOUPPER ${PROJECT_NAME} PROJECT_NAME_UPPERCASE )

include( cmake/CxxStandardSettings.cmake )
include( cmake/StandardSettings.cmake    )
include( cmake/UnixSpecials.cmake        )

# This variable is set by project() in CMake 3.21+
string(
  COMPARE EQUAL
          "${CMAKE_SOURCE_DIR}"
          "${PROJECT_SOURCE_DIR}"
          PROJECT_IS_TOP_LEVEL)
if(PROJECT_IS_TOP_LEVEL)
  # Consider the CTest module, which creates targets and options!
  # Only needed if you want to enable submissions to a CDash server.
  include(CTest)
endif()

set(GIT_SHA
    "Unknown"
    CACHE STRING "SHA this build was generated from")
string(
  SUBSTRING "${GIT_SHA}"
            0
            8
            GIT_SHORT_SHA)



get_property(BUILDING_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(BUILDING_MULTI_CONFIG)
  if(NOT CMAKE_BUILD_TYPE)
    # Make sure that all supported configuration types have their
    # associated conan packages available. You can reduce this
    # list to only the configuration types you use, but only if one
    # is not forced-set on the command line for VS
    message(TRACE "Setting up multi-config build types")
    set(CMAKE_CONFIGURATION_TYPES
        Debug
        Release
        RelWithDebInfo
        MinSizeRel
        CACHE STRING "Enabled build types" FORCE)
  else()
    message(TRACE "User chose a specific build type, so we are using that")
    set(CMAKE_CONFIGURATION_TYPES
        ${CMAKE_BUILD_TYPE}
        CACHE STRING "Enabled build types" FORCE)
  endif()
endif()

add_library(project_options INTERFACE)
target_compile_features(project_options INTERFACE cxx_std_${CMAKE_CXX_STANDARD})

add_library(project_warnings INTERFACE)
if(ENABLE_DEVELOPER_MODE)
  if(MSVC)
    target_compile_options(project_warnings INTERFACE /W4)
  else()
    target_compile_options(project_warnings INTERFACE -Wall -Wextra -Wpedantic)
  endif()
endif()

add_library(${PROJECT_NAME_UPPERCASE}::project_options ALIAS project_options)
add_library(${PROJECT_NAME_UPPERCASE}::project_warnings ALIAS project_warnings)

# configure files based on CMake configuration options
add_subdirectory(configured_files)

# Adding the src:
add_subdirectory(src)
# Adding the externals:
# add_subdirectory(external)
# Adding the util directory:
# add_subdirectory( utils )

# Adding the tests:
if(ENABLE_TESTING)
  enable_testing()
  add_subdirectory(test)
endif()

# set the startup project for the "play" button in MSVC
#set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT intro)
