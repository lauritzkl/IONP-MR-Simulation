set( CMAKE_C_STANDARD "11" CACHE STRING "C standard to use" )

# Not ideal to use this global variable, but necessary to make sure
# that tooling and projects use the same version
set(CMAKE_CXX_STANDARD "17" CACHE STRING "C++ standard to use")

# strongly encouraged to enable this globally to avoid conflicts between
# -Wpedantic being enabled and -std=c++20 and -std=gnu++20 for example
# when compiling with PCH enabled
option(CMAKE_CXX_EXTENSIONS "Specifying whether compiler specific extensions are requested" OFF)

if(MSVC)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /D_CRT_SECURE_NO_WARNINGS")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /D_CRT_SECURE_NO_WARNINGS")
endif()
