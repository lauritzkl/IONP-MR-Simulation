/**
 * @file csvparserwrapper.h
 *
 * Legacy include wrapper for the vendored CSV parser. Runtime configuration
 * and geometry parsing are implemented in Data.cpp.
 */

#pragma warning( push )
#pragma warning( disable : 4100 )
#pragma warning( disable : 4458 )
#pragma warning( disable : 4244 )
#pragma warning( disable : 4127 )
#include "csv.hpp"
#pragma warning( pop )
