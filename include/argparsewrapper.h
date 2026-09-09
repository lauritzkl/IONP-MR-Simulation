/**
 * @file argparsewrapper.h
 *
 * Legacy wrapper for an external argument parser. The simulator currently
 * parses its CLI in main.cpp; this wrapper is not a default build dependency.
 */

#pragma warning(push)
#pragma warning(disable: 4100)
#include <argparse/argparse.hpp>
#pragma warning(pop)
