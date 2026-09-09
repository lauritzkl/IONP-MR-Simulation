/**
 * @file profiler.cpp
 *
 * Implementation of the legacy timing utility. stop adds elapsed time since
 * the latest start; reset clears accumulated duration, not the start point.
 * Reported nanoseconds are storage units, not a guarantee of clock resolution.
 */

#include <iostream>
#include <vector>
#include <string>
#include <chrono> // Wall-clock measurement.

#include "mptmacros.h"

#include "profiler.h"

// Allocate unnamed timers and clear all accumulated durations.
void Profiler::init ( uint numTimers )
{
    duration = std::vector<long long>( numTimers, 0 );
    tStart = std::vector<std::chrono::high_resolution_clock::time_point>( numTimers );
    name = std::vector<std::string>( numTimers, "" );
}

// Allocate named timers; print() pairs these names with accumulated seconds.
void Profiler::init ( const std::vector<std::string>& names )
{
    duration = std::vector<long long>( names.size(), 0 );
    tStart = std::vector<std::chrono::high_resolution_clock::time_point>( names.size() );
    name = names;
}

// Clear one accumulated duration without changing its last start point.
void Profiler::reset ( uint i )
{
    checkSize( i );
    duration[i] = 0.;
}

// Clear accumulated durations for every timer.
void Profiler::resetAll ()
{
    for ( long long & dur: duration )
        dur = 0.;
}

// Record the start point; elapsed time is accumulated only by stop().
void Profiler::start ( uint i )
{
    checkSize( i );
    tStart[i] = std::chrono::high_resolution_clock::now();
}

// Add elapsed time since start(), allowing multiple measured intervals.
void Profiler::stop ( uint i )
{
    checkSize( i );
    auto tEnd = std::chrono::high_resolution_clock::now();
    duration[i] += std::chrono::duration_cast<std::chrono::nanoseconds>( tEnd - tStart[i] ).count();
}

// Number of initialized timer slots.
unsigned long int Profiler::getSize ()
{
    return duration.size();
}

// Accumulated duration expressed as integer nanoseconds.
long long Profiler::getTimeInNS ( uint i )
{
    checkSize( i );
    return duration[i];
}

// Convert the stored nanoseconds to floating-point seconds.
double Profiler::getTimeInS ( uint i )
{
    checkSize( i );
    return SC_D(duration[i]) * 1e-9;
}

// Return the configured timer label.
std::string Profiler::getName ( uint i )
{
    checkSize( i );
    return name[i];
}

// Validate an index before accessing storage; retains the legacy string throw.
void Profiler::checkSize ( uint i )
{
    if ( i >= SC_UI(duration.size()) || SC_I(i) < 0 )
    {
        std::cerr << "ERROR in Profiler: " << i << " is no valid timer index\nNumTimers = " << duration.size() << std::endl;
        throw( "Timer index out of range" );
    }
}

// Write all timer names and accumulated durations to the supplied stream.
void Profiler::print ( std::ostream& out )
{
    for ( uint i = 0; i < SC_UI(Profiler::getSize()); i++ )
        out << "# ----- " << Profiler::getName(i) << ": " << Profiler::getTimeInS(i) << " s -----\n";
    out << std::flush;
}
