/**
 * @file profiler.h
 *
 * Legacy cumulative wall-clock timers, not part of the default simulation
 * build. Call init before start/stop. This header relies on standard types and
 * uint supplied by its including context; its static storage is per translation
 * unit. Timer operations are not synchronized for concurrent use.
 */

#ifndef PROFILER_H
#define PROFILER_H

#include <chrono>

namespace Profiler {
    void init( uint numTimers );
    void init( const std::vector<std::string>& names );
    void reset( uint i );
    void resetAll();
    void start( uint i );
    void stop( uint i );
    unsigned long int getSize();
    double getTimeInS( uint i );
    long long getTimeInNS( uint i );
    std::string getName( uint i );

    void checkSize( uint i );
    void print( std::ostream& out = std::cout );

    typedef std::chrono::high_resolution_clock::time_point timePoint;
    static std::vector<long long>       duration;
    static std::vector<timePoint>       tStart;
    static std::vector<std::string>     name;
}

#endif // PROFILER_H
