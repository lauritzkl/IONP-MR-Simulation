/**
 * @file RandomList.hpp
 *
 * On-demand random stream (the legacy name does not imply a stored list).
 * For uniform sampling var0/var1 are lower/upper bounds; for Gaussian sampling
 * they are mean/standard deviation. Each object owns its engine and distributions.
 */

/*****************************************************************************
*Class to create lists with random numbers

*07-10-2022 (Lauritz Klünder)
******************************************************************************/

#ifndef RANDOMLIST_HPP_
#define RANDOMLIST_HPP_

#include <cstdint>
#include <random>

#include "RandNumberBetween.hpp"

class RandomList
{
public:
    RandomList( double var0, double var1, bool var2 );
    RandomList( double var0, double var1, bool var2, std::uint32_t seed );

    double getNumber(  );

private:
    bool gauss;
    std::mt19937 random_engine;
    std::uniform_real_distribution<double> uniform_distribution;
    std::normal_distribution<double> gaussian_distribution;

};

#endif
