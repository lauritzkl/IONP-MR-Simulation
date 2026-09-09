/**
 * @file RandNumberBetween.hpp
 *
 * Random-number helpers with optional run-wide seeding. Set the seed before
 * constructing simulation streams. Explicitly seeded constructors bypass the
 * global seed counter; changing construction order changes implicitly seeded
 * streams. Standard-library distributions can differ across toolchains.
 */

/*****************************************************************************
*Class to define a mathematical a random number (header .hpp)
*
* how to: create an object with creator and give the boundaries eg. RandNumberBetween r(0,100)
* draw a random number using () operator e.g r(), output is different each time eg. 0<=r()<100

*2-7-2021
******************************************************************************/
#ifndef RANDNUMBERBETWEEN_HPP_
#define RANDNUMBERBETWEEN_HPP_

#include <algorithm>
#include <iostream>
#include <optional>
#include <random>
#include <vector>
#include <numeric>
#include <chrono> //random generator
#include <cstdlib>
#include <ctime>
#include <atomic>

#include "mptmacros.h"

namespace rng_config {
inline std::optional<std::uint32_t>& base_seed() {
    static std::optional<std::uint32_t> seed;
    return seed;
}

inline std::atomic<std::uint32_t>& seed_counter() {
    static std::atomic<std::uint32_t> counter{0};
    return counter;
}

// Configure before parallel execution; resetting here restarts seed allocation.
inline void set_seed(std::uint32_t seed) {
    base_seed() = seed;
    seed_counter().store(0);
}

inline std::uint32_t next_seed() {
    auto& seed_opt = base_seed();
    if(seed_opt.has_value()) {
        return *seed_opt + seed_counter().fetch_add(1);
    }
    return SC_UI(std::chrono::system_clock::now().time_since_epoch().count());
}
} // namespace rng_config

class RandNumberBase {
public:
    virtual ~RandNumberBase() = default;
    virtual double operator()() = 0;
};

// mean/std parameterize the underlying normal distribution of log(value).
class RandNumberLognorm : public RandNumberBase {
public:
    RandNumberLognorm(double mean, double std)
        : random_engine( rng_config::next_seed() ),
        distribution{ mean, std } {}

    double operator() () override {
        return distribution(random_engine);
    }
private:
    std::mt19937 random_engine;
    std::lognormal_distribution<double> distribution;
};

// mean and standard deviation of the sampled normal distribution.
class RandNumberGaussian : public RandNumberBase {
public:
        RandNumberGaussian(double mean, double std)
            : random_engine( rng_config::next_seed() ),
            distribution{ mean, std } {}
        RandNumberGaussian(double mean, double std, std::uint32_t seed)
            : random_engine(seed),
            distribution{ mean, std } {}

    double operator() () override {
        return distribution(random_engine);
    }
private:
    std::mt19937 random_engine;
    std::normal_distribution<double> distribution;
};

// Uniform real samples use the standard distribution's [low, high) interval.
class RandNumberBetween {
public:
        RandNumberBetween(double low, double high)
            //: random_engine{ std::random_device{} ()}, //same sequence
            : random_engine( rng_config::next_seed() ),
            distribution{ low, high } {}
        RandNumberBetween(double low, double high, std::uint32_t seed)
            : random_engine(seed),
            distribution{ low, high } {}

    double operator() () {
        return distribution(random_engine);
    }

private:
    std::mt19937 random_engine;
    std::uniform_real_distribution<double> distribution;
};

// Discrete weighted bins, copied into the sampler. Callers must supply
// matching bin/frequency lengths, nonnegative counts and positive total count.
// The scalar method seeds per call; randomVector uses one engine for the vector.
class RandNumberHist {
public:

    RandNumberHist(const std::vector<double>& bins_var, const std::vector<int>& frequencies_var)
        : bins(bins_var),
        frequencies(frequencies_var) {}

        // Function to generate a random number based on the histogram
    double generateRandomNumber() {
        // Calculate the total frequency
        int totalFrequency = 0;
        for(int freq : frequencies) {
            totalFrequency += freq;
        }

        // Generate a random number between 1 and the total frequency
        std::mt19937 rng(rng_config::next_seed());
        std::uniform_int_distribution<int> dist(1, totalFrequency);
        int randomNumber = dist(rng);

        // Find the bin corresponding to the random number
        int cumulativeFrequency = 0;
        for(size_t i = 0; i < bins.size(); ++i) {
            cumulativeFrequency += frequencies[i];
            if(randomNumber <= cumulativeFrequency) {
                return bins[i];
            }
        }

        // This line should never be reached
        return -1.0;
    }

    std::vector<double> randomVector(const int lenVec){
        std::mt19937 rng(rng_config::next_seed());
        std::uniform_int_distribution<int> dist(1, std::accumulate(frequencies.begin(), frequencies.end(), 0));
        std::vector<double> Vec;
        for(int i = 0; i<lenVec; i++){
            int randomNumber = dist(rng);
            int cumulativeFrequency = 0;
            for(size_t j = 0; j < bins.size(); ++j) {
                cumulativeFrequency += frequencies[j];
                if(randomNumber <= cumulativeFrequency) {
                    Vec.push_back(bins[j]);
                    break;
                }
            }
        }

        return Vec;
    }

private:
    std::vector<double> bins;
    std::vector<int> frequencies;
};

#endif
