/**
 * @file RandomList.cpp
 *
 * Distribution dispatch and legacy seed allocation for RandomList. Preserve
 * both seed requests in the implicit-seed path to retain subsequent streams.
 */

#include "RandNumberBetween.hpp"
#include "RandomList.hpp"

namespace {

// Both requests are intentional, even though only one result is used.
// Removing either request changes the seed allocation of later streams.
std::uint32_t legacy_random_list_seed(bool gauss) {
    const std::uint32_t uniform_seed = rng_config::next_seed();
    const std::uint32_t gaussian_seed = rng_config::next_seed();
    return gauss ? gaussian_seed : uniform_seed;
}

} // namespace

RandomList::RandomList(double low, double high, bool use_gaussian)
    : RandomList(low, high, use_gaussian, legacy_random_list_seed(use_gaussian)) {
}

RandomList::RandomList(double low,
                       double high,
                       bool use_gaussian,
                       std::uint32_t seed)
    : gauss(use_gaussian),
      random_engine(seed),
      uniform_distribution(low, high),
      gaussian_distribution(low, high) {
}

double RandomList::getNumber(  ){
    return gauss ? gaussian_distribution(random_engine)
                 : uniform_distribution(random_engine);
}
