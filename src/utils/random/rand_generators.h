//
// Created by genshen on 2019-02-10.
//

#ifndef MISA_KMC_RAND_GENERATOR_H
#define MISA_KMC_RAND_GENERATOR_H

#include <cstdint>
#include <random>

#include <xoshiro_gen.h>
#include "legacy_random.hpp"
#include "md_building_config.h"

namespace md_rand {
#ifdef RAND_LCG
    typedef std::minstd_rand type_rng;
#endif

#ifdef RAND_MT
    typedef std::mt19937 type_rng;
#endif

#ifdef RAND_STC
    typedef std::ranlux24 type_rng;
#endif

#ifdef RAND_XOSHIRO
    // return type is uint32_t for 128;
    // return type is uint64_t for 256;
    typedef util::random::xoroshiro128_plus type_rng;
#endif

#ifdef RAND_LEGACY
    typedef LegacyRand type_rng;
#endif

#ifdef RAND_LINUX_REAL
#endif

    extern type_rng rng;

    /**
     * set seed for rng
     */
    void seed(const uint32_t seed);

    /**
     * \brief returns a uint32_t random type, between \param low to \param high.
     */
    uint32_t rand32();
}

#endif //MISA_KMC_RAND_GENERATOR_H
