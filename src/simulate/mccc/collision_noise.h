/**
 * @file collision_noise.h
 * @brief Counter-based standard normal numbers for collision steps
 *
 * The number for (seed, marker id, step, component) depends on nothing else,
 * so two codes that use this generator, the same seed and the same fixed
 * step apply identical collision kicks to the same marker. SIMPLE carries a
 * copy of this file (src/collision_noise.c); keep both identical.
 *
 * Mixing is splitmix64 (Steele, Lea, Flood 2014); the normal deviate is
 * Box-Muller from two 53-bit uniforms.
 */
#ifndef COLLISION_NOISE_H
#define COLLISION_NOISE_H

#include <math.h>
#include <stdint.h>

static inline uint64_t collision_noise_mix(uint64_t z) {
    z += 0x9e3779b97f4a7c15ULL;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static inline double collision_noise_normal(uint64_t seed, uint64_t id,
                                            uint64_t step, uint64_t component) {
    uint64_t key = collision_noise_mix(
        seed ^ collision_noise_mix(
            id ^ collision_noise_mix(step ^ collision_noise_mix(component))));
    uint64_t a = collision_noise_mix(key);
    uint64_t b = collision_noise_mix(key ^ 0xd1b54a32d192ed03ULL);
    double u1 = ((double)(a >> 11) + 0.5) * 0x1.0p-53;
    double u2 = (double)(b >> 11) * 0x1.0p-53;
    return sqrt(-2.0 * log(u1)) * cos(6.283185307179586 * u2);
}

#endif
