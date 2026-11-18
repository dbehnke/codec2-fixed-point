#ifndef CODEC2_MOD_UTIL_H
#define CODEC2_MOD_UTIL_H

#include <stdint.h>

/* The float codec's LCG (per decoder state), used for unvoiced harmonic phases. */
static inline int codec2_rand(uint32_t *prng_state)
{
	*prng_state = *prng_state * 1103515245U + 12345U;
	return ((unsigned)(*prng_state >> 16) & 0x7FFF);
}

#endif /* CODEC2_MOD_UTIL_H */
