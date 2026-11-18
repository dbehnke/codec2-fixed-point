#ifndef CODEC2_MOD_SYNTHESIS_H
#define CODEC2_MOD_SYNTHESIS_H

#include <stdint.h>
#include "codec2_internal.h"

/* Fixed-point synthesis of one 10 ms subframe (docs/06-fxp-synthesis.md).
 *   W: Wo as a uint32 angle (2^32 = 2 pi), L harmonics, voiced flag.
 *   lA[1..L]: log2 A[m] in Q16 (INT32_MIN/2 for A = 0), after the LPC correction.
 *   phA[1..L]: arg A(e^{jw}) at harmonic m's bin, uint32 angle (from fxp_aks_to_mag2). */
void synthesise_one_frame(codec2_decoder_t *c2, int16_t *speech, uint32_t W, int L, int voiced,
			  const int32_t *lA, const uint32_t *phA);

void synthesis_init(codec2_decoder_t *c2);

#ifdef FXP_SYNTH_DEBUG
/* Defined by a test harness: the parameters that drive the synthesis numerics (final phases th[m]),
 * plus the postfilter state: frame energy e and bg_est after the update (dB Q16). */
void fxp_synth_debug_hook(int L, uint32_t W, int voiced, int32_t e, int32_t bg, const int32_t *lA, const uint32_t *th);
#endif

#endif /* CODEC2_MOD_SYNTHESIS_H */
