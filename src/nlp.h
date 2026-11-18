#ifndef CODEC2_MOD_NLP_H
#define CODEC2_MOD_NLP_H

#include "codec2_internal.h"
#include "fxp_fft.h"

#define NLP_SCRATCH_WORDS 260 /* x (32 words) + re, im (2 x 114) */

/* pitch as the NLP DFT bin (pitch period 2560 / bin samples); *prev_bin is the tracking memory.
 * scratch: NLP_SCRATCH_WORDS int32, clobbered. */
int nlp(
    nlp_t *restrict nlp,
    const int16_t *restrict Sn, /* input speech vector */
    int *restrict prev_bin,     /* previous frame's bin */
    int32_t *restrict scratch
);

void nlp_init(nlp_t *nlp);

extern const int16_t fxp_nlp_fir_q15[NLP_NTAP]; /* fxp_nlp_tables.c */
extern const int16_t fxp_nlp_hann_q15[NDEC];
extern const uint32_t fxp_nlp_cs512[512]; /* cos | sin << 16, Q14 */

#endif /* CODEC2_MOD_NLP_H */
