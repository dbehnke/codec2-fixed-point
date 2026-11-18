#ifndef CODEC2_MOD_ANALYSIS_H
#define CODEC2_MOD_ANALYSIS_H

#include "codec2_internal.h"

/* Result of one 10 ms analysis frame (docs/09-fxp-analysis.md). The refined pitch period is the exact
 * rational p = pn / pd samples (Wo = 2 pi / p); L = number of harmonics. */
typedef struct
{
	uint32_t pn, pd;
	int L;
	int voiced;
} enc_frame_t;

void analyse_one_frame(codec2_encoder_t *c2, enc_frame_t *f, const int16_t *speech);

/* Wo index (WO_BITS) of Wo = 2 pi pd / pn, as the float encode_Wo but exact. */
int encode_Wo_p(uint32_t pn, uint32_t pd);

/* generated tables (tools/gen/gen-analysis-tables.py) */
extern const uint16_t fxp_hann_q15[NW];
extern const int16_t fxp_winspec_q15[257];
extern const uint32_t fxp_ana_floor_q24;

#endif /* CODEC2_MOD_ANALYSIS_H */
