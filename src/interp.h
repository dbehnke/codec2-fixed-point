#ifndef CODEC2_MOD_INTERP_H
#define CODEC2_MOD_INTERP_H

#include <stdint.h>
#include "codec2_internal.h"

/* Returns the interpolated subframe's wq (see WO_Q_MIN); may clear *voiced (isolated voiced subframe). */
int interp_Wo_q(int *voiced, int prev_voiced, int prev_wq, int next_voiced, int next_wq);
int32_t interp_energy_log2_q16(int32_t prev_log2e, int32_t next_log2e);
void interpolate_lsp_hz2(int16_t *interp, const int16_t *prev, const int16_t *next);

#endif /* CODEC2_MOD_INTERP_H */
