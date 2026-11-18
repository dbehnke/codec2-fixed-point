#ifndef CODEC2_MOD_LPC_H
#define CODEC2_MOD_LPC_H

#include "codec2_internal.h"

/* log2 E returned for the zero-energy trap (the float code's E = 0) */
#define LOG2E_ZERO (INT32_MIN / 2)

/* Decoder spectral amplitudes (docs/05, docs/07): log2 A[m] in Q16 and arg A at the harmonics
 * (uint32 angle) from the LSPs in half-Hz, log2 E (Q16) and Wo as a uint32 angle W, with the LPC
 * correction applied. Uses c2->scratch. */
void aks_to_mag2_lsp(codec2_decoder_t *c2, const int16_t *hz2, int32_t log2E, uint32_t W, int L,
                     int32_t *lA, uint32_t *phA);

/* Encoder LPC/LSP analysis of c2->Sn (docs/11): the 10 LSPs in Hz (Q16) and the return value log2 E (Q16;
 * LOG2E_ZERO for zero energy). Uses c2->scratch. */
int32_t speech_to_lsps_q16(codec2_encoder_t *c2, int32_t *restrict lsp_hz);

extern const int32_t fxp_bw_gamma_q30[LPC_ORD + 1];
extern const int32_t fxp_acos_poly_q29[6];
extern const int32_t fxp_lpc_log2_wscale_q16;
extern const uint32_t fxp_lpc_floor_r0;
extern const int32_t fxp_lpc_hz_q29;
extern const int32_t fxp_lpc_eidx_q29;

#endif /* CODEC2_MOD_LPC_H */
