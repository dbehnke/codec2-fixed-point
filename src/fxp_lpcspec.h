/*
 * fxp_lpcspec.h — fixed-point LPC spectral amplitudes for the decoder (replaces float aks_to_mag2 +
 * lsp_to_lpc). Works directly from the decoded LSPs, which are integers in half-Hz units
 * (decode_lspds_hz2() sums integer Hz codebook deltas; interpolate_lsp_hz2() averages two frames).
 */
#ifndef CODEC2_FXP_LPCSPEC_H
#define CODEC2_FXP_LPCSPEC_H

#include <stdint.h>
#include "fxp.h"
#include "fxp_fft.h"
#include "fxp_math.h"

#define FXP_LPC_ORD 10
#define FXP_SPEC_BINS 256 /* bins 0..255 of the 512-pt grid, w_k = 2 pi k / 512 */

/* cos(pi * hz2 / 8000) in Q30 for an LSP of hz2 half-Hz (any hz2 >= 0). Error < 1e-8. */
int32_t fxp_lsp_cos_q30(int32_t hz2);

/* LPC polynomial from LSP cosines (Q30): ak[0..10] in Q19 (ak[0] = 1 << 19). Same P/Q convention as
 * float lsp_to_lpc(): even-indexed LSPs build P(z) (times 1 + z^-1), odd-indexed build Q(z) (times
 * 1 - z^-1), A = (P + Q) / 2. Computed in Q22 (binomial bound, see .c), output rounded to Q19. */
void fxp_lsp_to_lpc_q19(const int32_t *c_q30, int32_t *ak_q19);

/* Per-bin evaluation of A(e^jw_k) from LSP cosines in product form:
 *   A = e^{-j 11 w / 2} (Pr + j Qr) / 2,  Pr = 2 cos(w/2) prod_even 2 (cos w - c_i),
 *                                         Qr = 2 sin(w/2) prod_odd  2 (cos w - c_i).
 * Returns log2|A|^2 (+ 1e-6 floor, as the float code) in Q16, and if ang != NULL, arg A as a uint32
 * angle (2^32 = 2 pi) in the same orientation as the float DFT of ak. */
int32_t fxp_lsp_spec_bin(const int32_t *c_q30, int k, uint32_t *ang);

/* Round x / 2^23 (a frequency in 512-point bins, from a uint32 angle) half up, as the float code's
 * (int)(f + 0.5) does with exact arithmetic. Codec2 3200's Wo lies on a grid of 1/80 bin per harmonic
 * (W0_MIN + k (W0_MAX - W0_MIN) / 128 = (128 + 7k) / 40 bins, or the mean of two), so m Wo and (m -+ 0.5) Wo
 * often fall exactly on a half bin. W carries a few units of error from float Wo, which would
 * resolve such ties at random. FXP_BIN_TIE (2^12 units = 0.0005 bin) resolves them upward. It cannot move a
 * non-tie, since those are at least 1/160 bin = 52429 units from a half bin; |m * W error| < 2^12 needs
 * |W error| < 51 units for m <= 80. */
#define FXP_BIN_TIE (1u << 12)
static inline uint32_t fxp_bin_round(uint32_t x)
{
	return (x + (1u << 22) + FXP_BIN_TIE) >> 23;
}

/* Decoder spectral amplitudes (fixed-point aks_to_mag2).
 *   lsp_hz2[10]: LSPs in half-Hz.  log2E_q16: log2 of the frame energy E.
 *   W: fundamental Wo as a uint32 angle (2^32 = 2 pi).  L: number of harmonics (L W <= 2^31).
 *   Outputs: log2 A[m] in Q16 (A_log2_q16[1..L], INT32_MIN/2 for an empty band);
 *            if phA != NULL, arg A at harmonic m's bin round(m Wo 512 / 2 pi) (clamped to 255), as
 *            uint32 angles phA[1..L] (synthesis uses -arg A).
 *   scratch: 514 int32 words (may alias the decoder fft_buffer). */
void fxp_aks_to_mag2(const int16_t *lsp_hz2, int32_t log2E_q16, uint32_t W, int L, int32_t *A_log2_q16,
		     uint32_t *phA, int32_t *scratch);

#endif
