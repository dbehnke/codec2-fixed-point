/*
 * fxp_fft.h — fixed-point 512-point real FFT (forward and inverse) for Codec2 3200.
 *
 * Data: Q15-style int16 with a block exponent. The functions return an exponent e such that
 *   true_value = stored_value * 2^e   (relative to the units of the input values).
 * Both transforms are unnormalised, matching kiss_fftr / kiss_fftri:
 *   forward  X[k] = sum_n x[n] e^{-j 2 pi k n / 512},  k = 0..256
 *   inverse  x[n] = sum_k X[k] e^{+j 2 pi k n / 512},  n = 0..511 (X Hermitian, X[0], X[256] real)
 * Precision: the input is normalised to use the 16-bit range, then every stage applies a conditional
 * shift (block floating point), so small signals keep their resolution. See docs/04-fxp-fft.md.
 */
#ifndef CODEC2_FXP_FFT_H
#define CODEC2_FXP_FFT_H

#include <stdint.h>
#include "fxp.h"

typedef struct { int16_t r, i; } cq15_t;

#define FXP_RFFT_N 512
#define FXP_RFFT_BUF (FXP_RFFT_N / 2 + 1) /* cq15_t elements needed by the in-place transforms */

extern const cq15_t fxp_w512[384];
extern const uint8_t fxp_digrev256[256];
extern const cq15_t fxp_tw256[243]; /* stage twiddles in reading order (copies of fxp_w512) */
extern const int16_t fxp_cos512_q14[512];

typedef struct { int32_t r, i; } cq31_t;

/* Sparse real DFT: X[k] = sum_{n<len} a[n] e^{-j 2 pi k n / 512}, k = 0..256, for short sequences
 * (LPC polynomials, len = 11). Exact int32 sums of int16 x Q14 products, no intermediate rounding:
 * true X = stored * 2^-14 (in the units of a). Precondition (checked by FXP_ASSERT): sum|a[n]| < 2^17,
 * which bounds every partial sum by 2^31 - see fxp_rdft.c. Use fxp_rdft_taps_shift() to scale taps.
 * X needs 257 entries; it may be the same memory as a float complex buffer of 257 entries. */
void fxp_rdft512_sparse(const int16_t *a, int len, cq31_t *X);

/* Left shift s (may be negative = right shift) that brings int32 taps to the most precise int16 form
 * allowed by fxp_rdft512_sparse: max|a<<s| < 2^15 and sum|a<<s| < 2^17. maxabs/l1 are the max and the
 * L1 norm of the taps (l1 as uint32). Returns 0 for all-zero taps. */
int fxp_rdft_taps_shift(uint32_t maxabs, uint32_t l1);

/* In-place 256-point complex forward FFT, natural order in and out. Any int16 input is accepted
 * (including -32768). Returns the number of right shifts applied (output = DFT * 2^-return). */
int fxp_cfft256(cq15_t *x);

/* Forward real FFT, in place.
 * In:  buf[n].r = x[2n], buf[n].i = x[2n+1], n = 0..255 (buf[256] is scratch).
 * Out: buf[k] = X[k], k = 0..256 (buf[0].i = buf[256].i = 0). Returns e (X_true = buf * 2^e). */
int fxp_rfft512(cq15_t buf[FXP_RFFT_BUF]);

/* Inverse real FFT, in place.
 * In:  buf[k] = X[k], k = 0..256 (imaginary parts of X[0] and X[256] are ignored).
 * Out: buf[n].r = x[2n], buf[n].i = x[2n+1], n = 0..255. Returns e (x_true = buf * 2^e). */
int fxp_irfft512(cq15_t buf[FXP_RFFT_BUF]);

#endif /* CODEC2_FXP_FFT_H */
