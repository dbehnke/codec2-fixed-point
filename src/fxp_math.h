/*
 * fxp_math.h — fixed-point log2 / exp2 / cosine helpers (conventions in fxp.h).
 * Log-domain values are int32 Q16 (log2 units): range +-32768, resolution 1.5e-5 (9e-5 dB).
 */
#ifndef CODEC2_FXP_MATH_H
#define CODEC2_FXP_MATH_H

#include <stdint.h>
#include "fxp.h"

extern const int32_t fxp_cosq_q30[513];
extern const int32_t fxp_log2_tab_q16[129];
extern const uint32_t fxp_exp2_tab_q30[129];
extern const uint32_t fxp_atan_tab_u32[257];

#define FXP_Q16_ONE 65536

/* log2(v) in Q16 for v > 0 (integer part included; v = 1 -> 0). Max error ~1.2e-5 (table + linear
 * interpolation, h = 1/128). v = 0 returns INT32_MIN / 2 (callers treat as "minus infinity"). */
int32_t fxp_log2_q16(uint32_t v);

/* The same function inline, for hot loops (fxp_log2_q16 is this body). */
static inline int32_t fxp_log2_q16_inl(uint32_t v)
{
	if (v == 0)
		return INT32_MIN / 2;
	int n = clz32(v);
	uint32_t u = v << n;                 /* bit 31 set: u = 2^31 (1 + f) */
	uint32_t idx = (u >> 24) & 127;      /* top 7 bits of f */
	uint32_t fr = (u >> 8) & 0xffff;     /* next 16 bits */
	int32_t lo = fxp_log2_tab_q16[idx], d = fxp_log2_tab_q16[idx + 1] - lo; /* d <= 739 */
	int32_t y = lo + (int32_t)(((uint32_t)d * fr + 0x8000u) >> 16); /* d * fr < 2^26 */
	return ((31 - n) << 16) + y;
}

/* log2(m * 2^e) in Q16, m > 0. */
static inline int32_t fxp_log2_pf(uint32_t m, int e)
{
	return fxp_log2_q16(m) + e * FXP_Q16_ONE;
}

/* 2^x for x <= 0 (Q16), result Q30 (1.0 = 2^30). Relative error < 2e-5. x < -31 returns 0. */
uint32_t fxp_exp2_neg_q30(int32_t x);

/* The same function inline, for hot loops (fxp_exp2_neg_q30 is this body). */
static inline uint32_t fxp_exp2_neg_q30_inl(int32_t x)
{
	FXP_ASSERT(x <= 0);
	int ip = -((-x + 0xffff) >> 16);      /* floor(x / 65536) without shifting a negative value */
	if (ip < -31)
		return 0;
	uint32_t f = (uint32_t)(x - ip * 65536); /* [0, 65536) */
	uint32_t idx = f >> 9, fr = f & 511;
	uint32_t lo = fxp_exp2_tab_q30[idx], d = fxp_exp2_tab_q30[idx + 1] - lo; /* d < 2^24 */
	uint32_t y = lo + ((((d >> 4) * fr) + 16) >> 5);                        /* (d>>4)*fr < 2^29 */
	return (y + ((1u << (-ip)) >> 1)) >> (-ip);
}

/* 2^x for any Q16 x as mantissa in [2^30, 2^31) (Q30) and integer exponent: 2^x = *m * 2^(*e - 30). */
void fxp_exp2_split(int32_t x, uint32_t *m, int *e);

/* cos(pi j / 1024) in Q30 for any integer j (period 2048). */
int32_t fxp_cos2048_q30(int32_t j);

/* (a * b) >> 30 without 64-bit math, a and b int32. Requires |a| < 2^29 and |b| <= 2^30 (bounds the
 * partial products; see fxp_math.c). Rounded, error <= 1 LSB. */
int32_t fxp_mul_q30(int32_t a, int32_t b);

/* Angles as uint32 "turns": 2^32 = 2 pi, so addition and integer multiples wrap exactly (mod 2 pi). */

/* atan2(y, x) as a uint32 angle in [0, 2^32). Any int32 inputs; (0, 0) -> 0. Error < 3e-5 rad (16-bit min/max ratio). */
uint32_t fxp_atan2_u32(int32_t y, int32_t x);

/* cos and sin of a uint32 angle, Q15 (clamped to +-32767). Error < 4e-5. */
void fxp_sincos_q15(uint32_t a, int16_t *c, int16_t *s);

#endif
