/* fxp_math.c — see fxp_math.h. */
#include "fxp_math.h"

int32_t fxp_log2_q16(uint32_t v)
{
	return fxp_log2_q16_inl(v);
}

uint32_t fxp_exp2_neg_q30(int32_t x)
{
	return fxp_exp2_neg_q30_inl(x);
}

void fxp_exp2_split(int32_t x, uint32_t *m, int *e)
{
	int ip = x >= 0 ? (x >> 16) : -((-x + 0xffff) >> 16);
	uint32_t f = (uint32_t)(x - ip * 65536);
	uint32_t idx = f >> 9, fr = f & 511;
	uint32_t lo = fxp_exp2_tab_q30[idx], d = fxp_exp2_tab_q30[idx + 1] - lo;
	*m = lo + ((((d >> 4) * fr) + 16) >> 5); /* [2^30, 2^31] */
	*e = ip;
}

int32_t fxp_cos2048_q30(int32_t j)
{
	uint32_t q = (uint32_t)j & 2047;
	if (q > 1024)
		q = 2048 - q;                    /* cos is even */
	return q > 512 ? -fxp_cosq_q30[1024 - q] : fxp_cosq_q30[q];
}

int32_t fxp_mul_q30(int32_t a, int32_t b)
{
	/* a = ah 2^16 + al, b = bh 2^16 + bl, al, bl in [0, 65536).
	 * |ah| <= 2^13, |bh| <= 2^14: ah*bh*4 <= 2^29; |ah*bl| < 2^29; |bh*al| < 2^30; the low product
	 * al*bl / 2^16 (via (al>>2)(bl>>2) >> 12, < 2^16) and the rounding constant keep the sum < 2^31.
	 * Rounded to nearest: error <= 1 LSB (no systematic bias, which matters when products chain). */
	FXP_ASSERT(a < (1 << 29) && a > -(1 << 29) && b <= (1 << 30) && b >= -(1 << 30));
	int32_t ah = a >> 16, bh = b >> 16;
	int32_t al = a & 0xffff, bl = b & 0xffff;
	int32_t lo = (int32_t)((((uint32_t)al >> 2) * ((uint32_t)bl >> 2)) >> 12);
	return ah * bh * 4 + ((ah * bl + bh * al + lo + 8192) >> 14);
}

uint32_t fxp_atan2_u32(int32_t y, int32_t x)
{
	uint32_t ax = (uint32_t)(x < 0 ? -(uint32_t)x : (uint32_t)x), ay = (uint32_t)(y < 0 ? -(uint32_t)y : (uint32_t)y);
	if ((ax | ay) == 0)
		return 0;
	/* bring max to [2^15, 2^16] (16-bit ratio precision), then t = min/max in Q16, rounded */
	uint32_t mx = ax > ay ? ax : ay;
	int s = 16 - clz32(mx);
	if (s > 0)
	{
		ax = (ax + (1u << (s - 1))) >> s;
		ay = (ay + (1u << (s - 1))) >> s;
	}
	else
	{
		ax <<= -s;
		ay <<= -s;
	}
	int swap = ay > ax;
	uint32_t num = swap ? ax : ay, den = swap ? ay : ax; /* den in [2^15, 2^16], num <= den */
	uint32_t t = num >= den ? 1u << 16 : ((num << 16) + (den >> 1)) / den; /* num < 2^16 here */
	uint32_t idx = t >> 8, fr = t & 255;
	uint32_t lo = fxp_atan_tab_u32[idx];
	uint32_t th = idx < 256 ? lo + (((fxp_atan_tab_u32[idx + 1] - lo) * fr + 128) >> 8) : lo; /* d * fr < 2^30 */
	if (swap)
		th = (1u << 30) - th; /* pi/2 - th */
	if (x < 0)
		th = (1u << 31) - th; /* pi - th */
	if (y < 0)
		th = 0u - th;         /* -th (mod 2 pi) */
	return th;
}

void fxp_sincos_q15(uint32_t a, int16_t *c, int16_t *s)
{
	/* 2048 steps per turn from the Q30 cosine table, linear interpolation on the low 21 bits
	 * (interpolation error <= (pi/1024)^2 / 8 = 1.2e-6), then round to Q15. */
	int32_t j = (int32_t)(a >> 21);
	int32_t f = (int32_t)((a >> 10) & 2047);
	int32_t c0 = fxp_cos2048_q30(j), c1 = fxp_cos2048_q30(j + 1);
	int32_t s0 = fxp_cos2048_q30(j - 512), s1 = fxp_cos2048_q30(j - 511);
	int32_t cv = c0 + ((((c1 - c0) >> 5) * f) >> 6); /* |c1 - c0| < 2^22: (>> 5) * f < 2^28 */
	int32_t sv = s0 + ((((s1 - s0) >> 5) * f) >> 6);
	cv = (cv + (1 << 14)) >> 15;
	sv = (sv + (1 << 14)) >> 15;
	*c = (int16_t)(cv > 32767 ? 32767 : (cv < -32767 ? -32767 : cv));
	*s = (int16_t)(sv > 32767 ? 32767 : (sv < -32767 ? -32767 : sv));
}
