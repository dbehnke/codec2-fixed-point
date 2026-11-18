/*
 * fxp.h — fixed-point primitives for the Codec2 3200 port.
 *
 * Conventions (apply to every fixed-point file in this tree):
 *  - Qn means a signed integer x representing x * 2^-n. Q15 lives in int16_t, range [-1, 1 - 2^-15].
 *  - Products of two int16_t are formed in int32_t (exact; CK803S `mulsh`, 1 instruction).
 *    |a*b| <= 2^30, so the sum of two such products also fits in int32_t.
 *  - Right shifts of negative signed values are arithmetic (floor). C leaves this
 *    implementation-defined; GCC and clang define it as arithmetic on every target we use.
 *    Left shifts are only applied to non-negative values or done in uint32_t (avoids UB).
 *  - Rounding: "_R" helpers round half up (add 2^(s-1) before shifting). Plain shifts truncate
 *    toward -inf. Each call site states which it uses.
 *  - Saturation is explicit (sat16). Code that relies on a proven bound instead of saturating says
 *    so in a comment and checks the bound with FXP_ASSERT in host test builds (-DFXP_CHECK).
 *  - No 64-bit arithmetic: CK803S has no 32x32->64 multiply (each one is a __muldi3 call).
 */
#ifndef CODEC2_FXP_H
#define CODEC2_FXP_H

#include <stdint.h>

#ifdef FXP_CHECK
#include <assert.h>
#define FXP_ASSERT(c) assert(c)
#else
#define FXP_ASSERT(c) ((void)0)
#endif

typedef int16_t q15_t;

/* Saturate int32 to int16 range. */
static inline int16_t sat16(int32_t x)
{
	return (int16_t)(x > 32767 ? 32767 : (x < -32768 ? -32768 : x));
}

/* Arithmetic right shift with round-half-up; s in [1, 31]. Caller guarantees x + 2^(s-1) does not
 * overflow (true whenever |x| < 2^31 - 2^30). */
static inline int32_t shr_r(int32_t x, int s)
{
	return (x + (1 << (s - 1))) >> s;
}

/* Left shift without UB for negative x; s in [0, 31]. Caller guarantees the result fits. */
static inline int32_t shl(int32_t x, int s)
{
	return (int32_t)((uint32_t)x << s);
}

/* Q15 x Q15 -> Q15, rounded. Exact product in int32, then round. Result fits int16 unless both
 * inputs are -32768 (product 2^30 -> 32768); callers must not pass that pair (twiddles are never -32768). */
static inline int16_t mul_q15_r(int16_t a, int16_t b)
{
	return (int16_t)(((int32_t)a * b + (1 << 14)) >> 15);
}

/* int32 x Q15 -> int32 scaled by 2^-15 (Speex MULT16_32_Q15). Uses two 32-bit-low products, no
 * 64-bit math. a = hi*2^16 + lo with lo in [0, 65535]; lo*b fits int32 (|b| <= 2^15, lo < 2^16).
 * Truncates (floor of each partial); error < 2 LSB. Not valid for (a < -2^30 && b == -32768). */
static inline int32_t mul32x16_q15(int32_t a, int16_t b)
{
	return (a >> 16) * (int32_t)b * 2 + (((int32_t)(a & 0xffff) * (int32_t)b) >> 15);
}

/* Saturating int32 add. */
static inline int32_t add_sat32(int32_t a, int32_t b)
{
	int32_t r = (int32_t)((uint32_t)a + (uint32_t)b);
	if (((a ^ r) & (b ^ r)) < 0)
		r = a < 0 ? INT32_MIN : INT32_MAX;
	return r;
}

/* Saturating left shift of an int32 by s in [0, 31]. */
static inline int32_t wsat(int32_t x, int s)
{
	if (s == 0)
		return x;
	if (x > (INT32_MAX >> s))
		return INT32_MAX;
	if (x < (INT32_MIN >> s))
		return INT32_MIN;
	return shl(x, s);
}

/* Count leading zeros of a nonzero uint32 (CK803S: ff1). Returns 32 for 0. */
static inline int clz32(uint32_t x)
{
	return x ? __builtin_clz(x) : 32;
}

/* Number of left shifts that bring max|x| (given as an OR-accumulated bound m = OR of (x ^ (x>>31)))
 * to the top of a 16-bit word without overflow: result r such that (m << r) < 2^15. m == 0 -> 15. */
static inline int norm16_shift(uint32_t m)
{
	return m ? clz32(m) - 17 : 15;
}

#endif /* CODEC2_FXP_H */
