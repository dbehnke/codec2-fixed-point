/*
 * fxp_fft.c — fixed-point 512-point real FFT with block floating point.
 *
 * Structure: 256-point complex radix-4 DIT FFT (base-4 digit-reversed reorder, 4 stages) plus the standard
 * even/odd split that turns it into a 512-point real transform. The inverse uses
 * IFFT(Z) = conj(FFT(conj(Z))), so there is a single complex kernel.
 *
 * Overflow analysis (M = max |component| of the stage input, w a Q15 twiddle with |w| <= 1):
 *  - t = c*w rounded: |t_r|, |t_i| <= |c||w| + 0.5 <= sqrt(2) M + 0.5   (Cauchy-Schwarz).
 *  - Radix-4 outputs are sums of 4 such terms: <= 5.243 M + 1.5, then >> s with rounding.
 *  - s is chosen from B = OR of (x ^ (x >> 31)) over the stage input (B >= max|x| - 1); the rule in
 *    r4_shift() keeps every stage output below 21477 < 2^15 (worked bounds next to it).
 *  - Twiddle products: |cr*wr - ci*wi| <= 2 * 32768 * 32767 + 2^14 < 2^31, so int32 never overflows.
 * Real split (forward): with Z the complex FFT output, f1 = Z[k] + conj Z[N-k], f2 = Z[k] - conj Z[N-k],
 *  components of f1, f2 <= 2M; tw = f2 * (-j W512^k): <= 2 sqrt(2) M + 0.5. X = (f1 + tw) >> (1 + s')
 *  with s' = 0 for B < 2^13 (X <= 4.83 * 8192 / 2 < 19780) and s' = 1 for B < 2^14
 *  (X <= 4.83 * 16384 / 4 < 19780). Z (|Z| < 21477 from the FFT) is pre-shifted by one if B >= 2^14,
 *  keeping f2 <= 32768 so the twiddle product sum stays < 2^31.
 * The inverse pre-split has the same bounds (it computes the same expressions with conjugated
 * twiddles), so it uses the same shift rule.
 *
 * All rounding is round-half-up. FXP_ASSERT checks the bounds in host test builds (-DFXP_CHECK).
 */
#include "fxp_fft.h"

static inline uint32_t mag_bound(int32_t x)
{
	return (uint32_t)(x ^ (x >> 31)); /* |x| for x >= 0, |x| - 1 for x < 0 */
}

static uint32_t block_bound(const cq15_t *x, int n)
{
	uint32_t b = 0;
	for (int i = 0; i < n; i++)
		b |= mag_bound(x[i].r) | mag_bound(x[i].i);
	return b;
}

/* Scale n elements by 2^sh (sh > 0: left, exact; sh < 0: right, rounded). */
static void block_scale(cq15_t *x, int n, int sh)
{
	if (sh > 0) {
		for (int i = 0; i < n; i++) {
			x[i].r = (int16_t)shl(x[i].r, sh);
			x[i].i = (int16_t)shl(x[i].i, sh);
		}
	} else if (sh < 0) {
		for (int i = 0; i < n; i++) {
			x[i].r = (int16_t)shr_r(x[i].r, -sh);
			x[i].i = (int16_t)shr_r(x[i].i, -sh);
		}
	}
}

/* Radix-4 stage shift from the OR bound B of the stage input (outputs <= 5.243 M + 1.5 before >> s):
 *   B < 2^12 -> s = 0: <= 5.243 * 4096 + 1.5 < 21477
 *   B < 2^13 -> s = 1: <= (5.243 * 8192 + 2.5) / 2 < 21477
 *   B < 2^14 -> s = 2: <= (5.243 * 16384 + 4.5) / 4 < 21477
 *   else     -> s = 3: <= (5.243 * 32768 + 8.5) / 8 < 21477
 * Every stage output stays below 2^15, so s <= 3 always suffices. */
static inline int r4_shift(uint32_t b)
{
	return b < (1u << 12) ? 0 : b < (1u << 13) ? 1 : b < (1u << 14) ? 2 : 3;
}

#define FFT_INLINE static inline __attribute__((always_inline))

/* One radix-4 DIT butterfly on p[0], p[q], p[2q], p[3q]; a0 = p[0], a1..a3 the (twiddled) other inputs.
 * Bound: |a_k| components <= |x_k| + 0.5 <= sqrt(2) M + 0.5, outputs <= M + 3 (sqrt(2) M + 0.5)
 * = 5.243 M + 1.5, then >> s (see r4_shift). Returns the OR magnitude bound of the 8 outputs when want,
 * else 0. neg stores the conjugate (imaginary parts negated; |y| < 21477, so the negation is exact). */
FFT_INLINE uint32_t r4(cq15_t *p, const int q, const int s, const int32_t rnd, int32_t a1r, int32_t a1i,
		       int32_t a2r, int32_t a2i, int32_t a3r, int32_t a3i, const int want, const int neg)
{
	int32_t a0r = p[0].r, a0i = p[0].i;
	int32_t b0r = a0r + a2r + rnd, b0i = a0i + a2i + rnd;
	int32_t b1r = a0r - a2r + rnd, b1i = a0i - a2i + rnd;
	int32_t b2r = a1r + a3r, b2i = a1i + a3i;
	int32_t b3r = a1r - a3r, b3i = a1i - a3i;
	int32_t y0r = (b0r + b2r) >> s, y0i = (b0i + b2i) >> s;
	int32_t y2r = (b0r - b2r) >> s, y2i = (b0i - b2i) >> s;
	int32_t y1r = (b1r + b3i) >> s, y1i = (b1i - b3r) >> s; /* b1 - j b3 */
	int32_t y3r = (b1r - b3i) >> s, y3i = (b1i + b3r) >> s; /* b1 + j b3 */
	FXP_ASSERT(y0r >= -32768 && y0r <= 32767 && y0i >= -32768 && y0i <= 32767);
	FXP_ASSERT(y1r >= -32768 && y1r <= 32767 && y1i >= -32768 && y1i <= 32767);
	FXP_ASSERT(y2r >= -32768 && y2r <= 32767 && y2i >= -32768 && y2i <= 32767);
	FXP_ASSERT(y3r >= -32768 && y3r <= 32767 && y3i >= -32768 && y3i <= 32767);
	p[0].r = (int16_t)y0r; p[0].i = (int16_t)(neg ? -y0i : y0i);
	p[q].r = (int16_t)y1r; p[q].i = (int16_t)(neg ? -y1i : y1i);
	p[2 * q].r = (int16_t)y2r; p[2 * q].i = (int16_t)(neg ? -y2i : y2i);
	p[3 * q].r = (int16_t)y3r; p[3 * q].i = (int16_t)(neg ? -y3i : y3i);
	if (!want)
		return 0;
	return mag_bound(y0r) | mag_bound(y0i) | mag_bound(y1r) | mag_bound(y1i)
	     | mag_bound(y2r) | mag_bound(y2i) | mag_bound(y3r) | mag_bound(y3i);
}

/* Rounded Q15 product c * w (|sum of products| < 2^31, see the file header). */
#define CMUL_R(c, w) (((c##r) * (w).r - (c##i) * (w).i + (1 << 14)) >> 15)
#define CMUL_I(c, w) (((c##r) * (w).i + (c##i) * (w).r + (1 << 14)) >> 15)

/* One radix-4 stage with span q (a compile-time constant at every call). Butterfly j of each group of 4q
 * uses W512^(m j 128 / q), m = 1..3, read in order from tw; j = 0 has unit twiddles and no multiply.
 * The butterflies of a stage touch disjoint elements, so the group-major order gives the same values as any
 * other order. Returns the OR magnitude bound of the stage output (0 unless want). */
FFT_INLINE uint32_t r4_stage(cq15_t *x, const int q, const int s, const cq15_t *tw, const int want, const int neg)
{
	const int32_t rnd = (1 << s) >> 1;
	uint32_t nb = 0;
	for (cq15_t *g = x; g < x + 256; g += 4 * q) {
		nb |= r4(g, q, s, rnd, g[q].r, g[q].i, g[2 * q].r, g[2 * q].i, g[3 * q].r, g[3 * q].i, want, neg);
		const cq15_t *t = tw;
		for (cq15_t *p = g + 1; p < g + q; p++, t += 3) {
			int32_t c1r = p[q].r, c1i = p[q].i, c2r = p[2 * q].r, c2i = p[2 * q].i, c3r = p[3 * q].r, c3i = p[3 * q].i;
			nb |= r4(p, q, s, rnd, CMUL_R(c1, t[0]), CMUL_I(c1, t[0]), CMUL_R(c2, t[1]), CMUL_I(c2, t[1]),
				 CMUL_R(c3, t[2]), CMUL_I(c3, t[2]), want, neg);
		}
	}
	return nb;
}

/* The first three stages on digit-reversed input with OR bound b; returns the summed shift, *b the bound. */
static int r4_stages123(cq15_t *x, uint32_t *b)
{
	int s1 = r4_shift(*b);
	uint32_t nb = r4_stage(x, 1, s1, 0, 1, 0);
	int s2 = r4_shift(nb);
	nb = r4_stage(x, 4, s2, fxp_tw256, 1, 0);
	int s3 = r4_shift(nb);
	*b = r4_stage(x, 16, s3, fxp_tw256 + 9, 1, 0);
	return s1 + s2 + s3;
}

/* Forward complex FFT of digit-reversed x with OR bound b; *bout = OR bound of the output. */
static int cfft_fwd(cq15_t *x, uint32_t b, uint32_t *bout)
{
	int total = r4_stages123(x, &b);
	int s = r4_shift(b);
	*bout = r4_stage(x, 64, s, fxp_tw256 + 54, 1, 0);
	return total + s;
}

/* As cfft_fwd, but stores the conjugate of the output and computes no bound. */
static int cfft_fwd_conj(cq15_t *x, uint32_t b)
{
	int total = r4_stages123(x, &b);
	int s = r4_shift(b);
	r4_stage(x, 64, s, fxp_tw256 + 54, 0, 1);
	return total + s;
}

/* base-4 digit-reversal permutation */
static void digrev(cq15_t *x)
{
	for (int i = 0; i < 256; i++) {
		int j = fxp_digrev256[i];
		if (j > i) {
			cq15_t t = x[i];
			x[i] = x[j];
			x[j] = t;
		}
	}
}

/* The digit reversal fused with the normalising block scale (block_scale's arithmetic: 2^sh, sh > 0 exact,
 * sh < 0 rounded): every element is scaled exactly once. */
#define DIGREV_MAP(OP)                                                                  \
	for (int i = 0; i < 256; i++) {                                                     \
		int j = fxp_digrev256[i];                                                       \
		if (j > i) {                                                                    \
			cq15_t t = x[i];                                                            \
			x[i].r = (int16_t)OP(x[j].r); x[i].i = (int16_t)OP(x[j].i);                 \
			x[j].r = (int16_t)OP(t.r); x[j].i = (int16_t)OP(t.i);                       \
		} else if (j == i) {                                                            \
			x[i].r = (int16_t)OP(x[i].r); x[i].i = (int16_t)OP(x[i].i);                 \
		}                                                                               \
	}
static void digrev_scale(cq15_t *x, int sh)
{
	if (sh > 0) {
#define OP_SHL(v) shl((v), sh)
		DIGREV_MAP(OP_SHL)
#undef OP_SHL
	} else if (sh < 0) {
#define OP_SHR(v) shr_r((v), -sh)
		DIGREV_MAP(OP_SHR)
#undef OP_SHR
	} else {
		digrev(x);
	}
}
#undef DIGREV_MAP

/* The conjugation before the inverse's complex FFT fused with the digit reversal; returns the OR bound of
 * the stored (conjugated, permuted) elements, as block_bound would. */
static uint32_t digrev_conj_bound(cq15_t *x)
{
	uint32_t b = 0;
	for (int i = 0; i < 256; i++) {
		int j = fxp_digrev256[i];
		if (j > i) {
			cq15_t t = x[i];
			x[i].r = x[j].r; x[i].i = (int16_t)-x[j].i;
			x[j].r = t.r; x[j].i = (int16_t)-t.i;
			b |= mag_bound(x[i].r) | mag_bound(x[i].i) | mag_bound(x[j].r) | mag_bound(x[j].i);
		} else if (j == i) {
			x[i].i = (int16_t)-x[i].i;
			b |= mag_bound(x[i].r) | mag_bound(x[i].i);
		}
	}
	return b;
}

int fxp_cfft256(cq15_t *x)
{
	uint32_t bout;
	digrev(x);
	return cfft_fwd(x, block_bound(x, 256), &bout);
}

/* Shared split for both directions. For the forward transform (inv = 0):
 *   X[k] = (f1 + f2 * (-j W^k)) / 2,  X[N-k] = conj(f1 - f2 * (-j W^k)) / 2
 * For the inverse pre-split (inv = 1), with input X and conjugated twiddles:
 *   Z[k] = f1 + f2 * (j conj W^k),  Z[N-k] = conj(f1 - f2 * (j conj W^k))
 * Every output is stored as (f1 +- tw) >> sh, where sh (the return value) = pre-shift + 1 + s'.
 * n is the number of valid elements the bound b covers (256 forward, 257 inverse).
 * buf[0] and buf[256] are handled by the callers. */
static int split(cq15_t *buf, int n, int inv, uint32_t b)
{
	int pre = 0;
	if (b >= (1u << 14)) { /* keep f2 <= 32768 so twiddle products fit int32 */
		block_scale(buf, n, -1);
		b = block_bound(buf, n);
		pre = 1;
	}
	const int s = 1 + (b < (1u << 13) ? 0 : 1);
	const int32_t rnd = 1 << (s - 1);
	for (int k = 1; k <= 128; k++) {
		const int32_t wr = fxp_w512[k].r;
		const int32_t wi = inv ? -fxp_w512[k].i : fxp_w512[k].i; /* conj for the inverse */
		cq15_t zk = buf[k], znk = buf[256 - k];
		int32_t f1r = (int32_t)zk.r + znk.r, f1i = (int32_t)zk.i - znk.i;
		int32_t f2r = (int32_t)zk.r - znk.r, f2i = (int32_t)zk.i + znk.i;
		/* p = f2 * w */
		int32_t pr = (f2r * wr - f2i * wi + (1 << 14)) >> 15;
		int32_t pi = (f2r * wi + f2i * wr + (1 << 14)) >> 15;
		/* forward: tw = -j p = (pi, -pr); inverse: tw = +j p = (-pi, pr) */
		int32_t twr = inv ? -pi : pi, twi = inv ? pr : -pr;
		int32_t xr = (f1r + twr + rnd) >> s, xi = (f1i + twi + rnd) >> s;
		int32_t yr = (f1r - twr + rnd) >> s, yi = (-(f1i - twi) + rnd) >> s;
		FXP_ASSERT(xr >= -32768 && xr <= 32767 && xi >= -32768 && xi <= 32767);
		FXP_ASSERT(yr >= -32768 && yr <= 32767 && yi >= -32768 && yi <= 32767);
		buf[k].r = (int16_t)xr; buf[k].i = (int16_t)xi;
		buf[256 - k].r = (int16_t)yr; buf[256 - k].i = (int16_t)yi; /* k = 128 writes the same slot twice, consistently */
	}
	return pre + s;
}

int fxp_rfft512(cq15_t buf[FXP_RFFT_BUF])
{
	/* block_normalise and the digit reversal in one pass. After a left (exact) normalising shift the OR bound
	 * is in [2^13, 2^14): (b << sh) | (2^sh - 1) for b << sh in [2^13, 2^14) with sh <= 13. Only r4_shift reads
	 * it, so 2^13 stands for it; a rounded right shift needs the bound measured. */
	uint32_t b = block_bound(buf, 256);
	int sh = b ? clz32(b) - 18 : 0;
	digrev_scale(buf, sh);
	b = b == 0 ? 0 : sh >= 0 ? 1u << 13 : block_bound(buf, 256);
	FXP_ASSERT(r4_shift(b) == r4_shift(block_bound(buf, 256)));
	int e = -sh + cfft_fwd(buf, b, &b);
	/* DC and Nyquist: X[0] = Zr + Zi, X[256] = Zr - Zi (real); shift with the rest below */
	int32_t z0r = buf[0].r, z0i = buf[0].i;
	int ssh = split(buf, 256, 0, b);
	/* split() applied (pre + s) shifts where the /2 is part of s. X[0] has no /2 in the formula, so
	 * its value is (z0r + z0i) scaled by 2^-(ssh-1) relative to the others' 2^-ssh: use >> (ssh - 1). */
	int s0 = ssh - 1;
	int32_t x0 = z0r + z0i, xn = z0r - z0i;
	if (s0 > 0) {
		x0 = shr_r(x0, s0);
		xn = shr_r(xn, s0);
	}
	FXP_ASSERT(x0 >= -32768 && x0 <= 32767 && xn >= -32768 && xn <= 32767);
	buf[0].r = (int16_t)x0; buf[0].i = 0;
	buf[256].r = (int16_t)xn; buf[256].i = 0;
	return e + ssh - 1; /* split's /2 is part of the X definition: undo it in the exponent */
}

int fxp_irfft512(cq15_t buf[FXP_RFFT_BUF])
{
	buf[0].i = 0;
	buf[256].i = 0;
	/* normalise as fxp_rfft512; split() only compares the bound with 2^13 and 2^14 */
	uint32_t b = block_bound(buf, 257);
	int sh = b ? clz32(b) - 18 : 0;
	block_scale(buf, 257, sh);
	b = b == 0 ? 0 : sh >= 0 ? 1u << 13 : block_bound(buf, 257);
	FXP_ASSERT(b == 0 || sh < 0 || (block_bound(buf, 257) >= (1u << 13) && block_bound(buf, 257) < (1u << 14)));
	int e = -sh;
	int32_t x0 = buf[0].r, xn = buf[256].r;
	int ssh = split(buf, 257, 1, b);
	/* Z[0] = (X0 + XN) + j (X0 - XN). The other Z[k] = f1 + tw are stored >> ssh with no /2 in their
	 * definition, so Z[0] takes the same shift (x0, xn were captured before any pre-shift). */
	int32_t z0r = shr_r(x0 + xn, ssh), z0i = shr_r(x0 - xn, ssh);
	FXP_ASSERT(z0r >= -32768 && z0r <= 32767 && z0i >= -32768 && z0i <= 32767);
	buf[0].r = (int16_t)z0r; buf[0].i = (int16_t)z0i;
	/* IFFT(Z) = conj(FFT(conj Z)): the first conj is folded into the digit reversal, the second into the
	 * last radix-4 stage (|i| <= 19780 after the split and < 21477 after the FFT, so both negations are exact) */
	e += cfft_fwd_conj(buf, digrev_conj_bound(buf));
	return e + ssh;
}
