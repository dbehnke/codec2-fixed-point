/*
 * lpc.c — encoder LPC/LSP analysis in fixed point (stage E3, docs/11-fxp-lpc.md) and the decoder's LPC
 * spectrum entry point.
 *
 * Same algorithm as the float Codec2 encoder (speech_to_uq_lsps): Hann-windowed speech, autocorrelation,
 * order-10 LPC, energy E = sum a R, 15 Hz bandwidth expansion, LSP root search on a 0.01 grid of x = cos w
 * with 6 bisections, then w = fast_acosf(x).
 *
 * Number formats:
 *   x[279]    int16 windowed speech: round(Sn h_q15 / 2^s), s >= 0 the smallest giving |x| <= 32767 (as
 *             dft_speech). Float: Wn = x 2^(s - 15) scale, scale = 1 / sqrt(512 m).
 *   R[11]     int64 exact sum x[i] x[i + j]; then int32 r = R 2^-sh with r[0] in [2^29, 2^30)
 *   Schur     reflection coefficients rc (Q30) by the Schur recursion: every intermediate is bounded by the
 *             current prediction error, and both generator rows are renormalised (common shift S) after each
 *             order, so each order keeps ~29 bits relative to its own error. Prediction error E = C1[0] 2^-S.
 *   a[11]     Q23 direct-form LPC by the step-up recursion from rc. |a_j| <= C(10, j) <= 252 < 256 for any
 *             minimum-phase A(z), so Q23 cannot overflow.
 *   P, Q      Q23 (|coefficient| <= 252: roots on the unit circle), then block-normalised per polynomial
 *   x grid    integer u in units of 1/6400: the float grid 1 - 0.01 n is u = 6400 - 64 n, and 6 bisections
 *             of a 64-wide interval end on integers. x = u / 6400 in Q29.
 *   T_k, sums Q29 Chebyshev recursion; only the signs of the sums are used
 *   LSPs      Hz in Q16 (fast_acosf in Q29 with an exact 64-bit integer square root)
 *   log2 E    Q16
 */
#include "codec2_internal.h"
#include "lpc.h"
#include "fxp.h"
#include "fxp_math.h"
#include "fxp_lpcspec.h"
#include "analysis.h"

#define W_START (M_PITCH / 2 - NW / 2) /* first windowed sample of Sn */
#define U_ONE 6400                     /* x = 1 in grid units */
#define U_STEP 64                      /* LSP_DELTA1 = 0.01 */

_Static_assert(NW % 2 == 1 && W_START == 21, "LPC window span");

static inline uint32_t uabs32(int32_t v) { return v < 0 ? 0u - (uint32_t)v : (uint32_t)v; }
static inline int bitlen(uint32_t v) { return v ? 32 - clz32(v) : 0; }
static inline int bitlen64(uint64_t v) { return v >> 32 ? 32 + bitlen((uint32_t)(v >> 32)) : bitlen((uint32_t)v); }

/* (a b) >> 29, rounded down per partial product (error < 4 LSB). Requires |a b| < 2^59.9 and |a|, |b| < 2^31
 * with |a >> 16| |b >> 16| < 2^28 (checked under FXP_CHECK). No 64-bit multiply. */
static inline int32_t mul_q29(int32_t a, int32_t b)
{
	int32_t ah = a >> 16, bh = b >> 16;
	uint32_t al = (uint32_t)a & 0xffffu, bl = (uint32_t)b & 0xffffu;
	FXP_ASSERT(ah * bh < (1 << 28) && ah * bh > -(1 << 28));
	return ah * bh * 8 + ((ah * (int32_t)bl) >> 13) + (((int32_t)al * bh) >> 13) + (int32_t)((al * bl) >> 29);
}

/* (a b) >> 30, same construction; |a >> 16| |b >> 16| < 2^29 */
static inline int32_t mul_q30(int32_t a, int32_t b)
{
	int32_t ah = a >> 16, bh = b >> 16;
	uint32_t al = (uint32_t)a & 0xffffu, bl = (uint32_t)b & 0xffffu;
	FXP_ASSERT(ah * bh < (1 << 29) && ah * bh > -(1 << 29));
	return ah * bh * 4 + ((ah * (int32_t)bl) >> 14) + (((int32_t)al * bh) >> 14) + (int32_t)((al * bl) >> 30);
}

/* v 2^sh for any sign of sh (|sh| < 32), rounding half up on right shifts; the caller guarantees the result fits */
static inline int32_t shift_r(int32_t v, int sh)
{
	if (sh >= 0)
		return shl(v, sh);
	return ((v >> (-sh - 1)) + 1) >> 1;
}

/* round(n 2^30 / d) for 0 <= n <= d < 2^30, d > 0 (restoring division, no 64-bit arithmetic) */
static uint32_t div_q30(uint32_t n, uint32_t d)
{
	uint32_t q = 0;
	for (int i = 0; i < 31; i++)
	{
		n <<= 1;
		q <<= 1;
		if (n >= d)
		{
			n -= d;
			q |= 1;
		}
	}
	return (q + 1) >> 1; /* q = floor(n 2^31 / d) <= 2^31 - 1 */
}

/* floor(sqrt(n)), bit by bit (no multiply) */
static uint32_t isqrt64(uint64_t n)
{
	uint64_t res = 0, bit = (uint64_t)1 << 62;
	while (bit > n)
		bit >>= 2;
	while (bit)
	{
		if (n >= res + bit)
		{
			n -= res + bit;
			res = (res >> 1) + bit;
		}
		else
			res >>= 1;
		bit >>= 2;
	}
	return (uint32_t)res;
}

/* Reflection coefficients rc[0..9] (Q30) of the normalised autocorrelation r (r[0] in [2^29, 2^30)) by the
 * Schur recursion (sign convention of the float Levinson: a_i = rc). As the float code, |rc| > 1 is replaced
 * by 0. Returns the prediction error mantissa; *esh is its exponent: E / r[0]-units = ret 2^-esh. */
static int32_t schur(const int32_t *restrict r, int32_t *restrict rc, int *restrict esh)
{
	int32_t C0[LPC_ORD + 1], C1[LPC_ORD + 1];
	int S = 0;
	for (int n = 0; n <= LPC_ORD; n++)
		C0[n] = C1[n] = r[n];
	for (int k = 0; k < LPC_ORD; k++)
	{
		int32_t num = C0[k + 1], den = C1[0];
		uint32_t an = uabs32(num);
		int32_t q = 0;
		if (den > 0 && an <= (uint32_t)den)
		{
			q = (int32_t)div_q30(an, (uint32_t)den);
			q = num > 0 ? -q : q; /* rc = -num / den */
		}
		rc[k] = q;
		/* |C| < 2^29 before the update, so |t + rc u| < 2^30 */
		uint32_t o = 0;
		for (int n = 0; n < LPC_ORD - k; n++)
		{
			int32_t t1 = C0[n + k + 1], t2 = C1[n];
			C0[n + k + 1] = t1 + mul_q30(q, t2);
			C1[n] = t2 + mul_q30(q, t1);
			o |= uabs32(C0[n + k + 1]) | uabs32(C1[n]);
		}
		/* renormalise the rows still in use to max < 2^29 */
		int sh = o ? 29 - bitlen(o) : 0;
		if (sh)
		{
			for (int n = 0; n < LPC_ORD - k; n++)
			{
				C0[n + k + 1] = shift_r(C0[n + k + 1], sh);
				C1[n] = shift_r(C1[n], sh);
			}
			S += sh;
		}
	}
	*esh = S;
	return C1[0];
}

/* Chebyshev sum c[0] T5 + c[1] T4 + ... + c[5] T0 at x (Q29); |c| < 2^27 */
static int32_t cheb(const int32_t *c, int32_t x)
{
	int32_t T0 = 1 << 29, T1 = x;
	int32_t s = c[5] + mul_q29(c[4], T1);
	int32_t T2 = 2 * mul_q29(x, T1) - T0;
	s += mul_q29(c[3], T2);
	int32_t T3 = 2 * mul_q29(x, T2) - T1;
	s += mul_q29(c[2], T3);
	int32_t T4 = 2 * mul_q29(x, T3) - T2;
	s += mul_q29(c[1], T4);
	int32_t T5 = 2 * mul_q29(x, T4) - T3;
	s += mul_q29(c[0], T5);
	return s;
}

/* x = u / 6400 in Q29: u 2^29 / 6400 = u 83886 + u 0.08 */
static inline int32_t u_to_q29(int32_t u)
{
	return u * 83886 + (u * 2 + (u < 0 ? -12 : 12)) / 25;
}

static inline int same_sign(int32_t a, int32_t b) { return (a > 0 && b > 0) || (a < 0 && b < 0); }

/* LSP roots of the Q23 LPC polynomial a (bandwidth expanded), alternating P' and Q', as the float lpc_to_lsp.
 * u[j] = cos(w_j) in 1/6400 units. Returns the number of roots found. */
static int lpc_to_lsp_u(const int32_t *a, int32_t *u)
{
	int32_t P[LPC_ORD / 2 + 1], Q[LPC_ORD / 2 + 1], cp[6], cq[6];
	P[0] = Q[0] = 1 << 23;
	for (int i = 1; i <= LPC_ORD / 2; i++)
	{
		/* modulo 2^32: a[i] +- a[11 - i] can exceed int32, the results (|P|, |Q| <= 252 2^23) cannot */
		P[i] = (int32_t)((uint32_t)a[i] + (uint32_t)a[LPC_ORD + 1 - i] - (uint32_t)P[i - 1]);
		Q[i] = (int32_t)((uint32_t)a[i] - (uint32_t)a[LPC_ORD + 1 - i] + (uint32_t)Q[i - 1]);
	}
	/* Chebyshev coefficients {2 P0, ..., 2 P4, P5}, halved (only signs matter): {P0, ..., P4, P5 / 2},
	 * normalised so that each |c| < 2^27 (the sum of six stays < 2^30) */
	for (int pq = 0; pq < 2; pq++)
	{
		const int32_t *src = pq ? Q : P;
		int32_t *c = pq ? cq : cp;
		uint32_t o = uabs32(src[5]) >> 1;
		for (int i = 0; i < 5; i++)
			o |= uabs32(src[i]);
		int sh = 27 - bitlen(o);
		for (int i = 0; i < 5; i++)
			c[i] = shift_r(src[i], sh);
		c[5] = shift_r(src[5], sh - 1);
	}

	int roots = 0;
	int32_t xl = U_ONE;
	for (int j = 0; j < LPC_ORD; j++)
	{
		const int32_t *c = (j & 1) ? cq : cp;
		int32_t xr = xl;
		int32_t sl = cheb(c, u_to_q29(xl));
		while (xr >= -U_ONE)
		{
			xr = xl - U_STEP;
			int32_t sr = cheb(c, u_to_q29(xr));
			if ((sr <= 0 && sl >= 0) || (sr >= 0 && sl <= 0))
			{
				int32_t xm = xl;
				roots++;
				for (int b = 0; b < 6; b++)
				{
					xm = (xl + xr) >> 1; /* xl - xr is even until the last step */
					int32_t sm = cheb(c, u_to_q29(xm));
					if (same_sign(sm, sl))
					{
						xl = xm;
						sl = sm;
					}
					else
						xr = xm;
				}
				u[j] = xm;
				xl = xm;
				break;
			}
			sl = sr;
			xl = xr;
		}
	}
	return roots;
}

/* fast_acosf(u / 6400) in Hz, Q16: (4000 / pi) sqrt(1 - |x|) p(|x|), mirrored for x < 0 */
static int32_t u_to_hz_q16(int32_t u)
{
	int neg = u < 0;
	int32_t au = neg ? -u : u;
	if (au > U_ONE)
		au = U_ONE;
	int32_t x = u_to_q29(au);
	int32_t p = fxp_acos_poly_q29[0];
	for (int k = 1; k < 6; k++)
		p = mul_q29(p, x) + fxp_acos_poly_q29[k];
	/* sqrt(1 - x) = sqrt(d) / 80, d = 6400 - |u|; in Q29: sqrt(d 2^44) 2^7 / 80 = isqrt(d 2^44) 8 / 5 */
	uint32_t sq = isqrt64((uint64_t)(U_ONE - au) << 44); /* <= 80 2^22 */
	int32_t s = (int32_t)(sq / 5 * 8 + sq % 5 * 8 / 5);
	int32_t y = mul_q29(s, p); /* acos in Q29 */
	int32_t hz = mul_q29(y, fxp_lpc_hz_q29);
	return neg ? (4000 << 16) - hz : hz;
}

/* Benign LSPs i pi / 10 of the float code (used for zero energy and failed root searches), Hz Q16 */
static void benign_lsps(int32_t *lsp_hz)
{
	for (int i = 0; i < LPC_ORD; i++)
		lsp_hz[i] = (400 * i) << 16;
}

int32_t speech_to_lsps_q16(codec2_encoder_t *c2, int32_t *restrict lsp_hz)
{
	const int16_t *Sn = c2->Sn;
	int16_t *x = (int16_t *)c2->scratch;

	/* window: |Sn h| <= 2^30 */
	uint32_t mx = 0;
	int32_t *x32 = c2->scratch + (NW + 1) / 2; /* products, after the int16 x */
	_Static_assert(sizeof(c2->scratch) >= ((NW + 1) / 2 + NW) * sizeof(int32_t), "LPC window scratch");
	for (int i = 0; i < NW; i++)
	{
		int32_t v = (int32_t)Sn[W_START + i] * (int32_t)fxp_hann_q15[i];
		uint32_t a = uabs32(v);
		x32[i] = v;
		if (a > mx)
			mx = a;
	}
	int s = 0;
	if (mx > 32767)
	{
		s = 17 - clz32(mx);
		if (((mx >> (s - 1)) + 1) >> 1 > 32767)
			s++;
	}
	for (int i = 0; i < NW; i++)
		x[i] = (int16_t)(s ? ((x32[i] >> (s - 1)) + 1) >> 1 : x32[i]);

	/* autocorrelation, exact: |x| <= 32767 so each product <= 32767^2 < 2^30, and a pair of products
	 * (<= 2 32767^2 < 2^31 - 1) is summed in int32 before the one 64-bit add */
	_Static_assert(2LL * 32767 * 32767 <= INT32_MAX, "autocorrelation pair sum");
	int64_t R[LPC_ORD + 1];
	for (int j = 0; j <= LPC_ORD; j++)
	{
		int64_t acc = 0;
		const int n = NW - j;
		int i = 0;
		for (; i + 1 < n; i += 2)
			acc += (int32_t)x[i] * x[i + j] + (int32_t)x[i + 1] * x[i + 1 + j];
		if (i < n)
			acc += (int32_t)x[i] * x[i + j];
		R[j] = acc;
	}

	/* zero-energy trap of the float code: sum (Sn w)^2 < 1e-6. With s >= 1, R[0] >= 2^28 is far above. */
	if (s == 0 && (uint64_t)R[0] < fxp_lpc_floor_r0)
	{
		benign_lsps(lsp_hz);
		return LOG2E_ZERO;
	}

	/* normalise: r[0] in [2^29, 2^30) */
	int sh = bitlen64((uint64_t)R[0]) - 30;
	int32_t r[LPC_ORD + 1];
	for (int j = 0; j <= LPC_ORD; j++)
	{
		int64_t v = R[j];
		if (sh > 0)
			v = ((v >> (sh - 1)) + 1) >> 1;
		else
			v = (int64_t)((uint64_t)v << -sh); /* |R[j]| <= R[0] < 2^30, no overflow; a shift, not a 64-bit multiply */
		r[j] = (int32_t)v;
	}

	int32_t rc[LPC_ORD];
	int esh;
	int32_t em = schur(r, rc, &esh);
	if (em < 1)
		em = 1;
	/* log2 E = log2(em 2^-esh 2^sh 2^(2 s - 30) scale^2) */
	int32_t log2E = fxp_log2_q16((uint32_t)em) + (sh - esh + 2 * s - 30) * 65536 + fxp_lpc_log2_wscale_q16;

	/* step-up to direct form (Q23), as the float Levinson's a[j] = a_prev[j] + k a_prev[i - j] */
	int32_t a[LPC_ORD + 1], ap[LPC_ORD + 1];
	a[0] = 1 << 23;
	for (int i = 1; i <= LPC_ORD; i++)
	{
		for (int j = 1; j < i; j++)
			ap[j] = a[j];
		int32_t k = rc[i - 1];
		a[i] = shift_r(k, -7);
		for (int j = 1; j < i; j++)
			a[j] = ap[j] + mul_q30(k, ap[i - j]);
	}
	/* bandwidth expansion */
	for (int i = 1; i <= LPC_ORD; i++)
		a[i] = mul_q30(a[i], fxp_bw_gamma_q30[i]);

	int32_t u[LPC_ORD];
	if (lpc_to_lsp_u(a, u) != LPC_ORD)
	{
		benign_lsps(lsp_hz);
		return log2E;
	}
	for (int i = 0; i < LPC_ORD; i++)
		lsp_hz[i] = u_to_hz_q16(u[i]);
	return log2E;
}

void aks_to_mag2_lsp(codec2_decoder_t *c2, const int16_t *hz2, int32_t log2E, uint32_t W, int L,
                     int32_t *lA, uint32_t *phA)
{
    _Static_assert(sizeof(c2->scratch) >= 514 * sizeof(int32_t), "fxp_aks_to_mag2 needs 514 scratch words");
    fxp_aks_to_mag2(hz2, log2E, W, L, lA, phA, c2->scratch);

    /* apply_lpc_correction in the log domain: A[1] *= 0.032 when Wo < pi 150 / 4000 */
    if (W < 80530637u && lA[1] > INT32_MIN / 4)
        lA[1] += -325438; /* log2(0.032) in Q16 */
}
