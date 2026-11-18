/*
 * fxp_lpcspec.c — fixed-point aks_to_mag2. Formats and error budget: docs/05-fxp-aks-to-mag2.md.
 *
 * Float reference (lpc.c aks_to_mag2), per bin k = 0..255:
 *   A2 = |A|^2 + 1e-6,  A2g = |A(z/g)|^2 + 1e-6  (g = 0.5),  R = sqrt(A2g / A2)
 *   Pw = (1/A2) (R + 1e-5)^(2 beta)  (beta = 0.2),  e_before = sum 1/A2,  e_after = sum Pw
 *   Pw *= 1.96 for k < 64;  gain = E e_before / e_after;  A[m] = sqrt(gain sum_{am<=k<bm} Pw)
 * Fixed point does all of this in the log2 domain (Q16), where powers become multiplies:
 *   log2 Pw = -1.2 log2 A2 + 0.2 log2 A2g  (the 1e-5 inside the power is dropped: R >= ~0.03 in practice,
 *   so it changes Pw by < 0.03%), sums are taken relative to their maximum term (exp2 in Q30, accumulated
 *   in int32 with headroom), and A[m] = 2^(0.5 (log2 E + log2 e_before - log2 e_after + log2 Em)).
 *
 * |A| comes from the LSP product form (no cancellation, so formant peaks keep full relative precision
 * even for near-unstable filters). A(z/g) is well conditioned for g = 0.5 and uses the sparse DFT of
 * ak g^i, with ak from a fixed-point lsp_to_lpc.
 */
#include "fxp_lpcspec.h"

#ifdef FXP_LPCSPEC_DEBUG
int32_t fxp_dbg_la[256], fxp_dbg_lag[256], fxp_dbg_leb, fxp_dbg_lea;
#endif

int32_t fxp_lsp_cos_q30(int32_t hz2)
{
	/* angle = pi hz2 / 8000 = (pi / 1024) * (hz2 * 16 / 125) */
	FXP_ASSERT(hz2 >= 0 && hz2 < (1 << 26));
	int32_t t = hz2 * 16;
	int32_t j = t / 125, rem = t - j * 125;        /* delta = rem/125 * pi/1024 rad < 0.0031 */
	int32_t c = fxp_cos2048_q30(j);
	if (rem == 0)
		return c;
	int32_t s = fxp_cos2048_q30(j - 512);          /* sin(pi j / 1024) */
	int32_t d = (rem * 843315 + 16) >> 5;          /* delta in Q30: pi / 128000 * 2^35 = 843314.85; < 2^22 */
	int32_t d8 = d >> 8;                           /* Q22, < 2^14 */
	int32_t half_d2 = (d8 * d8) >> 15;             /* delta^2 / 2 in Q30 (Q44 product < 2^28) */
	/* second-order Taylor: cos(a + d) = cos a - d sin a - d^2/2 cos a  (remainder <= d^3/6 < 5e-9);
	 * the products are full Q30 (fxp_mul_q30: |d| < 2^22, |s|, |c| <= 2^30), error <= 4 LSB */
	return c - fxp_mul_q30(d, s) - fxp_mul_q30(half_d2, c);
}

static int32_t lsp_c_q29(int32_t c_q30)
{
	return (c_q30 + 1) >> 1;
}

void fxp_lsp_to_lpc_q19(const int32_t *c_q30, int32_t *ak_q19)
{
	/* P, Q built by convolving second-order sections 1 - 2c z^-1 + z^-2, in Q22. Bounds: with |c| <= 1
	 * each product is dominated coefficient-wise by (1 + z^-1)^(2s), so after s sections every
	 * coefficient is <= C(2s, s): 70 before the last section (70 * 2^22 < 2^29, as fxp_mul_q30 needs),
	 * 252 after it (< 2^30), and P[j] + P[j-1] <= C(11, 5) = 462 (< 2^31). */
	int32_t p[FXP_LPC_ORD + 2] = {1 << 22}, q[FXP_LPC_ORD + 2] = {1 << 22};
	int n = 1;
	for (int s = 0; s < FXP_LPC_ORD / 2; s++, n += 2)
	{
		int32_t c0 = c_q30[2 * s], c1 = c_q30[2 * s + 1];
		for (int j = n + 1; j >= 0; j--) /* in place, high to low */
		{
			int32_t p1 = j >= 1 && j - 1 < n ? p[j - 1] : 0, p2 = j >= 2 ? p[j - 2] : 0;
			int32_t q1 = j >= 1 && j - 1 < n ? q[j - 1] : 0, q2 = j >= 2 ? q[j - 2] : 0;
			int32_t p0 = j < n ? p[j] : 0, q0 = j < n ? q[j] : 0;
			p[j] = p0 - 2 * fxp_mul_q30(p1, c0) + p2;
			q[j] = q0 - 2 * fxp_mul_q30(q1, c1) + q2;
		}
	}
	/* A[j] = (P[j] + P[j-1] + Q[j] - Q[j-1]) / 2, j = 0..10, then Q22 -> Q19 with rounding:
	 * (pj + qj) / 16 computed as pj/16 + qj/16 to stay in int32 (|pj|, |qj| <= 462 * 2^22) */
	for (int j = 0; j <= FXP_LPC_ORD; j++)
	{
		int32_t pj = p[j] + (j ? p[j - 1] : 0), qj = q[j] - (j ? q[j - 1] : 0);
		ak_q19[j] = ((pj >> 1) + (qj >> 1) + 4) >> 3;
	}
}

/* ---- magnitude pseudo-float: value = m * 2^e, m unsigned in [2^14, 2^15] or 0 ---- */
typedef struct
{
	uint32_t m;
	int e;
} pf_t;

/* Round d >= 0 to 15 significant bits: d ~= d16 * 2^(17 - n), d16 in [2^14, 2^15] (2^15 when rounding carries).
 * (d << n) has bit 31 set; its top 16 bits are rounded to 15. Rounding error <= 2^-15 relative. d = 0 gives
 * d16 = 0 (clz of d | 1 keeps the shift defined), and a zero mantissa stays zero through PF_MULN. */
#define PF_NORM(d, n, d16)                                                                                 \
	do                                                                                                 \
	{                                                                                                  \
		(n) = __builtin_clz((d) | 1u);                                                             \
		(d16) = ((((d) << (n)) >> 16) + 1) >> 1;                                                   \
	} while (0)

/* (M, E) *= d16 * 2^(17 - n) up to the caller's constants: M in [2^14, 2^15], d16 in [2^14, 2^15], so the
 * product p is in [2^28, 2^30] and a rounding shift by 14 + t, t = p >> 29 (0, 1 or 2), brings it back to
 * [2^14, 2^15] with no fix-up. E gains t - n; the per-factor constant 14 + 17 + de is added once by the
 * caller. M = 0 or d16 = 0 gives M = 0. */
#define PF_MULN(M, E, d16, n)                                                                              \
	do                                                                                                 \
	{                                                                                                  \
		uint32_t p_ = (M) * (d16);                                                                 \
		uint32_t t_ = p_ >> 29;                                                                    \
		(M) = ((p_ >> (13 + t_)) + 1) >> 1;                                                        \
		(E) += (int)t_ - (n);                                                                      \
	} while (0)

/* (a^2 + b^2) as m * 2^e with m < 2^31 (uint32). */
static inline void pf_sumsq(pf_t a, pf_t b, uint32_t *m, int *e)
{
	if (a.m == 0 || b.m == 0)
	{
		pf_t x = a.m ? a : b;
		*m = x.m * x.m;
		*e = 2 * x.e;
		return;
	}
	int ea = 2 * a.e, eb = 2 * b.e, E = ea > eb ? ea : eb;
	uint32_t sa = a.m * a.m, sb = b.m * b.m; /* each < 2^30 */
	sa = E - ea >= 32 ? 0 : sa >> (E - ea);
	sb = E - eb >= 32 ? 0 : sb >> (E - eb);
	*m = sa + sb;
	*e = E;
}

/* m * 2^e + 1e-6, then log2 in Q16. 1e-6 = 1125899907 * 2^-50. */
static inline __attribute__((always_inline)) int32_t log2_plus_floor(uint32_t m, int e)
{
	const uint32_t fm = 1125899907u;
	const int fe = -50;
	if (m == 0)
		return fxp_log2_q16_inl(fm) + fe * FXP_Q16_ONE;
	/* normalise m to [2^30, 2^31) so both operands have the same scale before aligning */
	int n = clz32(m) - 1;
	m = n >= 0 ? m << n : m >> -n;
	e -= n;
	uint32_t s;
	int E;
	if (e >= fe)
	{
		E = e;
		s = m + (e - fe >= 32 ? 0 : fm >> (e - fe));
	}
	else
	{
		E = fe;
		s = fm + (fe - e >= 32 ? 0 : m >> (fe - e));
	}
	return fxp_log2_q16_inl(s) + E * FXP_Q16_ONE; /* s < 2^32; fxp_log2_pf inline */
}

/* Product form at bin k with the LSP cosines in Q29 (c29[i] = round(c_q30[i] / 2)). Hot loop of the decoder
 * (256 bins x 2 subframes): straight-line, two 15-bit mantissa multiplies per LSP pair, signs tracked as masks. */
static inline __attribute__((always_inline)) int32_t lsp_spec_bin_q29(const int32_t *c29, int k, uint32_t *ang)
{
	FXP_ASSERT(k >= 0 && k < FXP_SPEC_BINS);
	/* fxp_cos2048_q30 values read straight from the quarter-wave table (0 <= k < 256):
	 * cos w_k = cos(pi 4k / 1024), cos(w/2) = cos(pi 2k / 1024), sin(w/2) = cos(pi (512 - 2k) / 1024) */
	int32_t c4 = 4 * k <= 512 ? fxp_cosq_q30[4 * k] : -fxp_cosq_q30[1024 - 4 * k];
	int32_t x = lsp_c_q29(c4);                       /* cos w_k, Q29 */
	uint32_t ch = (uint32_t)fxp_cosq_q30[2 * k];     /* cos(w/2) > 0, Q30 */
	int32_t shs = fxp_cosq_q30[512 - 2 * k];         /* sin(w/2) >= 0, Q30; 0 only at k = 0 */
	int n;
	uint32_t Pm, Qm;
	/* P = 2 cos(w/2) = ch 2^-29 ~= Pm 2^(17 - n - 29); Q likewise (Qm = 0 at k = 0) */
	PF_NORM(ch, n, Pm);
	int Pe = -12 - n;
	PF_NORM((uint32_t)shs, n, Qm);
	int Qe = -12 - n;
	int32_t sgnP = 0, sgnQ = 0; /* 0 or -1 */
#pragma GCC unroll 5
	for (int i = 0; i < FXP_LPC_ORD; i += 2)
	{
		/* factors 2 (x - c): d in Q29 (|d| <= 2^30), value d16 2^(17 - n) 2^-28 = d16 2^(-11 - n) */
		int32_t dp = x - c29[i], dq = x - c29[i + 1];
		int32_t mp = dp >> 31, mq = dq >> 31;
		uint32_t up = (uint32_t)((dp ^ mp) - mp), uq = (uint32_t)((dq ^ mq) - mq);
		sgnP ^= mp;
		sgnQ ^= mq;
		uint32_t d16;
		PF_NORM(up, n, d16);
		PF_MULN(Pm, Pe, d16, n);
		PF_NORM(uq, n, d16);
		PF_MULN(Qm, Qe, d16, n);
	}
	/* per factor 14 (PF_MULN shift) + 17 (PF_NORM) - 28 (2 (x - c) of a Q29 value) = 3 */
	pf_t P = {Pm, Pe + 3 * (FXP_LPC_ORD / 2)}, Q = {Qm, Qe + 3 * (FXP_LPC_ORD / 2)};
	uint32_t m;
	int e;
	pf_sumsq(P, Q, &m, &e);
	int32_t la = log2_plus_floor(m, e - 2); /* |A|^2 = (Pr^2 + Qr^2) / 4 */
	if (ang)
	{
		/* arg A = arg(Pr + j Qr) - 11 w / 2; (Pr, Qr) on a common exponent (<= 2^15) for atan2 */
		int32_t a = 0, b = 0;
		int E = P.m == 0 ? Q.e : (Q.m == 0 ? P.e : (P.e > Q.e ? P.e : Q.e));
		if (P.m)
			a = E - P.e >= 16 ? 0 : (int32_t)(P.m >> (E - P.e));
		if (Q.m)
			b = E - Q.e >= 16 ? 0 : (int32_t)(Q.m >> (E - Q.e));
		a = (a ^ sgnP) - sgnP;
		b = (b ^ sgnQ) - sgnQ;
		/* 11 w / 2 = pi 22k / 1024 = 22k * 2^21 in uint32 turns (wraps mod 2 pi) */
		*ang = fxp_atan2_u32(b, a) - (uint32_t)k * (22u << 21);
	}
	return la;
}

int32_t fxp_lsp_spec_bin(const int32_t *c_q30, int k, uint32_t *ang)
{
	int32_t c29[FXP_LPC_ORD];
	for (int i = 0; i < FXP_LPC_ORD; i++)
		c29[i] = lsp_c_q29(c_q30[i]);
	return lsp_spec_bin_q29(c29, k, ang);
}

static void taps_q(const int32_t *a, int n, int16_t *t, int *sh_out)
{
	uint32_t mx = 0, l1 = 0;
	for (int i = 0; i < n; i++)
	{
		uint32_t v = (uint32_t)(a[i] < 0 ? -a[i] : a[i]);
		mx = v > mx ? v : mx;
		l1 += v;
	}
	int s = fxp_rdft_taps_shift(mx, l1);
	for (;;) /* rounding may exceed the limits by a few LSB: drop a bit and retry */
	{
		uint32_t tl1 = 0;
		int ok = 1;
		for (int i = 0; i < n; i++)
		{
			int32_t v = s >= 0 ? shl(a[i], s) : shr_r(a[i], -s);
			if (v > 32767 || v < -32767)
				ok = 0;
			t[i] = (int16_t)(ok ? v : 0);
			tl1 += (uint32_t)(v < 0 ? -v : v);
		}
		if (ok && tl1 < (1u << 17))
			break;
		s--;
	}
	*sh_out = s;
}

void fxp_aks_to_mag2(const int16_t *lsp_hz2, int32_t log2E_q16, uint32_t W, int L, int32_t *A_log2_q16,
		     uint32_t *phA, int32_t *scratch)
{
	int32_t c[FXP_LPC_ORD], c29[FXP_LPC_ORD];
	for (int i = 0; i < FXP_LPC_ORD; i++)
	{
		c[i] = fxp_lsp_cos_q30(lsp_hz2[i]);
		c29[i] = lsp_c_q29(c[i]);
	}

	/* Scratch layout (514 words): first the DFT of A(z/g) as cq31_t X[257]; each X[k] is converted to
	 * log2|A(z/g)|^2 and written to word k (k <= 2k, X[k] already read, words < k only hold consumed
	 * bins). Then log2|A|^2 goes to words 256..511 (bins 128..256 of X, consumed), and log2 Pw replaces
	 * the A(z/g) values in words 0..255. */
	cq31_t *X = (cq31_t *)scratch;
	int32_t *pw = scratch, *la = scratch + FXP_SPEC_BINS;

	/* A(z/g), g = 1/2: taps ak_i 2^-i, sparse DFT */
	int32_t ak[FXP_LPC_ORD + 1], ag[FXP_LPC_ORD + 1];
	fxp_lsp_to_lpc_q19(c, ak);
	ag[0] = ak[0];
	for (int i = 1; i <= FXP_LPC_ORD; i++)
		ag[i] = shr_r(ak[i], i);
	int16_t t[FXP_LPC_ORD + 1];
	int tsh;
	taps_q(ag, FXP_LPC_ORD + 1, t, &tsh);
	fxp_rdft512_sparse(t, FXP_LPC_ORD + 1, X);
	const int xe = -14 - tsh - 19; /* A(z/g) = X * 2^xe */
	for (int k = 0; k < FXP_SPEC_BINS; k++)
	{
		int32_t re = X[k].r, im = X[k].i;
		uint32_t b = (uint32_t)(re ^ (re >> 31)) | (uint32_t)(im ^ (im >> 31));
		int s = b ? 17 - clz32(b) : 0; /* |re|, |im| -> <= 2^15 */
		if (s > 0)
		{
			re = shr_r(re, s);
			im = shr_r(im, s);
		}
		else
		{
			re = shl(re, -s);
			im = shl(im, -s);
		}
		uint32_t m2 = (uint32_t)(re * re) + (uint32_t)(im * im); /* <= 2^31 */
		scratch[k] = log2_plus_floor(m2, 2 * (s + xe));
	}

	/* log2 A2 per bin, A's phase at the harmonic bins, then log2 Pw = -la - (la - lag) / 5.
	 * Harmonic bin b_m = round(m Wo / (2 pi / 512)) = fxp_bin_round(m W), clamped to 255 (as sample_phase);
	 * m W <= L W <= 2^31 since L Wo <= pi. b_m is non-decreasing in m. */
	int32_t max_ib = INT32_MIN, max_pw = INT32_MIN;
	int mh = 1;
	for (int k = 0; k < FXP_SPEC_BINS; k++)
	{
		int32_t l;
		uint32_t bmh = mh <= L ? fxp_bin_round((uint32_t)mh * W) : 0;
		if (bmh > FXP_SPEC_BINS - 1)
			bmh = FXP_SPEC_BINS - 1;
		/* one call site, so the kernel inlines; the phase only at harmonic bins */
		uint32_t a;
		const int hit = phA && mh <= L && (int)bmh == k;
		l = lsp_spec_bin_q29(c29, k, hit ? &a : 0);
		if (hit)
		{
			do /* several harmonics can share a bin (clamping at 255) */
			{
				phA[mh++] = a;
				bmh = mh <= L ? fxp_bin_round((uint32_t)mh * W) : 0;
				if (bmh > FXP_SPEC_BINS - 1)
					bmh = FXP_SPEC_BINS - 1;
			} while (mh <= L && (int)bmh == k);
		}
		la[k] = l;
#ifdef FXP_LPCSPEC_DEBUG
		fxp_dbg_la[k] = l;
		fxp_dbg_lag[k] = pw[k];
#endif
		int32_t v = -l - (l - pw[k]) / 5;
		pw[k] = v;
		max_pw = v > max_pw ? v : max_pw;
		max_ib = -l > max_ib ? -l : max_ib;
	}
	const int32_t boost = 63625; /* log2(1.96) in Q16 */

	/* e_before = sum 2^-la, e_after = sum 2^pw (unboosted): terms Q22 relative to the max, 256 terms < 2^30 */
	uint32_t sb = 0, sa = 0;
	for (int k = 0; k < FXP_SPEC_BINS; k++)
	{
		sb += fxp_exp2_neg_q30_inl(-la[k] - max_ib) >> 8;
		sa += fxp_exp2_neg_q30_inl(pw[k] - max_pw) >> 8;
	}
	int32_t l_eb = fxp_log2_pf(sb, -22) + max_ib;
	int32_t l_ea = fxp_log2_pf(sa, -22) + max_pw;
	int32_t l_gain = log2E_q16 + l_eb - l_ea;
#ifdef FXP_LPCSPEC_DEBUG
	fxp_dbg_leb = l_eb;
	fxp_dbg_lea = l_ea;
#endif

	/* bands: sum of boosted Pw relative to the band maximum (keeps weak bands precise) */
	for (int m = 1; m <= L; m++)
	{
		/* band edges as the float code: (int)((m -+ 0.5) Wo / (2 pi / 512) + 0.5) */
		int a = (int)fxp_bin_round((uint32_t)m * W - (W >> 1));
		int bnd = (int)fxp_bin_round((uint32_t)m * W + (W >> 1));
		if (bnd > FXP_SPEC_BINS)
			bnd = FXP_SPEC_BINS;
		FXP_ASSERT(a >= 0 && bnd - a <= 64);
		if (bnd <= a)
		{
			A_log2_q16[m] = INT32_MIN / 2;
			continue;
		}
		int32_t mx = INT32_MIN;
		for (int k = a; k < bnd; k++)
		{
			int32_t vb = k < 64 ? pw[k] + boost : pw[k];
			mx = vb > mx ? vb : mx;
		}
		uint32_t s = 0;
		for (int k = a; k < bnd; k++)
		{
			int32_t vb = k < 64 ? pw[k] + boost : pw[k];
			s += fxp_exp2_neg_q30_inl(vb - mx) >> 6; /* Q24, <= 64 terms < 2^30 */
		}
		int32_t l_em = fxp_log2_pf(s, -24) + mx;
		A_log2_q16[m] = (l_gain + l_em) >> 1; /* log2 sqrt(gain * Em) */
	}
}
