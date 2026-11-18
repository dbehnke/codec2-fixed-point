/*
 * nlp.c — fixed-point NLP pitch estimator (encoder stage E2), see docs/10-fxp-nlp.md.
 *
 * Same algorithm as the float Codec2 NLP: square the speech, DC notch, 48-tap low-pass decimated by 5, 64-sample
 * Hann window, 512-point DFT, global peak over bins 16..128, then the sub-multiple search with pitch tracking.
 *
 * Number formats:
 *   sq[]      int32, the notch output y + 1 in units of 2^-k (|y| <= 2^30), k = nk the per-frame notch
 *             scale: 0 for loud frames, up to NLP_KMAX near silence (see nlp())
 *   sq_fir[]  int32 block-float: each frame's 16 low-pass outputs share an exponent fir_e = k + se (value =
 *             sq_fir 2^-fir_e), se <= 15 the largest that keeps them below 2^30. Exact accumulation (Q15
 *             taps), one rounding.
 *   x[]       int16 DFT input: sq_fir * w normalised per frame (max < 2^15, sum |x| < 2^17)
 *   P[]       uint32 power of the direct DFT at bins 15..128, normalised per frame
 * The decisions compare powers within one frame only, so no absolute scale is carried.
 *
 * Why a direct DFT and not fxp_rfft512: the decisions pick maxima of a smooth (8x zero-padded) spectrum and
 * compare neighbouring bins, so near-ties are common. With the 16-bit block-float FFT the bin differed from a
 * double-precision NLP on 1.7% of frames even for an exact input; the int16 input alone (double DFT) agreed on
 * 99.98-100%. The direct DFT accumulates exactly in int32 (as fxp_rdft.c) and only needs bins 15..128.
 */
#include "codec2_internal.h"
#include "nlp.h"
#include "fxp.h"
#include "fxp_fft.h"
#include <string.h>

#define START_POS (M_PITCH - N_SAMP)     /* first new sample in Sn */
#define BIN_MIN (PE_FFT_SIZE * DEC / P_MAX) /* 16 */
#define BIN_MAX (PE_FFT_SIZE * DEC / P_MIN) /* 128 */
#define BIN_LO (BIN_MIN - 1)                 /* lowest bin read: left neighbour in the sub-multiple search */
#define NBINS (BIN_MAX - BIN_LO + 1)         /* 114 */
#define NOUT (N_SAMP / DEC)                  /* 16 low-pass outputs per frame */
#ifndef NLP_KMAX
#define NLP_KMAX 28 /* largest notch scale exponent */
#endif

_Static_assert(N_SAMP % DEC == 0, "NLP decimation");
_Static_assert(BIN_MAX < PE_FFT_SIZE / 2, "NLP search range");
_Static_assert(NDEC == NLP_NBLK * NOUT, "NLP block exponents");
_Static_assert(40961LL * 45096 < (1LL << 31) && NLP_KMAX + 15 <= 127, "NLP FIR bounds");
_Static_assert(NLP_SCRATCH_WORDS >= NDEC / 2 + 2 * NBINS && NLP_SCRATCH_WORDS >= NDEC + NDEC / 2, "NLP scratch");

/* floor(3 g / 10) and floor(3 g / 20) without overflow: the float thresholds CNLP gmax and CNLP gmax / 2 as
 * "P > thresh" on integers (P > 0.3 g <=> 10 P > 3 g <=> P > floor(3 g / 10)) */
static inline uint32_t frac3_10(uint32_t g) { return g / 10 * 3 + g % 10 * 3 / 10; }
static inline uint32_t frac3_20(uint32_t g) { return g / 20 * 3 + g % 20 * 3 / 20; }

static inline uint32_t uabs32(int32_t v) { return v < 0 ? 0u - (uint32_t)v : (uint32_t)v; }
static inline int bitlen(uint32_t v) { return v ? 32 - clz32(v) : 0; }

/* v 2^sh for any sign of sh, rounding half up on right shifts; the caller guarantees the result fits */
static inline int32_t shift_r(int32_t v, int sh)
{
	if (sh >= 0)
		return shl(v, sh);
	int n = sh < -32 ? 32 : -sh;
	return ((v >> (n - 1)) + 1) >> 1;
}

/* P[k - BIN_LO] is the power at bin k */
static int post_process_sub_multiples(const uint32_t *P, uint32_t gmax, int gmax_bin, int prev_bin)
{
	int cmax_bin = gmax_bin;
	uint32_t th = frac3_10(gmax), th_track = frac3_20(gmax);

	for (int mult = 2; gmax_bin / mult >= BIN_MIN; mult++)
	{
		int b = gmax_bin / mult; /* search interval [0.8 b, 1.2 b]: the float (int)(0.8f b), (int)(1.2f b) */
		int bmin = 4 * b / 5, bmax = 6 * b / 5;
		if (bmin < BIN_MIN)
			bmin = BIN_MIN;

		/* lower threshold to favour the previous frame's pitch (pitch tracking) */
		uint32_t thresh = prev_bin > bmin && prev_bin < bmax ? th_track : th;

		uint32_t lmax = 0;
		int lmax_bin = bmin;
		for (int k = bmin; k <= bmax; k++)
			if (P[k - BIN_LO] > lmax)
			{
				lmax = P[k - BIN_LO];
				lmax_bin = k;
			}

		/* BIN_LO <= lmax_bin - 1 and lmax_bin + 1 <= 6 * 64 / 5 + 1 = 77 <= BIN_MAX */
		if (lmax > thresh && lmax > P[lmax_bin - 1 - BIN_LO] && lmax > P[lmax_bin + 1 - BIN_LO])
			cmax_bin = lmax_bin;
	}
	return cmax_bin;
}

/* 48-tap low-pass, decimated by 5, of nlp->sq: appends one block of 16 outputs to sq_fir with its exponent.
 * frac: 16 int32 of scratch. */
static void lowpass(nlp_t *restrict nlp, int32_t *restrict frac)
{
	const int32_t *restrict sq = nlp->sq;
	int32_t *restrict sq_fir = nlp->sq_fir;

	/* out[i] = sum_k h[k] sq[47 + 5 i - k] (the float sq[240 + 5 i - k]).
	 * Exact: v = hi 2^15 + lo with lo in [0, 2^15); both partial sums fit int32 since sum |h| = 45096 (Q15),
	 * |v| <= 2^30 + 2^NLP_KMAX + 8 so |hi| <= 40961, and lo < 2^15. sum h v / 2^15 = ah + al / 2^15 = ah' + al' / 2^15 with al' in [0, 2^15). */
	memmove(sq_fir, sq_fir + NOUT, (NDEC - NOUT) * sizeof(int32_t));
	memmove(nlp->fir_e, nlp->fir_e + 1, NLP_NBLK - 1);
	uint32_t ora = 0;
	for (int i = 0; i < NOUT; i++)
	{
		const int32_t *v = &sq[NLP_NTAP - 1 + DEC * i];
		/* linear phase, h[k] = h[47 - k] (checked in test_nlp): the same two exact sums, taps paired */
		int32_t ah = 0, al = 0;
		for (int k = 0; k < NLP_NTAP / 2; k++)
		{
			int32_t h = fxp_nlp_fir_q15[k];
			int32_t p = v[-k], q = v[k - (NLP_NTAP - 1)];
			ah += ((p >> 15) + (q >> 15)) * h;
			al += ((p & 0x7fff) + (q & 0x7fff)) * h;
		}
		ah += al >> 15;
		sq_fir[NDEC - NOUT + i] = ah;
		frac[i] = al & 0x7fff;
		ora |= uabs32(ah);
	}
	/* block exponent: the largest se <= 15 with (|ah| + 1) 2^se <= 2^30 */
	int se = 30 - bitlen(ora);
	se = se < 0 ? 0 : se > 15 ? 15 : se;
	for (int i = 0; i < NOUT; i++)
	{
		int32_t *o = &sq_fir[NDEC - NOUT + i];
		*o = shl(*o, se) + ((frac[i] + ((1 << (15 - se)) >> 1)) >> (15 - se));
	}
	nlp->fir_e[NLP_NBLK - 1] = (int8_t)(se + nlp->nk); /* value = sq_fir 2^-(se + k) */
}

/* Returns the pitch as the NLP DFT bin: f0 = 3.125 bin Hz, pitch period 2560 / bin samples. *prev_bin is the
 * pitch tracking memory. scratch: NLP_SCRATCH_WORDS int32, used only during the call. */
int nlp(nlp_t *restrict nlp, const int16_t *restrict Sn, int *restrict prev_bin, int32_t *restrict scratch)
{
	int32_t *restrict sq = nlp->sq;
	const int32_t *restrict sq_fir = nlp->sq_fir;

	/* Square and DC notch y[n] = x[n] - x[n-1] + 0.95 y[n-1], x = Sn^2 2^k, in units of 2^-k.
	 * Scale: with z = y - x, z[n] = 0.95 z[n-1] - 0.05 x[n-1], so for 0 <= x <= X, |z| <= max(|z[0]|, X) and
	 * |y| <= max(|z[0]|, X) + X. k (0..NLP_KMAX) is the largest with max(|z[0]|, X) 2^k <= 2^29 (X includes
	 * mem_x), giving |y| <= 2^30. Loud frames run at k = 0; near-silent ones keep up to 28 fractional bits,
	 * where integer rounding would otherwise decide the pitch.
	 * 0.95 y is formed as floor((19 y + e) / 20) with the remainder e carried to the next sample (error
	 * feedback): no dead band, and the error stays below one unit. y = 20 q + r, 0 <= r < 20, keeps 19 y
	 * from overflowing: (19 y + e) / 20 = 19 q + (19 r + e) / 20. */
	int32_t mx = nlp->mem_x, my = nlp->mem_y, e = nlp->mem_e;
	int kp = nlp->nk;
	{
		uint32_t X = 0;
		for (int i = 0; i < N_SAMP; i++)
		{
			uint32_t x = (uint32_t)(Sn[START_POS + i] * Sn[START_POS + i]);
			if (x > X)
				X = x;
		}
		/* state in real units, rounded up */
		uint32_t mr = ((uint32_t)mx + ((1u << kp) - 1)) >> kp;
		uint32_t zr = (uabs32(my - mx) + ((1u << kp) - 1)) >> kp;
		if (mr > X)
			X = mr;
		if (zr > X)
			X = zr;
		/* the FIR history (sq[0..46], previous frame) is rescaled too: keep it within 2^30 as well */
		uint32_t oh = 0;
		for (int j = 0; j < NLP_NTAP - 1; j++)
			oh |= uabs32(sq[j]);
		uint32_t hr = (oh + ((1u << kp) - 1)) >> kp;
		if (hr > X)
			X = hr;
		int k = 29 - bitlen(X);
		k = k < 0 ? 0 : k > NLP_KMAX ? NLP_KMAX : k;
		if (k != kp)
		{
			/* rescale the notch state and the FIR history (the remainder e is dropped: < 1 old unit) */
			int d = k - kp;
			mx = shift_r(mx, d);
			my = shift_r(my, d);
			e = 0;
			for (int j = 0; j < NLP_NTAP - 1; j++)
			{
				FXP_ASSERT(d <= 0 || uabs32(sq[j]) <= (0x7fffffffu >> d));
				sq[j] = shift_r(sq[j], d);
			}
			nlp->nk = kp = k;
		}
	}
	for (int i = 0; i < N_SAMP; i++)
	{
		int32_t s = Sn[START_POS + i];
		int32_t x = shl(s * s, kp);
		int32_t q = my / 20, r = my - 20 * q;
		if (r < 0)
		{
			r += 20;
			q--;
		}
		int32_t t = 19 * r + e;
		int32_t y = x - mx + (19 * q + t / 20);
		e = t % 20;
		mx = x;
		my = y;
		FXP_ASSERT(y >= -(1 << 30) - 8 && y <= (1 << 30) + 8);
		sq[NLP_NTAP - 1 + i] = y + (1 << kp); /* the float code's + 1.0 (keeps the silent-input spectrum non-zero) */
	}
	nlp->mem_x = mx;
	nlp->mem_y = my;
	nlp->mem_e = e;

	lowpass(nlp, scratch);
	memmove(sq, sq + N_SAMP, (NLP_NTAP - 1) * sizeof(int32_t));

	/* Window and normalise. u = sq_fir aligned across the four blocks and scaled so max |u| < 2^31:
	 * block b holds values below 2^(L_b - e_b) (L_b = bit length of its largest magnitude), G = max of those,
	 * u = sq_fir 2^(31 - G - e_b). y = u w (Q15). x = round(y / 2^s) with the smallest s giving max |x| <= 32767
	 * and sum |x| < 2^17 (the exact-int32 DFT precondition, as fxp_rdft.c). */
	int Lb[NLP_NBLK], G = -64;
	for (int b = 0; b < NLP_NBLK; b++)
	{
		uint32_t o = 0;
		for (int i = 0; i < NOUT; i++)
			o |= uabs32(sq_fir[b * NOUT + i]);
		Lb[b] = o ? bitlen(o) - nlp->fir_e[b] : -64;
		if (Lb[b] > G)
			G = Lb[b];
	}
	int32_t *yw = scratch + NDEC / 2; /* y[64], after the int16 x[64] */
	uint32_t ymax = 0, yl1 = 0;
	for (int i = 0; i < NDEC; i++)
	{
		int b = i / NOUT;
		int32_t u = G > -64 ? shift_r(sq_fir[i], 31 - G - nlp->fir_e[b]) : 0;
		int32_t y = mul32x16_q15(u, fxp_nlp_hann_q15[i]);
		uint32_t a = uabs32(y);
		yw[i] = y;
		if (a > ymax)
			ymax = a;
		yl1 += a >> 6;
	}
	/* max: round(ymax / 2^s) <= 32767. l1: sum round(|y| / 2^s) <= sum |y| / 2^s + 32 < 64 (yl1 + 64) / 2^s + 32,
	 * which is < 2^17 when yl1 + 64 < 2^(10 + s). */
	int s = 0;
	if (ymax > 32767)
	{
		s = 17 - clz32(ymax);
		if (((ymax >> (s - 1)) + 1) >> 1 > 32767)
			s++;
	}
	int s2 = bitlen(yl1 + 64) - 10;
	if (s2 > s)
		s = s2;
	int16_t *x = (int16_t *)scratch;
	for (int i = 0; i < NDEC; i++)
		x[i] = (int16_t)(s ? ((yw[i] >> (s - 1)) + 1) >> 1 : yw[i]);
#ifdef FXP_CHECK
	{
		uint32_t l1 = 0;
		for (int i = 0; i < NDEC; i++)
			l1 += uabs32(x[i]);
		FXP_ASSERT(l1 < (1u << 17));
	}
#endif

	/* direct DFT of the 64 windowed samples (zero-padded 512 grid), bins BIN_LO..BIN_MAX:
	 * re[k] = sum_i x[i] C(i k), im[k] = sum_i x[i] S(i k), C, S the Q14 table cos and sin of 2 pi m / 512.
	 * x[0] = x[63] = 0 (w = 0). Every partial sum is a signed subset of terms with sum |x C| <= sum |x| 2^14 <
	 * 2^31, so all sums are exact int32 and any order of summation gives the same integers.
	 * Bin pairing (k, 128 - k): i (128 - k) = 128 i - i k, and 128 i is a multiple of pi/2, so by the table's exact
	 * symmetries (C(-m) = C(m), C(m + 256) = -C(m), S(m) = C(m - 128); checked in test_nlp) the twiddles of bin
	 * 128 - k are those of bin k with signs and roles set by i mod 4:
	 *   i = 0: ( C, -S)   i = 1: ( S,  C)   i = 2: (-C,  S)   i = 3: (-S, -C)    (mod 4)
	 * so the pair shares every table read and product. Bins 15..63 pair with 113..65; 64 and 114..128 are single. */
	int32_t *re = scratch + NDEC / 2, *im = re + NBINS;
	uint32_t orx = 0;
	const uint32_t *cs = fxp_nlp_cs512;
#define CS_TERM(j)                                         \
	w = cs[m];                                             \
	a = x[i + (j)] * (int32_t)(int16_t)w;                  \
	b = x[i + (j)] * (int32_t)(int16_t)(w >> 16);          \
	m = (m + k) & 511;
	for (int k = BIN_LO; k < 64; k++)
	{
		int32_t c = 0, s = 0, c2 = 0, s2 = 0, a, b;
		uint32_t w;
		unsigned m = 0;
		for (int i = 0; i < NDEC; i += 4)
		{
			CS_TERM(0) c += a; s += b; c2 += a; s2 -= b;
			CS_TERM(1) c += a; s += b; c2 += b; s2 += a;
			CS_TERM(2) c += a; s += b; c2 -= a; s2 += b;
			CS_TERM(3) c += a; s += b; c2 -= b; s2 -= a;
		}
		re[k - BIN_LO] = c;
		im[k - BIN_LO] = s;
		re[128 - k - BIN_LO] = c2;
		im[128 - k - BIN_LO] = s2;
		orx |= uabs32(c) | uabs32(s) | uabs32(c2) | uabs32(s2);
	}
	for (int k = 64; k <= BIN_MAX; k = k == 64 ? 128 - BIN_LO + 1 : k + 1)
	{
		int32_t c = 0, s = 0, a, b;
		uint32_t w;
		unsigned m = 0;
		for (int i = 0; i < NDEC; i += 4)
		{
			CS_TERM(0) c += a; s += b;
			CS_TERM(1) c += a; s += b;
			CS_TERM(2) c += a; s += b;
			CS_TERM(3) c += a; s += b;
		}
		re[k - BIN_LO] = c;
		im[k - BIN_LO] = s;
		orx |= uabs32(c) | uabs32(s);
	}
#undef CS_TERM
	/* power, normalised so |re|, |im| < 2^15 (P < 2^31); written over re */
	int sp = bitlen(orx) - 15;
	uint32_t *P = (uint32_t *)re;
	for (int j = 0; j < NBINS; j++)
	{
		int32_t r = re[j], i = im[j];
		if (sp > 0)
		{
			r = (r + (1 << (sp - 1))) >> sp;
			i = (i + (1 << (sp - 1))) >> sp;
		}
		P[j] = (uint32_t)(r * r) + (uint32_t)(i * i);
	}

	/* global peak */
	uint32_t gmax = 0;
	int gmax_bin = BIN_MIN;
	for (int k = BIN_MIN; k <= BIN_MAX; k++)
		if (P[k - BIN_LO] > gmax)
		{
			gmax = P[k - BIN_LO];
			gmax_bin = k;
		}

	int best_bin = post_process_sub_multiples(P, gmax, gmax_bin, *prev_bin);
	*prev_bin = best_bin;
	return best_bin;
}

void nlp_init(nlp_t *nlp)
{
	memset(nlp, 0, sizeof(*nlp));
}
