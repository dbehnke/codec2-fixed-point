/*
 * analysis.c — fixed-point Codec2 encoder analysis front end (stage E1), see docs/09-fxp-analysis.md.
 *
 * Replaces the float dft_speech / two_stage_pitch_refinement / estimate_amplitudes / est_voicing_mbe.
 * The NLP pitch estimator (nlp.c) is integer since stage E2 (docs/10), the LPC/LSP analysis (lpc.c) since E3 (docs/11).
 *
 * Number formats:
 *   Sn[]      int16 input history (the float code's Sn held the same integers)
 *   Sw[]      cq15 spectrum of the Hann-windowed speech, block exponent: the float encoder's
 *             Sw_float = Sw * 2^g * scale, g = e_fft + s - 15, scale = 1 / sqrt(512 sum h^2)
 *   pw[]      uint32 power spectrum, pw = round(|Sw|^2 / 256), mirrored to 512 bins
 *   pitch     exact rational p = pn / pd samples (Wo = 2 pi pd / pn)
 *   energies  64-bit sums in |Sw|^2 units; compared as log2 Q16 with the float code's 1e-4 floor
 *             carried into these units (fxp_ana_floor_q24 * 2^(-24 - 2g))
 */
#include "codec2_internal.h"
#include "analysis.h"
#include "nlp.h"
#include "fxp.h"
#include "fxp_fft.h"
#include "fxp_math.h"
#include <string.h>

/* log2 thresholds in Q16 */
#define LG_V_THRESH 130624 /* log2(10^0.6), V_THRESH = 6 dB */
#define LG_10DB 217706     /* log2(10) */
#define LG_4DB 87082       /* -log2(0.39810717) */

/* layout in the encoder scratch: Sw (257 cq15) then pw (512 uint32) */
#define PW_OFFSET (FXP_RFFT_BUF * sizeof(cq15_t))
_Static_assert(PW_OFFSET % 4 == 0, "pw alignment");
_Static_assert(PW_OFFSET + FFT_ENC * sizeof(uint32_t) <= sizeof(((codec2_encoder_t *)0)->scratch),
	       "encoder scratch too small for Sw + pw");

/* Windowed, centred 512-point real FFT of Sn into Sw (bins 0..256). x32 is int32[512] scratch (it may alias
 * pw: the products are consumed before the power spectrum is written). Returns g: Sw_float = Sw 2^g scale. */
static int dft_speech(cq15_t *Sw, int32_t *x32, const int16_t *Sn)
{
	static const int nw2 = NW / 2;
	uint32_t mx = 0;

	/* x[i] = Sn[160 + i] h[139 + i], x[373 + i] = Sn[21 + i] h[i] (i < 139), zero elsewhere.
	 * |Sn h| <= 2^15 * 2^15 = 2^30. */
	for (int i = 0; i < nw2; i++)
	{
		int32_t a = (int32_t)Sn[M_PITCH / 2 + i] * (int32_t)fxp_hann_q15[nw2 + i];
		int32_t b = (int32_t)Sn[M_PITCH / 2 - nw2 + i] * (int32_t)fxp_hann_q15[i];
		x32[i] = a;
		x32[nw2 + i] = b;
		uint32_t ua = a < 0 ? 0u - (uint32_t)a : (uint32_t)a;
		uint32_t ub = b < 0 ? 0u - (uint32_t)b : (uint32_t)b;
		if (ua > mx)
			mx = ua;
		if (ub > mx)
			mx = ub;
	}

	/* smallest right shift s >= 0 with round(max / 2^s) <= 32767 */
	int s = 0;
	if (mx > 32767)
	{
		s = 17 - clz32(mx); /* mx >> s < 2^15 */
		if (((mx >> (s - 1)) + 1) >> 1 > 32767)
			s++;
	}

	int16_t *x = (int16_t *)Sw; /* buf[n].r = x[2n], buf[n].i = x[2n + 1] */
	memset(x, 0, FXP_RFFT_BUF * sizeof(cq15_t));
	for (int i = 0; i < nw2; i++)
	{
		int32_t a = x32[i], b = x32[nw2 + i];
		if (s)
		{
			a = ((a >> (s - 1)) + 1) >> 1; /* round half up */
			b = ((b >> (s - 1)) + 1) >> 1;
		}
		x[i] = (int16_t)a;
		x[FFT_ENC - nw2 + i] = (int16_t)b;
	}
	return fxp_rfft512(Sw) + s - 15;
}

/* pw[k] = round(|Sw[k]|^2 / 256) for k = 0..256, mirrored: pw[512 - k] = pw[k]. |Sw|^2 <= 2^31. The float
 * pitch refinement deliberately reads past Nyquist (up to bin ~341). */
static void power_spectrum(const cq15_t *Sw, uint32_t *pw)
{
	for (int k = 0; k <= FFT_ENC / 2; k++)
	{
		int32_t r = Sw[k].r, i = Sw[k].i;
		pw[k] = ((uint32_t)(r * r) + (uint32_t)(i * i) + 128u) >> 8;
	}
	for (int k = 1; k < FFT_ENC / 2; k++)
		pw[FFT_ENC - k] = pw[k];
}

/* Harmonic sum search over the candidates p = (n0 + i dn) / d, i = 0..cnt-1, as hs_pitch_refinement:
 * E(p) = sum_{m=1..L} pw[round(m 512 / p)], ties up. The bin b_m = floor((1024 d m + n) / (2 n)) is stepped
 * exactly with a quotient/remainder pair. The pitch with the largest E (first on ties) replaces (*bn, *bd);
 * if every E is 0 the input pitch is kept. E <= 82 * 2^23 < 2^30. */
static void hs_pitch_refinement(const uint32_t *pw, uint32_t n0, uint32_t dn, uint32_t d, int cnt, int L,
				uint32_t *bn, uint32_t *bd)
{
	uint32_t Em = 0;
	for (int i = 0; i < cnt; i++)
	{
		uint32_t n = n0 + (uint32_t)i * dn;
		uint32_t n2 = 2 * n, num = 1024 * d;
		uint32_t qi = num / n2, ri = num - qi * n2;
		uint32_t q = 0, r = n, E = 0;
		for (int m = 1; m <= L; m++)
		{
			q += qi;
			r += ri;
			if (r >= n2)
			{
				r -= n2;
				q++;
			}
			FXP_ASSERT(q < FFT_ENC);
			E += pw[q];
		}
		if (E > Em)
		{
			Em = E;
			*bn = n;
			*bd = d;
		}
	}
}

/* two_stage_pitch_refinement from the NLP bin (pitch 2560 / bin): coarse p0 - 5 .. p0 + 5 step 1, then fine
 * pc - 1 .. pc + 1 step 1/4; limit to [P_MIN, P_MAX]; L = floor(p / 2), less one if Wo L >= 0.95 pi. */
static void two_stage_pitch_refinement(const uint32_t *pw, int bin, enc_frame_t *f)
{
	FXP_ASSERT(bin >= PE_FFT_SIZE * DEC / P_MAX && bin <= PE_FFT_SIZE * DEC / P_MIN);
	const uint32_t b = (uint32_t)bin;
	uint32_t n = PE_FFT_SIZE * DEC, d = b; /* p0 = 2560 / bin */

	hs_pitch_refinement(pw, n - 5 * b, b, b, 11, (int)(n / (2 * b)), &n, &d);
	const uint32_t nc = n; /* d == b here */
	hs_pitch_refinement(pw, 4 * nc - 4 * b, b, 4 * b, 9, (int)(nc / (2 * b)), &n, &d);

	if (n > P_MAX * d)
	{
		n = P_MAX;
		d = 1;
	}
	else if (n < P_MIN * d)
	{
		n = P_MIN;
		d = 1;
	}
	int L = (int)(n / (2 * d));
	if (40 * (uint32_t)L * d >= 19 * n) /* 2 L >= 0.95 p */
		L--;
	f->pn = n;
	f->pd = d;
	f->L = L;
}

/* log2(v 2^ev + F) in Q16, F = fxp_ana_floor_q24 2^ef. Both terms as 31-bit mantissas, aligned and added. */
static int32_t log2_plus_floor(uint64_t v, int ev, int ef)
{
	uint32_t mf = fxp_ana_floor_q24;
	int nf = clz32(mf) - 1;
	mf <<= nf; /* [2^30, 2^31) */
	ef -= nf;
	if (v == 0)
		return fxp_log2_pf(mf, ef);
	int bits = 64 - __builtin_clzll(v);
	uint32_t mv;
	if (bits > 31)
	{
		mv = (uint32_t)(v >> (bits - 31));
		ev += bits - 31;
	}
	else
	{
		mv = (uint32_t)v << (31 - bits);
		ev -= 31 - bits;
	}
	if (ev >= ef)
	{
		int dd = ev - ef;
		if (dd < 32)
			mv += mf >> dd; /* < 2^32 */
		return fxp_log2_pf(mv, ev);
	}
	int dd = ef - ev;
	if (dd < 32)
		mf += mv >> dd;
	return fxp_log2_pf(mf, ef);
}

/* Residual of one spectrum component against the fit a w 2^-sh (round half up), fit clamped to +-2^24. */
static inline int32_t fit_residual(int32_t x, int32_t a, int32_t w, int sh)
{
	int32_t p = a * w; /* |a| < 2^16, |w| < 2^15 */
	int32_t fit;
	FXP_ASSERT(sh >= -31);
	if (sh > 0) /* sh >= 33: |p| < 2^31 so round(p 2^-sh) = 0, which the 31-bit shift also gives */
		fit = ((p >> (sh > 32 ? 31 : sh - 1)) + 1) >> 1;
	else
		fit = wsat(p, -sh > 31 ? 31 : -sh);
	if (fit > (1 << 24))
		fit = 1 << 24;
	else if (fit < -(1 << 24))
		fit = -(1 << 24);
	return x - fit;
}

static inline uint64_t sq64(int32_t dr)
{
	uint32_t a = dr < 0 ? 0u - (uint32_t)dr : (uint32_t)dr;
	if (a < 65536u)
		return a * a;
	/* rare (|residual| >= 2^16, never observed on the corpus or the review's probes): three 32-bit products, a <= 2^31 */
	uint32_t ah = a >> 16, al = a & 0xffffu;
	return ((uint64_t)(ah * ah) << 32) + ((uint64_t)(2u * ah * al) << 16) + al * al;
}

/* Amplitude a of x / den as a 16-bit mantissa: x / den = a 2^-k (in the units of d16 = den 2^-sd). */
static inline int32_t am_mantissa(int32_t x, uint32_t d16, int *k)
{
	if (x == 0)
	{
		*k = 0;
		return 0;
	}
	uint32_t ax = x < 0 ? 0u - (uint32_t)x : (uint32_t)x;
	*k = clz32(ax) - 1;                 /* |x| << k in [2^30, 2^31) */
	return shl(x, *k) / (int32_t)d16; /* |.| < 2^31 / 2^15 = 2^16 */
}

/* estimate_amplitudes + est_voicing_mbe. Band powers A[m]^2 feed only the energy ratios (3200 mode does not
 * transmit A); the voicing fit uses the window spectrum fxp_winspec_q15 (any scale: only W ratios matter). */
static int est_voicing_mbe(const cq15_t *Sw, const uint32_t *pw, uint32_t n, uint32_t d, int L, int g)
{
	const int l_1000hz = L / 4, l_2000hz = L / 2, l_4000hz = L;
	const uint32_t n2 = 2 * n, step = 1024 * d;
	const uint32_t qi = step / n2, ri = step - qi * n2;

	/* band edges floor(((2m + 1) 512 d + n) / (2 n)) = round((m + 0.5) Wo_bin), stepped exactly */
	uint32_t q = (512 * d + n) / n2, r = (512 * d + n) - q * n2;
	uint64_t sig = 0, elow = 0, ehigh = 0; /* pw units */
	for (int m = 1; m <= l_4000hz; m++)
	{
		uint32_t am = q;
		q += qi;
		r += ri;
		if (r >= n2)
		{
			r -= n2;
			q++;
		}
		uint32_t bm = q > FFT_ENC / 2 ? FFT_ENC / 2 : q;
		uint32_t t = 0; /* <= 27 * 2^23 */
		for (uint32_t i = am; i < bm; i++)
			t += pw[i];
		if (m <= l_2000hz)
		{
			elow += t;
			if (m <= l_1000hz)
				sig += t;
		}
		if (m >= l_2000hz)
			ehigh += t;
	}

	uint64_t err = 0; /* |Sw|^2 units */
	for (int l = 1; l <= l_1000hz; l++)
	{
		/* al = ceil((l - 0.5) Wo_bin), bl = ceil((l + 0.5) Wo_bin), Wo_bin = 512 d / n;
		 * offset = floor(256.5 - l Wo_bin) centres the window spectrum on harmonic l */
		int32_t al = (int32_t)(((uint32_t)(2 * l - 1) * 256 * d + n - 1) / n);
		int32_t bl = (int32_t)(((uint32_t)(2 * l + 1) * 256 * d + n - 1) / n);
		FXP_ASSERT(513 * n > 1024 * (uint32_t)l * d);
		int32_t off = (int32_t)((513 * n - 1024 * (uint32_t)l * d) / n2);
		if (off < -al)
			off = -al;
		if (off + bl >= FFT_ENC)
			off = FFT_ENC - bl - 1;
		FXP_ASSERT(al >= 0 && bl <= FFT_ENC / 2 && bl - al <= 32);

		/* least-squares amplitude Am = sum Sw W / sum W^2; products >> 5 keep 32 terms in int32 */
		int32_t nr = 0, ni = 0, den = 0;
		for (int m = al; m < bl; m++)
		{
			int idx = off + m - FFT_ENC / 2;
			int32_t w = fxp_winspec_q15[idx < 0 ? -idx : idx];
			nr += (Sw[m].r * w) >> 5;
			ni += (Sw[m].i * w) >> 5;
			den += (w * w) >> 5;
		}
		if (den <= 0)
			continue;

		/* den = d16 2^sd with d16 in [2^15, 2^16); Am = a 2^-(k + sd) */
		int sd = 16 - clz32((uint32_t)den);
		uint32_t d16 = sd >= 0 ? (uint32_t)den >> sd : (uint32_t)den << -sd;
		int kr, ki;
		int32_t ar = am_mantissa(nr, d16, &kr);
		int32_t ai = am_mantissa(ni, d16, &ki);

		for (int m = al; m < bl; m++)
		{
			int idx = off + m - FFT_ENC / 2;
			int32_t w = fxp_winspec_q15[idx < 0 ? -idx : idx];
			err += sq64(fit_residual(Sw[m].r, ar, w, kr + sd));
			err += sq64(fit_residual(Sw[m].i, ai, w, ki + sd));
		}
	}

	/* compare 1e-4 + sum in log2: the floor 1e-4 (float units) is fxp_ana_floor_q24 2^(-24 - 2 g) in |Sw|^2 */
	const int ef = -24 - 2 * g;
	const int32_t lsig = log2_plus_floor(sig, 8, ef);
	const int32_t lelow = log2_plus_floor(elow, 8, ef);
	const int32_t lehigh = log2_plus_floor(ehigh, 8, ef);
	const int32_t lerr = log2_plus_floor(err, 0, ef);

	int voiced = lsig > lerr + LG_V_THRESH;
	/* Type 1 errors: strongly voiced speech declared unvoiced */
	if (!voiced && lelow > lehigh + LG_10DB)
		voiced = 1;
	/* Type 2 errors: strongly unvoiced speech declared voiced; low (<= 60 Hz, p >= 400 / 3) pitch on noise */
	if (voiced)
	{
		if (lelow < lehigh - LG_10DB)
			voiced = 0;
		if (lelow < lehigh - LG_4DB && 3 * n >= 400 * d)
			voiced = 0;
	}
	return voiced;
}

void analyse_one_frame(codec2_encoder_t *c2, enc_frame_t *f, const int16_t *speech)
{
	cq15_t *Sw = (cq15_t *)c2->scratch;
	uint32_t *pw = (uint32_t *)((uint8_t *)c2->scratch + PW_OFFSET);

	memmove(c2->Sn, &c2->Sn[N_SAMP], (M_PITCH - N_SAMP) * sizeof(int16_t));
	memcpy(&c2->Sn[M_PITCH - N_SAMP], speech, N_SAMP * sizeof(int16_t));

	/* the NLP uses the scratch first; Sw and pw then overwrite it */
	int bin = nlp(&c2->nlp, c2->Sn, &c2->prev_bin_enc, c2->scratch);

	int g = dft_speech(Sw, (int32_t *)pw, c2->Sn);
	power_spectrum(Sw, pw);
	two_stage_pitch_refinement(pw, bin, f);
	f->voiced = est_voicing_mbe(Sw, pw, f->pn, f->pd, f->L, g);
}

int encode_Wo_p(uint32_t pn, uint32_t pd)
{
	/* index = floor(128 (Wo - W0_MIN) / (W0_MAX - W0_MIN) + 0.5) = floor((256 (160 pd - pn) + 7 pn) / (14 pn)) */
	if (pn >= P_MAX * pd)
		return 0;
	int idx = (int)((256 * (P_MAX * pd - pn) + 7 * pn) / (14 * pn));
	return idx > (1 << WO_BITS) - 1 ? (1 << WO_BITS) - 1 : idx;
}
