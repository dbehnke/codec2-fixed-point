/*
 * synthesis.c — fixed-point Codec2 sinusoidal synthesis (decoder), see docs/06-fxp-synthesis.md.
 *
 * Replaces the float sample_phase / phase_synth_zero_order / postfilter / synthesise / ear_protection.
 * Number formats:
 *   angles      uint32, 2^32 = 2 pi (wrap-around is the mod 2 pi of the float code)
 *   lA[m]       log2 A[m] in Q16 (INT32_MIN/2 = empty band, A = 0)
 *   bg_est      background level in dB, Q16
 *   Sn_[]       synthesised signal, int32 Q4 (sample value * 16)
 *   Pn[]        Parzen window, int16 Q14 (1.0 = 16384)
 *   Sw[] / sw[] irfft buffer (cq15 in the decoder scratch), block exponent from fxp_irfft512
 */
#include "codec2_internal.h"
#include "synthesis.h"
#include "util.h"
#include "fxp.h"
#include "fxp_fft.h"
#include "fxp_math.h"
#include "fxp_lpcspec.h"
#include <string.h>

#define LA_EMPTY_LIM (INT32_MIN / 4)    /* lA at or below this: empty band (A = 0) */
#define Q16(x) ((int32_t)((x) * 65536.0 + ((x) < 0 ? -0.5 : 0.5)))
#define RAND_TO_ANGLE 131076u           /* 2^32 / CODEC2_RAND_MAX = 131076.00003 */
#define SN_Q 4                          /* Sn_ fraction bits */

static void make_synthesis_window(int16_t *Pn)
{
	/* Parzen window as the float code with n0 - TW = 0 and n1 + TW = 2 N_SAMP (TW = N_SAMP / 2):
	 * a linear ramp up over 2 TW samples then down, Pn[i] = i / (2 TW) and (2 N_SAMP - i) / (2 TW). */
	_Static_assert(TW == N_SAMP / 2, "window formula assumes TW = N_SAMP / 2");
	for (int i = 0; i < 2 * N_SAMP; i++)
	{
		int n = i < N_SAMP ? i : 2 * N_SAMP - i;
		Pn[i] = (int16_t)((n * 16384 * 2 + 2 * TW) / (4 * TW)); /* round(n / (2 TW) * 2^14) */
	}
}

/* 10 log10(sum A^2 / L) in Q16 dB. Returns -120 dB (the float code's 1e-12 floor) if all A = 0. */
static int32_t frame_energy_db(const int32_t *lA, int L)
{
	int32_t lmax = INT32_MIN;
	for (int m = 1; m <= L; m++)
		if (lA[m] > lmax)
			lmax = lA[m];
	if (lmax <= LA_EMPTY_LIM)
		return Q16(-120.0);

	/* sum 2^(2 (lA - lmax)) in Q22: each term <= 2^22, L <= 80 -> < 2^29 */
	uint32_t s = 0;
	for (int m = 1; m <= L; m++)
		if (lA[m] > LA_EMPTY_LIM && lA[m] - lmax > -(32 << 15))
			s += fxp_exp2_neg_q30(2 * (lA[m] - lmax)) >> 8;
	int32_t l2 = fxp_log2_pf(s, -22) + 2 * lmax - fxp_log2_q16((uint32_t)L); /* log2(sum / L) */
	/* 10 log10(2) = 3 + 0.0103: |l2| < 2^29 needed by fxp_mul_q30 */
	int32_t e = 3 * l2 + fxp_mul_q30(l2, 11059494);
	return e < Q16(-120.0) ? Q16(-120.0) : e;
}

void synthesise_one_frame(codec2_decoder_t *c2, int16_t *speech, uint32_t W, int L, int voiced,
			  const int32_t *lA, const uint32_t *phA)
{
	_Static_assert(sizeof(c2->scratch) >= FXP_RFFT_BUF * sizeof(cq15_t), "scratch too small for the irfft");
	uint32_t th[MAX_AMP + 1]; /* final harmonic phases */

	/* --- phase_synth_zero_order: excitation phase track, H = conj(A) so arg P = ph - arg A --- */
	c2->ex_phase += W * N_SAMP;
	const uint32_t phi0 = c2->ex_phase;
	for (int m = 1; m <= L; m++)
	{
		uint32_t ph = voiced ? phi0 * (uint32_t)m
					    : (uint32_t)codec2_rand(&c2->next_rn) * RAND_TO_ANGLE;
		th[m] = ph - phA[m];
	}

	/* --- postfilter: background estimate, random phase for voiced harmonics below it --- */
	int32_t e = frame_energy_db(lA, L);
	if (!voiced)
	{
		if (e < Q16(BG_THRESH))
			c2->bg_est += fxp_mul_q30(e - c2->bg_est, 107374182); /* BG_BETA = 0.1 in Q30 */
	}
	else
	{
		/* A < 10^((bg + margin) / 20)  <=>  log2 A < (bg + margin) log2(10) / 20 */
		int32_t lthresh = fxp_mul_q30(c2->bg_est + Q16(BG_MARGIN), 178344657);
		for (int m = 1; m <= L; m++)
			if (lA[m] < lthresh)
				th[m] = (uint32_t)codec2_rand(&c2->next_rn) * RAND_TO_ANGLE;
	}

	/* --- synthesise: Sw[b_m] = A[m] e^{j th[m]}, scaled so max A -> 2^14 --- */
	cq15_t *Sw = (cq15_t *)c2->scratch;
	memset(Sw, 0, FXP_RFFT_BUF * sizeof(cq15_t));
	int32_t lmax = INT32_MIN;
	for (int m = 1; m <= L; m++)
		if (lA[m] > lmax)
			lmax = lA[m];
	int ip = lmax > LA_EMPTY_LIM ? (lmax + 65535) >> 16 : 0; /* ceil(log2 max A) */
	for (int m = 1; m <= L; m++)
	{
		uint32_t b = fxp_bin_round((uint32_t)m * W); /* round(m Wo / (2 pi / 512)), ties up */
		if (b > FFT_DEC / 2 - 1)
			b = FFT_DEC / 2 - 1;
		int32_t a14 = 0;
		int32_t d = lA[m] - ip * 65536;
		if (lA[m] > LA_EMPTY_LIM && d > -(31 << 16))
			a14 = (int32_t)((fxp_exp2_neg_q30(d) + (1u << 15)) >> 16); /* Q14, <= 16384 */
		int16_t c, s;
		fxp_sincos_q15(th[m], &c, &s);
		Sw[b].r = (int16_t)((a14 * c + (1 << 14)) >> 15); /* |a14 c| < 2^29 */
		Sw[b].i = (int16_t)((a14 * s + (1 << 14)) >> 15);
	}
	/* true sw[n] = buf[n] 2^(ex), sw in sample units */
	int ex = fxp_irfft512(Sw) + ip - 14;
	const int16_t *sw = (const int16_t *)Sw;

	/* --- overlap-add with the Parzen window into Sn_ (Q4): sw Pn 2^(ex - 14 + 4) --- */
	int32_t *Sn_ = c2->Sn_;
	const int16_t *Pn = c2->Pn;
	int sh = 14 - SN_Q - ex; /* right shift of the int32 product sw * Pn (|.| <= 2^29) */
	FXP_ASSERT(sh >= -31); /* wsat(., -sh) needs a left shift < 32 (would need A ~ 2^55) */
#define WIN(x, w) (sh >= 31 ? 0 : sh > 0 ? shr_r((int32_t)(x) * (w), sh) : wsat((int32_t)(x) * (w), -sh))
	memmove(Sn_, &Sn_[N_SAMP], (N_SAMP - 1) * sizeof(int32_t));
	Sn_[N_SAMP - 1] = 0;
	for (int i = 0; i < N_SAMP - 1; i++)
		Sn_[i] = add_sat32(Sn_[i], WIN(sw[FFT_DEC - N_SAMP + 1 + i], Pn[i]));
	for (int j = 0; j < N_SAMP + 1; j++)
		Sn_[N_SAMP - 1 + j] = WIN(sw[j], Pn[N_SAMP - 1 + j]);
#undef WIN

	/* --- ear_protection: above 30000, gain (30000 / max)^2 --- */
	uint32_t mx = 0;
	for (int i = 0; i < N_SAMP; i++)
	{
		uint32_t a = Sn_[i] < 0 ? 0u - (uint32_t)Sn_[i] : (uint32_t)Sn_[i];
		if (a > mx)
			mx = a;
	}
	if (mx > (30000u << SN_Q))
	{
		uint32_t q = (mx + (1u << (SN_Q - 1))) >> SN_Q; /* max in samples, > 30000 */
		uint32_t g = (30000u << 15) / q;                /* Q15, < 32768 */
		int16_t g2 = (int16_t)((g * g + (1u << 14)) >> 15);
		for (int i = 0; i < N_SAMP; i++)
			Sn_[i] = mul32x16_q15(Sn_[i], g2);
	}

	/* --- output: Q4 to samples truncating toward zero and clamping to +-32767, as the float code's
	 * (int16_t) cast (rounding instead would raise quiet-signal levels by ~0.01 dB vs the reference) --- */
	for (int i = 0; i < N_SAMP; i++)
	{
		int32_t v = Sn_[i] >= 0 ? Sn_[i] >> SN_Q : -(int32_t)((0u - (uint32_t)Sn_[i]) >> SN_Q);
		speech[i] = (int16_t)(v > 32767 ? 32767 : (v < -32767 ? -32767 : v));
	}
#ifdef FXP_SYNTH_DEBUG
	fxp_synth_debug_hook(L, W, voiced, e, c2->bg_est, lA, th);
#endif
}

void synthesis_init(codec2_decoder_t *c2)
{
	make_synthesis_window(c2->Pn);
}
