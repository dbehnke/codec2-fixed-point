/*
 * fxp_rdft.c — sparse real DFT for short sequences (LPC A(z), 11 taps) at 257 bins of a 512-pt grid.
 *
 * Why not the FFT: |A(e^jw)| has a 40-60 dB dynamic range; its minima are the formant peaks of 1/|A|^2,
 * where the decoder's harmonic amplitudes and phases come from. The 16-bit block-float FFT rounds at
 * every stage relative to the spectral peak, giving up to ~2.4 dB harmonic amplitude error there. A
 * direct DFT with exact int32 accumulation only carries tap (2^-15 of max) and twiddle (Q14) quantisation:
 * 0.005-0.009 dB rms, <= 0.22 dB max (host measurement, see docs/04-fxp-fft.md).
 *
 * Symmetry: with theta = 2 pi n k / 512,  W^{n(256-k)} = (-1)^n e^{+j theta}. Splitting the sum into even
 * (E) and odd (O) n with C = sum a cos(theta), S = sum a sin(theta):
 *   X[k]       = (Ec + Oc) - j (Es + Os)
 *   X[256 - k] = (Ec - Oc) + j (Es - Os)
 * so k = 0..128 yields all 257 bins (k = 128 is produced by both formulas and they agree).
 *
 * Overflow: |a[n] * twiddle| <= |a[n]| * 2^14, so every partial or combined sum is bounded by
 * sum|a| * 2^14 < 2^17 * 2^14 = 2^31 (precondition). No saturation needed.
 */
#include "fxp_fft.h"

int fxp_rdft_taps_shift(uint32_t maxabs, uint32_t l1)
{
	if (maxabs == 0)
		return 0;
	/* want maxabs << s <= 2^15 - 1 and l1 << s <= 2^17 - 1 */
	int s1 = clz32(maxabs) - 17; /* (maxabs << s1) < 2^15 */
	int s2 = clz32(l1) - 15;     /* (l1 << s2) < 2^17 */
	return s1 < s2 ? s1 : s2;
}

#define COS(m) ((int32_t)fxp_cos512_q14[(m)])
#define SIN(m) ((int32_t)fxp_cos512_q14[((m) + 384) & 511]) /* sin t = cos(t - pi/2) */

void fxp_rdft512_sparse(const int16_t *a, int len, cq31_t *X)
{
#ifdef FXP_CHECK
	{
		uint32_t l1 = 0;
		for (int n = 0; n < len; n++)
			l1 += (uint32_t)(a[n] < 0 ? -a[n] : a[n]);
		FXP_ASSERT(l1 < (1u << 17));
	}
#endif
	for (int k = 0; k <= 128; k++)
	{
		int32_t ec = 0, es = 0, oc = 0, os = 0;
		unsigned m = 0;
		int n = 0;
		for (; n + 1 < len; n += 2)
		{
			int32_t a0 = a[n], a1 = a[n + 1];
			ec += a0 * COS(m);
			es += a0 * SIN(m);
			m = (m + k) & 511;
			oc += a1 * COS(m);
			os += a1 * SIN(m);
			m = (m + k) & 511;
		}
		if (n < len)
		{
			int32_t a0 = a[n];
			ec += a0 * COS(m);
			es += a0 * SIN(m);
		}
		X[k].r = ec + oc;
		X[k].i = -(es + os);
		X[256 - k].r = ec - oc;
		X[256 - k].i = es - os;
	}
}
