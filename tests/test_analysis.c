/*
 * test_analysis.c — host checks of the fixed-point encoder analysis (src/analysis.c, docs/09).
 * Includes analysis.c to reach its static helpers. Run by scripts/test-fxp.sh (FXP_CHECK + UBSan + ASan).
 *  1. dft_speech: Sw 2^g against a double DFT of Sn * hann (scale and accuracy)
 *  2. harmonic-bin stepping: every candidate of every NLP bin against round(m 512 d / n) in double
 *  3. log2_plus_floor against double
 *  4. encode_Wo_p against the float encode_Wo formula on every reachable pitch rational
 */
#include "c2_float_consts.h"
#include "../src/analysis.c"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int fails;
#define CHECK(c, ...) do { if (!(c)) { fails++; printf("FAIL: " __VA_ARGS__); printf("\n"); } } while (0)

static uint32_t rng = 12345;
static int32_t rnd(int32_t lo, int32_t hi)
{
	rng = rng * 1103515245u + 12345u;
	return (int32_t)((uint32_t)lo + (rng >> 1) % ((uint32_t)hi - (uint32_t)lo + 1u));
}

static void test_dft(void)
{
	static int16_t Sn[M_PITCH];
	static cq15_t Sw[FXP_RFFT_BUF];
	static int32_t x32[FFT_ENC];
	double worst = 1e9;
	const int amps[] = {1, 30, 1000, 32767};
	for (int trial = 0; trial < 64; trial++)
	{
		int a = amps[trial % 4];
		double f = 0.01 + 0.4 * (trial / 64.0);
		for (int i = 0; i < M_PITCH; i++)
			Sn[i] = (int16_t)fmax(-32768, fmin(32767, lrint(a * (0.7 * sin(2 * M_PI * f * i) + 0.3 * (rnd(-1000, 1000) / 1000.0)))));
		int g = dft_speech(Sw, x32, Sn);
		double es = 0, ee = 0;
		for (int k = 0; k <= 256; k++)
		{
			double re = 0, im = 0;
			for (int i = 0; i < NW / 2; i++)
			{
				double v0 = Sn[160 + i] * (double)fxp_hann_q15[139 + i] / 32768.0;
				double v1 = Sn[21 + i] * (double)fxp_hann_q15[i] / 32768.0;
				int n1 = FFT_ENC - NW / 2 + i;
				re += v0 * cos(2 * M_PI * k * i / 512) + v1 * cos(2 * M_PI * k * n1 / 512);
				im -= v0 * sin(2 * M_PI * k * i / 512) + v1 * sin(2 * M_PI * k * n1 / 512);
			}
			double sr = ldexp(Sw[k].r, g), si = ldexp(Sw[k].i, g);
			es += re * re + im * im;
			ee += (re - sr) * (re - sr) + (im - si) * (im - si);
		}
		double snr = 10 * log10(es / (ee + 1e-30));
		if (snr < worst)
			worst = snr;
	}
	printf("dft_speech vs double: worst spectrum SNR %.1f dB (64 frames, amplitude 1..32767)\n", worst);
	CHECK(worst > 55, "dft_speech SNR %.1f dB", worst);
}

static void test_bins(void)
{
	long n_checked = 0;
	for (uint32_t b = 16; b <= 128; b++)
	{
		/* coarse and fine candidate sets (fine around every coarse candidate) */
		for (int j = -5; j <= 5; j++)
		{
			uint32_t nc = 2560 + (uint32_t)(j + 5) * b - 5 * b;
			for (int fine = 0; fine <= 9; fine++)
			{
				uint32_t n = fine ? 4 * nc - 4 * b + (uint32_t)(fine - 1) * b : nc, d = fine ? 4 * b : b;
				int L = fine ? (int)(nc / (2 * b)) : (int)(2560 / (2 * b));
				uint32_t n2 = 2 * n, num = 1024 * d, qi = num / n2, ri = num - qi * n2, q = 0, r = n;
				for (int m = 1; m <= L; m++)
				{
					q += qi;
					r += ri;
					if (r >= n2)
					{
						r -= n2;
						q++;
					}
					uint32_t ref = (uint32_t)floor((double)m * 512.0 * d / n + 0.5);
					CHECK(q == ref && q < FFT_ENC, "bin n %u d %u m %d: %u vs %u", n, d, m, q, ref);
					n_checked++;
				}
			}
		}
	}
	printf("harmonic bins: %ld stepped values equal round(m 512 d / n)\n", n_checked);
}

static void test_log2(void)
{
	double worst = 0;
	for (int t = 0; t < 200000; t++)
	{
		int ef = rnd(-60, 20), ev = rnd(0, 8);
		uint64_t v = t % 7 == 0 ? 0 : ((uint64_t)(uint32_t)rnd(0, 0x7fffffff) << rnd(0, 30)) >> rnd(0, 40);
		double ref = log2(ldexp((double)v, ev) + ldexp((double)fxp_ana_floor_q24, ef));
		double got = log2_plus_floor(v, ev, ef) / 65536.0;
		if (fabs(got - ref) > worst)
			worst = fabs(got - ref);
	}
	printf("log2_plus_floor: max error %.2e (log2 units)\n", worst);
	CHECK(worst < 1e-4, "log2_plus_floor error %.2e", worst);
}

static int encode_Wo_float(float Wo)
{
	float norm = (Wo - W0_MIN) / (W0_MAX - W0_MIN);
	int index = (int)(128 * norm + 0.5f);
	return index < 0 ? 0 : index > 127 ? 127 : index;
}

static void test_wo(void)
{
	int n_checked = 0, diff = 0;
	for (uint32_t d = 1; d <= 512; d++)
		for (uint32_t n = 20 * d; n <= 160 * d; n++)
		{
			double x = 128.0 * (160.0 * d - n) / (7.0 * n) + 0.5; /* exact index argument */
			int ref = (int)floor(x);
			if (ref > 127)
				ref = 127;
			int got = encode_Wo_p(n, d);
			CHECK(got == ref, "encode_Wo_p %u/%u: %d vs %d", n, d, got, ref);
			if (fabs(x - floor(x + 1e-9)) > 1e-6) /* away from exact ties: the float formula must agree */
				diff += encode_Wo_float(TWO_PI * d / n) != got;
			n_checked++;
		}
	printf("encode_Wo_p: %d rationals exact; float encode_Wo differs on %d (near-ties in float)\n", n_checked, diff);
}

int main(void)
{
	test_dft();
	test_bins();
	test_log2();
	test_wo();
	printf(fails ? "test_analysis: %d FAILURES\n" : "test_analysis: all passed\n", fails);
	return fails != 0;
}
