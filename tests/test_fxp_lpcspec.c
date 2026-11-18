/* Host tests for fixed-point aks_to_mag2 (fxp_lpcspec.c) and fxp_math.c against double references. */
#include "c2_float_consts.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "fxp_lpcspec.h"

extern const uint16_t delta_lsp_cb[10][32];
#ifdef FXP_LPCSPEC_DEBUG
extern int32_t fxp_dbg_la[256], fxp_dbg_lag[256], fxp_dbg_leb, fxp_dbg_lea;
#endif

static int fails;
#define CHECK(c, ...) do { if (!(c)) { fails++; printf("FAIL: " __VA_ARGS__); printf("\n"); } } while (0)

static double urand(void) { return rand() / (RAND_MAX + 1.0); }

/* --- double references (same algorithm as float lsp_to_lpc / aks_to_mag2, exact math) --- */
static void ref_lsp_to_lpc(const double *w, double *ak)
{
	double p[12] = {1}, q[12] = {1};
	int n = 1;
	for (int s = 0; s < 5; s++, n += 2)
		for (int j = n + 1; j >= 0; j--)
		{
			double c0 = cos(w[2 * s]), c1 = cos(w[2 * s + 1]);
			double p1 = j >= 1 && j - 1 < n ? p[j - 1] : 0, p2 = j >= 2 ? p[j - 2] : 0;
			double q1 = j >= 1 && j - 1 < n ? q[j - 1] : 0, q2 = j >= 2 ? q[j - 2] : 0;
			p[j] = (j < n ? p[j] : 0) - 2 * c0 * p1 + p2;
			q[j] = (j < n ? q[j] : 0) - 2 * c1 * q1 + q2;
		}
	for (int j = 0; j <= 10; j++)
		ak[j] = 0.5 * (p[j] + (j ? p[j - 1] : 0) + q[j] - (j ? q[j - 1] : 0));
}

/* float lsp_to_lpc() exactly as codec2 (lpc.c) builds it, in double: independent check of the P/Q
 * convention used by ref_lsp_to_lpc. */
static void ref_codec2_lsp_to_lpc(const double *freq, double *ak)
{
	double Wp[4 * 5 + 2] = {0}, xf[10];
	for (int i = 0; i < 10; i++)
		xf[i] = cos(freq[i]);
	double xin1 = 1.0, xin2 = 1.0;
	for (int j = 0; j <= 10; j++)
	{
		double *pw = Wp, *n1 = 0, *n2, *n3, *n4 = 0;
		double xout1 = 0, xout2 = 0;
		for (int i = 0; i < 5; i++)
		{
			n1 = pw + (i * 4);
			n2 = n1 + 1;
			n3 = n2 + 1;
			n4 = n3 + 1;
			xout1 = xin1 - 2 * xf[2 * i] * *n1 + *n2;
			xout2 = xin2 - 2 * xf[2 * i + 1] * *n3 + *n4;
			*n2 = *n1;
			*n4 = *n3;
			*n1 = xin1;
			*n3 = xin2;
			xin1 = xout1;
			xin2 = xout2;
		}
		xout1 = xin1 + *(n4 + 1);
		xout2 = xin2 - *(n4 + 2);
		ak[j] = (xout1 + xout2) * 0.5;
		*(n4 + 1) = xin1;
		*(n4 + 2) = xin2;
		xin1 = 0.0;
		xin2 = 0.0;
	}
}

static void ref_dft(const double *a, int k, double *re, double *im)
{
	double r = 0, i = 0;
	for (int n = 0; n <= 10; n++)
	{
		r += a[n] * cos(2 * M_PI * k * n / 512);
		i -= a[n] * sin(2 * M_PI * k * n / 512);
	}
	*re = r;
	*im = i;
}

/* random decoded LSP frame in half-Hz: a codebook frame, or the average of two (interpolated frame) */
static void rand_lsp_hz2(int16_t *hz2)
{
	int a[10] = {0}, b[10] = {0}, ha = 0, hb = 0, interp = rand() & 1;
	for (int i = 0; i < 10; i++)
	{
		ha += delta_lsp_cb[i][rand() & 31];
		hb += delta_lsp_cb[i][rand() & 31];
		a[i] = ha;
		b[i] = hb;
	}
	for (int i = 0; i < 10; i++)
		hz2[i] = (int16_t)(interp ? a[i] + b[i] : 2 * a[i]);
}

/* near-unstable frames: minimum codebook deltas (25 Hz) packed into pairs */
static void tight_lsp_hz2(int16_t *hz2)
{
	int h = 0;
	for (int i = 0; i < 10; i++)
	{
		h += (i & 1) ? 25 : delta_lsp_cb[i][rand() & 31];
		hz2[i] = (int16_t)(2 * h);
	}
}

int main(void)
{
	srand(1234);

	/* ---- log2 / exp2 ---- */
	double w_log = 0, w_exp = 0, w_spl = 0;
	for (int t = 0; t < 2000000; t++)
	{
		uint32_t v = (uint32_t)(urand() * 4294967295.0) >> (rand() & 31);
		if (!v)
			continue;
		double e = fabs(fxp_log2_q16(v) / 65536.0 - log2((double)v));
		w_log = e > w_log ? e : w_log;
		int32_t x = -(int32_t)(urand() * 31 * 65536);
		double ex = exp2(x / 65536.0), got = fxp_exp2_neg_q30(x) / 1073741824.0;
		double rel = fabs(got - ex) / ex;
		if (ex > 1e-4) /* below that the Q30 output quantum dominates */
			w_exp = rel > w_exp ? rel : w_exp;
		int32_t y = (int32_t)((urand() - 0.5) * 80 * 65536);
		uint32_t m;
		int ee;
		fxp_exp2_split(y, &m, &ee);
		double rs = fabs(ldexp((double)m, ee - 30) / exp2(y / 65536.0) - 1);
		w_spl = rs > w_spl ? rs : w_spl;
	}
	printf("log2: max abs err %.2e   exp2_neg: max rel err %.2e   exp2_split: %.2e\n", w_log, w_exp, w_spl);
	CHECK(w_log < 3e-5, "log2");
	CHECK(w_exp < 3e-5, "exp2_neg");
	CHECK(w_spl < 3e-5, "exp2_split");

	/* ---- cos tables / mul_q30 ---- */
	double w_cos = 0;
	for (int h = 0; h < 40000; h++)
	{
		double e = fabs(fxp_lsp_cos_q30(h) / 1073741824.0 - cos(M_PI * h / 8000.0));
		w_cos = e > w_cos ? e : w_cos;
	}
	printf("lsp cos: max abs err %.2e\n", w_cos);
	CHECK(w_cos < 1e-8, "lsp cos");
	double w_mul = 0;
	for (int t = 0; t < 1000000; t++)
	{
		int32_t a = (int32_t)((urand() * 2 - 1) * ((1 << 29) - 1)), b = (int32_t)((urand() * 2 - 1) * (1 << 30));
		double e = fabs(fxp_mul_q30(a, b) - (double)a * b / 1073741824.0);
		w_mul = e > w_mul ? e : w_mul;
	}
	printf("mul_q30: max err %.2f LSB\n", w_mul);
	CHECK(w_mul <= 1.01, "mul_q30");

	/* ---- atan2 / sincos (uint32 angles, 2^32 = 2 pi) ---- */
	const double U2R = 2 * M_PI / 4294967296.0;
	double w_at = 0, w_sc = 0;
	for (int t = 0; t < 2000000; t++)
	{
		int32_t y, x;
		if (t < 16) /* axes, diagonals, extremes */
		{
			static const int32_t v[4] = {0, 1, INT32_MAX, -INT32_MAX};
			y = v[t & 3];
			x = v[t >> 2];
			if (!x && !y)
				continue;
		}
		else
		{
			int sy = rand() & 31, sx = (rand() & 1) ? sy : rand() & 31;
			y = (int32_t)((urand() * 2 - 1) * 2147483647.0) >> sy;
			x = (int32_t)((urand() * 2 - 1) * 2147483647.0) >> sx;
			if ((x == 0 && y == 0) || (abs(x) < 64 && abs(y) < 64)) /* tiny vectors: arg ill defined */
				continue;
		}
		double e = fabs(remainder(fxp_atan2_u32(y, x) * U2R - atan2((double)y, (double)x), 2 * M_PI));
		w_at = fmax(w_at, e);
		uint32_t a = (uint32_t)(urand() * 4294967295.0);
		int16_t c, sn;
		fxp_sincos_q15(a, &c, &sn);
		w_sc = fmax(w_sc, fmax(fabs(c / 32768.0 - cos(a * U2R)), fabs(sn / 32768.0 - sin(a * U2R))));
	}
	printf("atan2_u32: max err %.2e rad   sincos_q15: max err %.2e (%.2f LSB)\n", w_at, w_sc, w_sc * 32768);
	CHECK(w_at < 3e-5, "atan2"); /* 16-bit ratio: 2^-16 rad-equivalent quantisation of min/max */
	CHECK(w_sc < 1.6 / 32768, "sincos");

	/* ---- harmonic bins / band edges on the decoder's Wo grid (exact half-up rounding) ----
	 * Wo = W0_MIN + k step (decode_Wo, float as the decoder) or the mean of two (interp_Wo); W as
	 * codec2_decode. Exact bins per harmonic: (128 + 7k) / 40, or (256 + 7(k1 + k2)) / 80. */
	long n_bin = 0, n_bin_bad = 0;
	{
		const float w0min = (2.0f * (float)M_PI) / 160.0f, w0max = (2.0f * (float)M_PI) / 20.0f;
		const float step = (w0max - w0min) / 128;
		for (int k1 = 0; k1 < 128; k1++)
			for (int k2 = k1; k2 < 128; k2++)
			{
				float a = w0min + step * k1, b = w0min + step * k2;
				float Wo = k1 == k2 ? a : a + 0.5f * (b - a); /* interp_Wo */
				uint32_t W = (uint32_t)lrintf(Wo * 683565275.6f);
				long num = k1 == k2 ? 2 * (128 + 7 * k1) : 256 + 7 * (k1 + k2); /* bins/harmonic * 80 */
				int L = (int)((float)M_PI / Wo);
				for (int m = 1; m <= L; m++)
					for (int h = -1; h <= 1; h++) /* (m - 0.5), m, (m + 0.5) */
					{
						long x2 = (2 * m + h) * num;          /* 160 * frequency in bins */
						long ref = (x2 + 80) / 160;           /* floor(f + 0.5), f >= 0 */
						uint32_t ang = (uint32_t)m * W + (h < 0 ? -(W >> 1) : h > 0 ? (W >> 1) : 0);
						n_bin++;
						n_bin_bad += (long)fxp_bin_round(ang) != ref;
					}
			}
	}
	printf("grid Wo bins (incl. half-bin ties): %ld of %ld differ from exact half-up rounding\n", n_bin_bad, n_bin);
	CHECK(n_bin_bad == 0, "bin rounding on the Wo grid");

	/* ---- per-frame: lsp_to_lpc, |A|^2, phase, A[m] ---- */
	long n_edge_diff = 0;
	double w_ak = 0, w_conv = 0, w_la = 0, w_ph = 0, w_Am[2] = {0}, s_Am[2] = {0}, s_la = 0;
	long n_Am[2] = {0}, n_la = 0, n_Am_big[2] = {0};
	int32_t scratch[514];
	uint32_t phA[81];
	double w_phm = 0;
	long n_phm = 0;
	const int NF = 20000;
	for (int f = 0; f < NF; f++)
	{
		int16_t hz2[10];
		if (f % 4 == 3)
			tight_lsp_hz2(hz2);
		else
			rand_lsp_hz2(hz2);
		double w[10], ak[11], ak2[11];
		int32_t c[10], aq[11];
		for (int i = 0; i < 10; i++)
		{
			w[i] = M_PI * hz2[i] / 8000.0;
			c[i] = fxp_lsp_cos_q30(hz2[i]);
		}
		ref_lsp_to_lpc(w, ak);
		if (f < 1000)
		{
			ref_codec2_lsp_to_lpc(w, ak2);
			for (int j = 0; j <= 10; j++)
				w_conv = fmax(w_conv, fabs(ak[j] - ak2[j]));
		}
		fxp_lsp_to_lpc_q19(c, aq);
		for (int j = 0; j <= 10; j++)
			w_ak = fmax(w_ak, fabs(aq[j] / 524288.0 - ak[j]));

		/* exact spectrum + reference aks_to_mag2 */
		double A2[256], A2gv[256], Pw[256], eb = 0, ea = 0;
		double ag[11];
		for (int j = 0; j <= 10; j++)
			ag[j] = ak[j] * pow(0.5, j);
		for (int k = 0; k < 256; k++)
		{
			double re, im, gr, gi;
			ref_dft(ak, k, &re, &im);
			ref_dft(ag, k, &gr, &gi);
			A2[k] = re * re + im * im + 1e-6;
			double A2g = gr * gr + gi * gi + 1e-6;
			A2gv[k] = A2g;
			Pw[k] = pow(sqrt(A2g / A2[k]) + 1e-5, 0.4) / A2[k];
			eb += 1 / A2[k];
			ea += Pw[k];
			uint32_t p;
			int32_t la = fxp_lsp_spec_bin(c, k, &p);
			double el = fabs(la / 65536.0 - log2(A2[k])) * 10 * log10(2.0); /* dB */
			w_la = fmax(w_la, el);
			s_la += el * el;
			n_la++;
			double dph = p * U2R - atan2(im, re);
			dph = fabs(remainder(dph, 2 * M_PI));
			if (re * re + im * im > 1e-10)
				w_ph = fmax(w_ph, dph);
		}
		for (int k = 0; k < 64; k++)
			Pw[k] *= 1.96;
		double E = exp2(urand() * 20 - 4); /* energy range of decode_energy */
		double gain = E * eb / ea;
		double Wo = 2 * M_PI / (20 + urand() * 140);
		int L = (int)(M_PI / Wo);
		uint32_t W = (uint32_t)lrint(Wo / U2R);
		int16_t am[81], bm[81];
		double Aref[81];
		for (int m = 1; m <= L; m++)
		{
			/* band edges: the kernel's integer form of (int)((m -+ 0.5) Wo / (2 pi / 512) + 0.5) */
			am[m] = (int16_t)(((uint64_t)m * W - W / 2 + (1u << 22) + FXP_BIN_TIE) >> 23);
			int b = (int)(((uint64_t)m * W + W / 2 + (1u << 22) + FXP_BIN_TIE) >> 23);
			bm[m] = (int16_t)(b > 256 ? 256 : b);
			int bf = (int)((m + 0.5) * Wo / (2 * M_PI / 512) + 0.5), af = (int)((m - 0.5) * Wo / (2 * M_PI / 512) + 0.5);
			n_edge_diff += (af != am[m]) + ((bf > 256 ? 256 : bf) != bm[m]);
			double Em = 0;
			for (int i = am[m]; i < bm[m]; i++)
				Em += Pw[i];
			Aref[m] = sqrt(Em * gain);
		}
		int32_t Al[81];
		fxp_aks_to_mag2(hz2, (int32_t)lrint(log2(E) * 65536), W, L, Al, phA, scratch);
		for (int m = 1; m <= L; m++)
		{
			int bh = (int)(((uint64_t)m * W + (1u << 22) + FXP_BIN_TIE) >> 23);
			bh = bh > 255 ? 255 : bh;
			double re, im;
			ref_dft(ak, bh, &re, &im);
			if (re * re + im * im > 1e-10)
			{
				w_phm = fmax(w_phm, fabs(remainder(phA[m] * U2R - atan2(im, re), 2 * M_PI)));
				n_phm++;
			}
			if (Aref[m] == 0)
				continue;
			double e = fabs(Al[m] / 65536.0 - log2(Aref[m])) * 20 * log10(2.0);
			if (getenv("DIAG") && e > 0.1)
			{
				double Em = 0, mx = 0, smin = 1e300;
				for (int i = am[m]; i < bm[m]; i++) { Em += Pw[i]; mx = fmax(mx, Pw[i]); }
				for (int k = 0; k < 256; k++) smin = fmin(smin, A2[k]);
#ifdef FXP_LPCSPEC_DEBUG
				{
					double wla = 0, wlag = 0; int kla = 0, klag = 0;
					for (int k = 0; k < 256; k++)
					{
						double a = fabs(fxp_dbg_la[k] / 65536.0 - log2(A2[k])), g = fabs(fxp_dbg_lag[k] / 65536.0 - log2(A2gv[k]));
						if (a > wla) { wla = a; kla = k; }
						if (g > wlag) { wlag = g; klag = k; }
					}
					printf("  dbg: la worst %.4f (k=%d A2=%.2e)  lag worst %.4f (k=%d A2g=%.2e)  leb err %.4f  lea err %.4f (log2)\n",
					       wla, kla, A2[kla], wlag, klag, A2gv[klag], fxp_dbg_leb / 65536.0 - log2(eb), fxp_dbg_lea / 65536.0 - log2(ea / 1.0));
				}
#endif
				printf("f=%d tight=%d m=%d L=%d band=[%d,%d) err=%.3f dB  Aref=%.3e  Em/mx=%.3f  minA2=%.2e eb=%.3e ea=%.3e\n",
				       f, f % 4 == 3, m, L, am[m], bm[m], e, Aref[m], Em / mx, smin, eb, ea);
			}
			int tg = f % 4 == 3;
			w_Am[tg] = fmax(w_Am[tg], e);
			s_Am[tg] += e * e;
			n_Am[tg]++;
			n_Am_big[tg] += e > 0.1;
		}
	}
	printf("P/Q convention vs codec2 lsp_to_lpc: max |diff| %.2e\n", w_conv);
	printf("lsp_to_lpc Q19: max abs err %.2e (LSB %.2e)\n", w_ak, 1 / 524288.0);
	printf("|A|^2 (product form): rms %.5f dB  max %.5f dB  (%ld bins)\n", sqrt(s_la / n_la), w_la, n_la);
	printf("arg A: max err %.2e rad (all bins)   %.2e rad (harmonic bins, %ld)\n", w_ph, w_phm, n_phm);
	printf("band edges differing from the float formula: %ld\n", n_edge_diff);
	for (int tg = 0; tg < 2; tg++)
		printf("A[m] %-28s rms %.5f dB  max %.5f dB  >0.1 dB: %ld of %ld\n",
		       tg ? "(synthetic 5x25 Hz pairs):" : "(codebook frames):", sqrt(s_Am[tg] / n_Am[tg]), w_Am[tg], n_Am_big[tg], n_Am[tg]);
	CHECK(w_conv < 1e-9, "P/Q convention");
	CHECK(w_ak < 4e-5, "lsp_to_lpc");
	CHECK(w_la < 0.01, "|A|^2");
	CHECK(w_ph < 1e-3, "phase");
	CHECK(w_phm < 1e-3, "phase at harmonics");
	CHECK(w_Am[0] < 0.05, "A[m] codebook frames");
	/* All five LSP pairs at the 25 Hz minimum cluster the zeros of A(z/0.5) too: |A(z/g)|^2 reaches
	 * ~1e-5 against an L1 of ~30, beyond 16-bit taps / Q14 twiddles (docs/05). Bounded, not exact. */
	CHECK(w_Am[1] < 1.0, "A[m] synthetic near-unstable frames");

	printf(fails ? "FAILED (%d)\n" : "all lpcspec tests passed\n", fails);
	return fails != 0;
}
