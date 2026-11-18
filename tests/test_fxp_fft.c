/* Host accuracy test for fxp_rfft512 / fxp_irfft512 against a double-precision DFT and kiss_fftr(i).
 * Build: see scripts/test-fxp.sh. Prints SNR (dB) per signal class; exits non-zero on failure. */
#include "c2_float_consts.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fxp_fft.h"
#include "kiss_fftr.h"

#define N 512
static double rnd(unsigned *s) { *s = *s * 1103515245u + 12345u; return ((*s >> 8) & 0xffffff) / 16777216.0 - 0.5; }

static void dft(const double *x, double *Xr, double *Xi)
{
	for (int k = 0; k <= N / 2; k++) {
		double r = 0, i = 0;
		for (int n = 0; n < N; n++) { double a = -2 * M_PI * k * n / N; r += x[n] * cos(a); i += x[n] * sin(a); }
		Xr[k] = r; Xi[k] = i;
	}
}
static void idft(const double *Xr, const double *Xi, double *x)
{
	for (int n = 0; n < N; n++) {
		double s = Xr[0] + Xr[N / 2] * ((n & 1) ? -1 : 1);
		for (int k = 1; k < N / 2; k++) { double a = 2 * M_PI * k * n / N; s += 2 * (Xr[k] * cos(a) - Xi[k] * sin(a)); }
		x[n] = s;
	}
}

static double worst_fwd = 1e9, worst_inv = 1e9, worst_bin_db = 0;

/* LPC case: worst per-bin |A|^2 error in dB over bins 0..255 (formant peaks have small |A|). */
static void lpc_bin_error(const char *name, const double *a, int n)
{
	int16_t x[N] = {0};
	double amax = 0;
	for (int i = 0; i < n; i++) if (fabs(a[i]) > amax) amax = fabs(a[i]);
	for (int i = 0; i < n; i++) x[i] = (int16_t)lrint(a[i] * 32767 / amax);
	cq15_t buf[FXP_RFFT_BUF];
	for (int i = 0; i < N / 2; i++) { buf[i].r = x[2 * i]; buf[i].i = x[2 * i + 1]; }
	int e = fxp_rfft512(buf);
	double worst = 0, pmin = 1e300, pmax = 0;
	for (int k = 0; k < N / 2; k++) {
		double r = 0, im = 0;
		for (int i = 0; i < n; i++) { r += x[i] * cos(-2 * M_PI * k * i / N); im += x[i] * sin(-2 * M_PI * k * i / N); }
		double pr = r * r + im * im, pf = pow(ldexp(buf[k].r, e), 2) + pow(ldexp(buf[k].i, e), 2);
		double d = fabs(10 * log10((pf + 1e-9) / (pr + 1e-9)));
		if (d > worst) worst = d;
		if (pr < pmin) pmin = pr;
		if (pr > pmax) pmax = pr;
	}
	printf("%-34s worst bin |A|^2 error %.4f dB (|A|^2 range %.1f dB)\n", name, worst, 10 * log10(pmax / pmin));
	if (worst > worst_bin_db) worst_bin_db = worst;
}

/* Sparse DFT (fxp_rdft512_sparse) */
static int rdft_fail;
/* exactness: int64 reference with the same Q14 table must match bit for bit (no overflow, no rounding) */
static void rdft_exact(const char *name, const int16_t *a, int len)
{
	cq31_t X[257];
	fxp_rdft512_sparse(a, len, X);
	for (int k = 0; k <= 256; k++) {
		long long r = 0, im = 0;
		for (int n = 0; n < len; n++) {
			unsigned m = (unsigned)(n * k) & 511;
			r += (long long)a[n] * fxp_cos512_q14[m];
			im -= (long long)a[n] * fxp_cos512_q14[(m + 384) & 511];
		}
		if (r != X[k].r || im != X[k].i) { printf("%-34s MISMATCH at bin %d\n", name, k); rdft_fail = 1; return; }
	}
	printf("%-34s bit-exact vs int64 reference\n", name);
}
static double worst_rdft_db, worst_rdft_ph;
/* accuracy on float LPC taps, scaled exactly as a caller would (fxp_rdft_taps_shift) */
static void rdft_lpc(const char *name, const double *a, int len)
{
	double amax = 0, l1 = 0;
	for (int i = 0; i < len; i++) { if (fabs(a[i]) > amax) amax = fabs(a[i]); l1 += fabs(a[i]); }
	/* taps as Q20 int32 first (stand-in for a fixed-point lsp_to_lpc output), then normalised */
	int32_t q[32]; uint32_t qm = 0, ql1 = 0;
	for (int i = 0; i < len; i++) { q[i] = (int32_t)lrint(a[i] * (1 << 20)); uint32_t v = (uint32_t)(q[i] < 0 ? -q[i] : q[i]); if (v > qm) qm = v; ql1 += v; }
	int sh = fxp_rdft_taps_shift(qm, ql1);
	int16_t t[32];
	for (int i = 0; i < len; i++) t[i] = (int16_t)(sh >= 0 ? shl(q[i], sh) : shr_r(q[i], -sh));
	cq31_t X[257];
	fxp_rdft512_sparse(t, len, X);
	double worst = 0, wph = 0, pmin = 1e300, pmax = 0;
	for (int k = 0; k < 256; k++) {
		double r = 0, im = 0;
		for (int i = 0; i < len; i++) { r += a[i] * cos(-2 * M_PI * k * i / N); im += a[i] * sin(-2 * M_PI * k * i / N); }
		double fr = ldexp(X[k].r, -14 - sh - 20), fi = ldexp(X[k].i, -14 - sh - 20);
		double pr = r * r + im * im, pf = fr * fr + fi * fi;
		double d = fabs(10 * log10(pf / pr));
		double ph = fabs(remainder(atan2(fi, fr) - atan2(im, r), 2 * M_PI));
		if (d > worst) worst = d;
		if (ph > wph) wph = ph;
		if (pr < pmin) pmin = pr;
		if (pr > pmax) pmax = pr;
	}
	printf("%-34s worst bin |A|^2 error %.4f dB, phase %.4f rad (range %.1f dB, tap shift %d)\n", name, worst, wph, 10 * log10(pmax / pmin), sh);
	if (worst > worst_rdft_db) worst_rdft_db = worst;
	if (wph > worst_rdft_ph) worst_rdft_ph = wph;
}

/* Build a stable LPC polynomial from formant poles: radius r, frequencies f (Hz at 8 kHz). */
static int lpc_from_poles(double *a, const double *f, const double *r, int np)
{
	int n = 1; a[0] = 1;
	for (int p = 0; p < np; p++) {
		double c1 = -2 * r[p] * cos(2 * M_PI * f[p] / 8000), c2 = r[p] * r[p];
		double b[16] = {0};
		for (int i = 0; i < n; i++) { b[i] += a[i]; b[i + 1] += c1 * a[i]; b[i + 2] += c2 * a[i]; }
		n += 2; memcpy(a, b, n * sizeof(double));
	}
	return n;
}

/* forward: x given as int16 */
static double test_fwd(const int16_t *x16)
{
	double x[N], Xr[N / 2 + 1], Xi[N / 2 + 1];
	cq15_t buf[FXP_RFFT_BUF];
	for (int n = 0; n < N; n++) x[n] = x16[n];
	for (int n = 0; n < N / 2; n++) { buf[n].r = x16[2 * n]; buf[n].i = x16[2 * n + 1]; }
	int e = fxp_rfft512(buf);
	dft(x, Xr, Xi);
	double sig = 0, err = 0;
	for (int k = 0; k <= N / 2; k++) {
		double r = ldexp(buf[k].r, e), i = ldexp(buf[k].i, e);
		sig += Xr[k] * Xr[k] + Xi[k] * Xi[k];
		err += (r - Xr[k]) * (r - Xr[k]) + (i - Xi[k]) * (i - Xi[k]);
	}
	return sig == 0 ? (err == 0 ? 999 : -999) : 10 * log10(sig / (err + 1e-300));
}

/* inverse: spectrum given as int16 pairs (Hermitian half) */
static double test_inv(const int16_t *Xr16, const int16_t *Xi16)
{
	double Xr[N / 2 + 1], Xi[N / 2 + 1], x[N];
	cq15_t buf[FXP_RFFT_BUF];
	for (int k = 0; k <= N / 2; k++) { Xr[k] = Xr16[k]; Xi[k] = (k == 0 || k == N / 2) ? 0 : Xi16[k]; buf[k].r = Xr16[k]; buf[k].i = Xi16[k]; }
	int e = fxp_irfft512(buf);
	idft(Xr, Xi, x);
	double sig = 0, err = 0;
	for (int n = 0; n < N; n++) {
		double y = ldexp(n & 1 ? buf[n / 2].i : buf[n / 2].r, e);
		sig += x[n] * x[n]; err += (y - x[n]) * (y - x[n]);
	}
	return sig == 0 ? (err == 0 ? 999 : -999) : 10 * log10(sig / (err + 1e-300));
}

static void report(const char *name, double snr, int inv)
{
	printf("%-34s %s SNR %7.2f dB\n", name, inv ? "inv" : "fwd", snr);
	if (inv) { if (snr < worst_inv) worst_inv = snr; } else if (snr < worst_fwd) worst_fwd = snr;
}

/* fxp_tw256 (stage twiddles in reading order) must be exact copies of fxp_w512 entries */
static int test_tw256(void)
{
	int n = 0, bad = 0;
	for (int q = 4; q <= 64; q *= 4)
		for (int j = 1; j < q; j++)
			for (int m = 1; m <= 3; m++, n++) {
				cq15_t a = fxp_tw256[n], b = fxp_w512[m * j * (128 / q)];
				bad += a.r != b.r || a.i != b.i;
			}
	printf("stage twiddle table: %d entries, %d differ from fxp_w512\n", n, bad);
	return bad != 0 || n != 243;
}

int main(void)
{
	unsigned seed = 1;
	if (test_tw256())
		return 1;
	int16_t x[N], Xr[N / 2 + 1], Xi[N / 2 + 1];
	char name[64];

	/* kiss scaling check: forward and inverse must match kiss_fftr / kiss_fftri */
	{
		static char mf[16384], mi[16384];
		size_t lf = sizeof mf, li = sizeof mi;
		kiss_fftr_cfg cf = kiss_fftr_alloc(N, 0, mf, &lf), ci = kiss_fftr_alloc(N, 1, mi, &li);
		if (!cf || !ci) { printf("kiss alloc failed\n"); return 1; }
		float xf[N], yf[N]; kiss_fft_cpx Xk[N / 2 + 1];
		cq15_t buf[FXP_RFFT_BUF];
		for (int n = 0; n < N; n++) { x[n] = (int16_t)(8000 * rnd(&seed)); xf[n] = x[n]; }
		for (int n = 0; n < N / 2; n++) { buf[n].r = x[2 * n]; buf[n].i = x[2 * n + 1]; }
		kiss_fftr(cf, xf, Xk);
		int e = fxp_rfft512(buf);
		double se = 0, ss = 0;
		for (int k = 0; k <= N / 2; k++) {
			double dr = ldexp(buf[k].r, e) - Xk[k].r, di = ldexp(buf[k].i, e) - Xk[k].i;
			se += dr * dr + di * di; ss += Xk[k].r * Xk[k].r + Xk[k].i * Xk[k].i;
		}
		report("vs kiss_fftr (noise 8000)", 10 * log10(ss / se), 0);
		kiss_fftri(ci, Xk, yf);
		e += fxp_irfft512(buf); /* buf holds our forward output (exponent e): round trip */
		se = 0; ss = 0;
		for (int n = 0; n < N; n++) {
			double y = ldexp(n & 1 ? buf[n / 2].i : buf[n / 2].r, e);
			se += (y - yf[n]) * (y - yf[n]); ss += (double)yf[n] * yf[n];
		}
		report("round trip vs kiss_fftri", 10 * log10(ss / se), 1);
		printf("kiss_fftri(kiss_fftr(x))[0] / x[0] = %.3f (unnormalised inverse => %d)\n", yf[0] / x[0], N);
	}

	/* forward: white noise at many amplitudes */
	for (int amp = 1; amp <= 32767; amp = amp * 4 + 3) {
		for (int n = 0; n < N; n++) x[n] = (int16_t)lrint(2 * amp * rnd(&seed));
		snprintf(name, sizeof name, "white amp %d", amp); report(name, test_fwd(x), 0);
	}
	/* full-scale extremes */
	for (int n = 0; n < N; n++) x[n] = (n & 1) ? -32768 : 32767; report("alternating full scale", test_fwd(x), 0);
	for (int n = 0; n < N; n++) x[n] = -32768; report("DC -32768", test_fwd(x), 0);
	for (int n = 0; n < N; n++) x[n] = 32767; report("DC +32767", test_fwd(x), 0);
	memset(x, 0, sizeof x); x[0] = 32767; report("impulse", test_fwd(x), 0);
	memset(x, 0, sizeof x); report("zeros", test_fwd(x), 0);
	for (int n = 0; n < N; n++) x[n] = (int16_t)lrint(32767 * cos(2 * M_PI * 37 * n / N)); report("cos bin 37 full", test_fwd(x), 0);
	for (int n = 0; n < N; n++) x[n] = (int16_t)lrint(30000 * cos(2 * M_PI * 37.3 * n / N) + 2 * cos(2 * M_PI * 150 * n / N));
	report("cos 37.3 + -84 dB tone", test_fwd(x), 0);
	/* LPC-like polynomial: 11 taps, rest zeros (aks_to_mag2 case) */
	{
		const double a[11] = {1, -2.6, 3.9, -4.1, 3.5, -2.4, 1.4, -0.66, 0.25, -0.07, 0.01};
		memset(x, 0, sizeof x);
		for (int i = 0; i < 11; i++) x[i] = (int16_t)lrint(a[i] * 32767 / 4.1);
		report("LPC polynomial (11 taps)", test_fwd(x), 0);
		memset(x, 0, sizeof x);
		for (int i = 0; i < 11; i++) x[i] = (int16_t)lrint(a[i] * 32767 / 4.1 * pow(0.5, i));
		report("LPC*gamma^i polynomial", test_fwd(x), 0);
	}
	{
		double a[16];
		const double f[5] = {500, 1500, 2500, 3300, 3800};
		const double r1[5] = {0.95, 0.93, 0.9, 0.88, 0.85}, r2[5] = {0.99, 0.985, 0.98, 0.97, 0.96}, r3[5] = {0.995, 0.993, 0.99, 0.99, 0.98};
		int n = lpc_from_poles(a, f, r1, 5); lpc_bin_error("LPC formants r~0.9", a, n);
		n = lpc_from_poles(a, f, r2, 5); lpc_bin_error("LPC formants r~0.98", a, n);
		n = lpc_from_poles(a, f, r3, 5); lpc_bin_error("LPC formants r~0.99 (very peaky)", a, n);
	}
	{
		printf("-- sparse real DFT (LPC spectrum kernel) --\n");
		double a[16];
		const double f[5] = {500, 1500, 2500, 3300, 3800};
		const double r1[5] = {0.95, 0.93, 0.9, 0.88, 0.85}, r2[5] = {0.99, 0.985, 0.98, 0.97, 0.96}, r3[5] = {0.995, 0.993, 0.99, 0.99, 0.98};
		int n = lpc_from_poles(a, f, r1, 5); rdft_lpc("rdft LPC formants r~0.9", a, n);
		n = lpc_from_poles(a, f, r2, 5); rdft_lpc("rdft LPC formants r~0.98", a, n);
		n = lpc_from_poles(a, f, r3, 5); rdft_lpc("rdft LPC formants r~0.99 (very peaky)", a, n);
		int16_t t[11];
		for (int i = 0; i < 11; i++) t[i] = (int16_t)((i < 4) ? ((i & 1) ? -32767 : 32767) : 0); /* L1 = 131068 */
		rdft_exact("rdft worst case +-32767 x4", t, 11);
		for (int i = 0; i < 11; i++) t[i] = (int16_t)((i & 1) ? -11915 : 11915); /* L1 = 131065 */
		rdft_exact("rdft worst case alternating x11", t, 11);
		for (int trial = 0; trial < 2000; trial++) {
			double v[11], l1 = 0, vm = 0;
			for (int i = 0; i < 11; i++) { v[i] = rnd(&seed) * ((trial & 1) ? 1 : exp(-0.5 * i)); l1 += fabs(v[i]); vm = fmax(vm, fabs(v[i])); }
			double sc = fmin(131000 / l1, 32767 / vm); /* both preconditions: |a| < 2^15, L1 < 2^17 */
			for (int i = 0; i < 11; i++) t[i] = (int16_t)(v[i] * sc);
			cq31_t X[257]; fxp_rdft512_sparse(t, 11, X);
			if (trial % 500 == 0) rdft_exact("rdft random taps at L1 limit", t, 11);
		}
		for (int i = 0; i < 11; i++) t[i] = 0;
		rdft_exact("rdft zeros", t, 11);
		t[0] = 1; rdft_exact("rdft unit impulse", t, 1);
		int16_t t13[13]; for (int i = 0; i < 13; i++) t13[i] = (int16_t)(1000 - 150 * i);
		rdft_exact("rdft 13 taps (odd tail)", t13, 13);
		rdft_exact("rdft 12 taps (even)", t13, 12);
	}
	/* windowed speech-like: noise burst in the middle */
	for (int n = 0; n < N; n++) x[n] = (int16_t)lrint(20000 * rnd(&seed) * exp(-pow((n - 256) / 80.0, 2)));
	report("gaussian-windowed noise", test_fwd(x), 0);

	/* inverse: sparse harmonic spectra (synthesis case), at several levels */
	for (int L = 10; L <= 80; L += 35) for (int amp = 30; amp <= 30000; amp *= 10) {
		memset(Xr, 0, sizeof Xr); memset(Xi, 0, sizeof Xi);
		double wo = 255.0 / L;
		for (int l = 1; l <= L; l++) {
			int b = (int)(l * wo + 0.5); if (b > 255) b = 255;
			double ph = 2 * M_PI * (rnd(&seed) + 0.5), a = amp * (1.0 / (1 + l * 0.05));
			Xr[b] = (int16_t)lrint(a * cos(ph)); Xi[b] = (int16_t)lrint(a * sin(ph));
		}
		snprintf(name, sizeof name, "harmonics L=%d amp %d", L, amp); report(name, test_inv(Xr, Xi), 1);
	}
	for (int k = 0; k <= N / 2; k++) { Xr[k] = (int16_t)lrint(60000 * rnd(&seed)); Xi[k] = (int16_t)lrint(60000 * rnd(&seed)); }
	report("dense full-scale spectrum", test_inv(Xr, Xi), 1);
	for (int k = 0; k <= N / 2; k++) { Xr[k] = -32768; Xi[k] = -32768; }
	report("all -32768", test_inv(Xr, Xi), 1);
	memset(Xr, 0, sizeof Xr); memset(Xi, 0, sizeof Xi); Xr[0] = 1000; report("DC only", test_inv(Xr, Xi), 1);
	memset(Xr, 0, sizeof Xr); memset(Xi, 0, sizeof Xi); Xr[256] = 1000; report("Nyquist only", test_inv(Xr, Xi), 1);
	memset(Xr, 0, sizeof Xr); memset(Xi, 0, sizeof Xi); report("zeros", test_inv(Xr, Xi), 1);

	printf("worst forward SNR %.2f dB, worst inverse SNR %.2f dB, worst LPC bin error %.3f dB\n", worst_fwd, worst_inv, worst_bin_db);
	printf("sparse DFT: worst LPC bin error %.4f dB, worst phase %.4f rad\n", worst_rdft_db, worst_rdft_ph);
	int ok = worst_fwd > 55 && worst_inv > 55 && !rdft_fail && worst_rdft_db < 0.3 && worst_rdft_ph < 0.03;
	printf("%s (regression threshold: SNR > 55 dB vs double DFT; LPC bin error reported)\n", ok ? "PASS" : "FAIL");
	return ok ? 0 : 1;
}
