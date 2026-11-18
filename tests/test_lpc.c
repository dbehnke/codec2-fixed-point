/*
 * test_lpc.c — stage E3 (docs/11-fxp-lpc.md): the fixed-point encoder LPC/LSP analysis and quantisers
 * against the same algorithm in double (tests/lpc_ref.h), with the float implementation it replaced as the
 * yardstick. Component tests of the integer helpers, then (with .raw files as arguments) per-frame index
 * agreement on real speech.
 *
 * lpc.c is included so that its static helpers can be tested directly.
 *
 * With LPC_DUMP=<dir> set, writes <dir>/<n>.idx per input file (n = argument number): per frame 33 int16, the
 * energy and 10 LSP indices of the integer, float and double analyses (scripts/lpc-vs-encoders.py compares
 * them with the float encoders' bitstreams).
 */
#include "c2_float_consts.h"
#include "../src/lpc.c"
#include "quantise.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static int fails;
#define CHECK(c, ...) do { if (!(c)) { printf("FAIL: " __VA_ARGS__); printf("\n"); fails++; } } while (0)

#define REAL float
#define FN(n) flt_##n
#include "lpc_ref.h"
#undef REAL
#undef FN
#define REAL double
#define FN(n) dbl_##n
#include "lpc_ref.h"
#undef REAL
#undef FN

static uint64_t rng = 88172645463325252ull;
static uint32_t rnd(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return (uint32_t)(rng >> 16); }
static int32_t rnds(int bits) { return (int32_t)(rnd() & ((1u << bits) - 1)) - (1 << (bits - 1)); }

static void test_helpers(void)
{
	/* mul_q29 / mul_q30 vs exact floor */
	int64_t w29 = 0, w30 = 0;
	for (int i = 0; i < 1000000; i++)
	{
		int32_t a = rnds(31), b = rnds(28);
		if (i & 1) { int32_t t = a; a = b; b = t; }
		int64_t e = ((int64_t)a * b) >> 29, d = e - mul_q29(a, b);
		if (llabs(d) > w29) w29 = llabs(d);
		e = ((int64_t)a * b) >> 30; d = e - mul_q30(a, b);
		if (llabs(d) > w30) w30 = llabs(d);
	}
	printf("mul_q29 max err %lld LSB, mul_q30 max err %lld LSB\n", (long long)w29, (long long)w30);
	CHECK(w29 <= 3 && w30 <= 3, "mul_q29/q30 error");

	/* div_q30: round(n 2^30 / d) */
	int bad = 0;
	for (int i = 0; i < 1000000; i++)
	{
		uint32_t d = (rnd() & ((1u << 30) - 1)) | 1, n = rnd() % (d + 1);
		if (i < 4) { d = (1u << 30) - 1; n = i < 2 ? d : 0; }
		uint64_t e = (((uint64_t)n << 31) / d + 1) >> 1;
		bad += div_q30(n, d) != e;
	}
	CHECK(bad == 0, "div_q30: %d wrong", bad);

	/* isqrt64 exact */
	bad = 0;
	for (int i = 0; i < 1000000; i++)
	{
		uint64_t n = ((uint64_t)rnd() << 32 | rnd()) >> (rnd() % 20);
		if (i < 1000) n = (uint64_t)i * i - (i & 1);
		uint64_t r = isqrt64(n);
		bad += !(r * r <= n && (r + 1) * (r + 1) > n);
	}
	CHECK(bad == 0, "isqrt64: %d wrong", bad);

	/* acos to Hz on the whole grid vs double fast_acosf */
	double wh = 0;
	for (int32_t u = -6464; u <= 6464; u++)
	{
		double ref = 4000 / M_PI * dbl_acos(u / 6400.0), d = fabs(u_to_hz_q16(u) / 65536.0 - ref);
		if (d > wh) wh = d;
	}
	printf("u_to_hz_q16: max err %.5f Hz\n", wh);
	CHECK(wh < 0.005, "u_to_hz_q16 error %.5f Hz", wh);

	/* u_to_q29 is round(u 2^29 / 6400) */
	bad = 0;
	for (int32_t u = -6464; u <= 6464; u++)
		bad += u_to_q29(u) != (int32_t)llround(u * 536870912.0 / 6400);
	CHECK(bad == 0, "u_to_q29: %d wrong", bad);

	/* energy quantiser: agrees with the double formula except within 2^-13 index units of a decision level */
	bad = 0;
	int near = 0;
	for (int32_t l = -12 * 65536; l <= 24 * 65536; l += 7)
	{
		double lg = l / 65536.0, x = 0.64 * 10 * log10(2) * lg + 6.9;
		int ref = dbl_encode_energy(pow(2, lg));
		if (encode_energy_log2(l) != ref)
		{
			if (fabs(x - floor(x + 0.5)) * 65536 < 8) near++;
			else bad++;
		}
	}
	printf("encode_energy_log2: %d off-boundary mismatches, %d within 2^-13 of a level\n", bad, near);
	CHECK(bad == 0, "encode_energy_log2");
	CHECK(encode_energy_log2(LOG2E_ZERO) == 0, "LOG2E_ZERO index");

	/* LSP quantiser: identical to the scalar reference on Q16 inputs */
	bad = 0;
	for (int i = 0; i < 200000; i++)
	{
		int32_t q[LPC_ORD];
		double h[LPC_ORD];
		int32_t acc = 0;
		for (int k = 0; k < LPC_ORD; k++)
		{
			acc += (int32_t)(rnd() % (400 << 16));
			if (i & 1) /* hit codebook midpoints exactly */
				acc = (acc >> 15) << 15;
			q[k] = acc;
			h[k] = acc / 65536.0;
		}
		int a[LPC_ORD], b[LPC_ORD];
		encode_lspds_q16(a, q);
		dbl_encode_lspds(b, h);
		for (int k = 0; k < LPC_ORD; k++)
			bad += a[k] != b[k];
	}
	CHECK(bad == 0, "encode_lspds_q16: %d index mismatches", bad);

	/* benign LSPs (zero energy / failed search): the integer 400 i Hz quantise as the float pi i / 10 */
	{
		int32_t q[LPC_ORD];
		float f[LPC_ORD];
		int a[LPC_ORD], b[LPC_ORD];
		benign_lsps(q);
		for (int i = 0; i < LPC_ORD; i++)
			f[i] = 4000.0f / (float)M_PI * (((float)M_PI / LPC_ORD) * (float)i);
		encode_lspds_q16(a, q);
		flt_encode_lspds(b, f);
		bad = 0;
		for (int i = 0; i < LPC_ORD; i++)
			bad += a[i] != b[i];
		CHECK(bad == 0, "benign LSP indices differ from the float");
	}
}

/* Schur vs double Levinson on autocorrelations of windowed AR(2) noise (some of it clipped) at levels
 * 2^0 .. 2^14, with float32 Levinson on the same r as the yardstick: the integer rc must be at least as close
 * to double as float's, in max and rms. (Absolute rc errors grow with the conditioning of r: ill-conditioned
 * cases reach ~2e-5 for both.) */
static void test_schur(void)
{
	double wa = 0, wf = 0, sa = 0, sf = 0;
	long nk = 0;
	for (int t = 0; t < 20000; t++)
	{
		int16_t x[NW];
		double st[2] = {0};
		double c2 = -(0.5 + 0.49 * (rnd() & 0xffff) / 65536.0), c1 = ((rnd() & 0xffff) / 32768.0 - 1) * 0.98 * (1 - c2); /* stable */
		double amp = pow(2, rnd() % 15);
		for (int i = 0; i < NW; i++)
		{
			double v = rnds(16) / 32768.0 + c1 * st[0] + c2 * st[1];
			st[1] = st[0]; st[0] = v;
			double s = v * amp * (0.5 - 0.5 * cos(2 * M_PI * i / (NW - 1)));
			x[i] = (int16_t)(s > 32767 ? 32767 : s < -32767 ? -32767 : s);
		}
		int64_t R[LPC_ORD + 1];
		double Rd[LPC_ORD + 1];
		for (int j = 0; j <= LPC_ORD; j++)
		{
			R[j] = 0;
			for (int i = 0; i < NW - j; i++) R[j] += x[i] * x[i + j];
			Rd[j] = (double)R[j];
		}
		if (R[0] < 1000) continue;
		int sh = bitlen64((uint64_t)R[0]) - 30;
		int32_t r[LPC_ORD + 1], rc[LPC_ORD];
		for (int j = 0; j <= LPC_ORD; j++)
			r[j] = (int32_t)(sh > 0 ? ((R[j] >> (sh - 1)) + 1) >> 1 : R[j] * ((int64_t)1 << -sh));
		int esh;
		int32_t em = schur(r, rc, &esh);
		/* double Levinson on the same normalised r */
		float af[LPC_ORD + 1], apf[LPC_ORD + 1], errf = (float)r[0];
		double a[LPC_ORD + 1], ap[LPC_ORD + 1], err = r[0];
		for (int i = 1; i <= LPC_ORD; i++)
		{
			double sum = 0;
			for (int j = 1; j < i; j++) sum += ap[j] * r[i - j];
			double k = -(r[i] + sum) / err;
			if (fabs(k) > 1) k = 0;
			a[i] = k;
			for (int j = 1; j < i; j++) a[j] = ap[j] + k * ap[i - j];
			err *= 1 - k * k;
			for (int j = 1; j <= i; j++) ap[j] = a[j];
			float sumf = 0;
			for (int j = 1; j < i; j++) sumf += apf[j] * (float)r[i - j];
			float kf = -((float)r[i] + sumf) / errf;
			if (fabsf(kf) > 1) kf = 0;
			af[i] = kf;
			for (int j = 1; j < i; j++) af[j] = apf[j] + kf * apf[i - j];
			errf *= 1 - kf * kf;
			for (int j = 1; j <= i; j++) apf[j] = af[j];
			double d = fabs(rc[i - 1] / 1073741824.0 - k), df = fabs(kf - k);
			if (d > wa) wa = d;
			if (df > wf) wf = df;
			sa += d * d;
			sf += df * df;
			nk++;
		}
		(void)em; (void)esh; (void)Rd;
	}
	printf("schur: |rc - k| max %.2e rms %.2e; float Levinson max %.2e rms %.2e\n", wa, sqrt(sa / nk), wf, sqrt(sf / nk));
	CHECK(wa <= wf && sa <= sf, "schur less accurate than float Levinson");
}

typedef struct
{
	long n, e_id, e_fd, lsp_n, lsp_id, lsp_fd, frm_id, frm_fd, roots_d;
	double hz_max, hz_sq, le_max, le_gain, hz_gain;
	int hz_k, hz_rd;
	long hz_n;
} stats_t;

static FILE *dumpf;

static void run_file(const int16_t *x, long len, stats_t *s)
{
	static codec2_encoder_t c2;
	c2fxp_encoder_init(&c2);
	static float wf[NW];
	static double wd[NW];
	static int init;
	if (!init) { flt_window(wf); dbl_window(wd); init = 1; }
	for (long f = 0; f + 2 * N_SAMP <= len; f += 2 * N_SAMP)
	{
		enc_frame_t fr;
		analyse_one_frame(&c2, &fr, x + f);
		analyse_one_frame(&c2, &fr, x + f + N_SAMP);
		int32_t hq[LPC_ORD];
		int32_t lE = speech_to_lsps_q16(&c2, hq);
		float hf[LPC_ORD], Ef;
		double hd[LPC_ORD], Ed;
		flt_analyse(c2.Sn, wf, hf, &Ef);
		int rd = dbl_analyse(c2.Sn, wd, hd, &Ed);
		int ii[LPC_ORD], fi[LPC_ORD], di[LPC_ORD];
		encode_lspds_q16(ii, hq);
		flt_encode_lspds(fi, hf);
		dbl_encode_lspds(di, hd);
		int ei = encode_energy_log2(lE), ef = flt_encode_energy(Ef), ed = dbl_encode_energy(Ed);
		if (dumpf)
		{
			int16_t rec[33];
			rec[0] = (int16_t)ei, rec[11] = (int16_t)ef, rec[22] = (int16_t)ed;
			for (int k = 0; k < LPC_ORD; k++)
				rec[1 + k] = (int16_t)ii[k], rec[12 + k] = (int16_t)fi[k], rec[23 + k] = (int16_t)di[k];
			fwrite(rec, 2, 33, dumpf);
		}
		s->n++;
		s->e_id += ei != ed;
		s->e_fd += ef != ed;
		int fid = ei != ed, ffd = ef != ed;
		for (int k = 0; k < LPC_ORD; k++)
		{
			s->lsp_n++;
			s->lsp_id += ii[k] != di[k];
			s->lsp_fd += fi[k] != di[k];
			fid |= ii[k] != di[k];
			ffd |= fi[k] != di[k];
		}
		s->frm_id += fid;
		s->frm_fd += ffd;
		if (rd != LPC_ORD)
			s->roots_d++;
		else
		{
			for (int k = 0; k < LPC_ORD; k++)
			{
				double d = fabs(hq[k] / 65536.0 - hd[k]);
				if (d > s->hz_max) { s->hz_max = d; s->hz_k = k; s->hz_gain = dbl_last_r0 / Ed; }
				s->hz_sq += d * d;
				s->hz_n++;
			}
			double le = fabs(lE / 65536.0 - log2(Ed));
			if (le > s->le_max) { s->le_max = le; s->le_gain = dbl_last_r0 / Ed; }
		}
	}
}

int main(int argc, char **argv)
{
	test_helpers();
	test_schur();
	if (argc > 1)
	{
		stats_t s = {0};
		for (int a = 1; a < argc; a++)
		{
			FILE *fp = fopen(argv[a], "rb");
			if (!fp) continue;
			fseek(fp, 0, SEEK_END);
			long len = ftell(fp) / 2;
			fseek(fp, 0, SEEK_SET);
			int16_t *x = malloc(len * 2);
			if (fread(x, 2, len, fp) != (size_t)len) len = 0;
			fclose(fp);
			if (getenv("LPC_DUMP"))
			{
				char path[1024];
				snprintf(path, sizeof path, "%s/%d.idx", getenv("LPC_DUMP"), a);
				dumpf = fopen(path, "wb");
			}
			run_file(x, len, &s);
			if (dumpf)
				fclose(dumpf), dumpf = NULL;
			free(x);
		}
		printf("frames %ld (double trap/failed search %ld)\n", s.n, s.roots_d);
		printf("  energy index != double:  int %.3f%%  float %.3f%%\n", 100.0 * s.e_id / s.n, 100.0 * s.e_fd / s.n);
		printf("  LSP index != double:     int %.3f%%  float %.3f%% (of %ld indices)\n", 100.0 * s.lsp_id / s.lsp_n,
		       100.0 * s.lsp_fd / s.lsp_n, s.lsp_n);
		printf("  frame packet != double:  int %.3f%%  float %.3f%%\n", 100.0 * s.frm_id / s.n, 100.0 * s.frm_fd / s.n);
		printf("  LSP |int - double|: max %.4f Hz, rms %.5f Hz; log2 E max err %.5f\n", s.hz_max,
		       sqrt(s.hz_sq / (s.hz_n ? s.hz_n : 1)), s.le_max);
		printf("  (worst LSP: index %d, R0 / E %.3g; worst log2 E: R0 / E %.3g)\n", s.hz_k, s.hz_gain, s.le_gain);
		CHECK(s.lsp_id <= s.lsp_fd && s.e_id <= s.e_fd, "integer agrees with double less often than float");
	}
	printf(fails ? "test_lpc: %d FAILURES\n" : "test_lpc: all passed\n", fails);
	return fails != 0;
}
