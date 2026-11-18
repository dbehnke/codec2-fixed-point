/*
 * test_nlp.c — host checks of the fixed-point NLP pitch estimator (src/nlp.c, docs/10).
 * Includes nlp.c to reach its static helpers. Run by scripts/test-fxp.sh (FXP_CHECK + UBSan + ASan).
 *  1. threshold helpers: frac3_10 / frac3_20 equal floor(3g/10), floor(3g/20); 4b/5, 6b/5 equal the float
 *     (int)(0.8f b), (int)(1.2f b)
 *  2. the decimating FIR is exact: equals round(sum h sq / 2^15) in 64-bit for random and extreme inputs
 *  3. the whole estimator against a double-precision reference (same algorithm, exact DFT) and against the
 *     float reference (the E1 float code: float arithmetic, kiss_fftr): bin agreement on synthetic speech-like
 *     signals over a 60 dB level range. With file arguments (raw s16 8 kHz), the same comparison on those files.
 */
#include "c2_float_consts.h"
#include "../src/nlp.c"
#include "kiss_fftr.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int fails;
#define CHECK(c, ...) do { if (!(c)) { fails++; printf("FAIL: " __VA_ARGS__); printf("\n"); } } while (0)

static uint32_t rng = 777;
static uint32_t rnd32(void) { rng = rng * 1103515245u + 12345u; return rng >> 1; }

/* float polyphase taps of the E1 code, time order h[p + 5 t] */
static const double fir_ref[50] = { /* h[48] = h[49] = 0 pad the 5 x 10 bank */
	-0.00108181240, -0.00110083440, -0.00092768838, -0.00042289438, 0.00055034190, 0.00200298490,
	0.00370585090, 0.00514494150, 0.00559246660, 0.00430367540, 0.00080284511, -0.00482046100,
	-0.01170581000, -0.01819927500, -0.02206528200, -0.02092061000, -0.01280883100, 0.00322047750,
	0.02668381100, 0.05552062400, 0.08630594400, 0.11480192000, 0.13674206000, 0.14867556000,
	0.14867556000, 0.13674206000, 0.11480192000, 0.08630594400, 0.05552062400, 0.02668381100,
	0.00322047750, -0.01280883100, -0.02092061000, -0.02206528200, -0.01819927500, -0.01170581000,
	-0.00482046100, 0.00080284511, 0.00430367540, 0.00559246660, 0.00514494150, 0.00370585090,
	0.00200298490, 0.00055034190, -0.00042289438, -0.00092768838, -0.00110083440, -0.00108181240,
	0.0, 0.0};

/* reference NLP in double (exact DFT) or float (kiss_fftr), same structure as the E1 float code */
typedef struct
{
	double sq[M_PITCH], sq_fir[NDEC], mem_x, mem_y;
	int prev;
	int isfloat;
	float fsq[M_PITCH], fsq_fir[NDEC], fmem_x, fmem_y;
} ref_t;

static int sub_multiples_ref(const double *Fw, double gmax, int gmax_bin, int prev, int isfloat)
{
	int cmax_bin = gmax_bin;
	for (int mult = 2; gmax_bin / mult >= 16; mult++)
	{
		int b = gmax_bin / mult;
		int bmin = (int)(0.8f * b), bmax = (int)(1.2f * b);
		if (bmin < 16)
			bmin = 16;
		double thresh = isfloat ? (prev > bmin && prev < bmax ? (double)(0.3f * 0.5f * (float)gmax) : (double)(0.3f * (float)gmax))
					: (prev > bmin && prev < bmax ? 0.15 * gmax : 0.3 * gmax);
		double lmax = 0;
		int lmax_bin = bmin;
		for (int k = bmin; k <= bmax; k++)
			if (Fw[k] > lmax)
			{
				lmax = Fw[k];
				lmax_bin = k;
			}
		if (lmax > thresh && lmax > Fw[lmax_bin - 1] && lmax > Fw[lmax_bin + 1])
			cmax_bin = lmax_bin;
	}
	return cmax_bin;
}

static double whann(int i) { return 0.5 - 0.5 * cos(2 * M_PI * i / 63); }

static int nlp_ref(ref_t *r, const int16_t *Sn, double *Fw_out)
{
	double Fw[PE_FFT_SIZE / 2 + 1];
	if (!r->isfloat)
	{
		for (int i = START_POS; i < M_PITCH; i++)
		{
			double x = (double)Sn[i] * Sn[i];
			double y = x - r->mem_x + 0.95 * r->mem_y;
			r->mem_x = x;
			r->mem_y = y;
			r->sq[i] = y + 1.0;
		}
		memmove(r->sq_fir, r->sq_fir + 16, 48 * sizeof(double));
		for (int i = 0; i < 16; i++)
		{
			double acc = 0;
			for (int k = 0; k < 48; k++)
				acc += fir_ref[k] * r->sq[START_POS + 5 * i - k];
			r->sq_fir[48 + i] = acc;
		}
		double xw[NDEC];
		for (int i = 0; i < NDEC; i++)
			xw[i] = r->sq_fir[i] * whann(i);
		for (int k = 0; k <= 129; k++)
		{
			double re = 0, im = 0;
			for (int i = 0; i < NDEC; i++)
			{
				re += xw[i] * cos(2 * M_PI * k * i / 512);
				im -= xw[i] * sin(2 * M_PI * k * i / 512);
			}
			Fw[k] = re * re + im * im;
		}
		memmove(r->sq, r->sq + N_SAMP, START_POS * sizeof(double));
	}
	else
	{
		static uint8_t mem[8192];
		static kiss_fftr_cfg cfg;
		if (!cfg)
		{
			size_t len = sizeof(mem);
			cfg = kiss_fftr_alloc(PE_FFT_SIZE, 0, mem, &len);
		}
		for (int i = START_POS; i < M_PITCH; i++)
		{
			float x = Sn[i];
			r->fsq[i] = x * x;
		}
		for (int i = START_POS; i < M_PITCH; i++)
		{
			float notch = r->fsq[i] - r->fmem_x;
			notch += 0.95f * r->fmem_y;
			r->fmem_x = r->fsq[i];
			r->fmem_y = notch;
			r->fsq[i] = notch + 1.0f;
		}
		memmove(r->fsq_fir, r->fsq_fir + 16, 48 * sizeof(float));
		for (int i = 0; i < 16; i++)
		{
			float acc = 0;
			for (int p = 0; p < 5; p++) /* same summation order as the polyphase float code */
				for (int t = 0; t < 10; t++)
					acc += (float)fir_ref[p + 5 * t] * r->fsq[START_POS + 5 * i - p - 5 * t];
			r->fsq_fir[48 + i] = acc;
		}
		float buf[PE_FFT_SIZE] = {0};
		kiss_fft_cpx X[PE_FFT_SIZE / 2 + 1];
		for (int i = 0; i < NDEC; i++)
			buf[i] = r->fsq_fir[i] * (0.5f - 0.5f * cosf(2 * (float)M_PI * i / (NDEC - 1)));
		kiss_fftr(cfg, buf, X);
		for (int k = 0; k <= 129; k++)
			Fw[k] = X[k].r * X[k].r + X[k].i * X[k].i;
		for (int i = 0; i < START_POS; i++)
			r->fsq[i] = r->fsq[i + N_SAMP];
	}
	double gmax = 0;
	int gmax_bin = 16;
	for (int k = 16; k <= 128; k++)
		if (Fw[k] > gmax)
		{
			gmax = Fw[k];
			gmax_bin = k;
		}
	int b = sub_multiples_ref(Fw, gmax, gmax_bin, r->prev, r->isfloat);
	r->prev = b;
	if (Fw_out)
		memcpy(Fw_out, Fw, sizeof(Fw));
	return b;
}

static void test_helpers(void)
{
	for (uint32_t t = 0; t < 2000000; t++)
	{
		uint32_t g = t < 1000 ? t : t < 2000 ? 0xffffffffu - t : rnd32() * 2 + (t & 1);
		uint64_t a = 3ull * g / 10, b = 3ull * g / 20;
		CHECK(frac3_10(g) == a && frac3_20(g) == b, "frac3 %u", g);
		if (fails > 5)
			return;
	}
	for (int b = 1; b <= 64; b++)
		CHECK(4 * b / 5 == (int)(0.8f * b) && 6 * b / 5 == (int)(1.2f * b), "bmin/bmax b=%d", b);
	printf("threshold helpers: floor(3g/10), floor(3g/20) exact on 2M values; 4b/5, 6b/5 equal the float casts\n");
}

static void test_fir(void)
{
	static nlp_t st;
	static int32_t scratch[NLP_SCRATCH_WORDS];
	long n = 0, nse = 0;
	for (int trial = 0; trial < 6000; trial++)
	{
		/* random history in the full notch range, the extremes, and tiny values (block exponent > 0) */
		int kind = trial % 5;
		for (int j = 0; j < NLP_HIST; j++)
		{
			uint32_t u = rnd32();
			int32_t v = kind == 0 ? (int32_t)(u % 2 ? (1 << 30) + 3 : -(1 << 30) - 3)
				  : kind == 1 ? (int32_t)(u % 7) - 3
				  : kind == 2 ? (int32_t)(u % 4001) - 2000
				  : kind == 3 ? (int32_t)(u % (1u << (u % 31 + 1))) - (int32_t)(1u << (u % 31)) / 2
				  : (int32_t)(u % ((1u << 31) + 7u)) - (1 << 30) - 3;
			st.sq[j] = v;
		}
		int64_t Tq[16];
		uint32_t ora = 0;
		for (int i = 0; i < 16; i++)
		{
			Tq[i] = 0;
			for (int k = 0; k < NLP_NTAP; k++)
				Tq[i] += (int64_t)fxp_nlp_fir_q15[k] * st.sq[NLP_NTAP - 1 + 5 * i - k];
			int64_t fl = Tq[i] >> 15;
			ora |= (uint32_t)(fl < 0 ? -fl : fl);
		}
		int se = 30 - (ora ? 32 - __builtin_clz(ora) : 0);
		se = se < 0 ? 0 : se > 15 ? 15 : se;
		lowpass(&st, scratch);
		CHECK(st.fir_e[NLP_NBLK - 1] == se, "fir trial %d: exponent %d vs %d", trial, st.fir_e[NLP_NBLK - 1], se);
		for (int i = 0; i < 16; i++)
		{
			int64_t ref = ((Tq[i] * ((int64_t)1 << se)) + (se < 15 ? 1 << 14 : 0)) >> 15;
			CHECK(st.sq_fir[NDEC - 16 + i] == ref, "fir trial %d i %d: %d vs %lld", trial, i, st.sq_fir[NDEC - 16 + i], (long long)ref);
			n++;
		}
		nse += se > 0;
		if (fails > 5)
			return;
	}
	printf("decimating FIR: %ld outputs equal round(sum h sq 2^(e - 15)) (64-bit), block exponent e > 0 in %ld of 6000 blocks\n", n, nse);
}

typedef struct
{
	long n, agree_int, agree_flt, agree_int_flt, oct_int, oct_flt;
} stats_t;

static void run_signal(const int16_t *x, long len, stats_t *s)
{
	static nlp_t st;
	static ref_t rd, rf;
	nlp_init(&st);
	memset(&rd, 0, sizeof(rd));
	memset(&rf, 0, sizeof(rf));
	rd.prev = rf.prev = 16;
	rf.isfloat = 1;
	int prev = 16;
	int16_t Sn[M_PITCH];
	for (int i = 0; i < M_PITCH; i++)
		Sn[i] = 1;
	static int32_t scratch[NLP_SCRATCH_WORDS];
	for (long f = 0; f + N_SAMP <= len; f += N_SAMP)
	{
		memmove(Sn, Sn + N_SAMP, START_POS * sizeof(int16_t));
		memcpy(Sn + START_POS, x + f, N_SAMP * sizeof(int16_t));
		int bi = nlp(&st, Sn, &prev, scratch);
		int bd = nlp_ref(&rd, Sn, NULL);
		int bf = nlp_ref(&rf, Sn, NULL);
		s->n++;
		s->agree_int += bi == bd;
		s->agree_flt += bf == bd;
		s->agree_int_flt += bi == bf;
		s->oct_int += abs(bi - bd) > 2;
		s->oct_flt += abs(bf - bd) > 2;
	}
}

static void show(const char *lbl, const stats_t *s)
{
	printf("  %-22s frames %6ld  int=dbl %6.2f%%  flt=dbl %6.2f%%  int=flt %6.2f%%  |int-dbl|>2 %5.2f%%  |flt-dbl|>2 %5.2f%%\n",
	       lbl, s->n, 100.0 * s->agree_int / s->n, 100.0 * s->agree_flt / s->n, 100.0 * s->agree_int_flt / s->n,
	       100.0 * s->oct_int / s->n, 100.0 * s->oct_flt / s->n);
}

static void test_synthetic(void)
{
	enum { LEN = 8000 * 6 };
	static int16_t x[LEN];
	const double levels[] = {30000, 3000, 300, 30, 3, -1}; /* -1: steps between 30000 and 3 every 1930 samples */
	printf("NLP bin vs double reference, synthetic voiced/unvoiced signal (6 s per level; steps exercise the notch rescale):\n");
	for (int li = 0; li < 6; li++)
	{
		stats_t s = {0};
		double ph = 0, f0 = 90;
		for (long n = 0; n < LEN; n++)
		{
			/* gliding f0 90..260 Hz, 3 harmonics with decaying amplitude, voiced/noise alternating every 0.5 s */
			f0 = 90 + 85 * (1 + sin(2 * M_PI * n / 16000.0));
			ph += 2 * M_PI * f0 / 8000;
			int voiced = (n / 4000) % 2 == 0;
			double v = voiced ? sin(ph) + 0.5 * sin(2 * ph + 0.3) + 0.25 * sin(3 * ph + 1.1)
					  : ((double)(rnd32() % 2001) - 1000) / 1000.0;
			v += 0.02 * (((double)(rnd32() % 2001) - 1000) / 1000.0);
			double lv = levels[li] > 0 ? levels[li] : (n / 1930) % 2 ? 3 : 30000;
			if (levels[li] < 0 && lv > 3 && voiced)
				v = pow(v / 1.75, 7) * 1.8; /* peaky (high crest factor), so the FIR history dominates the rescale */
			double a = lv * v / 1.8;
			x[n] = (int16_t)fmax(-32768, fmin(32767, lrint(a)));
		}
		run_signal(x, LEN, &s);
		char lbl[32];
		if (levels[li] > 0)
			snprintf(lbl, sizeof(lbl), "peak %5.0f", levels[li]);
		else
			snprintf(lbl, sizeof(lbl), "steps 30000 <-> 3");
		show(lbl, &s);
		CHECK(s.agree_int >= s.n * (levels[li] < 0 ? 995 : 980) / 1000, "synthetic level %.0f: int agrees with double on %ld of %ld", levels[li], s.agree_int, s.n);
	}
}

/* preconditions of the exact rearrangements in nlp(): bin pairing (k, 128 - k) and paired FIR taps */
static void test_tables(void)
{
	int bad = 0;
	for (int m = 0; m < 512; m++)
	{
		int c = fxp_cos512_q14[m];
		bad += fxp_cos512_q14[(512 - m) & 511] != c;  /* C(-m) = C(m) */
		bad += fxp_cos512_q14[(m + 256) & 511] != -c; /* C(m + 256) = -C(m) */
		bad += (int16_t)fxp_nlp_cs512[m] != c;
		bad += (int16_t)(fxp_nlp_cs512[m] >> 16) != fxp_cos512_q14[(m + 384) & 511]; /* S(m) = C(m - 128) */
	}
	CHECK(bad == 0, "twiddle table symmetry / packing: %d violations", bad);
	for (int k = 0; k < NLP_NTAP; k++)
		CHECK(fxp_nlp_fir_q15[k] == fxp_nlp_fir_q15[NLP_NTAP - 1 - k], "FIR tap %d not symmetric", k);
	CHECK(fxp_nlp_hann_q15[0] == 0 && fxp_nlp_hann_q15[NDEC - 1] == 0, "NLP window end points not zero");
}

int main(int argc, char **argv)
{
	test_tables();
	test_helpers();
	test_fir();
	test_synthetic();
	if (argc > 1)
	{
		stats_t tot = {0};
		for (int a = 1; a < argc; a++)
		{
			FILE *fp = fopen(argv[a], "rb");
			if (!fp)
				continue;
			fseek(fp, 0, SEEK_END);
			long len = ftell(fp) / 2;
			fseek(fp, 0, SEEK_SET);
			int16_t *x = malloc(len * 2);
			if (fread(x, 2, len, fp) != (size_t)len)
				len = 0;
			fclose(fp);
			stats_t s = {0};
			run_signal(x, len, &s);
			free(x);
			tot.n += s.n; tot.agree_int += s.agree_int; tot.agree_flt += s.agree_flt;
			tot.agree_int_flt += s.agree_int_flt; tot.oct_int += s.oct_int; tot.oct_flt += s.oct_flt;
		}
		show("files", &tot);
	}
	printf(fails ? "test_nlp: %d FAILURES\n" : "test_nlp: all passed\n", fails);
	return fails != 0;
}
