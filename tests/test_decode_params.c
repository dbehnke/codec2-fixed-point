/* test_decode_params.c — integer decoder parameters (docs/07-fxp-decode-params.md) against the float
 * decoder's formulas (Codec2-mod split: decode_Wo, decode_energy, decode_lspds_scalar, interp_Wo,
 * interp_energy, interpolate_lsp), over every index and every interpolation case. */
#include "c2_float_consts.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "codec2_internal.h"
#include "quantise.h"
#include "interp.h"

static int fails;
#define CHECK(c, ...) do { if (!(c)) { fails++; if (fails < 20) { printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } } } while (0)

/* float reference, as the float decoder */
static float f_decode_Wo(int index) { return W0_MIN + ((W0_MAX - W0_MIN) / (1 << WO_BITS)) * index; }
static float f_decode_energy(int index) { return powf(10.0f, (E_MIN_DB + ((float)(E_MAX_DB - E_MIN_DB) / (1 << E_BITS)) * index) / 10.0f); }

int main(void)
{
	/* Wo: every decoded value, W0_MIN, and every interpolation (both orders), plus L = pi / Wo */
	int nw = 0, nl = 0;
	double wmax = 0, wbmax = 0;
	for (int a = -1; a < (1 << WO_BITS); a++)
		for (int b = 0; b < (1 << WO_BITS); b++)
			for (int k = 0; k < 3; k++)
			{
				float fa = a < 0 ? W0_MIN : f_decode_Wo(a), fb = f_decode_Wo(b);
				int qa = a < 0 ? WO_Q_MIN : decode_Wo_q(a), qb = decode_Wo_q(b);
				float fo = k == 0 ? fa : (k == 1 ? fa + 0.5f * (fb - fa) : fb);
				int v = 1, q = k == 0 ? qa : (k == 1 ? interp_Wo_q(&v, 1, qa, 1, qb) : qb);
				double exact = 2 * M_PI * q / 40960.0;
				uint32_t W = wo_q_to_angle(q);
				double dw = fabs(W - exact * 4294967296.0 / (2 * M_PI)); /* LSB of W */
				double db = fabs((double)fo - exact);                    /* float decoder vs exact */
				wmax = fmax(wmax, dw);
				wbmax = fmax(wbmax, db);
				CHECK(dw <= 0.5, "W rounding q %d", q);
				CHECK(db < 1e-6, "float Wo far from wq grid: q %d fo %.9g", q, fo);
				int Lf = (int)((float)M_PI / fo), Li = 20480 / q;
				CHECK(Lf == Li, "L q %d float %d int %d", q, Lf, Li);
				nl += Lf == Li;
				nw++;
			}
	printf("Wo: %d values, |W - exact| max %.3f LSB, float decoder Wo vs exact grid max %.2e rad, L equal %d/%d\n",
	       nw, wmax, wbmax, nl, nw);

	/* interp_Wo voicing rules */
	for (int pv = 0; pv < 2; pv++)
		for (int v = 0; v < 2; v++)
			for (int nv = 0; nv < 2; nv++)
			{
				int vi = v, q = interp_Wo_q(&vi, pv, 300, nv, 500);
				int ve = (v && !pv && !nv) ? 0 : v;
				int qe = !ve ? WO_Q_MIN : (pv && nv ? 400 : (pv ? 300 : 500));
				CHECK(vi == ve && q == qe, "interp voicing pv %d v %d nv %d", pv, v, nv);
			}

	/* energy: log2 E vs log2 of the float decoder's E, and interpolation vs sqrtf */
	double emax = 0, eimax = 0;
	for (int a = 0; a < (1 << E_BITS); a++)
	{
		int32_t la = decode_energy_log2_q16(a);
		double ex = (E_MIN_DB + 50.0 / 32 * a) / 10 * log2(10.0);
		emax = fmax(emax, fabs(la / 65536.0 - ex) * 65536);
		CHECK(fabs(la / 65536.0 - log2(f_decode_energy(a))) < 3e-5, "energy index %d", a);
		for (int b = 0; b < (1 << E_BITS); b++)
		{
			int32_t li = interp_energy_log2_q16(la, decode_energy_log2_q16(b));
			double fi = log2(sqrtf(f_decode_energy(a) * f_decode_energy(b)));
			eimax = fmax(eimax, fabs(li / 65536.0 - fi) * 65536);
		}
	}
	printf("energy: log2 E vs exact max %.3f LSB (Q16); interpolated vs float sqrt max %.3f LSB\n", emax, eimax);
	CHECK(emax <= 0.6 && eimax <= 4, "energy accuracy");

	/* LSPs: half-Hz vs the float decoder's radians, decoded and interpolated (random indexes) */
	srand(1);
	int16_t prev[LPC_ORD];
	float prevf[LPC_ORD];
	for (int i = 0; i < LPC_ORD; i++)
	{
		prev[i] = (int16_t)(2 * ((2 * 4000 * i + LPC_ORD + 1) / (2 * (LPC_ORD + 1)))); /* as codec2_decoder_init */
		prevf[i] = i * (float)M_PI / (LPC_ORD + 1);
	}
	int nbad = 0, ntot = 0;
	for (int f = 0; f < 200000; f++)
	{
		int idx[LPC_ORD];
		for (int i = 0; i < LPC_ORD; i++)
			idx[i] = rand() & 31;
		int16_t h[LPC_ORD], hi[LPC_ORD];
		float lf[LPC_ORD], lif[LPC_ORD];
		decode_lspds_hz2(h, idx);
		int hz = 0;
		for (int i = 0; i < LPC_ORD; i++)
		{
			hz += delta_lsp_cb[i][idx[i]];
			lf[i] = ((float)M_PI / 4000.0f) * (float)hz;
			lif[i] = prevf[i] + 0.5f * (lf[i] - prevf[i]);
		}
		interpolate_lsp_hz2(hi, prev, h);
		for (int i = 0; i < LPC_ORD; i++)
		{
			/* the previous bridge recovered half-Hz with lrintf(lsp 8000 / pi) */
			nbad += h[i] != (int16_t)lrintf(lf[i] * (8000.0f / (float)M_PI));
			if (hi[i] != (int16_t)lrintf(lif[i] * (8000.0f / (float)M_PI))) { nbad++; printf("frame %d i %d int %d float %.4f\n", f, i, hi[i], lif[i] * (8000.0 / M_PI)); }
			ntot += 2;
		}
		for (int i = 0; i < LPC_ORD; i++)
		{
			prev[i] = h[i];
			prevf[i] = lf[i];
		}
	}
	printf("LSPs: half-Hz differing from lrintf(float lsp 8000 / pi): %d of %d\n", nbad, ntot);
	CHECK(nbad == 0, "LSP half-Hz");

	if (fails)
	{
		printf("%d FAILURES\n", fails);
		return 1;
	}
	printf("all decode-parameter tests passed\n");
	return 0;
}
