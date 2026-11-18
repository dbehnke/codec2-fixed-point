/*
 * lpc_ref.h — reference copies of the float Codec2 encoder LPC/LSP path (speech_to_uq_lsps, encode_energy,
 * encode_lspds_scalar as in codec/src/lpc.c and quantise.c before stage E3), written once over a type.
 * Include with REAL = float (the E2 float implementation) or REAL = double (the same algorithm in double:
 * the reference for docs/11), and FN(name) giving distinct function names.
 */
#define W_END_R (W_START + NW)

/* Hann window times 1/sqrt(512 m), as the removed make_analysis_window */
static void FN(window)(REAL *w)
{
	REAL m = 0;
	for (int j = 0; j < NW; j++)
	{
		REAL v = (REAL)0.5 - (REAL)0.5 * (sizeof(REAL) == 4 ? (REAL)cosf((float)(TWO_PI * j / (NW - 1))) : (REAL)cos(2 * M_PI * j / (NW - 1)));
		w[j] = v;
		m += v * v;
	}
	REAL scale = (REAL)1 / (sizeof(REAL) == 4 ? (REAL)sqrtf((float)(m * FFT_ENC)) : (REAL)sqrt(m * FFT_ENC));
	for (int j = 0; j < NW; j++)
		w[j] *= scale;
}

static REAL FN(acos)(REAL x)
{
	if (x < -1) x = -1;
	if (x > 1) x = 1;
	int neg = x < 0;
	if (neg) x = -x;
	REAL p = (((((REAL)-0.004731162409361028 * x + (REAL)0.02013917916990988) * x - (REAL)0.04547268716508862) * x +
	           (REAL)0.08798969877894365) * x - (REAL)0.2145148619010314) * x + (REAL)1.57079461300486;
	REAL y = (sizeof(REAL) == 4 ? (REAL)sqrtf((float)(1 - x)) : (REAL)sqrt(1 - x)) * p;
	return neg ? (REAL)M_PI - y : y;
}

static REAL FN(cheb)(const REAL *c, REAL x)
{
	REAL T0 = 1, T1 = x, s = c[5] * T0 + c[4] * T1, tx = 2 * x;
	REAL T2 = tx * T1 - T0; s += c[3] * T2;
	REAL T3 = tx * T2 - T1; s += c[2] * T3;
	REAL T4 = tx * T3 - T2; s += c[1] * T4;
	REAL T5 = tx * T4 - T3; s += c[0] * T5;
	return s;
}

/* roots as x = cos w (not yet acos'ed) in xs; returns the root count */
static int FN(lpc_to_lsp)(const REAL *a, REAL *xs)
{
	REAL P[LPC_ORD + 1], Q[LPC_ORD + 1];
	P[0] = Q[0] = 1;
	for (int i = 1; i <= 5; i++)
	{
		P[i] = a[i] + a[LPC_ORD + 1 - i] - P[i - 1];
		Q[i] = a[i] - a[LPC_ORD + 1 - i] + Q[i - 1];
	}
	for (int i = 0; i < 5; i++) { P[i] *= 2; Q[i] *= 2; }
	int roots = 0;
	REAL xl = 1, xr, xm = 0, sl, sr, sm;
	for (int j = 0; j < LPC_ORD; j++)
	{
		const REAL *pt = (j & 1) ? Q : P;
		xr = xl;
		sl = FN(cheb)(pt, xl);
		while (xr >= -1)
		{
			xr = xl - (REAL)LSP_DELTA1;
			sr = FN(cheb)(pt, xr);
			if ((sr <= 0 && sl >= 0) || (sr >= 0 && sl <= 0))
			{
				roots++;
				for (int b = 0; b < 6; b++)
				{
					xm = (REAL)0.5 * (xl + xr);
					sm = FN(cheb)(pt, xm);
					if (sm * sl > 0) { xl = xm; sl = sm; } else xr = xm;
				}
				xs[j] = xm;
				xl = xm;
				break;
			}
			sl = sr;
			xl = xr;
		}
	}
	return roots;
}

static REAL FN(last_r0); /* R[0] of the last analyse() call */

/* LSPs in Hz and the energy E of Sn (M_PITCH samples); returns the root count (LPC_ORD + 1 for the trap) */
static int FN(analyse)(const int16_t *Sn, const REAL *w, REAL *lsp_hz, REAL *E)
{
	REAL Wn[NW], R[LPC_ORD + 1], a[LPC_ORD + 1], ap[LPC_ORD + 1], e = 0;
	for (int i = 0; i < NW; i++) { Wn[i] = Sn[W_START + i] * w[i]; e += Wn[i] * Wn[i]; }
	if (e < (REAL)LPC_ENERGY_FLOOR)
	{
		for (int i = 0; i < LPC_ORD; i++) lsp_hz[i] = (REAL)4000 / (REAL)M_PI * (((REAL)M_PI / LPC_ORD) * (REAL)i);
		*E = 0;
		return LPC_ORD + 1;
	}
	for (int j = 0; j <= LPC_ORD; j++)
	{
		REAL acc = 0;
		for (int i = 0; i < NW - j; i++) acc += Wn[i] * Wn[i + j];
		R[j] = acc;
	}
	FN(last_r0) = R[0];
	REAL err = R[0];
	a[0] = 1;
	for (int i = 1; i <= LPC_ORD; i++)
	{
		REAL sum = 0;
		for (int j = 1; j <= i - 1; j++) sum += ap[j] * R[i - j];
		REAL k = -(R[i] + sum) / err;
		if ((k < 0 ? -k : k) > 1) k = 0;
		a[i] = k;
		for (int j = 1; j <= i - 1; j++) a[j] = ap[j] + k * ap[i - j];
		err *= (1 - k * k);
		for (int j = 1; j <= i; j++) ap[j] = a[j];
	}
	REAL En = 0;
	for (int i = 0; i <= LPC_ORD; i++) En += a[i] * R[i];
	*E = En;
	static const double g[LPC_ORD + 1] = {1.0, 0.994, 0.988036, 0.982107784, 0.976215137296, 0.970357846472224,
	                                      0.964535699393391, 0.95874848519703, 0.952995994285848, 0.947278018320133, 0.941594350210212};
	for (int i = 0; i <= LPC_ORD; i++) a[i] *= (REAL)g[i];
	REAL xs[LPC_ORD];
	int roots = FN(lpc_to_lsp)(a, xs);
	for (int i = 0; i < LPC_ORD; i++)
		lsp_hz[i] = roots == LPC_ORD ? (REAL)4000 / (REAL)M_PI * FN(acos)(xs[i])
		                             : (REAL)4000 / (REAL)M_PI * (((REAL)M_PI / LPC_ORD) * (REAL)i);
	return roots;
}

static int FN(encode_energy)(REAL e)
{
	REAL m = e > (REAL)1e-12 ? e : (REAL)1e-12;
	REAL e_db = 10 * (sizeof(REAL) == 4 ? (REAL)log10f((float)m) : (REAL)log10(m));
	int index = (int)(32 * ((e_db - (REAL)E_MIN_DB) / ((REAL)E_MAX_DB - (REAL)E_MIN_DB)) + (REAL)0.5);
	return index < 0 ? 0 : index > 31 ? 31 : index;
}

static void FN(encode_lspds)(int *idx, const REAL *lsp_hz)
{
	int last = 0;
	for (int i = 0; i < LPC_ORD; i++)
	{
		const uint16_t *cb = delta_lsp_cb[i];
		REAL t = lsp_hz[i] - (REAL)last;
		int best;
		if (t <= (REAL)cb[0]) best = 0;
		else if (t >= (REAL)cb[31]) best = 31;
		else
		{
			int lo = 0, hi = 31;
			while (lo + 1 < hi) { int mid = (lo + hi) >> 1; if ((REAL)cb[mid] <= t) lo = mid; else hi = mid; }
			REAL el = t - (REAL)cb[lo], eh = (REAL)cb[hi] - t;
			best = el <= eh ? lo : hi;
		}
		idx[i] = best;
		last += cb[best];
	}
}
#undef W_END_R
