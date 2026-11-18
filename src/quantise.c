#include "codec2_internal.h"
#include "quantise.h"
#include "lpc.h"

static inline unsigned int binary_to_gray(unsigned int x)
{
	return x ^ (x >> 1);
}

static inline unsigned int gray_to_binary(unsigned int g)
{
	g ^= g >> 16;
	g ^= g >> 8;
	g ^= g >> 4;
	g ^= g >> 2;
	g ^= g >> 1;
	return g;
}

void pack(
	uint8_t *bitArray,		/* The output bit array. */
	unsigned int *bitIndex, /* Index into the string in BITS, not bytes. */
	unsigned int field,		/* The bit field to be packed. */
	unsigned int fieldWidth /* Width of the field in BITS, not bytes. */
)
{
	/* Convert to Gray code if multi-bit */
	if (fieldWidth > 1)
		field = binary_to_gray(field);

	do
	{
		unsigned int bI = *bitIndex;
		unsigned int bitsLeft = 8 - (bI & 0x7);
		unsigned int sliceWidth = bitsLeft < fieldWidth ? bitsLeft : fieldWidth;
		unsigned int wordIndex = bI >> 3;

		bitArray[wordIndex] |=
			(uint8_t)((field >> (fieldWidth - sliceWidth)) << (bitsLeft - sliceWidth));

		*bitIndex = bI + sliceWidth;
		fieldWidth -= sliceWidth;

	} while (fieldWidth != 0);
}

int unpack(
	const uint8_t *bitArray, /* The input bit string. */
	unsigned int *bitIndex,	 /* Index into the string in BITS, not bytes. */
	unsigned int fieldWidth	 /* Width of the field in BITS, not bytes. */
)
{
	static const uint8_t mask8[9] = {
		0x00, 0x01, 0x03, 0x07, 0x0F,
		0x1F, 0x3F, 0x7F, 0xFF};

	unsigned int field = 0;
	unsigned int origWidth = fieldWidth;

	do
	{
		unsigned int bI = *bitIndex;
		unsigned int bitsLeft = 8 - (bI & 0x7);
		unsigned int sliceWidth = bitsLeft < fieldWidth ? bitsLeft : fieldWidth;

		field |= (((bitArray[bI >> 3] >> (bitsLeft - sliceWidth)) &
				   mask8[sliceWidth])
				  << (fieldWidth - sliceWidth));

		*bitIndex = bI + sliceWidth;
		fieldWidth -= sliceWidth;
	} while (fieldWidth != 0);

	return (origWidth > 1) ? gray_to_binary(field) : field;
}

/* Decoder Wo as the integer wq (Wo = 2 pi wq / 40960, see WO_Q_*). Same value as the float
 * W0_MIN + step * index: step = (W0_MAX - W0_MIN) / 128 = 2 pi 14 / 40960. */
int decode_Wo_q(int index)
{
	return WO_Q_MIN + WO_Q_STEP * index;
}

/* Wo as a uint32 angle (2^32 = 2 pi) from wq: round(wq 2^32 / 40960) = round(wq 2^19 / 5). wq <= 2034
 * keeps wq 2^19 < 2^31, and x mod 5 is never 2.5, so (x + 2) / 5 rounds to nearest. */
uint32_t wo_q_to_angle(int wq)
{
	return ((uint32_t)wq * (1u << 19) + 2u) / 5u;
}

/* Energy index from log2 E (Q16). The float encode_energy: index = (int)(32 (10 log10 max(E, 1e-12) + 10) / 50 + 0.5),
 * clamped to 0..31, i.e. floor(0.64 10 log10(2) log2 E + 6.9) (truncation and floor agree after the clamp). */
int encode_energy_log2(int32_t log2E_q16)
{
	_Static_assert(E_BITS == 5 && E_MIN_DB == -10 && E_MAX_DB == 40, "energy grid assumed below");
	if (log2E_q16 < -8 * 65536) /* index < 0 well before here; also covers LOG2E_ZERO */
		return 0;
	if (log2E_q16 > 20 * 65536)
		return 31;
	/* v = (log2E fxp_lpc_eidx_q29) >> 29 (Q16), 16-bit split: |h| <= 20, K < 2^30 */
	int32_t h = log2E_q16 >> 16, kh = fxp_lpc_eidx_q29 >> 16;
	uint32_t l = (uint32_t)log2E_q16 & 0xffffu, kl = (uint32_t)fxp_lpc_eidx_q29 & 0xffffu;
	int32_t v = h * kh * 8 + ((h * (int32_t)kl) >> 13) + (int32_t)((l * (uint32_t)kh) >> 13) + (int32_t)((l * kl) >> 29);
	int index = (v + 452198) >> 16; /* 6.9 2^16 = 452198.4 */
	return index < 0 ? 0 : index > 31 ? 31 : index;
}

/* Decoder energy as log2 E in Q16. The float decoder used E = 10^(e_dB / 10) with
 * e_dB = E_MIN_DB + index (E_MAX_DB - E_MIN_DB) / 32 = (-640 + 100 index) / 64 dB, so
 * log2 E = d64 log2(10) / 640 = d64 * 340.16544 / 2^16 * 2^16: 340 + 10842 / 2^16 (error < 0.02 LSB per unit). */
int32_t decode_energy_log2_q16(int index)
{
	_Static_assert(E_BITS == 5 && E_MIN_DB == -10 && E_MAX_DB == 40, "energy grid assumed below");
	int32_t d64 = -640 + 100 * index; /* |d64| <= 2460 */
	return d64 * 340 + ((d64 * 10842 + 32768) >> 16);
}

/* LSP delta indices from the LSPs in Hz (Q16), as the float encode_lspds_scalar: each target (LSP minus the sum
 * of the previous quantised deltas) takes the nearest codebook entry, the lower index on ties. */
void encode_lspds_q16(int *indexes, const int32_t *lsp_hz_q16)
{
	int last_q_hz = 0;
	for (int i = 0; i < LPC_ORD; i++)
	{
		const uint16_t *cb = delta_lsp_cb[i];
		int32_t target = lsp_hz_q16[i] - (last_q_hz << 16); /* |target| < 2^28 */
		int best_j;
		if (target <= (int32_t)cb[0] << 16)
			best_j = 0;
		else if (target >= (int32_t)cb[31] << 16)
			best_j = 31;
		else
		{
			int lo = 0, hi = 31; /* first j with cb[j] > target */
			while (lo + 1 < hi)
			{
				int mid = (lo + hi) >> 1;
				if ((int32_t)cb[mid] << 16 <= target)
					lo = mid;
				else
					hi = mid;
			}
			int32_t e_lo = target - ((int32_t)cb[lo] << 16), e_hi = ((int32_t)cb[hi] << 16) - target;
			best_j = e_lo <= e_hi ? lo : hi;
		}
		indexes[i] = best_j;
		last_q_hz += cb[best_j];
	}
}

/* Decoded LSPs as integer half-Hz (lsp = pi hz2 / 8000): the codebook deltas are integer Hz.
 * Any index pattern is valid input (bit errors): the row maxima sum to 9800 Hz, so hz2 <= 19600 fits int16.
 * LSPs above 4000 Hz (> pi) alias in fxp_lsp_cos_q30's angle wrap exactly as cosf() aliases them in the float code. */
void decode_lspds_hz2(int16_t *hz2, const int *indexes)
{
	int lsp_hz = 0;
	for (int i = 0; i < LPC_ORD; i++)
	{
		lsp_hz += delta_lsp_cb[i][indexes[i]];
		hz2[i] = (int16_t)(2 * lsp_hz); /* <= 19600 */
	}
}
