#include "codec2_internal.h"
#include "codec2/version.h"
#include "analysis.h"
#include "synthesis.h"
#include "nlp.h"
#include "lpc.h"
#include "interp.h"
#include "quantise.h"

void c2fxp_encoder_init(codec2_encoder_t *c2)
{
	for (int i = 0; i < M_PITCH; i++)
		c2->Sn[i] = 1;

	c2->prev_bin_enc = PE_FFT_SIZE * DEC / P_MAX; /* f0 = 1 / P_MAX_S = 50 Hz, as the float code's prev_f0 */

	nlp_init(&c2->nlp);
}

void c2fxp_decoder_init(codec2_decoder_t *c2)
{
	c2->next_rn = 1;

	memset(c2->Sn_, 0, sizeof(c2->Sn_));

	/* the decoder is FFT-free in float: fxp_irfft512 and the LSP product form need no kiss configs */
	synthesis_init(c2);

	c2->bg_est = 0;
	c2->ex_phase = 0;

	c2->prev_voiced_dec = 0;
	c2->prev_wq_dec = WO_Q_MIN;  /* Wo = 2 pi / P_MAX */
	c2->prev_log2e_dec = 0;      /* E = 1 */
	/* lsp = i pi / (LPC_ORD + 1), i.e. 8000 i / 11 half-Hz. Stored as 2 round(4000 i / 11): even like every
	 * decoded LSP, so the first interpolation (prev + next) / 2 = next / 2 + round(4000 i / 11) is exactly
	 * the float decoder's lrintf of its float mean. */
	for (int i = 0; i < LPC_ORD; i++)
		c2->prev_hz2_dec[i] = (int16_t)(2 * ((2 * 4000 * i + LPC_ORD + 1) / (2 * (LPC_ORD + 1))));
}

void c2fxp_encode(codec2_encoder_t *c2, uint8_t bits[C2FXP_BYTES_PER_FRAME], const int16_t speech[C2FXP_SAMPLES_PER_FRAME])
{
	enc_frame_t f;
	int32_t lsp_hz[LPC_ORD]; /* Q16 */
	int32_t log2E;
	int Wo_index, e_index;
	int lspd_indexes[LPC_ORD];
	unsigned int nbit = 0;

	memset(bits, 0, BYTES_PER_FRAME); // 64 bits

	/* first 10ms analysis frame - we just want voicing */
	analyse_one_frame(c2, &f, speech);
	pack(bits, &nbit, f.voiced, 1);

	/* second 10ms analysis frame */
	analyse_one_frame(c2, &f, &speech[N_SAMP]);
	pack(bits, &nbit, f.voiced, 1);
	Wo_index = encode_Wo_p(f.pn, f.pd);
	pack(bits, &nbit, Wo_index, WO_BITS);

	log2E = speech_to_lsps_q16(c2, lsp_hz);
	e_index = encode_energy_log2(log2E);
	pack(bits, &nbit, e_index, E_BITS);

	encode_lspds_q16(lspd_indexes, lsp_hz);
	for (int i = 0; i < LSPD_SCALAR_INDEXES; i++)
	{
		pack(bits, &nbit, lspd_indexes[i], 5); // codebook indices are 5 bits
	}
}

void c2fxp_decode(codec2_decoder_t *c2, int16_t speech[C2FXP_SAMPLES_PER_FRAME], const uint8_t bits[C2FXP_BYTES_PER_FRAME])
{
	int voiced[2], wq[2];
	int32_t log2e[2];         /* log2 E, Q16 */
	int16_t hz2[2][LPC_ORD];  /* LSPs, half-Hz */
	int lspd_indexes[LPC_ORD];
	int32_t lA[MAX_AMP + 1];  /* log2 A[m], Q16 */
	uint32_t phA[MAX_AMP + 1]; /* arg A at the harmonics, 2^32 = 2 pi */

	/* unpack bits from channel: voicing of both 10 ms subframes, then the 20 ms parameters */
	unsigned int nbit = 0;
	voiced[0] = unpack(bits, &nbit, 1);
	voiced[1] = unpack(bits, &nbit, 1);
	wq[1] = decode_Wo_q(unpack(bits, &nbit, WO_BITS));
	log2e[1] = decode_energy_log2_q16(unpack(bits, &nbit, E_BITS));
	for (int i = 0; i < LSPD_SCALAR_INDEXES; i++)
		lspd_indexes[i] = unpack(bits, &nbit, 5); // codebook indices are 5 bits
	decode_lspds_hz2(hz2[1], lspd_indexes);

	/* interpolate the first 10 ms subframe between the previous and this 20 ms frame */
	wq[0] = interp_Wo_q(&voiced[0], c2->prev_voiced_dec, c2->prev_wq_dec, voiced[1], wq[1]);
	log2e[0] = interp_energy_log2_q16(c2->prev_log2e_dec, log2e[1]);
	interpolate_lsp_hz2(hz2[0], c2->prev_hz2_dec, hz2[1]);

	for (int i = 0; i < 2; i++)
	{
		uint32_t W = wo_q_to_angle(wq[i]);
		int L = 20480 / wq[i]; /* floor(pi / Wo), as the float (int) conversion */
		aks_to_mag2_lsp(c2, hz2[i], log2e[i], W, L, lA, phA);
		synthesise_one_frame(c2, &speech[N_SAMP * i], W, L, voiced[i], lA, phA);
	}

	/* update memories for next frame */
	c2->prev_voiced_dec = voiced[1];
	c2->prev_wq_dec = wq[1];
	c2->prev_log2e_dec = log2e[1];
	memcpy(c2->prev_hz2_dec, hz2[1], sizeof(c2->prev_hz2_dec));
}

const char *c2fxp_version(void) { return C2FXP_VERSION; }
