#ifndef CODEC2_MOD_INTERNAL_H
#define CODEC2_MOD_INTERNAL_H

/*
  This work uses code written by David Rowe VK5DGR et al.
  https://github.com/drowe67/codec2
*/

#include <stdint.h>
#include <string.h>



#define BYTES_PER_FRAME 8
#define SAMPLES_PER_FRAME (2 * N_SAMP)

#define MAX_AMP 80

#define FFT_ENC 512
#define FFT_DEC 512
#define V_THRESH 6.0f
#define V_THRESH_LIN (powf(10.0f, V_THRESH * 0.1f))
#define LPC_ORD 10


#define SAMP_RATE 8000
#define N_SAMP 80
#define M_PITCH 320
#define P_MIN 20
#define P_MAX 160
/* Decoder Wo as an integer wq: Wo = 2 pi wq / 40960. Exact for every decoded Wo (W0_MIN + step k),
 * their means (interp_Wo) and W0_MIN; pi / Wo = 20480 / wq. */
#define WO_Q_MIN 256  /* 2 pi / 160 = W0_MIN */
#define WO_Q_STEP 14  /* (W0_MAX - W0_MIN) / 2^WO_BITS = 2 pi 7 / 20480 */
#define NW 279
#define TW 40

#define PE_FFT_SIZE 512
#define DEC 5
#define NDEC (M_PITCH / DEC)
#define T 0.1f
#define CNLP 0.3f
#define NLP_NTAP 48
#define NLP_HIST (NLP_NTAP - 1 + N_SAMP) /* NLP squared-signal history: FIR memory + one frame */
#define NLP_NBLK (NDEC / (N_SAMP / DEC))  /* frames held in the decimated history (4) */

#define WO_BITS 7
#define WO_LEVELS (1 << WO_BITS)
#define E_BITS 5
#define E_LEVELS (1 << E_BITS)
#define E_MIN_DB (-10) /* energy grid, dB */
#define E_MAX_DB 40
#define BG_THRESH 40 /* postfilter background estimate, dB */
#define BG_MARGIN 6
#define LSPD_SCALAR_INDEXES 10



#define CODEC2_RAND_MAX 32767

/* asserts */
/* this is the exact Codec2 3200 bps configuration */
/* changing any of these parameters requires DSP math modification and buffer resizing */
_Static_assert(BYTES_PER_FRAME == 8, "Codec2 3200 only");
_Static_assert(LPC_ORD == 10, "Codec2 3200 assumes LPC_ORD == 10");
_Static_assert(MAX_AMP >= 80, "MAX_AMP must be >= 80 for Codec2 3200");
_Static_assert(FFT_ENC == 512, "Codec2 3200 assumes FFT_ENC == 512");
_Static_assert(FFT_DEC == FFT_ENC, "Encoder/decoder FFT sizes must match");
_Static_assert(PE_FFT_SIZE == 512, "Pitch estimator FFT size must be 512");
_Static_assert((M_PITCH % DEC) == 0, "M_PITCH must be divisible by DEC");
_Static_assert(NDEC == (M_PITCH / DEC), "NDEC must equal M_PITCH / DEC");
_Static_assert(NLP_NTAP == 48, "NLP FIR tables (fxp_nlp_tables.c) have 48 taps");
_Static_assert(WO_BITS <= 8, "encode_Wo_p() assumes <= 8-bit fields");
_Static_assert(E_BITS <= 8, "encode_energy_log2() assumes <= 8-bit fields");
_Static_assert(TW <= N_SAMP / 2, "Parzen window ramp too long");

/* consts */
extern const uint16_t delta_lsp_cb[LPC_ORD][32];

/* Codec state. The structs are defined with literal sizes in the public header <codec2/codec2_fxp.h> so a caller
 * can allocate them statically; these checks tie every array to the constants used here. Fields:
 *   nlp_t:            sq[NLP_HIST]  DC-notched Sn^2 + 1 in units of 2^-nk; sq_fir[NDEC] decimated low-pass output,
 *                     block-float (fir_e); nk notch scale exponent; mem_x / mem_y notch previous input / output;
 *                     mem_e 0.95 mem_y rounding remainder (1/20 units); fir_e[NLP_NBLK] block exponents
 *   codec2_encoder_t: Sn[M_PITCH] input history; prev_bin_enc NLP tracking memory (as the NLP bin); scratch reused
 *                     by the NLP (257 cq15_t), the analysis (Sw 257 cq15_t + pw 512 uint32) and the LPC window
 *   codec2_decoder_t: Sn_ synthesised speech Q4; Pn Parzen window Q14; ex_phase (2^32 = 2 pi); bg_est dB Q16;
 *                     prev_* previous frame voicing, Wo (see WO_Q_MIN), log2 E Q16, LSPs half-Hz; scratch for
 *                     fxp_aks_to_mag2, then the irfft buffer (cq15_t) */
#include "codec2/codec2_fxp.h"
typedef c2fxp_nlp_t nlp_t;
typedef c2fxp_encoder_t codec2_encoder_t;
typedef c2fxp_decoder_t codec2_decoder_t;
#define C2_COUNT(s, f) (sizeof(((s *)0)->f) / sizeof(((s *)0)->f[0]))
_Static_assert(C2_COUNT(nlp_t, sq) == NLP_HIST, "nlp sq");
_Static_assert(C2_COUNT(nlp_t, sq_fir) == NDEC, "nlp sq_fir");
_Static_assert(C2_COUNT(nlp_t, fir_e) == NLP_NBLK, "nlp fir_e");
_Static_assert(C2_COUNT(codec2_encoder_t, Sn) == M_PITCH, "encoder Sn");
_Static_assert(C2_COUNT(codec2_encoder_t, scratch) == FFT_ENC / 2 + 1 + FFT_ENC, "encoder scratch");
_Static_assert(C2_COUNT(codec2_decoder_t, Sn_) == 2 * N_SAMP, "decoder Sn_");
_Static_assert(C2_COUNT(codec2_decoder_t, Pn) == 2 * N_SAMP, "decoder Pn");
_Static_assert(C2_COUNT(codec2_decoder_t, prev_hz2_dec) == LPC_ORD, "decoder prev_hz2_dec");
_Static_assert(C2_COUNT(codec2_decoder_t, scratch) == 514, "decoder scratch");
_Static_assert(C2FXP_BYTES_PER_FRAME == BYTES_PER_FRAME && C2FXP_SAMPLES_PER_FRAME == SAMPLES_PER_FRAME, "frame");


#endif /* CODEC2_MOD_INTERNAL_H */
