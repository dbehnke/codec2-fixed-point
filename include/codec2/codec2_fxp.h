/*
 * codec2_fxp.h: native API of codec2-fixed-point, an integer-only Codec2 3200 encoder and decoder.
 *
 * No floating point, no libm, no heap: the caller owns the state (allocate it statically or on the heap).
 * Bitstream: standard Codec2 3200, 64 bits (8 bytes) per 20 ms frame of 160 samples at 8 kHz, as carried by
 * M17 voice streams. For the upstream libcodec2 API (codec2_create, codec2_encode, ...) include <codec2/codec2.h>.
 *
 * Copyright (C) 2010 David Rowe and the Codec2 contributors; Codec2-mod by the M17 Project;
 * fixed-point conversion by dbehnke. LGPL-2.1, see LICENSE.
 */
#ifndef C2FXP_PUBLIC_CODEC2_FXP_H
#define C2FXP_PUBLIC_CODEC2_FXP_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#if defined(C2FXP_BUILDING_SHARED) && (defined(__GNUC__) || defined(__clang__))
#define C2FXP_API __attribute__((visibility("default")))
#else
#define C2FXP_API
#endif

#define C2FXP_BITS_PER_FRAME 64
#define C2FXP_BYTES_PER_FRAME 8
#define C2FXP_SAMPLES_PER_FRAME 160
#define C2FXP_SAMPLE_RATE 8000

/* Codec state. The layouts are private: only their sizes are part of the API, so they can be allocated
 * without a heap. Initialise with c2fxp_encoder_init / c2fxp_decoder_init; do not read or write fields. */
typedef struct c2fxp_nlp {
    int32_t sq[127];
    int32_t sq_fir[64];
    int32_t nk, mem_x, mem_y, mem_e;
    int8_t fir_e[4];
} c2fxp_nlp_t;

typedef struct c2fxp_encoder {
    int16_t Sn[320];
    c2fxp_nlp_t nlp;
    int prev_bin_enc;
    int32_t scratch[769];
} c2fxp_encoder_t;

typedef struct c2fxp_decoder {
    uint32_t next_rn;
    int32_t Sn_[160];
    int16_t Pn[160];
    uint32_t ex_phase;
    int32_t bg_est;
    int prev_voiced_dec;
    int prev_wq_dec;
    int32_t prev_log2e_dec;
    int16_t prev_hz2_dec[10];
    int32_t scratch[514];
} c2fxp_decoder_t;

C2FXP_API void c2fxp_encoder_init(c2fxp_encoder_t *enc);
C2FXP_API void c2fxp_decoder_init(c2fxp_decoder_t *dec);

/* Encode one 20 ms frame: speech[160] 16-bit PCM at 8 kHz -> bits[8]. */
C2FXP_API void c2fxp_encode(c2fxp_encoder_t *enc, uint8_t bits[C2FXP_BYTES_PER_FRAME],
                            const int16_t speech[C2FXP_SAMPLES_PER_FRAME]);

/* Decode one frame: bits[8] -> speech[160]. */
C2FXP_API void c2fxp_decode(c2fxp_decoder_t *dec, int16_t speech[C2FXP_SAMPLES_PER_FRAME],
                            const uint8_t bits[C2FXP_BYTES_PER_FRAME]);

/* Library version string, e.g. "0.1.0". */
C2FXP_API const char *c2fxp_version(void);

#ifdef __cplusplus
}
#endif

#endif /* C2FXP_PUBLIC_CODEC2_FXP_H */
