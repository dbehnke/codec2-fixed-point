/*
 * codec2.h: libcodec2-compatible API of codec2-fixed-point.
 *
 * Same declarations as upstream Codec2's codec2.h (https://github.com/drowe67/codec2), so programs written for
 * libcodec2 build and link against this library unchanged. Only CODEC2_MODE_3200 is implemented:
 * codec2_create() returns NULL for every other mode, as upstream does for a mode disabled at build time.
 * The remaining entry points exist so existing programs link; see README.md for how each behaves in 3200.
 *
 * Copyright (C) 2010 David Rowe and the Codec2 contributors; Codec2-mod by the M17 Project;
 * fixed-point conversion by dbehnke. LGPL-2.1, see LICENSE.
 */
#ifndef __CODEC2__
#define __CODEC2__
#include <codec2/version.h>
#include <stdbool.h>
#include <codec2/codec2_fxp.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CODEC2_MODE_3200 0
#define CODEC2_MODE_2400 1
#define CODEC2_MODE_1600 2
#define CODEC2_MODE_1400 3
#define CODEC2_MODE_1300 4
#define CODEC2_MODE_1200 5
#define CODEC2_MODE_700C 8

/* Mode availability, as upstream's compile-time switches: only 3200 is built. */
#define CODEC2_MODE_EN_DEFAULT 0
#define CODEC2_MODE_3200_EN 1
#define CODEC2_MODE_2400_EN 0
#define CODEC2_MODE_1600_EN 0
#define CODEC2_MODE_1400_EN 0
#define CODEC2_MODE_1300_EN 0
#define CODEC2_MODE_1200_EN 0
#define CODEC2_MODE_700C_EN 0

#define CODEC2_MODE_ACTIVE(mode_name, var) \
  ((mode_name##_EN) == 0 ? 0 : (var) == mode_name)

struct CODEC2;

C2FXP_API struct CODEC2 *codec2_create(int mode);
C2FXP_API void codec2_destroy(struct CODEC2 *codec2_state);
C2FXP_API void codec2_encode(struct CODEC2 *codec2_state, unsigned char bytes[],
                             short speech_in[]);
C2FXP_API void codec2_decode(struct CODEC2 *codec2_state, short speech_out[],
                             const unsigned char bytes[]);
C2FXP_API void codec2_decode_ber(struct CODEC2 *codec2_state, short speech_out[],
                                 const unsigned char *bytes, float ber_est);
C2FXP_API int codec2_samples_per_frame(struct CODEC2 *codec2_state);
C2FXP_API int codec2_bits_per_frame(struct CODEC2 *codec2_state);
C2FXP_API int codec2_bytes_per_frame(struct CODEC2 *codec2_state);

C2FXP_API void codec2_set_lpc_post_filter(struct CODEC2 *codec2_state, int enable,
                                          int bass_boost, float beta, float gamma);
C2FXP_API int codec2_get_spare_bit_index(struct CODEC2 *codec2_state);
C2FXP_API int codec2_rebuild_spare_bit(struct CODEC2 *codec2_state, char unpacked_bits[]);
C2FXP_API void codec2_set_natural_or_gray(struct CODEC2 *codec2_state, int gray);
C2FXP_API void codec2_set_softdec(struct CODEC2 *c2, float *softdec);
C2FXP_API float codec2_get_energy(struct CODEC2 *codec2_state, const unsigned char *bits);

/* ML and VQ experiments: not supported (no-ops; codec2_get_var returns 0, codec2_enable_user_ratek NULL) */
C2FXP_API void codec2_open_mlfeat(struct CODEC2 *codec2_state, char *feat_filename,
                                  char *model_filename);
C2FXP_API void codec2_load_codebook(struct CODEC2 *codec2_state, int num, char *filename);
C2FXP_API float codec2_get_var(struct CODEC2 *codec2_state);
C2FXP_API float *codec2_enable_user_ratek(struct CODEC2 *codec2_state, int *K);

/* 700C post filter and equaliser: no-ops (700C is not built) */
C2FXP_API void codec2_700c_post_filter(struct CODEC2 *codec2_state, bool en);
C2FXP_API void codec2_700c_eq(struct CODEC2 *codec2_state, bool en);

#ifdef __cplusplus
}
#endif

#endif
