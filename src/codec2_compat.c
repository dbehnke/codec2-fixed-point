/*
 * codec2_compat.c: the upstream libcodec2 API (codec2.h) on top of the fixed-point 3200 codec.
 *
 * The codec itself is integer-only. This file adds the heap-allocated struct CODEC2 that upstream programs
 * expect, and the float arguments and returns of the upstream signatures (no float arithmetic beyond a table
 * read). Firmware that wants no heap and no float at all should use <codec2/codec2_fxp.h> directly.
 */
#include <stdlib.h>
#include "codec2/codec2.h"
#include "codec2_internal.h"
#include "quantise.h"

struct CODEC2 {
    int mode;
    c2fxp_encoder_t enc;
    c2fxp_decoder_t dec;
};

struct CODEC2 *codec2_create(int mode)
{
    if (!CODEC2_MODE_ACTIVE(CODEC2_MODE_3200, mode))
        return NULL;
    struct CODEC2 *c2 = malloc(sizeof *c2);
    if (!c2)
        return NULL;
    c2->mode = mode;
    c2fxp_encoder_init(&c2->enc);
    c2fxp_decoder_init(&c2->dec);
    return c2;
}

void codec2_destroy(struct CODEC2 *c2)
{
    free(c2);
}

/* short is 16 bits on every supported target, so int16_t and short buffers are interchangeable */
_Static_assert(sizeof(short) == sizeof(int16_t), "short must be 16 bits");

void codec2_encode(struct CODEC2 *c2, unsigned char bytes[], short speech_in[])
{
    c2fxp_encode(&c2->enc, bytes, (const int16_t *)speech_in);
}

void codec2_decode(struct CODEC2 *c2, short speech_out[], const unsigned char bytes[])
{
    c2fxp_decode(&c2->dec, (int16_t *)speech_out, bytes);
}

/* upstream ignores ber_est in 3200 (it is used only by 1300) */
void codec2_decode_ber(struct CODEC2 *c2, short speech_out[], const unsigned char *bytes, float ber_est)
{
    (void)ber_est;
    codec2_decode(c2, speech_out, bytes);
}

int codec2_samples_per_frame(struct CODEC2 *c2) { (void)c2; return C2FXP_SAMPLES_PER_FRAME; }
int codec2_bits_per_frame(struct CODEC2 *c2) { (void)c2; return C2FXP_BITS_PER_FRAME; }
int codec2_bytes_per_frame(struct CODEC2 *c2) { (void)c2; return C2FXP_BYTES_PER_FRAME; }

/* The decoder's LPC post filter is fixed at upstream's defaults (enabled, bass boost on, beta 0.2,
 * gamma 0.5). Other settings are not supported and are ignored. */
void codec2_set_lpc_post_filter(struct CODEC2 *c2, int enable, int bass_boost, float beta, float gamma)
{
    (void)c2; (void)enable; (void)bass_boost; (void)beta; (void)gamma;
}

/* 3200 has no spare bit: upstream returns -1 for this mode from both */
int codec2_get_spare_bit_index(struct CODEC2 *c2) { (void)c2; return -1; }
int codec2_rebuild_spare_bit(struct CODEC2 *c2, char unpacked_bits[]) { (void)c2; (void)unpacked_bits; return -1; }

/* Upstream's 3200 mode always packs Gray-coded fields regardless of this setting, and so does this codec. */
void codec2_set_natural_or_gray(struct CODEC2 *c2, int gray) { (void)c2; (void)gray; }

/* Soft decisions are used only by 700C upstream. */
void codec2_set_softdec(struct CODEC2 *c2, float *softdec) { (void)c2; (void)softdec; }

/* Frame energy from the 5-bit energy index (bits 9..13, Gray-coded), as upstream decode_energy():
 * 10^((-10 + 50 i / 32) / 10), correctly rounded to float. */
static const float energy_tab[32] = {
    0.100000001f, 0.143301263f, 0.2053525f, 0.294272721f, 0.421696514f, 0.604296386f, 0.865964353f,
    1.24093771f, 1.77827942f, 2.54829669f, 3.65174127f, 5.23299122f, 7.4989419f, 10.7460785f, 15.3992653f,
    22.0673409f, 31.622776f, 45.3158379f, 64.9381638f, 93.0572052f, 133.352142f, 191.095291f, 273.841949f,
    392.418976f, 562.341309f, 805.842163f, 1154.78198f, 1654.81714f, 2371.37378f, 3398.20825f, 4869.67529f,
    6978.30566f};

float codec2_get_energy(struct CODEC2 *c2, const unsigned char *bits)
{
    (void)c2;
    unsigned int nbit = 1 + 1 + WO_BITS;
    return energy_tab[unpack(bits, &nbit, E_BITS)];
}

void codec2_open_mlfeat(struct CODEC2 *c2, char *feat_filename, char *model_filename)
{
    (void)c2; (void)feat_filename; (void)model_filename;
}
void codec2_load_codebook(struct CODEC2 *c2, int num, char *filename) { (void)c2; (void)num; (void)filename; }
float codec2_get_var(struct CODEC2 *c2) { (void)c2; return 0.0f; }
float *codec2_enable_user_ratek(struct CODEC2 *c2, int *K) { (void)c2; (void)K; return NULL; }
void codec2_700c_post_filter(struct CODEC2 *c2, bool en) { (void)c2; (void)en; }
void codec2_700c_eq(struct CODEC2 *c2, bool en) { (void)c2; (void)en; }
