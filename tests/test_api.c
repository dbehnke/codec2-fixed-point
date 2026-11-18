/* test_api.c: public API checks for both headers, and a bit-exact regression of the whole codec.
 *
 * The expected hashes come from the research tree (codec2research 75b25c2) encoding and decoding the
 * deterministic signal in signal.h. The codec is integer-only, so the hashes are the same on every platform. */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <codec2/codec2.h>
#include "signal.h"

#define EXP_SIGNAL 0xcd8c456au
#define EXP_BITS 0x9959853fu
#define EXP_PCM 0xf1391accu

static int fails;
#define CHECK(c) do { if (!(c)) { printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #c); fails++; } } while (0)

static c2fxp_encoder_t enc;
static c2fxp_decoder_t dec;

int main(void)
{
    /* modes: only 3200 is built */
    int modes[] = {CODEC2_MODE_2400, CODEC2_MODE_1600, CODEC2_MODE_1400, CODEC2_MODE_1300,
                   CODEC2_MODE_1200, CODEC2_MODE_700C, -1, 99};
    for (unsigned i = 0; i < sizeof modes / sizeof modes[0]; i++)
        CHECK(codec2_create(modes[i]) == NULL);

    struct CODEC2 *c2 = codec2_create(CODEC2_MODE_3200);
    CHECK(c2 != NULL);
    if (!c2) return 1;
    CHECK(codec2_samples_per_frame(c2) == 160);
    CHECK(codec2_bits_per_frame(c2) == 64);
    CHECK(codec2_bytes_per_frame(c2) == 8);
    CHECK(codec2_get_spare_bit_index(c2) == -1);
    CHECK(CODEC2_VERSION_MAJOR == 1 && CODEC2_VERSION_MINOR == 2);
    CHECK(strcmp(c2fxp_version(), C2FXP_VERSION) == 0);
    codec2_set_natural_or_gray(c2, 0);           /* no effect in 3200, as upstream */
    codec2_set_lpc_post_filter(c2, 1, 1, 0.2f, 0.5f);

    /* the same signal through the upstream API and the native API */
    c2fxp_encoder_init(&enc);
    c2fxp_decoder_init(&dec);
    sig_t s;
    sig_init(&s);
    uint32_t hs = FNV_INIT, hb = FNV_INIT, hp = FNV_INIT;
    int same_bits = 1, same_pcm = 1, energy_ok = 1;
    for (int f = 0; f < SIG_FRAMES; f++) {
        short pcm[160], out[160];
        int16_t pcm2[160], out2[160];
        unsigned char b[8];
        uint8_t b2[8];
        sig_frame(&s, f, pcm2);
        memcpy(pcm, pcm2, sizeof pcm);
        hs = fnv1a(hs, pcm, sizeof pcm);
        codec2_encode(c2, b, pcm);
        c2fxp_encode(&enc, b2, pcm2);
        same_bits &= !memcmp(b, b2, 8);
        hb = fnv1a(hb, b, 8);
        if (f & 1)
            codec2_decode(c2, out, b);
        else
            codec2_decode_ber(c2, out, b, 0.0f);
        c2fxp_decode(&dec, out2, b2);
        same_pcm &= !memcmp(out, out2, sizeof out);
        hp = fnv1a(hp, out, sizeof out);
        /* energy index: bits 9..13, Gray-coded */
        unsigned g = ((unsigned)b[1] >> 2) & 31, idx = g;
        for (unsigned sh = 1; sh < 5; sh <<= 1) idx ^= idx >> sh;
        double ref = pow(10.0, (-10.0 + 50.0 * idx / 32.0) / 10.0); /* upstream decode_energy() */
        energy_ok &= fabs(codec2_get_energy(c2, b) - ref) <= 1e-6 * ref;
    }
    codec2_destroy(c2);
    printf("signal %08x bits %08x pcm %08x\n", hs, hb, hp);
    CHECK(hs == EXP_SIGNAL);
    CHECK(hb == EXP_BITS);
    CHECK(hp == EXP_PCM);
    CHECK(same_bits);
    CHECK(same_pcm);
    CHECK(energy_ok);
    printf(fails ? "test_api: %d FAILED\n" : "test_api: all passed\n", fails);
    return fails != 0;
}
