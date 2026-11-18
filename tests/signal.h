/* signal.h: deterministic integer-only test speech for regression hashes (identical on every platform).
 * 300 frames of 160 samples: voiced pulse train through a two-pole resonator with a slowly moving pitch and
 * formant, plus segments of silence, noise only, clipped loud speech, very quiet speech and a fast pitch sweep. */
#ifndef C2FXP_TEST_SIGNAL_H
#define C2FXP_TEST_SIGNAL_H
#include <stdint.h>

#define SIG_FRAMES 300

typedef struct { uint32_t rng; int32_t y1, y2; int phase; } sig_t;

static void sig_init(sig_t *s) { s->rng = 12345u; s->y1 = s->y2 = 0; s->phase = 0; }

static int32_t sig_noise(sig_t *s) { s->rng = s->rng * 1664525u + 1013904223u; return (int32_t)(s->rng >> 20) - 2048; }

static void sig_frame(sig_t *s, int f, int16_t out[160])
{
    for (int i = 0; i < 160; i++) {
        int n = f * 160 + i;
        int seg = f < 50 ? 0 : f < 80 ? 1 : f < 130 ? 2 : f < 200 ? 3 : f < 250 ? 4 : 5;
        int p = 40 + ((n / 400) % 81);                 /* pitch period 40..120 samples */
        if (seg == 5) p = 30 + ((n / 37) % 110);        /* fast sweep */
        int32_t e = 0;
        if (seg != 1 && seg != 2 && ++s->phase >= p) { s->phase = 0; e = 4000; }
        e += sig_noise(s) >> (seg == 2 ? 0 : 3);
        if (seg == 1) e = 0;
        int32_t a1 = 16000 + ((n / 800) % 5) * 2500;    /* 2 r cos(theta), Q14; r^2 = 0.9: poles at radius 0.95 */
        int32_t y = (int32_t)(((int64_t)a1 * s->y1 - (int64_t)14746 * s->y2) >> 14) + e;
        if (y > 1 << 20) y = 1 << 20;
        if (y < -(1 << 20)) y = -(1 << 20);
        s->y2 = s->y1; s->y1 = y;
        int32_t v = y;
        if (seg == 3) v *= 16;
        if (seg == 4) v >>= 6;
        out[i] = (int16_t)(v > 32767 ? 32767 : v < -32768 ? -32768 : v);
    }
}

static uint32_t fnv1a(uint32_t h, const void *p, unsigned n)
{
    const uint8_t *b = p;
    while (n--) { h ^= *b++; h *= 16777619u; }
    return h;
}
#define FNV_INIT 2166136261u
#endif
