/* c2fxp: encode or decode raw files with codec2-fixed-point (Codec2 3200).
 * usage: c2fxp enc in.raw out.bit | c2fxp dec in.bit out.raw
 * PCM is 8 kHz mono signed 16-bit native-endian; bits are 8 bytes per 20 ms frame, as upstream c2enc/c2dec 3200.
 * A trailing partial PCM frame is zero padded. "-" means stdin or stdout. */
#include <stdio.h>
#include <string.h>
#include <codec2/codec2_fxp.h>

static c2fxp_encoder_t enc;
static c2fxp_decoder_t dec;

int main(int argc, char **argv)
{
    if (argc != 4 || (strcmp(argv[1], "enc") && strcmp(argv[1], "dec"))) {
        fprintf(stderr, "usage: %s enc in.raw out.bit | dec in.bit out.raw   (codec2-fixed-point %s)\n",
                argv[0], c2fxp_version());
        return 2;
    }
    FILE *in = strcmp(argv[2], "-") ? fopen(argv[2], "rb") : stdin;
    FILE *out = strcmp(argv[3], "-") ? fopen(argv[3], "wb") : stdout;
    if (!in || !out) { perror("open"); return 1; }
    int16_t pcm[C2FXP_SAMPLES_PER_FRAME];
    uint8_t bits[C2FXP_BYTES_PER_FRAME];
    if (!strcmp(argv[1], "enc")) {
        c2fxp_encoder_init(&enc);
        size_t n;
        while ((n = fread(pcm, sizeof pcm[0], C2FXP_SAMPLES_PER_FRAME, in)) > 0) {
            if (n < C2FXP_SAMPLES_PER_FRAME) memset(pcm + n, 0, (C2FXP_SAMPLES_PER_FRAME - n) * sizeof pcm[0]);
            c2fxp_encode(&enc, bits, pcm);
            fwrite(bits, 1, sizeof bits, out);
        }
    } else {
        c2fxp_decoder_init(&dec);
        while (fread(bits, 1, sizeof bits, in) == sizeof bits) {
            c2fxp_decode(&dec, pcm, bits);
            fwrite(pcm, sizeof pcm[0], C2FXP_SAMPLES_PER_FRAME, out);
        }
    }
    if (in != stdin) fclose(in);
    if (out != stdout) fclose(out);
    return 0;
}
