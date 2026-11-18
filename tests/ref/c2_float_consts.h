/* Float constants of the float Codec2-mod encoder/decoder, for the host reference code in tests/ only.
 * The codec itself (codec/) has no float constants; its integer equivalents are documented where used. */
#ifndef C2_FLOAT_CONSTS_H
#define C2_FLOAT_CONSTS_H
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define TWO_PI (2.0f * (float)M_PI)
#define P_MAX_S 0.0200f
#define W0_MIN (TWO_PI / 160.0f)
#define W0_MAX (TWO_PI / 20.0f)
#define LPCPF_GAMMA 0.5f
#define LPCPF_BETA 0.2f
#define LPCPF_TWO_BETA (2.0f * LPCPF_BETA)
#define LSP_DELTA1 0.01f
#define BG_BETA 0.1f
#define LPC_ENERGY_FLOOR 1e-6f
#define FFT_R (TWO_PI / FFT_ENC)
#define FFT_1_R (FFT_ENC / TWO_PI)
#endif
