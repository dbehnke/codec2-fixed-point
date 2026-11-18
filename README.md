# codec2-fixed-point

An integer-only Codec2 3200 encoder and decoder, packaged as a drop-in `libcodec2` for M17.

- **Standard bitstream:** Codec2 3200, 64 bits (8 bytes) per 20 ms frame. That is what M17 voice streams carry, and any Codec2 3200 decoder plays it.
- **No floating point, no libm, no heap** in the codec. It is meant for radios without an FPU, such as the DM32's C-SKY CK803S.
- **libcodec2 compatible:** it installs as `libcodec2` with upstream's `codec2.h`, `codec2.pc` and `find_package(codec2)` config, and upstream's soname (`libcodec2.so.1.2`). M17 software builds against it unchanged; you only point the build at a different install prefix.
- **Native API** for firmware: `<codec2/codec2_fxp.h>`, with caller-owned state that can be allocated statically.

Status: version 0.1.0. Verified on the host (macOS and Linux, clang and GCC) and in a CK803S instruction-count emulator. **Not yet timed on real radio hardware.**

## Using it with M17 software

Build and install into a prefix of its own, so it sits beside any system Codec2 without replacing it:

```bash
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=/opt/codec2-fxp
```

```bash
cmake --build build -j && cmake --install build
```

Then point the M17 program's build at that prefix. For CMake projects (`find_package(codec2)`):

```bash
cmake -S m17-tools -B m17-build -DCMAKE_PREFIX_PATH=/opt/codec2-fxp
```

For pkg-config projects:

```bash
export PKG_CONFIG_PATH=/opt/codec2-fxp/lib/pkgconfig
```

A program that is already built against upstream `libcodec2.so.1.2` and uses only 3200 can also be switched at run time. Use `LD_LIBRARY_PATH=/opt/codec2-fxp/lib` on Linux, or `DYLD_LIBRARY_PATH` on macOS.

## Using it in firmware

Compile `src/*.c` (except `codec2_compat.c` if you don't want the upstream API) with `include/` and `src/` on the include path:

```c
#include <codec2/codec2_fxp.h>

static c2fxp_encoder_t enc;   /* 4,504 B */
static c2fxp_decoder_t dec;   /* 3,060 B */

c2fxp_encoder_init(&enc);
c2fxp_decoder_init(&dec);
c2fxp_encode(&enc, bits, pcm);    /* int16_t pcm[160] -> uint8_t bits[8] */
c2fxp_decode(&dec, pcm, bits);    /* uint8_t bits[8]  -> int16_t pcm[160] */
```

GCC or clang is required: the code uses `__builtin_clz` and `always_inline`, and relies on arithmetic right shifts of negative values.

## The upstream API: what each call does

| Call | Behaviour |
|---|---|
| `codec2_create(CODEC2_MODE_3200)` | Allocates the encoder and decoder state. |
| `codec2_create(any other mode)` | Returns `NULL`, as upstream does for a mode disabled at build time. |
| `codec2_encode`, `codec2_decode`, `codec2_destroy` | As upstream. |
| `codec2_decode_ber` | Same as `codec2_decode`; upstream also ignores `ber_est` in 3200. |
| `codec2_samples_per_frame` / `bits_per_frame` / `bytes_per_frame` | 160 / 64 / 8. |
| `codec2_get_energy` | Upstream's `decode_energy()` value for the frame's energy index. |
| `codec2_set_natural_or_gray` | No effect. Upstream's 3200 mode always uses Gray coding, and so does this codec. |
| `codec2_get_spare_bit_index`, `codec2_rebuild_spare_bit` | Return −1; 3200 has no spare bit, as upstream. |
| `codec2_set_lpc_post_filter` | Ignored. The post filter is fixed at upstream's defaults (on, bass boost on, β 0.2, γ 0.5). |
| `codec2_set_softdec`, the ML/VQ calls, `codec2_700c_*` | No-ops, so existing programs link. |

Only the public API is exported from the shared library; the codec's internal functions are hidden.

## Limitations

- **3200 only.** M17 voice streams use 3200, but M17 voice-plus-data streams use Codec2 1600, which this library does not implement.
- **The encoder is not bit-identical to upstream's.** The integer analysis resolves near-ties differently on a few percent of frames. Every frame is still a valid 3200 frame for any decoder, and the measured quality is the same (see below). A test that compares packets byte for byte against upstream's encoder will report differences.
- **The decoder's post-filter settings are fixed** at upstream's defaults.
- **No hardware timing yet.** The CK803S figures below are emulator instruction counts, not cycles.

## Quality

Measured on a 492-item corpus derived from LibriSpeech test-clean: nominal, loud, quiet and 15 dB SNR noisy speech plus synthetic edge cases, with separate tune and eval speakers. On the held-out eval set:

- Fixed-point encoder and decoder: PESQ-NB 2.676.
- The floating-point Codec2-mod reference: PESQ-NB 2.676.
- Upstream Codec2 at both ends: PESQ-NB 2.677.
- This encoder into the upstream decoder: PESQ-NB 2.676.

The speed-ups since then are bit-identical, so these figures still hold. In the author's blind listening tests the fixed-point version was preferred more often than the floating-point Codec2-mod. See `docs/11-fxp-lpc.md` for the per-condition figures.

## Cost (CK803S emulator, instruction counts)

QEMU instruction counts per 20 ms frame at `-O2`, over 240 frames of speech. These are not cycles and not DM32 timing.

| | mean | max |
|---|---:|---:|
| encode | 355,728 | 380,001 |
| decode | 438,135 | 478,980 |

| Resource | Size |
|---|---:|
| Stack, encode | 680 B |
| Stack, decode (cumulative) | 1,444 B |
| Code (`.text` of the emulator image) | about 19 KB |
| State, encoder | 4.5 KB |
| State, decoder | 3.1 KB |

A frame at 168 MHz is 3.36 M cycles.

## Verification

All of the following pass for release 0.1.0:

- **Unit tests** (`ctest`): FFT, LPC spectrum, decoder parameters, encoder analysis, NLP and LPC/LSP kernels, with the `FXP_CHECK` bound asserts on.
  - `-DC2FXP_SANITIZE=ON` adds UBSan and ASan.
  - `test_api` checks both APIs and a bit-exact regression of the whole codec on a deterministic signal. The expected hashes are the same on every platform.
- **Research-tree parity:** all 492 corpus items encode and decode byte-identically to the research tree this library was extracted from (codec2research `75b25c2`).
- **Upstream tools:** upstream's own `c2enc` and `c2dec` give output byte-identical to this library on all 492 items, both rebuilt from source through `pkg-config` and as the unmodified upstream binaries with this library swapped in at run time.
- **Interop:** the genuine upstream decoder decodes this encoder's bitstream, and this decoder decodes upstream's.
- **M17:**
  - The M17 Project's `m17-tools` (`m17-mod`, `m17-demod`, commit `806ec78`) builds against this library on Debian 13 with GCC, with no source changes.
  - Speech went through `m17-mod`, then the M17 baseband, then `m17-demod`. All 525 codec frames received in the overlap are identical to this library's direct encoding.

## Build options

| Option | Default | |
|---|---|---|
| `BUILD_SHARED_LIBS` | ON | shared `libcodec2` (OFF: static) |
| `C2FXP_BUILD_TOOLS` | ON | `c2fxp enc|dec` raw-file tool |
| `C2FXP_BUILD_TESTS` | ON | `ctest` |
| `C2FXP_SANITIZE` | OFF | UBSan and ASan for the tests |

## Layout

- `include/codec2/`: public headers. `codec2.h` is the upstream API, `codec2_fxp.h` the native API, `version.h` the version.
- `src/`: the codec. `codec2_compat.c` is the upstream API layer.
- `tests/`: unit tests, the API regression, and the float reference FFT used only by tests (`tests/ref/`).
- `tools/gen/`: generators for the constant tables in `src/*_tables.c`.
- `docs/`: design notes for each stage of the conversion. They refer to paths in the research tree (`codec/src` there is `src` here).

## License and credits

LGPL-2.1, see `LICENSE`.

- Codec2 by David Rowe and contributors (https://github.com/drowe67/codec2).
- Codec2-mod, the 3200-only refactor this is built from, by the M17 Project (https://github.com/M17-Project/Codec2-mod).
- Fixed-point conversion by dbehnke.
