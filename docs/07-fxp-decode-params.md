# 07 — Fixed-point decoder parameters (Wo, E, LSPs): the decoder is float-free

**Status:** `codec2_decode` now has no float arithmetic. The last float bridges, which were Wo → `W`, E → `log2 E` and LSP radians → half-Hz, are replaced by integer decoding and interpolation of the bitstream fields. With docs 04–06, the decoder runs in fixed point end to end, from the 64 packed bits to the int16 PCM. The encoder is unchanged and still float.

**Measurement types used here:**
- **host:** macOS arm64, clang. Used for accuracy and quality.
- **emulator:** XUANTIE QEMU instruction counts, `-icount shift=0`. These are instructions, **not cycles**.

Nothing here is DM32 hardware timing.

## Why these can be exact

All three parameters come from small integer codebooks, so the float decoder's values lie on exact rational grids:

| field | float decoder | integer form (this stage) |
|---|---|---|
| Wo (7 bits) | `W0_MIN + k·(W0_MAX − W0_MIN)/128` | `wq = 256 + 14k`, with Wo = 2π·wq/40960 |
| Wo, interpolated subframe | `prev + 0.5·(next − prev)`, or W0_MIN when unvoiced | `(wq_prev + wq_next) >> 1`, exact because every wq is even; 256 when unvoiced |
| L | `(int)(π / Wo)` | `20480 / wq` |
| W (Wo as a uint32 angle, doc 06) | `lrintf(Wo · 2^32/2π)` (the bridge) | `wo_q_to_angle(wq) = (wq·2^19 + 2) / 5` = round(wq·2^32/40960) |
| E (5 bits) | `10^(e_dB/10)`, e_dB = −10 + 1.5625·i | `log2 E` Q16 = `d64·340 + ((d64·10842 + 2^15) >> 16)`, d64 = −640 + 100·i (e_dB in 1/64 dB) |
| E, interpolated | `sqrtf(prev · next)` | `(log2 E_prev + log2 E_next) >> 1` |
| LSPs | `π/4000 · Σ` integer-Hz codebook deltas | half-Hz `hz2 = 2·Σ deltas` (int16) |
| LSPs, interpolated | `prev + 0.5·(next − prev)` | `(hz2_prev + hz2_next) >> 1`, exact because every hz2 is even |

- **Initial LSP state:** the float initial LSPs (i·π/11) are off the half-Hz grid. They are stored as `2·round(4000·i/11)`, which is even. That makes the first interpolation equal to `lrintf` of the float mean for every i. The obvious `round(8000·i/11)` gave a 1 half-Hz difference for some i, with either rounding of the mean.
- **LPC correction:** the threshold Wo < π·150/4000 becomes `W < 80530637`. No grid value of `wq` falls between the float and integer thresholds.
- **Decoder state** keeps the previous frame as `prev_voiced_dec`, `prev_wq_dec`, `prev_log2e_dec` and `prev_hz2_dec[10]`. That replaces a float `model_t`, which carried a dead `A[81]`, plus float LSPs and E.

Functions: `decode_Wo_q`, `wo_q_to_angle`, `decode_energy_log2_q16` and `decode_lspds_hz2` (in `quantise.c`); `interp_Wo_q`, `interp_energy_log2_q16` and `interpolate_lsp_hz2` (in `interp.c`). They replace the float `decode_Wo`, `decode_energy`, `decode_lspds_scalar`, `interp_Wo`, `interp_energy` and `interpolate_lsp`. `aks_to_mag2_lsp` (in `lpc.c`) now takes the integers directly; only the LPC correction remains there.

## Cleanup in the same change

- Removed the dead float `aks_to_mag2` and the whole `fxp_bridge.c/.h` (the irfft, rfft and LPC-DFT float bridges). Nothing used them any more.
- The decoder's 4 KB `kiss_fft_cpx fft_buffer[512]` became `int32_t scratch[514]` (2,056 B). It is shared by the `fxp_aks_to_mag2` scratch (514 words) and the synthesis irfft buffer (257 `cq15_t`). Static asserts check both uses.
- `synthesise_one_frame` takes `(W, L, voiced)` instead of a `model_t`.

## Accuracy

### Unit tests (`tests/test_decode_params.c`, in `scripts/test-fxp.sh` with FXP_CHECK + UBSan + ASan)

The test compares against the float decoder's formulas over every index and every interpolation case:

- **Wo:** 49,536 decoded and interpolated values (both orders, including W0_MIN).
  - `W` is within 0.400 LSB of exact. The float decoder's Wo lies within 3.3e-8 rad of the `wq` grid.
  - L is equal in 49,536 of 49,536 cases.
  - The old float bridge's W was up to 32 LSB from exact (host).
- **Voicing rules of `interp_Wo`:** all 8 cases match.
- **Energy:** `log2 E` is within 0.47 LSB (Q16) of exact for all 32 indices, and the interpolated value is within 0.93 LSB of `log2(sqrtf(prev·next))` in float.
- **LSPs:** 0 of 4,000,000 half-Hz values (decoded and interpolated, random indexes, starting from the initial state) differ from the old bridge's `lrintf(lsp·8000/π)`.

### End to end (host)

- **Against the previous commit's decoder (doc 06) on the same bitstreams, eval speech (241 items):** waveform SNR median 60.4 dB, p10 56.2, min 53.5.
  - The only input changes are W (≤ 32 LSB, now exact) and log2 E (< 1 LSB), and they move the harmonic phases slightly.
  - No postfilter random-phase divergence appears (the minimum stays above 50 dB).
- **Synthesis vs double re-synthesis (`scripts/synth-ref.sh`):**
  - eval: median 58.74 dB, min 52.04, postfilter flips 4 (arith) and 2 (bg) of 5.9M;
  - tune: median 58.73 dB, flips 8 and 3 of 5.3M.
  - Same as doc 06.

## Quality (host, `results/fxp-params-suite/{tune,eval}`)

**Eval set (held out), all speech conditions:**

| metric | SS | FF (doc 06) | **FF (this stage)** | target | met |
|---|---:|---:|---:|---|---|
| PESQ-NB | 2.676 | 2.674 | **2.673** | FF within 0.05 of SS | yes (−0.003) |
| LSD (dB) | 9.011 | 9.002 | **9.002** | ≤ SS + 0.2 | yes (−0.009) |
| level error (dB) | 0.646 | 0.646 | **0.646** | within ±0.5 of SS | yes (0.000) |

**Tune set:** PESQ 2.708 vs SS 2.705, LSD 9.096 vs 9.103, level error 0.528 vs 0.528 (unchanged from doc 06).

- **Per condition, FF − SS PESQ:**
  - eval: loud −0.012, noisy +0.002, nominal 0.000, quiet +0.001, transitions −0.002;
  - tune: −0.003 … +0.005.
  - Eval loud remains the largest delta (doc 06: −0.010). It is decoder-side and inside the target, and the listening samples should include loud items.
- **Non-speech output levels:** equal to SS (eval silence 18.1, −60 dBFS noise 29.0, −40 dBFS noise 49.2).
- **Encoder:** unchanged, FvsS frames 100% identical on every condition; emulator `packets.bit` md5 51aea44b… identical.
- **Waveform agreement with the float S decoder:** eval median 24.4 dB, as in doc 06. This is the half-bin tie effect explained there.
- **Listening:** not done yet. Listening samples are required before a quality claim is final (doc 03, target 3).

## Cost (emulator, `results/fxp-params-O2`, same 240-frame corpus as doc 02)

| | doc 06 stage | this stage | note |
|---|---:|---:|---|
| decode mean / frame | 679,709 | **643,448** | −5%; softfloat removed |
| decode p95 / max per frame | 710,889 / 716,459 | 676,116 / 681,336 | |
| softfloat in decode | ≈ 30k / frame | **0** | no `__*sf*`, `lrintf`, `frexpf` or `ldexpf` in the decode profile |
| decoder state (`codec2_decoder_t`, CK803S) | 5,448 B | **3,060 B** | scratch 4,096 → 2,056 B; no float `model_t` |
| decoder init | 7.2k instr | **1.9k instr** | |
| stack high water through decode | 2,116 B | ≤ 1,476 B (measured; equals the encoder's high water) | static estimate from `.su`: about 1.4 KB (`codec2_decode` 828 B + synthesis 396 B + irfft/cfft 172 B) |
| `.text` (whole emulator image, encoder included) | 30,388 B | 27,336 B | −3,052 B |
| encode mean / frame | 10,925,126 | 10,925,126 | unchanged |

- **Where the decode goes:** `fxp_lsp_spec_bin` 336.6k (52%), the irfft 109k (17%), the rest of `fxp_aks_to_mag2` 169k, and the synthesis loops 28k. The parameter decode itself is a few hundred instructions.
- **Code size (O2):** `codec2_decode` 340 B, `codec2_decoder_init` 90 B, the integer decode/interp helpers 6–64 B each (200 B in total), `aks_to_mag2_lsp` 64 B.
- **Stack:** `codec2_decode` 828 B, mostly `lA[81]` and `phA[81]`.
- **Sanitizers:** the `-DFXP_CHECK -fsanitize=undefined,address -fno-sanitize-recover=all` build encoded and decoded all 492 corpus items without a report.

## What remains

> **Update (doc 08):** the spectrum kernel is now about half its cost. Decode is 479,584 instructions per frame on the emulator (was 643,448), with output bit-identical on 489 of 492 items.

- **Encoder:** still float (analysis FFT, NLP pitch estimator, LPC/LSP analysis, quantisers). That is the next large stage for a fixed-point M17 transmit path.
- **Decode cost:** half of it is `fxp_lsp_spec_bin` (256 bins × 2 subframes). Evaluating only the bins that the band sums and harmonic phases use, or a cheaper product-form loop, is the obvious profile-led optimisation.
- **Not yet measured:** decode stack in isolation (the harness reports the cumulative high water).
- **Required before firmware integration:** listening samples and the independent review loop.

## Reproduce

```bash
scripts/test-fxp.sh
```

```bash
scripts/host-build.sh
```

```bash
scripts/synth-ref.sh
```

```bash
.venv/bin/python scripts/ref-suite.py eval results/fxp-params-suite/eval
```

```bash
.venv/bin/python scripts/suite-summary.py results/fxp-params-suite/eval
```

Inside the gcc container (repo mounted at /w):

```bash
scripts/emu-build.sh codec results/fxp-params-O2
```

Then on the host:

```bash
scripts/emu-run.sh results/fxp-params-O2 profile
```

```bash
.venv/bin/python scripts/emu-analyze.py results/fxp-params-O2
```

```bash
.venv/bin/python scripts/profile-report.py results/fxp-params-O2
```
