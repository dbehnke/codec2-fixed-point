# 06 — Fixed-point synthesis (phase, postfilter, irfft, overlap-add)

**Status:** the decoder's synthesis stage runs entirely in fixed point (`codec/src/synthesis.c`). It replaces the float `sample_phase`, `phase_synth_zero_order`, `postfilter`, `synthesise` and `ear_protection`, and it takes the fixed-point outputs of doc 05 directly: `log2 A[m]` (Q16) and `arg A` at each harmonic (uint32 angle). The float `A[m]`/`Aw` bridge outputs and the LPC correction in float are gone. Inputs still crossing a float bridge: Wo (converted once per subframe to the integer angle `W`), E and the LSPs (`aks_to_mag2_lsp()` in `codec/src/lpc.c`). The encoder is unchanged.

**Update (doc 07):** those bridges are gone, as are the dead float `aks_to_mag2` and `fxp_bridge_irfft`; `fft_buffer` became a 2 KB `scratch`, and `synthesise_one_frame` takes `(W, L, voiced)` instead of a `model_t`. The figures below are for the stage as committed in doc 06.

**Measurement types used here:**
- **host:** macOS arm64, clang. Used for accuracy and quality.
- **emulator:** XUANTIE QEMU instruction counts, `-icount shift=0`. These are instructions, **not cycles**.

Nothing here is DM32 hardware timing.

## Number formats

| quantity | format | note |
|---|---|---|
| angles (Wo, phases, excitation phase) | uint32, 2^32 = 2π | wrap-around is the float code's mod 2π |
| `W` | uint32 angle of Wo | bridge: `W = lrintf(Wo · 683565275.6f)` (2^32 / 2π) |
| `lA[m]` | log2 A[m], Q16 | `INT32_MIN/2` = empty band (A = 0); "empty" test is `<= INT32_MIN/4` |
| `bg_est` | background level, dB in Q16 | |
| `Sn_[]` | int32 Q4 (sample × 16) | overlap-add memory, 160 words |
| `Pn[]` | int16 Q14 | Parzen window, `Pn[i] = i/80` rising then falling |
| irfft buffer | `cq15_t` with block exponent | `fxp_irfft512` from doc 04, in the decoder `fft_buffer` |

Helpers in `fxp_math.c`: `fxp_sincos_q15` (uint32 angle → Q15 cos/sin, 1 LSB max error) and `fxp_atan2_u32` (max error 2.45e-5 rad), from the generated table `fxp_atan_tab_u32[257]` (`tools/gen/gen-math-tables.py`). `fxp.h` gains `add_sat32` and `wsat` (saturating left shift).

## Design

### Phase (`phase_synth_zero_order`)

- `ex_phase += W · 80` per subframe, in uint32 (the float code's explicit wrap is free).
- Per harmonic: `ph = ex_phase · m` when voiced, otherwise `rand · 131076` (2^32 / 32767 = 131076.00003).
- `th[m] = ph − phA[m]`: the synthesis filter is conj(A), so its phase is subtracted. The float code's normalisation of `H·Ex` to unit magnitude is not needed because amplitude and phase are carried separately.

### Postfilter

- Frame energy in the log domain: `l2 = log2(Σ A² / L)` from the `lA` (terms below 2^-32 of the maximum are skipped), then `e_dB = 3·l2 + mul_q30(l2, 11059494)` (10·log10 2 = 3 + 0.0103). The float code's 1e-12 floor becomes −120 dB.
- Unvoiced, `e < BG_THRESH`: `bg += mul_q30(e − bg, 107374182)` (β = 0.1).
- Voiced: random phase for harmonics with `lA[m] < lthresh`, where `lthresh = mul_q30(bg + 6 dB, 178344657)` converts dB to log2 (log2(10)/20 in Q30). The integer comparison replaces `A[m] < 10^((bg+6)/20)`.

### Spectrum and inverse DFT

- `ip = ceil(max lA)`; each harmonic gets `a14 = round(2^(lA − ip) · 2^14)` via `fxp_exp2_neg_q30`, so the loudest harmonic is 2^14.
- `Sw[b] = round(a14 · (cos, sin)(th) / 2^15)` at bin `b = fxp_bin_round(m · W)` (the tie rule below), clamped to 255.
- `ex = fxp_irfft512(Sw) + ip − 14` is the true exponent of the time signal.

### Overlap-add, ear protection, output

- `Sn_ += sw · Pn` with a single right shift `sh = 14 − 4 − ex` (round-half-up `shr_r`), or a saturating left shift (`wsat`) when `sh < 0`. `add_sat32` for the accumulation. `|sw · Pn| ≤ 2^29`.
- Ear protection above 30000 (`30000 << 4` in Q4): `g = (30000 << 15) / max`, `g2 = g² >> 15`, applied with `mul32x16_q15`. Same rule as float (gain = 1/over²).
- Output truncates toward zero and clamps to ±32767, matching the float code's `(int16_t)` cast. Rounding instead raised quiet-speech levels by 0.01 dB and silence by 0.4 dB against S, so truncation was kept.

### The harmonic-bin tie rule

Codec2 3200 quantises Wo so that bins per harmonic are `(128 + 7k)/40` (or the mean of two such values for interpolated subframes). Harmonics therefore often land **exactly** on a half bin, and `round(m · Wo / (2π/512))` is a tie.

- Float rounds these ties by single-precision noise. For example, at Wo = 0.120877683, harmonic 10 is 98.5 bins: `10 · Wo / FFT_R` comes out as 98.4999… and rounds down, while the harmonic 5 of another frame (19.5 bins) rounds up. The float code even uses two different expressions (`Wo / FFT_R` in `sample_phase`, `Wo · FFT_1_R` in `synthesise`), so it does not always round a given tie the same way in both places.
- No integer rule can reproduce that. We round every tie **up**, with a small bias that also absorbs the error in `W`:

  `fxp_bin_round(x) = (x + 2^22 + FXP_BIN_TIE) >> 23`, `FXP_BIN_TIE = 2^12` (`fxp_lpcspec.h`)

  It is used everywhere a harmonic bin is chosen: the phase sample in the doc 05 kernel, the band edges, and the synthesis bin. Unit test: on the full Wo grid, 0 of 493,458 harmonic bins differ from exact (double) rounding with ties up. Without the bias, 5,654 of the 10,532 ties would round the wrong way.
- Both choices move the harmonic by half a bin (7.8 Hz) relative to the true frequency, in opposite directions; neither is more correct.

**Consequence for waveform comparisons.** A tie that F and S round differently puts that harmonic in a different bin with its phase sampled at a different bin, which changes the waveform but not the quality metrics. On eval nominal file 4, the first large F-vs-S difference (subframe 69) is exactly one such harmonic (10, at 98.5 bins: float bin 98, ours 99; phase differs by 0.23 rad, every other harmonic agrees to 1e-4 rad and 1e-4 dB). The float S decoder has a differing tie in 218 of that file's 846 subframes.

To check that ties are the *only* systematic difference, the HEAD (doc 05) float decoder was rebuilt with `fxp_bin_round` in its three bin computations, and F compared against it:

| comparison (eval speech, 241 items, waveform SNR) | median | p10 | min |
|---|---:|---:|---:|
| F vs float split decoder S | 24.2 dB | 15.8 dB | 11.2 dB |
| F vs HEAD float + same tie rule | **63.7 dB** | 57.1 dB | 15.7 dB |

(Doc 05's F vs S median was 58.3 dB because there both used the float bin expressions.) The remaining minimum is the postfilter random-phase effect (doc 04): a borderline decision shifts the random sequence. On file 4 the random-number state first diverges at subframe 210.

## Accuracy (host)

### Synthesis arithmetic vs double (`scripts/synth-ref.sh`, `results/fxp-synth-ref/summary.md`)

A debug build (`-DFXP_SYNTH_DEBUG`) dumps, per subframe, the inputs synthesis actually used (L, W, voicing, lA, final phases, e, bg). `tests/synth/synth_ref.py` re-synthesises the same frames in double (irfft, window, overlap-add, ear protection, truncation) and compares. The script first checks that the debug build's output is byte-identical to the normal build's.

| set | items | SNR median | SNR min | seg SNR median | seg p1 | max abs err | e err max | bg err max | postfilter decisions | flips (arith) | flips (bg) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| eval speech | 241 | 58.78 dB | 52.02 dB | 57.77 dB | 49.90 dB | 196 | 0.0001 dB | 0.0001 dB | 5,936,591 | 4 | 3 |
| tune speech | 245 | 58.74 dB | 49.55 dB | 57.56 dB | 49.71 dB | 161 | 0.0001 dB | 0.0001 dB | 5,261,004 | 6 | 5 |

- SNR = synthesis arithmetic only (the Q15 irfft dominates; doc 04 measured it at 58–62 dB). Segmental SNR uses 10 ms segments above −50 dBFS.
- "flips (arith)": the double rule on the same `bg` disagrees with the integer rule. "flips (bg)": a float-tracked `bg` would decide differently. Both are a handful per million.
- Quiet speech has the lowest SNR (median 55 dB) because the error floor is the output truncation; its maximum error is 11 LSB.

### Unit tests (`scripts/test-fxp.sh`, FXP_CHECK + UBSan + ASan)

- `fxp_atan2_u32` < 3e-5 rad (measured 2.45e-5), `fxp_sincos_q15` < 1.6 LSB.
- Kernel `arg A` at harmonic bins: max 1.44e-4 rad.
- Wo-grid bin test (above): 0 differences.
- Band edges for random off-grid Wo: 915 differ from the float formula. That is expected: the bias shifts the rounding point by 0.0005 bin.

## Quality (host, `results/fxp-synth-suite/{tune,eval}`)

Candidate F = the `codec/` tree (`build/host/c2fx`).

**Eval set (held out), all speech conditions:**

| metric | SS | FF (doc 05) | **FF (this stage)** | UU | target | met |
|---|---:|---:|---:|---:|---|---|
| PESQ-NB | 2.676 | 2.678 | **2.674** | 2.677 | FF within 0.05 of SS | yes (−0.002) |
| LSD (dB) | 9.011 | 9.005 | **9.002** | 9.022 | ≤ SS + 0.2 | yes (−0.009) |
| level error (dB) | 0.646 | 0.647 | **0.646** | 0.650 | within ±0.5 of SS | yes (0.000) |

**Tune set:** PESQ 2.708 vs SS 2.705, LSD 9.096 vs 9.103, level error 0.528 vs 0.528.

- **Per condition, FF − SS PESQ:** eval loud −0.010, noisy +0.002, nominal +0.001, quiet +0.001, transitions (1 item) −0.002. Tune −0.003 … +0.006.
  - The eval-loud −0.010 is the same for every bitstream decoded by F (UF, MF, SF all 2.638), so it comes from the decoder, most likely the tie direction. Tune loud is +0.002. It is well inside the 0.05 target, but it is the largest per-condition delta so far; listening should include loud items.
- **LSD:** FF ≤ SS on 9 of 10 set/condition pairs (eval: loud −0.014, nominal −0.012, transitions −0.016); the exception is tune nominal, +0.003 dB.
- **Non-speech output levels (eval):** silence 18.1 (FF) vs 18.1 (SS); −60 dBFS noise 29.0 vs 29.0; −40 dBFS noise 49.2 vs 49.2.
- **Encoder:** unchanged, FvsS frames 100% identical on every condition; emulator `packets.bit` md5 51aea44b… identical to doc 05.
- **Listening:** not done yet. Listening samples remain required before a quality claim is final (doc 03, target 3).

## Cost (emulator, `results/fxp-synth-O2`, same 240-frame corpus as doc 02)

| | doc 05 stage | this stage | note |
|---|---:|---:|---|
| decode mean / frame | 3,063,526 | **679,709** | −78% |
| decode p95 / max per frame | 4,973,569 / 6,100,183 | 710,889 / 716,459 | |
| `synthesise_one_frame` inclusive | — | 68.4k per call | 2 calls per frame |
| of which `fxp_irfft512` | ≈ 56.6k per call | 54.6k per call | |
| `aks_to_mag2_lsp` inclusive | — | 258.4k per call | kernel 252.6k (doc 05) + 5.8k bridge |
| stack high water through decode | 4,268 B | **2,116 B** | |
| decoder state (`codec2_decoder_t`, CK803S) | 16,592 B | **5,448 B** | no kiss configs, int32/int16 Sn_/Pn |
| decoder init | 6.57M instr | **7.2k instr** | no `kiss_fftr_alloc` |
| encode mean / frame | 10,925,126 | 10,925,126 | unchanged |

- **Where the decode now goes:** `fxp_lsp_spec_bin` 336.6k (50%), the irfft 109k (16%), the rest of `fxp_aks_to_mag2` 169k. Remaining softfloat in decode is about 30k per frame (the Wo/E/LSP bridges); `__muldi3` is only reached from softfloat `__mulsf3`.
- **Code size (O2):** `synthesise_one_frame` 1312 B, `fxp_sincos_q15` 352 B, `fxp_atan2_u32` 264 B, `fxp_atan_tab_u32` 1028 B.
- **Stack:** `synthesise_one_frame` 396 B (`th[81]`), `fxp_sincos_q15` 12 B, `fxp_atan2_u32` 8 B.
- **Sanitizers:** the `-DFXP_CHECK -fsanitize=undefined,address -fno-sanitize-recover=all` build of `codec/` encoded and decoded all 492 corpus items (tune and eval) without a report.

## Known follow-ups

- Shrink `fft_buffer`: it keeps its float-era size, but synthesis now only needs `FXP_RFFT_BUF` `cq15_t` words and the doc 05 kernel's int32 scratch.
- Remove the dead float `aks_to_mag2` in `lpc.c` and the now-unused `fxp_bridge_irfft`.
- Port the remaining float bridges: Wo → `W`, E → `log2E`, LSPs → `hz2` (the decoder's dequantisation and interpolation).
- Profile-led optimisation of `fxp_lsp_spec_bin` (half of the decode).

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
.venv/bin/python scripts/ref-suite.py eval results/fxp-synth-suite/eval
```

```bash
.venv/bin/python scripts/suite-summary.py results/fxp-synth-suite/eval
```

Inside the gcc container (repo mounted at /w):

```bash
scripts/emu-build.sh codec results/fxp-synth-O2
```

Then on the host:

```bash
scripts/emu-run.sh results/fxp-synth-O2 profile
```

```bash
.venv/bin/python scripts/emu-analyze.py results/fxp-synth-O2
```

```bash
.venv/bin/python scripts/profile-report.py results/fxp-synth-O2
```
