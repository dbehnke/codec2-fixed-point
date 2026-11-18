# 09 — Fixed-point encoder analysis front end (stage E1)

**Status:** the first encoder stage is now integer code. That covers the windowed analysis FFT, the two-stage pitch refinement, the harmonic band energies and the MBE voicing decision. The NLP pitch estimator and the LPC/LSP analysis are still float; they are stages E2 and E3.

Results:
- 99.7% of frames are bit-identical to the float encoder.
- Voicing differs on 0.01–0.02% of frames and Wo on 0.25–0.28%.
- E and the LSPs never differ.
- Quality targets are met on the held-out eval set.
- Encode cost drops from 10.93M to 5.95M emulator instructions per frame (−45.6%).

**Measurement types used here:**
- **host:** macOS arm64, clang. Used for accuracy and quality.
- **emulator:** XUANTIE QEMU instruction counts, `-icount shift=0`. These are instructions, **not cycles**.

Nothing here is DM32 hardware timing, and no real-time claim is made.

## Encoder plan

The float encoder profile (doc 08 stage, emulator, 10.93M instructions per frame) is about 90% softfloat.

| stage | float parts | cost before (inclusive, per frame) | state |
|---|---|---:|---|
| **E1 (this doc)** | `dft_speech` (kiss FFT 512), `two_stage_pitch_refinement`, `estimate_amplitudes`, `est_voicing_mbe`, `encode_Wo` | ≈ 5.2M | **integer** |
| E2 | `nlp` (squaring, DC notch, decimating FIR, its own kiss FFT 512, peak search) | 3.79M | float |
| E3 | `speech_to_uq_lsps` (window, autocorrelation, Levinson, energy, bandwidth expansion, Chebyshev LSP search, `fast_acosf`), `encode_energy`, `encode_lspds_scalar` | 1.97M | float |

In 3200 mode the bitstream carries only voicing ×2, Wo, E and the LSPs. The harmonic amplitudes A[m] are never transmitted; in the encoder they feed only the voicing energy ratios. So E1 has no A[m] output, only band power sums.

## Number formats (`codec/src/analysis.c`)

| quantity | format |
|---|---|
| `Sn[320]` input history | int16. The float encoder held the same integers in float, and the initial value is 1 in both. |
| windowed speech | `Sn · hann_q15` as exact int32 (≤ 2^30), then a rounding right shift s ≥ 0 to the int16 range |
| `Sw[257]` spectrum | cq15 from `fxp_rfft512`, block exponent. The float encoder's Sw_float = Sw · 2^g · scale, with g = e_fft + s − 15 and scale = 1/sqrt(512 Σh²). The scale is the float window normalisation, dropped here. |
| `pw[512]` power spectrum | uint32 round(\|Sw\|²/256), at most 2^23, mirrored to 512 bins because the float refinement reads past Nyquist (up to bin 341) |
| pitch | exact rational p = `pn / pd` samples (Wo = 2π·pd/pn) |
| energies (sig, elow, ehigh, err) | uint64 sums in \|Sw\|² units, compared as log2 Q16 |
| 1e-4 floor of the float sums | `fxp_ana_floor_q24 · 2^(−24 − 2g)` in \|Sw\|² units (`fxp_ana_floor_q24` = 1e-4·512·Σh² in Q24) |

`tools/gen/gen-analysis-tables.py` generates `codec/src/fxp_analysis_tables.c`:
- `fxp_hann_q15[279]`.
- `fxp_winspec_q15[257]`, the 512-point DFT of the centred Hann window normalised to 32767. The voicing fit depends only on ratios of W, so its scale is free.
- `fxp_ana_floor_q24`.

## How each step is done

### Pitch refinement: exact candidates

The NLP returns its peak bin (doc E2 will make it integer). Its pitch is 2560/bin samples.

The float code converted that bin to f0 in Hz and back. Both conversions are exact in float: on host, 0 of 256 bins round-tripped wrongly, and the initial 1/P_MAX_S also gives bin 16. NLP therefore now returns the bin, and its tracking memory is the previous bin (`prev_bin_enc`, initial 16). This NLP change leaves its output bit-identical.

The candidates are exact rationals:

| search | candidates p = N/D | L used |
|---|---|---|
| coarse (p0 ± 5, step 1) | N = 2560 + j·bin, D = bin, j = −5..5 | floor(1280/bin) |
| fine (pc ± 1, step 1/4) | N = 4Nc − 4bin + i·bin, D = 4bin, i = 0..8 | floor(Nc/(2bin)) |

- The harmonic bin is b_m = round(m·512·D/N) with ties up, i.e. floor((1024·D·m + N)/(2N)). It is stepped exactly with a quotient/remainder pair (no division in the inner loop), and E(p) = Σ pw[b_m].
- The largest E wins, the first one on ties. If every E is 0 the input pitch is kept, as in the float code.
- Then the pitch is limited to [20, 160] (N > 160D → 160/1 and N < 20D → 20/1). L = floor(N/(2D)), less one when 2L ≥ 0.95p, tested as 40·L·D ≥ 19·N. The float code applies the same rule.

The float code accumulates `p += pstep` and `currentB += bFloat` in float. Its candidate count and its bins therefore depend on rounding near ties. The exact form removes that dependence, and it accounts for most of the Wo differences below.

### Band energies

The band edges are round((m ± 0.5)·Wo_bin) = floor(((2m ± 1)·512·D + N)/(2N)), stepped exactly and clamped to 256. Each band power Σ pw feeds the 1000/2000/4000 Hz sums (l_1000 = L/4, l_2000 = L/2, l_4000 = L) with the float code's inclusive boundaries. All sums stay in pw units until the log2 compare, which adds 8 to the exponent.

### Voicing fit (harmonics below 1 kHz)

- **Band limits.** al = ceil((2l − 1)·256·D/N) and bl = ceil((2l + 1)·256·D/N). The window offset is floor((513N − 1024·l·D)/(2N)), and the float code's clamps are kept.
- **Window.** The window value is `fxp_winspec_q15[|offset + m − 256|]`, which is the float W[] arrangement.
- **Least-squares amplitude.** Am = Σ Sw·W / Σ W².
  - The products are int32 (≤ 2^30). They are shifted right by 5 before summing, so up to 32 terms fit in int32 (a band has at most 27).
  - Σ W² is normalised to a 16-bit d16 · 2^sd.
  - Each component becomes a 16-bit mantissa a = (x << k) / d16, with one 32-bit hardware divide per component and harmonic.
- **Residual.** The residual is computed directly, Sw − a·W·2^−(k+sd) per bin, rounded with the fit clamped to ±2^24. The sum of squares goes into uint64; a 64-bit multiply is used only when |residual| ≥ 2^16, which is rare. The direct form avoids the cancellation in Σ\|Sw\|² − \|Σ Sw W\|²/Σ W².

### Decisions

`log2_plus_floor(v, ev, ef)` = log2(v·2^ev + floor) in Q16. It aligns two 31-bit mantissas, adds them, and then calls `fxp_log2_pf`; the error is at most 2.4e-5 log2.

The float tests become comparisons of Q16 logs:

| float | integer |
|---|---|
| sig > 10^0.6 · err | lg(sig) > lg(err) + 130624 |
| elow > 10 · ehigh (unvoiced → voiced) | lg(elow) > lg(ehigh) + 217706 |
| elow < 0.1 · ehigh (voiced → unvoiced) | lg(elow) < lg(ehigh) − 217706 |
| elow < 0.398 · ehigh and Wo ≤ 2π·60/8000 | lg(elow) < lg(ehigh) − 87082 and 3N ≥ 400D |

### Wo index

`encode_Wo_p(N, D)` = floor((256·(160D − N) + 7N)/(14N)), clamped to 0..127. This is the float formula floor(128·(Wo − W0_MIN)/(W0_MAX − W0_MIN) + 0.5), evaluated exactly. The float `encode_Wo` is removed.

## Memory

- **Encoder scratch.** `fft_buffer` (4 KB) holds Sw (1,028 B) followed by pw (2,048 B). The int32 windowed products temporarily occupy the pw area before the FFT. A static assert checks the layout. `speech_to_uq_lsps` reuses the same buffer afterwards, as before.
- **Removed from the encoder state:** the float window spectrum `W[512]` (2,048 B) and its kiss FFT at init. `Sn` went from float to int16 (1,280 → 640 B).
  - The float window `w[320]` stays only for the float LPC analysis (E3).
  - The kiss config stays only for the float NLP (E2).
- **Codec state on the CK803S** (`state` symbol, encoder + decoder): 20,028 → **17,340 B**. The encoder accounts for 14,280 B of that (was 16,968).

## Accuracy (host)

`tests/test_analysis.c` is in `scripts/test-fxp.sh` (FXP_CHECK + UBSan + ASan). It includes `analysis.c` to reach the static helpers.

| test | result |
|---|---|
| `dft_speech` Sw·2^g vs a double DFT of Sn·hann, 64 frames, amplitudes 1…32767 | worst spectrum SNR 60.4 dB |
| stepped harmonic bins vs round(m·512·D/N) in double, every coarse and fine candidate of every NLP bin 16..128 | 291,919 of 291,919 equal |
| `log2_plus_floor` vs double, 200,000 random cases | max error 2.44e-5 log2 |
| `encode_Wo_p` vs the exact index, every N/D with D ≤ 512 and 20 ≤ N/D ≤ 160 | 18,386,432 of 18,386,432 exact. The float `encode_Wo` differs on 14, at near-ties in float. |

## Bitstream agreement (host, ref suite, FvsS = this encoder vs the float Codec2-mod split encoder)

| set | frames identical | v1 | v2 | Wo | E | LSP |
|---|---:|---:|---:|---:|---:|---:|
| tune (101,356 frames) | 0.9973 | 1.0000 | 0.9999 | 0.9975 | 1.0000 | 1.0000 |
| eval (113,702 frames) | 0.9970 | 0.9999 | 0.9998 | 0.9972 | 1.0000 | 1.0000 |

- **Comparison with the float encoders.** SvsU (float split vs upstream) is 0.9875 / 0.9840, and FvsU is 0.9851 / 0.9813. The integer front end adds about 0.3% Wo disagreement, against the float encoders' existing 1.2–1.6% frame disagreement with upstream.
- **Tune set Wo differences:** 0.09% of frames are ±1 index and 0.16% are larger. A larger difference is a different winner between near-equal harmonic sums; the float accumulation of `p` and `currentB` makes those choices rounding-dependent.
- **Non-speech** (silence, −40/−60 dBFS white noise): identical.
- **Emulator vs host.** On the 240-frame emulator corpus, the integer analysis gives the same output on both. The one frame that differs (an LSP bit) comes from the float LPC path; the previous encoder showed the same single-frame host/emulator difference.

## Quality (host, `results/fxp-ana-suite/{tune,eval}`)

**Eval set (held out), all speech conditions:**

| metric | SS | FF (doc 08) | **FF (this stage)** | target | met |
|---|---:|---:|---:|---|---|
| PESQ-NB | 2.676 | 2.673 | **2.674** | FF within 0.05 of SS | yes (−0.002) |
| LSD (dB) | 9.011 | 9.002 | **9.004** | ≤ SS + 0.2 | yes (−0.007) |
| level error (dB) | 0.646 | 0.646 | **0.647** | within ±0.5 of SS | yes (+0.001) |

**Per condition, FF − SS PESQ:**

| condition | eval | tune |
|---|---:|---:|
| loud | −0.005 | −0.003 |
| noisy | +0.001 | +0.005 |
| nominal | +0.002 | −0.003 |
| quiet | −0.004 | −0.001 |
| transitions | +0.017 | 0.000 |

- **Tune set overall:** PESQ 2.705 (SS 2.705; doc 08 FF 2.708), LSD 9.093, level error 0.529.
- **Non-speech output levels:** equal to SS.
- **Listening:** not done yet. Listening samples, including loud items, are still required before any quality claim is final.

## Cost (emulator, `results/fxp-ana-O2`, same 240-frame corpus)

| | doc 08 | this stage | change |
|---|---:|---:|---:|
| encode mean / frame | 10,925,126 | **5,945,608** | −45.6% |
| encode p95 / max per frame | 12,887,233 / 12,959,784 | 5,970,811 / 5,982,241 | the spread is now mostly gone |
| integer analysis (2 calls: FFT, refinement, bands, voicing), inclusive minus `nlp` | ≈ 5.2M (float) | **158k** | |
| — of which `fxp_rfft512` | | 107.6k | |
| `nlp` (float, E2) | 3.78M | 3.79M | |
| `speech_to_uq_lsps` (float, E3) | 1.96M | 1.97M | |
| encoder init | 6.68M | 5.07M | the remaining cost is the kiss config and float windows for NLP/LPC |
| stack high water through encode | 1,476 B | **1,212 B** | |
| `.text` (whole image) | 27,540 B | 27,632 B | `analysis.o` 2,378 B; tables 1,076 B `.rodata` |
| decode mean / frame | 479,584 | 479,616 | decoder unchanged; it decodes the new packets |

- **Calls from `analysis.o`:** only `fxp_rfft512`, `fxp_log2_q16`, `nlp`, the mem* functions, and two `__muldi3` sites on the rare |residual| ≥ 2^16 path. The divisions use the CK803S hardware divider; there is no `__divsi3`. The float calls in the object are all in `make_analysis_window` (init only).
- **Sanitizers:** the `-DFXP_CHECK -fsanitize=undefined,address -fno-sanitize-recover=all` build encoded and decoded all 492 corpus items without a report.

## What remains

> **Update (doc 10):** stage E2 (NLP) is integer. Encode is 2.40M instructions per frame on the emulator (was 5.95M); the kiss FFT is no longer linked. Only the LPC/LSP analysis (E3) is still float.

- **E2:** integer NLP, now 64% of encode. Squaring, notch and FIR become int32; the FFT becomes `fxp_rfft512`; the peak search becomes integer. Its output is already an integer bin.
- **E3:** integer LPC/LSP analysis and the E/LSP quantisers. The energy uses the log domain: log2 E = log2 E_hann − log2(Σhann²·512).
- **After E2 and E3:** drop the kiss FFT, the float windows and the libm calls from the encoder.
- **Required before firmware integration:** listening samples, decode stack measured in isolation, and the independent review loop.

## Reproduce

```bash
python3 tools/gen/gen-analysis-tables.py
```

```bash
scripts/test-fxp.sh
```

```bash
scripts/host-build.sh
```

```bash
.venv/bin/python scripts/ref-suite.py eval results/fxp-ana-suite/eval
```

```bash
.venv/bin/python scripts/suite-summary.py results/fxp-ana-suite/eval
```

Emulator, from the repo root:

```bash
docker run --rm -v "$PWD:/w" -w /w p32-csky-gcc14:14.2.0-arm64 sh scripts/emu-build.sh codec results/fxp-ana-O2
```

```bash
scripts/emu-run.sh results/fxp-ana-O2 profile
```

```bash
.venv/bin/python scripts/emu-analyze.py results/fxp-ana-O2
```

```bash
.venv/bin/python scripts/profile-report.py results/fxp-ana-O2
```
