# 11 — Fixed-point LPC/LSP analysis and E/LSP quantisers (stage E3)

**Status:** the whole encoder is integer code. No float, softfloat, libm or 64-bit multiply is linked into the codec: the emulator image (encoder and decoder) contains no `__*sf*` helpers, `sqrtf`, `cosf`, `log10f` or `__muldi3`, and `scripts/emu-build.sh` now fails if any codec object references one (see "Review fixes").

Results:
- The energy and LSP indices agree with a double-precision run of the same algorithm far more often than the float code they replace did. Eval LSP index mismatches vs double: 0.120% integer, 0.599% float. Energy: 0.016% integer, 0.229% float.
- Every float encoder (upstream libcodec2, Codec2-mod main and split) sits about 0.6% from double on LSP indices. The integer encoder's distance from upstream (0.661%) is therefore about the same as double's own (0.618%); see "Agreement with the float encoders".
- Quality targets are met on the held-out eval set. FF PESQ 2.676 equals SS, LSD 9.007 dB, level error 0.647 dB. Upstream decoding of the integer bitstream (FU) scores PESQ 2.676 against UU 2.677.
- Encode cost drops from 2.40M to 0.50M emulator instructions per frame (−79.2%). The LPC/LSP analysis drops from 1.97M to 97k (−95.1%).
- Encoder init drops from 1.48M to 1,666 instructions; the float window generation is gone.
- `.text` drops from 23,764 to 15,944 B (softfloat and libm are no longer linked) and `.bss` from 8,028 to 6,748 B. The codec's own `.rodata` (tables) is 8,961 B.

**Measurement types used here:**
- **host:** macOS arm64, clang. Used for accuracy and quality.
- **emulator:** XUANTIE QEMU instruction counts, `-icount shift=0`. These are instructions, **not cycles**.

Nothing here is DM32 hardware timing, and no real-time claim is made.

## Algorithm

The algorithm is the float Codec2-mod encoder's `speech_to_uq_lsps`, `encode_energy` and `encode_lspds_scalar`, step for step (`codec/src/lpc.c`, `codec/src/quantise.c`):

1. Hann-window the 279 centre samples of `Sn` (samples 21..299). The float window was the Hann window times 1/sqrt(512·m), with m = Σ hann².
2. If the windowed energy Σ(Sn·w)² < 1e-6, use the benign LSPs π·i/10 and E = 0.
3. Autocorrelation R[0..10], then Levinson–Durbin to order 10. A reflection coefficient with |k| > 1 is replaced by 0.
4. E = Σ a_i·R_i, taken before bandwidth expansion.
5. 15 Hz bandwidth expansion, a_i·0.994^i.
6. LSP roots: P'(x) and Q'(x) as Chebyshev sums, scanned from x = 1 downwards in steps of 0.01. Each sign change is bisected 6 times, and the search alternates between P' and Q'. If fewer than 10 roots are found, the benign LSPs are used.
7. w = fast_acosf(x), the codec's own polynomial approximation (not exact acos).
8. Energy index: round(32·(10·log10 E + 10)/50), clamped to 0..31.
9. LSP index: scalar delta quantisation in Hz. Each target is the LSP minus the sum of the previous quantised deltas, and takes the nearest codebook entry; ties go to the lower index.

## Number formats

| quantity | format |
|---|---|
| windowed speech x[279] | int16, round(Sn·h_q15 / 2^s). s ≥ 0 is the smallest shift giving \|x\| ≤ 32767 (the `dft_speech` rule). Float equivalent: Wn = x·2^(s−15)·scale. The window is the existing `fxp_hann_q15`. |
| R[0..10] | int64, exact Σ x[i]·x[i+j]. Products are int16×int16 < 2^30; the accumulation uses adds only. |
| r[0..10] | int32, R·2^−sh rounded, with r[0] in [2^29, 2^30) |
| reflection coefficients rc | Q30, by the Schur recursion (below) |
| prediction error E | mantissa C1[0] (~2^28) with exponent −S, relative to r |
| log2 E | Q16: log2(C1[0]) − S + sh + 2s − 30 + log2(1/(512·m)) |
| LPC a[0..10] | Q23, by the step-up recursion from rc. For any minimum-phase A(z), \|a_j\| ≤ C(10, j) ≤ 252 < 256, so Q23 cannot overflow. |
| bandwidth expansion g^i | Q30, `fxp_bw_gamma_q30` (the float32 values of the float table) |
| P, Q | Q23. The sums a_i ± a_(11−i) are formed modulo 2^32: they can exceed int32, but the results (\|P_i\|, \|Q_i\| ≤ 252) cannot. |
| Chebyshev coefficients | {P0..P4, P5/2}: the float's {2·P0..2·P4, P5}, halved, which keeps the sign of the sum. Block-normalised per polynomial so that every \|c\| < 2^27. |
| x grid | integer u in units of 1/6400. The float grid 1 − 0.01·n is u = 6400 − 64·n, and the 6 bisections of a 64-wide interval land exactly on integers. |
| x, T_k(x), sums | Q29. x = round(u·2^29/6400). \|T_5(1.01)\| ≈ 1.26, and a sum of six terms stays < 2^30. Only the signs of the sums are used. |
| acos | Q29. sqrt(1 − x) = isqrt(d·2^44)·8/5 with d = 6400 − \|u\|, using an exact bitwise 64-bit integer square root; the polynomial p(x) is evaluated by Horner in Q29 (`fxp_acos_poly_q29`). |
| LSPs | Hz in Q16. For x < 0, hz = 4000 − hz(\|x\|), which avoids π in Q29. |
| benign LSPs | 400·i Hz exactly. They quantise to the same indices as the float π·i/10 (checked in the unit test). |

**Multiplies.**
- `mul_q29` and `mul_q30` form (a·b) >> 29 or >> 30 from four 16-bit partial products. There is no 64-bit multiply and no `__muldi3`. The error is at most 2 LSB (tested).
- The divide for each reflection coefficient is a 31-step restoring division giving round(n·2^30/d) exactly.
- The energy quantiser multiplies log2 E by the slope 0.64·10·log10(2) in the same 16-bit split form, then takes floor(v + 6.9).

**Schur recursion.** The reflection coefficients come from the Schur recursion (the form of Opus `silk_schur`), not from Levinson on the autocorrelation:
- Two generator rows C0, C1 start equal to r.
- At order k: rc = −C0[k+1]/C1[0]. Then C0[n+k+1] += rc·C1[n] and C1[n] += rc·C0[n+k+1] (old values) for n = 0..9−k.
- Every entry is bounded by the current prediction error C1[0]. After each order, the entries still in use are renormalised by a common shift so that the largest is just below 2^29, and the shift is accumulated in S.
- So each order keeps about 29 bits relative to its own prediction error, however large the prediction gain.
- In exact arithmetic the result equals Levinson's, including the |k| > 1 → 0 rule: rc = 0 leaves the rows unchanged, as k = 0 leaves Levinson's state unchanged.

The final prediction error is C1[0]·2^−S. In exact arithmetic this is the float code's E = Σ a_i·R_i.

**Zero-energy trap.** The float condition Σ(Sn·w)² < 1e-6 becomes R[0] < 57,312,040 at s = 0 (`fxp_lpc_floor_r0` = ceil(1e-6·2^30·512·m)). With s ≥ 1, R[0] ≥ 2^28, which is far above the threshold, so the test only runs at s = 0. The trapped frame returns `LOG2E_ZERO`, which quantises to energy index 0, like the float E = 0.

**Constants** come from `tools/gen/gen-lpc-tables.py`, which writes `codec/src/fxp_lpc_tables.c`. m is the float32 Hann energy summed in float32 (104.24999…), as the float encoder computed it.

## Interfaces

| old (float) | new |
|---|---|
| `speech_to_uq_lsps(c2, lsp, ak, &e, Sn, w)` | `int32_t speech_to_lsps_q16(c2, lsp_hz_q16)` returns log2 E in Q16 (or `LOG2E_ZERO`) |
| `encode_energy(float e, bits)` | `int encode_energy_log2(int32_t log2E_q16)` |
| `encode_lspds_scalar(idx, const float *lsp)` | `void encode_lspds_q16(idx, const int32_t *lsp_hz_q16)` |

**Removed:**
- `c2->w` (1,280 B of encoder state), `make_analysis_window` and `analysis_init`.
- The float `bw_gamma` table and the float `lsp_to_lpc` (it had no callers since the decoder port).
- `codec/src/util.c` (`fast_acosf`, `fast_cosf`, `fast_atan2f`; no callers remain). `util.h` keeps the inline helpers that synthesis still includes.

## Accuracy (host)

**Unit tests** (`tests/test_lpc.c`, run by `scripts/test-fxp.sh` with `-DFXP_CHECK`, UBSan and ASan):

| test | result |
|---|---|
| `mul_q29`, `mul_q30` vs exact floor (10^6 random pairs) | max error 2 LSB |
| `div_q30` vs round(n·2^30/d), `isqrt64` vs exact | exact on 10^6 cases each, plus edge cases |
| acos to Hz over the whole grid u = −6464..6464 vs double fast_acosf | max error 0.00009 Hz |
| `encode_energy_log2` vs the double formula (sweep of log2 E in 7/65536 steps) | 0 mismatches away from a decision level; 8 within 2^−13 index units of one |
| `encode_lspds_q16` vs the double scalar quantiser (200,000 random LSP sets, half on exact codebook midpoints) | identical |
| Schur vs double Levinson on the same r (20,000 windowed AR(2) signals at levels 2^0..2^14, some clipped) | \|rc − k\| rms 2.6e-7, max 1.7e-5; float32 Levinson on the same r: rms 2.6e-6, max 1.6e-4 |

Absolute rc errors grow with the conditioning of r. The worst cases are clipped, strongly resonant signals, and the integer code beats float on both max and rms.

**Per-frame agreement with double** (`test_lpc` with the corpus files as arguments). The double reference runs the same algorithm, including fast_acosf, the 0.01 grid and the 6 bisections, on the same `Sn`. It uses an exact double Hann window, an exact 0.01 step and an exact 1e-6 floor. The float column is the float E2 code in float32 (`tests/lpc_ref.h`).

| set | frames | E index ≠ double: int / float | LSP index ≠ double: int / float | any E/LSP field ≠ double: int / float |
|---|---:|---:|---:|---:|
| tune | 101,356 | 0.026% / 0.150% | 0.096% / 0.423% | 0.550% / 2.250% |
| eval | 113,702 | 0.016% / 0.229% | 0.120% / 0.599% | 0.686% / 3.120% |

- **LSP error vs double** (frames where double finds 10 roots): rms 0.13 Hz. The maximum is 16.5 Hz (tune) or 11.3 Hz (eval).
- **log2 E error:** max 0.015 (0.05 dB).
- **Where the maxima come from:** the root search itself, at the first LSP. When a root lies within a few thousandths of a grid unit of a bisection point, the sign of the sum there is decided below the analysis precision. The final midpoint then moves by 2 grid units, which is about 5 Hz near 100 Hz; flips at coarser steps give the larger maxima. The double reference flips just as readily under its own tiny perturbations, and so does the float.

**Where the remaining integer-vs-double differences come from.** The largest single source is the int16 rounding of the windowed speech, which adds noise at about the input's own LSB level. An experiment computed R exactly from the int32 products Sn·h_q15 (128-bit accumulation, host only). Mismatches vs double fell from 0.096% to 0.044% (LSP) and from 0.026% to 0.011% (E) on tune. The cost on the target would be about 3× the autocorrelation work (a split hi/lo product). It is not done: the integer path is already 4–5× closer to double than the float encoders. The option stays available if listening shows a need.

## Agreement with the float encoders

The ref-suite reports lower field agreement with upstream than in E2. For example, eval nominal FvsU LSP field agreement is 0.9908 (E2: 0.9968), and frame identity is 0.9461 (E2: 0.9756). Comparing the analyses' indices directly with each encoder's bitstream shows why (`scripts/lpc-vs-encoders.py`, eval, 113,702 frames):

| pair | E index differs | LSP index differs |
|---|---:|---:|
| upstream vs Codec2-mod main | 0.009% | 0.057% |
| upstream vs Codec2-mod split | 0.070% | 0.220% |
| upstream vs double | 0.236% | 0.618% |
| upstream vs **integer** | 0.244% | 0.661% |
| mod split vs double | 0.234% | 0.599% |
| **integer vs double** | **0.016%** | **0.120%** |

- **The float encoders cluster.** They round the same operations the same way, so they agree with each other on near-ties. All of them sit about 0.6% away from the exact algorithm.
- **The integer encoder sits next to double.** Its distance from upstream is double's distance plus about 0.04%.
- **Tune shows the same pattern:** upstream vs double 0.427%, upstream vs integer 0.471%, integer vs double 0.096%.
- **What a differing index is:** a near-tie resolved the other way. It is not an error in the bitstream format, and any standard decoder decodes it: FU PESQ 2.676 vs UU 2.677 on eval.
- **The test's float column is not bit-exact with the split build.** It differs from S on 0.064% of LSP indices, from compile and contraction differences. It is a float yardstick, not a replica.

## Bitstream change vs E2 (host)

| set | frames differing | E | LSP | voicing / Wo |
|---|---:|---:|---:|---:|
| tune | 2.49% | 0.15% | 2.43% | never |
| eval | 3.32% | 0.24% | 3.25% | never |

Silence and white-noise items are unchanged.

## Quality (host, `results/fxp-lpc-suite/{tune,eval}`)

**Eval set (held out), all speech conditions:**

| metric | SS | UU | FF (doc 10) | **FF (this stage)** | **FU (this stage)** | target | met |
|---|---:|---:|---:|---:|---:|---|---|
| PESQ-NB | 2.676 | 2.677 | 2.674 | **2.676** | **2.676** | FF within 0.05 of SS | yes (0.000) |
| LSD (dB) | 9.011 | 9.022 | 9.006 | **9.007** | 9.020 | ≤ SS + 0.2 | yes (−0.004) |
| level error (dB) | 0.646 | 0.650 | 0.647 | **0.647** | 0.649 | within ±0.5 of SS | yes (+0.001) |

**Per condition, FF − SS PESQ:**

| condition | eval (doc 10) | **eval** | tune (doc 10) | **tune** |
|---|---:|---:|---:|---:|
| loud | −0.010 | **−0.003** | −0.004 | **−0.006** |
| noisy | +0.002 | **+0.002** | +0.010 | **+0.009** |
| nominal | +0.004 | **+0.006** | +0.001 | **−0.002** |
| quiet | −0.004 | **−0.005** | +0.003 | **−0.007** |
| transitions | +0.017 | **+0.035** | 0.000 | **+0.007** |

- **Tune set overall:** PESQ 2.704 (SS 2.705; doc 10 FF 2.708), LSD 9.086, level error 0.527.
- **Eval loud** improved from −0.010 to −0.003.
- **Every per-condition change is within ±0.01** except transitions (+0.018 eval, +0.007 tune), where there are few items.
- **Listening (2026-11-17, one listener, informal):** the 13-item set from `scripts/listening-samples.py results/fxp-lpc-suite/eval build/listening` (seed 17) at commit 4f82b00: nine eval speech items, including each condition's largest FF − SS PESQ drop, and four generated items (100 Hz and 1 kHz tones, DTMF, speech with DC offset and hum). The listener heard "very subtle differences but very similar", judged the fixed-point codec good, and preferred ours over the float reference in more of the blind A/B pairs than not. Per-pair answers were not recorded, and one listener is not a formal test.

**Sanitizers.** All 492 corpus items (tune and eval) were encoded and decoded with `-DFXP_CHECK`, UBSan and ASan, with no recovery. No findings.

## Cost (emulator, `results/fxp-lpc-O2`, same 240-frame corpus)

| | doc 10 | this stage | change |
|---|---:|---:|---:|
| encode mean / frame | 2,401,292 | **498,792** | −79.2% |
| encode p95 / max per frame | 2,428,888 / 2,443,764 | 521,075 / 522,983 | |
| LPC/LSP analysis, inclusive | 1.97M (`speech_to_uq_lsps`) | **97,082** (`speech_to_lsps_q16`) | −95.1% |
| — of which Chebyshev sums | | 37,225 (266 calls, about 140 instr each) | |
| E/LSP quantisers | | 854 (`encode_lspds_q16`; energy inlined) | |
| encoder init | 1,481,464 | **1,666** | the float window is gone |
| stack high water through encode | 588 B | 672 B | `speech_to_lsps_q16` frame 496 B, with Schur and the root search inlined |
| `.text` (whole image) | 23,764 B | **15,944 B** | softfloat, libm and `__muldi3` are no longer linked |
| `.rodata`, codec objects only | | **8,961 B** | the image's 85,765 B includes the harness's 76,800 B test corpus |
| `.bss` (codec state) | 8,028 B | **6,748 B** | |
| encoder state, host `sizeof` | 5,784 B | 4,504 B | `w[320]` float removed |
| decode mean / frame | 479,616 | 479,617 | decoder unchanged |

**Encode now.** `nlp` 242.6k (48.6%), the two 512-point FFTs 107.6k, the rest of the analysis 47.5k, LPC/LSP 97.1k.

**Library calls and symbols:**
- `lpc.o` calls only `memcpy`, `fxp_log2_q16` and `fxp_aks_to_mag2` (decoder), plus const tables. `quantise.o` calls nothing outside the codec. Neither calls `__muldi3`, softfloat or a divide helper.
- `__muldi3` is no longer linked: the rare wide-residual square in `analysis.o` (`sq64`) is now three 32-bit products (exact on 5·10^7 random and edge inputs).
- No `__*sf*`, `sqrtf`, `cosf`, `log10f` or `__muldi3` symbol remains in the image. Doc 10 had all of them.

**Code size per object:**

| object | doc 10 | this stage |
|---|---:|---:|
| `lpc.o` | 2,922 B | 2,812 B |
| `quantise.o` | 743 B | 542 B |
| `analysis.o` | 2,366 B | 2,358 B |
| `util.o` | 564 B | removed |

**Packets:**
- `packets.bit` (md5 dadc5700…) is identical to a host `c2fx` encode of the same 240 frames, so host and emulator are bit-exact. It is unchanged by the review fixes.
- It differs from doc 10 by the E3 bitstream change.

## Scratch memory

`speech_to_lsps_q16` uses the shared encoder scratch (`int32_t[769]`) for the int16 window x[279] and the int32 products [279]. Everything else is on the stack. It runs after `analyse_one_frame`, when the NLP and analysis buffers in the scratch are dead.

## Independent review

A Fable reviewer (standing in for Astra, which is not available here) reviewed the whole fixed-point codec at commit b7afe9c, docs 01–11 and the evidence. Verdict: **ACCEPT**, no blockers or majors, seven minors. Its probes: 200,000 random 8-byte frames plus structured edge frames (E 0/31, Wo 0/127, all LSP indices 0/31, every voicing pattern) through the decoder, and silence, DC, full-scale square, Nyquist, tones, impulses, white noise and re-init through the encoder, all under `-DFXP_CHECK`, UBSan and ASan, with no report.

**Review fixes.** None changes the output: encode and decode of all 492 corpus items are byte-identical to b7afe9c, and so is the emulator `packets.bit`.

| finding | fix |
|---|---|
| `decode_lspds_hz2` comment claimed hz2 < 8000; bit-error frames reach 19600 | comment states the real bound (fits int16) and that LSPs above π alias like `cosf()` in the float code |
| `lpc.c` normalisation wrote a left shift as a 64-bit multiply (folded at -O2, not guaranteed) | written as a shift |
| `__muldi3` linked for `sq64`'s never-observed wide path | three 32-bit products; `__muldi3` is gone from the image |
| stale "LPC is still float" comment; `.rodata` not reported; `<math.h>`, float macros, `model_t`, float helpers in `util.h` still in the codec headers | comment fixed; codec-only `.rodata` reported; float constants moved to `tests/ref/c2_float_consts.h`, the rest deleted; `util.h` keeps only `codec2_rand` |
| kiss FFT sources in `codec/` (dead on target, removed only by gc-sections) | moved to `tests/ref/` (they are the FFT test reference) |
| shift ranges in synthesis overlap-add and `fit_residual` unasserted | `FXP_ASSERT(sh >= -31)` in both |
| (enforcement) | `scripts/emu-build.sh` fails if a codec object references softfloat, libm or `__muldi3` |

**Degenerate inputs differ from the float encoder by design.** For near-pure tones (≥ 60 dB SNR), pure DC and Nyquist alternation the LPC problem is ill-conditioned, and both encoders often fall back to the benign LSPs. Ours does so steadily; the float encoder jitters between benign and non-benign frames, which gives 30 dB level bursts on a clean 100 Hz tone. A full-scale DC input gives energy index 0 here (the Schur error reaches the floor) and 23 in the float code, whose E there is rounding noise. At 20–40 dB SNR, tones encode identically. None of this is worse than the float encoder, but test tones, DTMF, ADC offset and mains hum are realistic radio inputs, so the listening set includes them.

**Bit-error frames.** On invalid or corrupted frames (for example E 31, Wo 127, all LSP indices 0) the three decoders disagree by several dB, because the decoded LPC spectrum is near-singular; ours sits between upstream and the float split decoder, with no clipping or undefined behaviour. Matching the float decoder on such frames is not a goal.

## What remains

- **Required before firmware integration:** decode stack measured in isolation. The independent review has accepted, and informal listening found the fixed-point codec good.
- **Optional:**
  - Exact wide autocorrelation (see Accuracy). It roughly halves the remaining differences from double, for about 3× the autocorrelation cost.
  - NLP DFT savings (radix-2 split or fewer bins), now the largest encode cost.
  - Trimming the 496 B stack frame of `speech_to_lsps_q16`.

## Reproduce

```bash
python3 tools/gen/gen-lpc-tables.py
```

```bash
scripts/test-fxp.sh
```

```bash
scripts/host-build.sh
```

```bash
.venv/bin/python scripts/ref-suite.py eval results/fxp-lpc-suite/eval
```

```bash
.venv/bin/python scripts/suite-summary.py results/fxp-lpc-suite/eval
```

Per-frame agreement with double, and the comparison with the float encoders (the dump directory must exist):

```bash
cc -O2 -Icodec/inc -Itests/ref tests/test_lpc.c $(ls codec/src/*.c | grep -v '/lpc.c$') -lm -o build/test/tlpc && LPC_DUMP=build/test build/test/tlpc corpus/eval/*.raw
```

```bash
.venv/bin/python scripts/lpc-vs-encoders.py build/test corpus/eval/*.raw
```

Emulator, from the repo root:

```bash
docker run --rm -v "$PWD:/w" -w /w p32-csky-gcc14:14.2.0-arm64 sh scripts/emu-build.sh codec results/fxp-lpc-O2
```

```bash
scripts/emu-run.sh results/fxp-lpc-O2 profile
```

```bash
.venv/bin/python scripts/emu-analyze.py results/fxp-lpc-O2
```

```bash
.venv/bin/python scripts/profile-report.py results/fxp-lpc-O2
```
