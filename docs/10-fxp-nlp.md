# 10 — Fixed-point NLP pitch estimator (stage E2)

**Status:** the NLP pitch estimator is integer code. The whole encoder is now integer except the LPC/LSP analysis and the E/LSP quantisers (stage E3).

Results:
- The pitch bin agrees with a double-precision NLP on 99.74–99.89% of frames, in every corpus condition and at every level down to digital near-silence.
- Against E1, 0.12% (tune) and 0.14% (eval) of frames differ. Voicing differs on ≤ 0.03% of frames; E and the LSPs never differ.
- Quality targets are met on the held-out eval set.
- Encode cost drops from 5.95M to 2.40M emulator instructions per frame (−59.6%). The NLP itself drops from 3.79M to 0.24M.
- Codec RAM (`.bss`) drops from 19,596 to 8,028 B. The kiss FFT is no longer linked.

**Measurement types used here:**
- **host:** macOS arm64, clang. Used for accuracy and quality.
- **emulator:** XUANTIE QEMU instruction counts, `-icount shift=0`. These are instructions, **not cycles**.

Nothing here is DM32 hardware timing, and no real-time claim is made.

## Algorithm

The algorithm is the float Codec2 NLP, step for step (`codec/src/nlp.c`):

1. Square the 80 new samples.
2. DC notch y[n] = x[n] − x[n−1] + 0.95·y[n−1], then add 1.
3. 48-tap low-pass, decimated by 5, giving 16 outputs per frame. The history holds 64 outputs (4 frames).
4. 64-sample Hann window w[i] = 0.5 − 0.5·cos(2πi/63).
5. 512-point zero-padded DFT and its power.
6. Global peak over bins 16..128.
7. Sub-multiple search with pitch tracking (thresholds 0.3·gmax, or 0.15·gmax near the previous bin; local-peak test).

The output is the DFT bin (pitch period 2560/bin samples), as since E1.

## Number formats

| quantity | format |
|---|---|
| `sq[]` notch output | int32: y + 1 in units of 2^−k, \|y\| ≤ 2^30. k (`nk`, 0..28) is a per-frame notch scale; see below. |
| notch state `mem_x`, `mem_y` | int32, same units. `mem_e` is the remainder (0..19, in 1/20 units) of the 0.95 multiply. |
| low-pass taps | Q15 int16, `fxp_nlp_fir_q15[48]`, linear phase, sum 32768, Σ\|h\| = 45096 |
| `sq_fir[64]` low-pass output | int32 block-float. Each frame's 16 outputs share an exponent `fir_e` = k + se (value = sq_fir·2^−fir_e), with se ≤ 15 the largest that keeps them below 2^30. |
| window | Q15 int16, `fxp_nlp_hann_q15[64]` |
| DFT input x[64] | int16: sq_fir·w aligned across the four blocks and normalised per frame (max\|x\| ≤ 32767, Σ\|x\| < 2^17) |
| DFT output | exact int32 at bins 15..128, using `fxp_cos512_q14` |
| power P | uint32, normalised per frame (\|re\|, \|im\| < 2^15 before squaring) |
| thresholds | integer: P > 0.3·g ⇔ P > ⌊3g/10⌋, computed without overflow as g/10·3 + (g%10)·3/10 |

Every decision compares powers within one frame, so no absolute scale is carried past the notch.

### Notch and its scale

- **0.95·y** is formed as ⌊(19y + e)/20⌋, with the remainder e carried to the next sample (error feedback). There is no dead band, and the error stays below one unit. Writing y = 20q + r keeps 19y from overflowing: (19y + e)/20 = 19q + (19r + e)/20.
- **Range.** With z = y − x, z[n] = 0.95·z[n−1] − 0.05·x[n−1]. So for 0 ≤ x ≤ X, \|z\| ≤ max(\|z0\|, X) and \|y\| ≤ max(\|z0\|, X) + X.
- **Per-frame scale k.** k is the largest value (≤ 28) for which max(\|z0\|, X, mem_x, FIR history)·2^k ≤ 2^29. That keeps \|y\| ≤ 2^30. Loud frames run at k = 0. Near-silent frames keep up to 28 fractional bits.
- **Why the scale is needed.** Without it (k = 0 always), the notch's integer rounding decided the pitch on near-silent frames. 47% of frames with max\|Sn\| < 4 and 11% with max\|Sn\| < 16 disagreed with double.
- **Rescale on a change of k.** The notch state and the 47-sample FIR history are rescaled, and the remainder e is dropped (less than one old unit).
  - The FIR history is included in the range test. Without it, a peaky loud frame followed by a quiet one overflowed the history on the left shift. That was found when loud-set agreement fell from 99.77% to 99.38%. An `FXP_ASSERT` now guards the shift, and a mutation test confirms that the unit test's loud/silent step signal triggers it.

### Low-pass

- Exact: each input is split as v = hi·2^15 + lo with lo in [0, 2^15).
- Σ h·hi and Σ h·lo both fit int32: \|hi\| ≤ 40961 (from \|v\| ≤ 2^30 + 2^28 + 8), and 40961·45096 < 2^31.
- Combined, they give ⌊Σ h·v / 2^15⌋ and a 15-bit fraction. The block exponent se then keeps up to 15 of those fraction bits, with one rounding.
- `test_nlp` checks 96,000 outputs against a 64-bit reference: all equal.

### Why a direct DFT and not `fxp_rfft512`

The decisions pick maxima of a smooth (8× zero-padded) spectrum and compare neighbouring bins, so near-ties are common.

Diagnostic results (host, tune set):

| path | agreement with double |
|---|---:|
| double FIR + `fxp_rfft512` | 98.27% (1.7% of frames lost, even with an exact input) |
| int16 DFT input, then a double DFT | 99.98–100% |

The direct DFT needs only bins 15..128: 114 bins × 62 non-zero samples. With Σ\|x\| < 2^17 and Q14 cosines, every sum is below 2^31, so it is exact in int32. This is the same precondition as `fxp_rdft.c`.

## Accuracy (host)

`tests/test_nlp.c` is in `scripts/test-fxp.sh` (FXP_CHECK + UBSan + ASan). It holds a double-precision reference NLP (exact DFT) and the E1 float NLP (float arithmetic, kiss FFT). The float NLP agrees with double on 100% of all frames below.

| test | result |
|---|---|
| ⌊3g/10⌋, ⌊3g/20⌋ helpers, 2M values | exact |
| 4b/5, 6b/5 sub-multiple bounds | equal to the float casts for every b |
| decimating FIR vs 64-bit reference, 6,000 blocks | 96,000 of 96,000 outputs equal |
| synthetic voiced/unvoiced signal, peak 30000 / 3000 / 300 / 30 / 3 | 99.67% / 99.83% / 99.83% / 99.83% / 99.83% bin agreement with double |
| steps between peak 30000 (peaky) and 3, every 1930 samples | 99.67% |

**Corpus, NLP bin agreement with double (all frames; |int − dbl| > 2 bins in brackets):**

| condition | tune | eval |
|---|---:|---:|
| loud | 99.81% (0.02%) | 99.75% (0.03%) |
| nominal | 99.82% (0.02%) | 99.78% (0.02%) |
| quiet | 99.81% (0.01%) | 99.79% (0.01%) |
| noisy15 | 99.78% (0.02%) | 99.74% (0.03%) |
| transitions | 99.78% (0.00%) | 99.89% (0.00%) |

**Development history (tune set, agreement with double):**

| version | loud | nominal | quiet | noisy15 |
|---|---:|---:|---:|---:|
| `fxp_rfft512`, integer notch and FIR | 98.44% | 98.27% | 95.54% | 98.55% |
| + direct DFT, block-float FIR | 99.77% | 99.70% | 97.79% | 99.78% |
| + per-frame notch scale (history not in the range test) | 99.38% | 99.36% | 99.36% | 99.64% |
| **+ FIR history in the range test (this stage)** | **99.81%** | **99.82%** | **99.81%** | **99.78%** |

With this stage, the integer notch and FIR followed by a double DFT agree on 99.98–100%. The remaining ~0.2% comes from the int16 DFT input and the Q14 cosines, at near-ties between bins.

## Bitstream agreement (host)

**E2 vs E1:**

| set | frames differ | v1 | v2 | Wo | Wo ±1 | Wo > 1 | E | LSP |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| tune (101,356 frames) | 0.12% | 0.01% | 0.01% | 0.11% | 0.02% | 0.09% | 0 | 0 |
| eval (113,702 frames) | 0.14% | 0.01% | 0.02% | 0.13% | 0.03% | 0.10% | 0 | 0 |

**Reference suite, FvsS (this encoder vs the float Codec2-mod split encoder):**
- Frames identical: tune 0.9961 (E1 0.9973), eval 0.9956 (E1 0.9970).
- FvsU: 0.9839 / 0.9800, against SvsU 0.9875 / 0.9840.
- Silence and transitions match E1 exactly.
- On eval, one frame in each of the −40 and −60 dBFS white-noise items (150 frames each) now differs in Wo; one of them also differs in v2.

## Quality (host, `results/fxp-nlp-suite/{tune,eval}`)

**Eval set (held out), all speech conditions:**

| metric | SS | FF (doc 09) | **FF (this stage)** | target | met |
|---|---:|---:|---:|---|---|
| PESQ-NB | 2.676 | 2.674 | **2.674** | FF within 0.05 of SS | yes (−0.002) |
| LSD (dB) | 9.011 | 9.004 | **9.006** | ≤ SS + 0.2 | yes (−0.005) |
| level error (dB) | 0.646 | 0.647 | **0.647** | within ±0.5 of SS | yes (+0.001) |

**Per condition, FF − SS PESQ:**

| condition | eval | tune |
|---|---:|---:|
| loud | −0.010 | −0.004 |
| noisy | +0.002 | +0.010 |
| nominal | +0.004 | +0.001 |
| quiet | −0.004 | +0.003 |
| transitions | +0.017 | 0.000 |

- **Tune set overall:** PESQ 2.708 (SS 2.705; doc 09 FF 2.705), LSD 9.085, level error 0.528.
- **Eval loud** moved from −0.005 (doc 09) to −0.010. That is within the target and the same size as the per-condition spread of the float encoders (SU vs SS loud: −0.006). It is noted for the listening check.
- **Listening:** not done yet. Listening samples are still required before any quality claim is final.

## Cost (emulator, `results/fxp-nlp-O2`, same 240-frame corpus)

| | doc 09 | this stage | change |
|---|---:|---:|---:|
| encode mean / frame | 5,945,608 | **2,401,292** | −59.6% |
| encode p95 / max per frame | 5,970,811 / 5,982,241 | 2,428,888 / 2,443,764 | |
| `nlp`, inclusive (2 calls / frame) | 3.79M | **242,557** | −93.6% |
| `speech_to_uq_lsps` (float, E3) | 1.97M | 1.97M | now 82% of encode |
| encoder init | 5,068,720 | 1,481,464 | the kiss config is gone. The rest is float window generation in `analysis_init`. |
| stack high water through encode | 1,212 B | **588 B** | |
| `.text` (whole image) | 27,632 B | 23,764 B | kiss no longer linked; `nlp.o` 1,594 B, tables 224 B `.rodata` |
| `.bss` (codec state) | 19,596 B | **8,028 B** | kiss config and memory, and the float NLP buffers, are gone |
| encoder state, host `sizeof` | 17,352 B | 5,784 B | includes a 3,076 B scratch shared by NLP, analysis and the LPC window |
| decode mean / frame | 479,616 | 479,616 | decoder unchanged |

- **Emulator packets.** `packets.bit` is identical to doc 09 (md5 dd1c26db…) on the 240-frame emulator corpus.
- **Calls from `nlp.o`:** only `memcpy`, `memmove` and `memset`, plus the three const tables. There is no softfloat and no `__muldi3`, and nothing calls `__divsi3` (the divides use the hardware divider).
- **`nlp` stack frame:** 60 B.
- **Where the NLP's 121k instructions per call go:** mostly the direct DFT (114 bins × 62 samples × 2 multiply-accumulates). It could be reduced later, by a radix-2 split of the 64-point input or by fewer bins with an adaptive search. That is not needed for correctness and is not done.

## Scratch memory

The encoder struct now holds one `int32_t scratch[769]` (3,076 B), used in turn by three parts that never overlap in time:
- `nlp()`: x[64] int16, y[64] and re/im[2×114]; 260 words.
- `dft_speech`/`power_spectrum`: Sw and pw.
- The float LPC window, which needs `M_PITCH` floats until E3.

`analyse_one_frame` calls `nlp()` first, before the analysis FFT overwrites the scratch.

## What remains

> **Update (doc 12):** the NLP DFT (paired bins k and 128 − k) and the FIR (paired symmetric taps) are rearranged, bit-identical. Encode is 0.39M instructions per frame on the emulator, and `nlp` 134k (was 243k).

> **Update (doc 11):** stage E3 (LPC/LSP analysis and the E/LSP quantisers) is integer. The whole encoder is now integer, and no softfloat or libm is linked. Encode is 0.50M instructions per frame on the emulator (was 2.40M).

- **E3:** integer LPC/LSP analysis and the E/LSP quantisers. `speech_to_uq_lsps` is now 82% of encode (1.97M instructions per frame).
- **After E3:** remove the float window generation in `analysis_init`, so encoder init is free of libm.
- **Required before firmware integration:** listening samples, decode stack measured in isolation, and the independent review loop.

## Reproduce

```bash
python3 tools/gen/gen-nlp-tables.py
```

```bash
scripts/test-fxp.sh
```

```bash
scripts/host-build.sh
```

```bash
.venv/bin/python scripts/ref-suite.py eval results/fxp-nlp-suite/eval
```

```bash
.venv/bin/python scripts/suite-summary.py results/fxp-nlp-suite/eval
```

Per-corpus NLP agreement (pass raw files to the unit test):

```bash
cc -O2 -Icodec/inc tests/test_nlp.c codec/src/fxp_nlp_tables.c codec/src/fxp_fft.c codec/src/fxp_fft_tables.c codec/src/fxp_math.c codec/src/fxp_math_tables.c codec/src/kiss_fft.c codec/src/kiss_fftr.c -lm -o build/test/tnlp && build/test/tnlp corpus/eval/*_loud.raw
```

Emulator, from the repo root:

```bash
docker run --rm -v "$PWD:/w" -w /w p32-csky-gcc14:14.2.0-arm64 sh scripts/emu-build.sh codec results/fxp-nlp-O2
```

```bash
scripts/emu-run.sh results/fxp-nlp-O2 profile
```

```bash
.venv/bin/python scripts/emu-analyze.py results/fxp-nlp-O2
```

```bash
.venv/bin/python scripts/profile-report.py results/fxp-nlp-O2
```
