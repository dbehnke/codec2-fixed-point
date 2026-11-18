# 08 — Faster LSP spectrum kernel (`lsp_spec_bin_q29`)

**Status:** the product-form kernel that evaluates |A(e^jw)|² and arg A at each of the 256 bins (doc 05) is about half its previous cost. Decoder output is bit-identical to doc 07 on 489 of the 492 corpus items. The other three differ by 1–2 LSB in a few samples (explained below). Quality metrics are unchanged.

**Measurement types used here:**
- **host:** macOS arm64, clang. Used for accuracy and quality.
- **emulator:** XUANTIE QEMU instruction counts, `-icount shift=0`. These are instructions, **not cycles**.

Nothing here is DM32 hardware timing, and no real-time claim is made.

## What changed (`codec/src/fxp_lpcspec.c`)

The algorithm and the number format are the same as doc 05:

- P = 2cos(w/2)·Π 2(x − c_even) and Q = 2sin(w/2)·Π 2(x − c_odd), each held as a 15-bit pseudo-float (mantissa m, exponent e).
- |A|² = (P² + Q²)/4, then `log2_plus_floor`.
- arg A = atan2(Q, P) − 11w/2, computed at harmonic bins only.

What changed is how the loop runs:

1. **Branch-free mantissa multiply.**
   - `PF_NORM` rounds each factor d to 15 bits as `d16 = (((d << clz(d)) >> 16) + 1) >> 1`, in [2^14, 2^15]. This equals the old `(d + 2^(s−1)) >> s`.
   - `PF_MULN` forms p = m·d16 in [2^28, 2^30] and rounds it back with a shift of 14 + (p >> 29). The result lands in [2^14, 2^15] with no overflow fix-up branch.
   - The mantissa may now equal 2^15 (it was < 2^15). m² ≤ 2^30 still fits `pf_sumsq` and the atan2 alignment.
   - The constant part of each factor's exponent (14 + 17 − 28 = 3) is added once after the loop.
2. **No zero flags.** `clz(d | 1)` keeps the shift defined, and d = 0 gives d16 = 0. A zero mantissa then stays exactly zero through every later multiply. That is the same result as the old `T->m = 0` path, and it also covers sin(w/2) = 0 at k = 0.
3. **Signs as masks.** |d| comes from `(d ^ s) − s` with s = d >> 31, and the signs accumulate by XOR. They are applied to (P, Q) only when the phase is wanted.
4. **P and Q in one loop iteration,** unrolled 5×. There is no `T` pointer and no per-factor `m == 0` test.
5. **Direct table reads.** cos w, cos(w/2) and sin(w/2) are read from the quarter-wave table `fxp_cosq_q30` instead of 3 `fxp_cos2048_q30` calls per bin. For 0 ≤ k < 256 the values are identical.
6. **Q29 LSP cosines computed once per subframe** in `fxp_aks_to_mag2` (`c29[]`), not once per factor per bin.
   - The public `fxp_lsp_spec_bin(c_q30, k, ang)` used by the unit tests is now a wrapper: it builds `c29` and calls the static kernel.

### Why three items are not bit-identical

`pf_sumsq` aligns P² and Q² to the larger exponent and truncates the smaller term. A mantissa of exactly 2^15 at exponent e used to be stored as 2^14 at e + 1. The new form keeps the smaller term 2 bits more precisely when that exponent sets the alignment.

Test: renormalising m = 2^15 → 2^14 before `pf_sumsq` (scratch build) made all three items bit-identical to doc 07. The new code is kept, since it is the more precise form and needs no extra step.

## Accuracy (host)

`scripts/test-fxp.sh` (FXP_CHECK + UBSan + ASan), against exact double:

| test | doc 07 | this stage |
|---|---|---|
| \|A\|² product form, 5.12M bins | rms 0.00032 dB, max 0.00170 dB | same |
| arg A (all bins / harmonic bins) | max 1.44e-4 rad | same |
| A[m] codebook frames, 667,097 | rms 0.00043 dB, max 0.01959 dB | same |
| A[m] synthetic 5×25 Hz pairs | rms 0.00594 dB, max 0.41276 dB | same |

## Agreement and quality (host)

- **Same bitstreams (upstream encoder), new vs doc 07 decoder, all 492 items:**
  - 489 items are bit-identical.
  - 3 differ: 1221-135766-0000_loud (1 sample, ±1), 7021-79730-0001_noisy15 (28 samples, max ±2) and 7127-75946-0000_noisy15 (3 samples, ±1).
- **Suites (`results/fxp-spec-suite/{tune,eval}`):** `summary.md` is identical to `results/fxp-params-suite` at every printed digit.
  - Eval FF: PESQ 2.673 (SS 2.676), LSD 9.002 dB, level error 0.646 dB.
  - Tune FF: PESQ 2.708.
- **`scripts/synth-ref.sh`:** `results/fxp-synth-ref/summary.md` is unchanged. Eval median 58.74 dB; tune median 58.73 dB, flips 8/3.
- **Sanitizers:** the `-DFXP_CHECK -fsanitize=undefined,address -fno-sanitize-recover=all` build encoded and decoded all 492 corpus items without a report.
- **Encoder:** untouched. Emulator `packets.bit` md5 51aea44b… is identical.

## Cost (emulator, `results/fxp-spec-O2`, same 240-frame corpus)

| | doc 07 | this stage | change |
|---|---:|---:|---:|
| decode mean / frame | 643,448 | **479,584** | −25.5% |
| decode p95 / max per frame | 676,116 / 681,336 | 514,107 / 519,017 | |
| spectrum kernel, inclusive (512 calls / frame) | 336,550 | 172,602 | −48.7% |
| spectrum kernel, self | 290,941 (568 / bin) | 145,159 (284 / bin) | −50% |
| `fxp_cos2048_q30` calls / frame | 1,572 | 36 | |
| kernel stack frame (`.su`) | 68 B | 96 B | +28 B |
| `fxp_aks_to_mag2` stack frame | 228 B | 268 B | +40 B (`c29[10]`) |
| `.text` (whole image) | 27,336 B | 27,540 B | +204 B (unrolled loop) |
| stack high water through decode (cumulative) | 1,476 B | 1,476 B | still the encoder's |

Each stage was measured on the emulator:

| step | decode mean / frame |
|---|---:|
| doc 07 | 643,448 |
| branch-light `pf_mul` | 542,338 |
| without the zero flags, which removed register spills | 502,275 |
| folded rounding shift and 5× unroll | 496,734 |
| quarter-table reads in place of `fxp_cos2048_q30` | 479,584 |

**Where the decode now goes (inclusive, per frame):**

| part | instructions | share |
|---|---:|---:|
| spectrum kernel | 172.6k | 36% |
| irfft (`fxp_cfft256` + split + normalise) | 109k | 23% |
| `log2_plus_floor` + log2 | 47k | 10% |
| A(z/γ) sparse DFT | 42.6k | 9% |
| exp2 (all callers: band, energy and frame-energy sums) | 41.7k | 9% |
| `fxp_aks_to_mag2` own loops | 55k | 11% |
| synthesis loops | 16k | 3% |

## What remains

> **Update (doc 09):** encoder stage E1 (analysis front end) is integer. Encode is 5.95M instructions per frame on the emulator (was 10.93M); NLP and LPC/LSP are still float.

- **Further decode savings** (not done; each is a separate, measurable step):
  - A log2 that skips the `+1e-6` floor alignment when the value is far above it.
  - Computing band sums from the per-bin exp2 already needed for e_after.
  - A radix-4 `fxp_cfft256`.
  - None of these are needed for correctness.
- **The encoder** (float, about 10.9M instructions per frame on the emulator) dominates the total and is the next large stage.
- **Required before firmware integration:** listening samples, decode stack measured in isolation, and the independent review loop.

## Reproduce

```bash
scripts/test-fxp.sh
```

```bash
scripts/host-build.sh
```

```bash
.venv/bin/python scripts/ref-suite.py eval results/fxp-spec-suite/eval
```

```bash
.venv/bin/python scripts/suite-summary.py results/fxp-spec-suite/eval
```

```bash
scripts/synth-ref.sh
```

Inside the gcc container (repo mounted at /w):

```bash
scripts/emu-build.sh codec results/fxp-spec-O2
```

Then on the host:

```bash
scripts/emu-run.sh results/fxp-spec-O2 profile
```

```bash
.venv/bin/python scripts/emu-analyze.py results/fxp-spec-O2
```

```bash
.venv/bin/python scripts/profile-report.py results/fxp-spec-O2
```
