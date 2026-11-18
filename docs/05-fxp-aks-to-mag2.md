# 05 — Fixed-point spectral amplitudes (`aks_to_mag2`)

**Status:** the decoder's LPC-to-harmonic-amplitude stage runs in fixed point, working directly from the decoded LSPs. It replaces float `lsp_to_lpc` + `aks_to_mag2` (with its sparse DFT from doc 04) in `codec2_decode`. Inputs (float LSPs, E, Wo) and outputs (float `A[m]`, phase-only `Aw`) still cross a thin float bridge, `aks_to_mag2_lsp()` in `codec/src/lpc.c`. The encoder is unchanged.

**Update (doc 06):** the outputs are now fixed point end to end. `fxp_aks_to_mag2` takes the integer Wo angle `W` and returns `log2 A[m]` (Q16) and `arg A` at each harmonic (uint32 angle) straight into the fixed-point synthesis; the float `A[m]`/`Aw` outputs, the float LPC correction and `apply_lpc_correction` are gone. Harmonic bins and band edges now use the tie rule `fxp_bin_round` (doc 06). The figures below are for the stage as committed in doc 05.

**Measurement types used here:**
- **host:** macOS arm64, clang. Used for accuracy and quality.
- **emulator:** XUANTIE QEMU instruction counts, `-icount shift=0`. These are instructions, **not cycles**.

Nothing here is DM32 hardware timing.

## What the stage computes

Float reference (`aks_to_mag2`), per bin k = 0..255 of the 512-point grid (ω_k = 2πk/512):

- A2 = |A(e^jω)|² + 1e-6 and A2g = |A(e^jω/γ)|² + 1e-6, with γ = 0.5
- Pw = (1/A2)·(sqrt(A2g/A2) + 1e-5)^0.4
- e_before = Σ 1/A2 and e_after = Σ Pw, then Pw × 1.96 for k < 64
- gain = E·e_before/e_after
- A[m] = sqrt(gain·Σ_{am ≤ k < bm} Pw), band edges am, bm from Wo
- `Aw` = A(e^jω_k), of which synthesis uses **only the phase** (`sample_phase` → `phase_synth_zero_order`, which normalises)

## Design

### |A| and arg A from the LSPs in product form

The direct DFT of the 11 LPC taps loses precision at deep minima of |A| (doc 04, known limitation), which is exactly where harmonic amplitudes and phases are read. The LSPs give a form with no cancellation:

A(e^jω) = e^{−j11ω/2}·(P_r + jQ_r)/2

P_r = 2cos(ω/2)·Π_{i even} 2(cos ω − cos ω_i), Q_r = 2sin(ω/2)·Π_{i odd} 2(cos ω − cos ω_i)

|A|² = (P_r² + Q_r²)/4. The even/odd assignment matches float `lsp_to_lpc` (checked against a double transcription of it: 1.4e-14). Each factor has full relative precision, so the precision of |A| is the same at a −90 dB minimum as at a peak.

**Arithmetic (`fxp_lsp_spec_bin`):**
- cos ω_k, cos(ω_k/2), sin(ω_k/2) are table reads (`fxp_cos2048_q30`, Q30, quarter wave, 513 entries).
- P_r and Q_r are **pseudo-float**: a 15-bit unsigned mantissa in [2^14, 2^15), an integer exponent, and a separate sign. Each factor (cos ω − c_i in Q29, |·| ≤ 2^30) is normalised to 15 bits with rounding, multiplied as 15×15 → 30 bits, and renormalised. Six multiplies per polynomial; no overflow possible.
- |A|² aligns the two squared mantissas, adds the 1e-6 floor (1125899907·2^−50) and returns log2 in Q16.
- Phase: (P_r + jQ_r) on a common exponent, rotated by e^{−j11ω/2} (a Q15 cos/sin pair at index 22k of the same table). The output is an unnormalised `cq31_t`, |·| in about [2^29, 2^31).

### LSP cosines

The decoded LSPs are integer Hz (`decode_lspds_scalar` sums integer codebook deltas); interpolated frames average two, so all LSPs are integers in **half-Hz** (hz2). c_i = cos(π·hz2/8000):
- index j = ⌊16·hz2/125⌋ into the 2048-per-period cosine table, remainder δ = (16·hz2 mod 125)/125 · π/1024 < 0.0031 rad;
- second-order Taylor: cos(a+δ) = cos a − δ sin a − δ²/2 cos a, full-precision Q30 products (`fxp_mul_q30`).
- Works for any hz2 ≥ 0 (cumulative LSPs can exceed 4000 Hz, angle > π). Error 6.0e-9 over hz2 = 0..40000.

The bridge recovers hz2 from the float LSPs with `lrintf(lsp·8000/π)`; that is exact for every decoded frame except the decoder's very first, where `prev_lsps_dec` starts at iπ/11 (not integer Hz), so the first interpolated frame is rounded to the nearest half-Hz. Carrying integer LSPs through the decoder removes the bridge.

### A(z/γ): fixed-point `lsp_to_lpc` + sparse DFT

A(z/γ) has no LSP product form, but with γ = 0.5 its zeros sit at half the radius of A's zeros, so it is normally well conditioned and keeps the 16-bit sparse DFT of doc 04.
- `fxp_lsp_to_lpc_q19` convolves the second-order sections in **Q22**. With |c| ≤ 1 every partial product is dominated coefficient-wise by (1 + z^−1)^(2s), so coefficients are ≤ C(8,4) = 70 before the last section (keeps `fxp_mul_q30`'s |a| < 2^29), ≤ 252 after it, and P[j] + P[j−1] ≤ C(11,5) = 462 (< 2^31 in Q22). Output is rounded to Q19. Max error 1.6e-5 (8 LSB of Q19).
- ag_i = ak_i·2^−i (rounded shift), scaled to int16 taps with the L1 rule of doc 04 (one retry if rounding breaks it), sparse DFT, then |·|² from the int32 output normalised to 15-bit components, plus the floor, as log2 Q16.

### Weighting and band sums in the log2 domain

Everything after the spectra is done in log2 Q16, where the powers become multiply/adds:
- log2 Pw = −log2 A2 − (log2 A2 − log2 A2g)/5 (hardware `divs` on CK803S).
- The `+1e-5` inside the power is **dropped**. It matters only when sqrt(A2g/A2) is near 1e-5; for any realisable A the ratio is ≥ about 0.03, where it changes Pw by < 0.03%.
- e_before and e_after are sums of 2^(term − max) in Q22 (256 terms < 2^30); each band sum is relative to the band maximum in Q24 (band ≤ 64 bins; real bands are 3–26), with the 1.96 boost added as log2(1.96) = 63625 in Q16. Relative accumulation keeps weak bands as precise as strong ones.
- log2 A[m] = (log2 E + log2 e_before − log2 e_after + log2 Em)/2.
- `fxp_log2_q16`: clz + 128-entry table + linear interpolation, error 2.4e-5 (1.5e-4 dB). `fxp_exp2_neg_q30` / `fxp_exp2_split`: 128-entry table + interpolation, relative error < 1e-5.

**Hybrid parts (temporary):** log2 E comes from float E via `frexpf`; band edges use the float Wo expressions unchanged; A[m] is returned as float through `exp2_split` + `ldexpf`; phase is converted to float in place in `Aw` (`_Static_assert` on the layout).

### Memory

- Scratch: 514 words (2056 B) in the decoder's existing `fft_buffer`. The DFT of A(z/γ) is written there as `cq31_t[257]`; each bin is converted to log2|A(z/γ)|² in place at word k (safe: word k ≤ 2k, and words below k hold only consumed bins); log2 A2 then goes to words 256..511 (consumed bins), and log2 Pw overwrites the A(z/γ) values.
- The decoder no longer needs `ak[2][11]` or `A2[256]` on its stack (−1.1 KB).
- Tables (new, `codec/src/fxp_math_tables.c`, generated by `tools/gen/gen-math-tables.py`): `fxp_cosq_q30` 2052 B, `fxp_log2_tab_q16` 516 B, `fxp_exp2_tab_q30` 516 B. The cosine table could be halved later (Q15 sin/cos plus interpolation) if flash is tight.

## Accuracy

### Unit tests (host, `tests/test_fxp_lpcspec.c`, run by `scripts/test-fxp.sh` with FXP_CHECK + UBSan + ASan)

Frames are random draws from the real LSP codebooks (`delta_lsp_cb`), half of them interpolated; a quarter are deliberately pathological (all five LSP pairs at the 25 Hz codebook minimum). References are exact double.

| quantity | result |
|---|---|
| `fxp_log2_q16` | max abs error 2.4e-5 |
| `fxp_exp2_neg_q30` / `fxp_exp2_split` | max relative error 8.3e-6 / 3.7e-6 |
| LSP cosine | max abs error 6.0e-9 |
| `fxp_mul_q30` | max 0.50 LSB |
| `fxp_lsp_to_lpc_q19` | max abs error 1.6e-5 |
| \|A\|² (product form), 5.1M bins | rms 0.0003 dB, max 0.0018 dB |
| arg A | max 1.6e-4 rad |
| A[m], codebook frames (669k harmonics) | rms 0.0004 dB, **max 0.025 dB** |
| A[m], pathological 5×25 Hz frames (223k) | rms 0.006 dB, max 0.50 dB |

**Known limitation:** in the pathological frames all ten zeros of A cluster at low frequency, and γ = 0.5 does not move the zeros of A(z/γ) far enough from the unit circle: |A(z/γ)|² falls to about 1.5e-5 against a tap L1 in the tens, beyond what 16-bit taps and Q14 twiddles resolve (log2 A2g off by up to 0.7). The 0.2 power turns that into at most 0.5 dB on A[m]. These frames do not occur in the corpus (below). A higher-precision A(z/γ) evaluation would fix it if it ever matters.

### Real corpus, parameter level (host)

Each decoded subframe of the S-encoded speech items (nominal, quiet, loud, noisy) compared against an exact double `lsp_to_lpc` + DFT + `aks_to_mag2` from the same float LSPs; phase compared at the harmonic bins synthesis reads.

| method | set | harmonics | A[m] err rms | A[m] err max | phase err rms | phase err max |
|---|---|---:|---:|---:|---:|---:|
| sparse DFT (doc 04) | eval | 6.8M | 0.017 dB | 6.90 dB | 0.0015 rad | 2.40 rad |
| **LSP product form (this stage)** | tune | 8.4M | **0.0002 dB** | **0.012 dB** | **0.00003 rad** | 0.0013 rad |
| **LSP product form (this stage)** | eval | 9.3M | **0.0002 dB** | **0.012 dB** | **0.00003 rad** | 0.0011 rad |

No harmonic exceeds 0.1 dB. (The doc 04 row used nominal, quiet and loud only; this stage adds noisy.)

### Quality (host, `results/fxp-aks-suite/{tune,eval}`)

Candidate F = the `codec/` tree (`build/host/c2fx`), as in doc 04.

**Eval set (held out), all speech conditions:**

| metric | SS | FF (doc 04) | **FF (this stage)** | UF | UU | target | met |
|---|---:|---:|---:|---:|---:|---|---|
| PESQ-NB | 2.676 | 2.676 | **2.678** | 2.675 | 2.677 | FF within 0.05 of SS | yes (+0.002) |
| LSD (dB) | 9.011 | 9.001 | **9.005** | 9.012 | 9.022 | ≤ SS + 0.2 | yes (−0.006) |
| level error (dB) | 0.646 | 0.646 | **0.647** | 0.648 | 0.650 | within ±0.5 of SS | yes (+0.001) |

- **Per condition, FF − SS:** PESQ eval −0.001 … +0.008 (tune 0.000 … +0.002); LSD lower (better) on every condition of both sets (−0.001 … −0.038 dB); level error identical to 0.001 dB.
- **Non-speech output levels:** silence 18.0 (FF) vs 18.1 (SS); −60 dBFS noise 29.0 vs 29.0; −40 dBFS noise 49.2 vs 49.2 (eval; tune identical to 0.1 dB).
- **Encoder:** unchanged, FvsS frames 100% identical on every condition.
- **Waveform agreement with the S decoder** on the same bitstream: median **58.3 dB** (doc 04: 31–45 dB), min 11.8–14.5 dB. The minimum is the postfilter random-phase effect (doc 04): one borderline `A[m] < thresh` decision shifts the random sequence. With A[m] now within 0.012 dB of exact, those flips are rare.
- **Listening:** not done yet. Listening samples remain required before a quality claim is final (doc 03, target 3).

## Cost (emulator, `results/fxp-aks-O2`, same 240-frame corpus as doc 02)

| | doc 04 stage | this stage | note |
|---|---:|---:|---|
| decode mean / frame | 6,893,556 | 3,063,526 | −56% |
| decode p95 / max per frame | — / 9,978,500 | 4,973,569 / 6,100,183 | |
| `fxp_aks_to_mag2` inclusive | — | 262.8k per call | 2 calls per frame |
| of which `fxp_lsp_spec_bin` | — | 711 per bin (182k per call) | 256 bins |
| bridge `aks_to_mag2_lsp` glue | — | ≈ 107k per call | `lrintf`, `ldexpf`, int→float phase; goes when neighbours are ported |
| stack high water through decode | 5,384 B | 4,268 B | `ak`/`A2` gone from `codec2_decode` |
| encode mean / frame | 10,925,126 | 10,925,126 | unchanged; `packets.bit` md5 51aea44b… identical |

- **Code size (O2):** `fxp_lsp_spec_bin` 786 B, `fxp_aks_to_mag2` 766 B, `fxp_lsp_to_lpc_q19` 400 B, `log2_plus_floor` 102 B, `fxp_lsp_cos_q30` 78 B, log2/exp2/cos/mul helpers 64–72 B each, bridge `aks_to_mag2_lsp` 384 B. The float `aks_to_mag2` and `lsp_to_lpc` are no longer linked (`--gc-sections`).
- **Stack:** `fxp_aks_to_mag2` 220 B, `fxp_lsp_to_lpc_q19` 172 B, `fxp_lsp_spec_bin` 68 B, bridge 756 B (its `am`/`bm`/`Al` arrays; they go when Wo and A[m] are fixed point).
- **Remaining decode cost** is dominated by the softfloat in synthesis (`cosf`/`sinf` in phase synthesis, the irfft bridge) — the next candidates.

## Reproduce

```bash
scripts/test-fxp.sh
```

```bash
scripts/host-build.sh
```

```bash
.venv/bin/python scripts/ref-suite.py eval results/fxp-aks-suite/eval
```

```bash
.venv/bin/python scripts/suite-summary.py results/fxp-aks-suite/eval
```

Inside the gcc container (repo mounted at /w):

```bash
scripts/emu-build.sh codec results/fxp-aks-O2
```

Then on the host:

```bash
scripts/emu-run.sh results/fxp-aks-O2 profile
```

```bash
.venv/bin/python scripts/emu-analyze.py results/fxp-aks-O2
```

```bash
.venv/bin/python scripts/profile-report.py results/fxp-aks-O2
```

**Sanitizers on the codec:** the `-DFXP_CHECK -fsanitize=undefined,address` build of the `codec/` tree encoded and decoded all 492 corpus items (tune and eval) without a report.

**Parameter-level harness:** not committed (a scratch copy of `codec/` whose decoder runs the exact double reference next to `aks_to_mag2_lsp` and accumulates the errors above). `-DFXP_LPCSPEC_DEBUG` exports the kernel's intermediate log2 spectra and sums for this kind of diagnosis.
