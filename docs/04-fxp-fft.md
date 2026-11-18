# 04 — Fixed-point FFT stage (decoder side, hybrid)

**Status:** the decoder's three FFT uses run in fixed point. Everything around them is still float, and a temporary bridge converts at the boundaries (`codec/src/fxp_bridge.c`). The encoder is unchanged.

**Measurement types used here:**
- **host:** macOS arm64, clang. Used for accuracy and quality.
- **emulator:** XUANTIE QEMU instruction counts, `-icount shift=0`. These are instructions, **not cycles**.

Nothing here is DM32 hardware timing.

## Kernels

| kernel | file | used for | method |
|---|---|---|---|
| `fxp_irfft512` | `fxp_fft.c` | synthesis inverse DFT (`kiss_fftri`), 2 per frame | 256-pt complex radix-4 FFT plus real split; block floating point |
| `fxp_rdft512_sparse` | `fxp_rdft.c` | `aks_to_mag2` spectrum of A(z) and A(z/γ), 4 per frame | direct DFT of the 11 taps; exact int32 accumulation |
| `fxp_rfft512` | `fxp_fft.c` | not used yet (the encoder `dft_speech` and NLP are its future users) | forward version of the FFT above |

### Number formats

Conventions are in `codec/inc/fxp.h`: arithmetic right shifts, round half up, explicit saturation, no 64-bit math.

**FFT (`fxp_rfft512` / `fxp_irfft512`)**
- **Data:** int16 values with a block exponent `e`. True value = stored × 2^e.
- **Twiddles:** Q15 `fxp_w512[384]`, clamped to ±32767.
- **Scaling:**
  - The input is normalised so its peak lands in [2^13, 2^14).
  - Each radix-4 stage applies a data-dependent shift of 0–3 bits, chosen from an OR-accumulated magnitude bound. Stage outputs stay below 21477.
  - The overflow proofs are in the `fxp_fft.c` header. `FXP_ASSERT` checks them in test builds.
- **Normalisation:** unnormalised, like kiss: `irfft(rfft(x)) = 512·x`.

**Sparse DFT (`fxp_rdft512_sparse`)**
- **Inputs:** int16 taps and Q14 twiddles `fxp_cos512_q14[512]`. Sin is read at index `(m+384)&511`.
- **Output:** int32. True X = stored × 2^-14 (in tap units).
- **Overflow:** every partial sum is bounded by Σ|a|·2^14. The precondition Σ|a| < 2^17 keeps that below 2^31, so no saturation is needed.
- **Tap scaling:** `fxp_rdft_taps_shift(max, L1)` picks the largest shift with max < 2^15 and L1 < 2^17.
- **Symmetry:** the even/odd split gives X[k] and X[256−k] from one pass, so k runs over 0..128 only.

### Why the LPC spectrum doesn't use the FFT

The first integration put `fxp_rfft512` in `aks_to_mag2`. That gave correct metrics on average, but it had large errors at formant peaks:
- **Cause:** |A(e^jω)| spans 40–60 dB, and its minima are where the harmonic amplitudes and phases are read. The FFT rounds at every stage relative to the spectral *peak*.
- **Measurement:** on the host, per-harmonic `A[m]` and `Aw` phase were compared against an exact double DFT of the same `ak`, over the speech items of each set.

| LPC spectrum method | set | A[m] error rms | A[m] error max | Aw phase error rms | phase error max |
|---|---|---:|---:|---:|---:|
| float kiss (S, reference build) | 2 files | 0.0000 dB | 0.001 dB | 0.00000 rad | 0.0001 rad |
| 16-bit FFT `fxp_rfft512` | tune | 0.200 dB | 17.4 dB | 0.0113 rad | 3.07 rad |
| 16-bit FFT `fxp_rfft512` | eval | 0.275 dB | 20.0 dB | 0.0137 rad | 3.09 rad |
| **sparse DFT `fxp_rdft512_sparse`** | tune | **0.012 dB** | 1.10 dB | **0.00085 rad** | 0.25 rad |
| **sparse DFT `fxp_rdft512_sparse`** | eval | **0.017 dB** | 6.90 dB | **0.0015 rad** | 2.40 rad |

- **Speed:** on the emulator the sparse DFT is also 2.5× cheaper than the FFT path: 21.3k vs 53.3k instructions per call.
- **Unit test, peakiest synthetic LPC** (poles at r ≈ 0.99, 51 dB range): worst bin |A|² error 0.024 dB, against 0.62 dB for the FFT.
- **Exactness:** bit-exact against an int64 reference at the overflow limits (±32767 taps with L1 = 131068, random taps at the L1 limit, 12 and 13 taps).

**Known limitation: residual sparse-DFT outliers.**
- **Scale:** 25 of 6.8M eval harmonics exceed 1 dB.
- **Where they occur:** near-unstable interpolated LPCs, with taps up to ±30, L1 = 161, and an |A| minimum around 0.005, about −90 dB relative to L1.
- **Why the kernel can't fix them:** resolving those minima needs roughly 20-bit taps *and* twiddles. No 16×16 kernel can deliver that, and a 32×32 multiply has no cheap form on CK803S.
- **Planned fix:** when `aks_to_mag2` itself is ported, evaluate |A| and arg A from the LSPs in product form:

  A = e^{−j11ω/2}(P_r + jQ_r)/2

  with

  P_r = 2cos(ω/2)·Π(2cos ω − 2cos ω_i) and Q_r = 2sin(ω/2)·Π(…).

  This form has no cancellation, so it is well conditioned exactly where the direct sum isn't. The fixed-point decoder has the LSPs available anyway. A(z/γ) with γ = 0.5 is well conditioned and keeps the sparse DFT.

### Waveform SNR is not a precision metric for this decoder

Decoder-vs-decoder waveform SNR between the candidate and S looked poor: 16 dB median with the FFT. Swapping in quantised direct DFTs showed it is *non-monotonic* in precision; 16-bit taps scored worse than 15-bit on one file.

- **Cause:** `postfilter()` draws a random phase only for harmonics with `A[m] < thresh`. A single borderline flip shifts the `codec2_rand` sequence, and every later random phase differs.
- **Same effect in the float codecs:** it explains why the float S and U decoders also agree only to about 16 dB (doc 03).
- **Consequence:** precision is judged at the parameter level (table above) and by the perceptual and spectral metrics, not by waveform SNR.

## Integration (hybrid bridge)

`fxp_bridge.c` converts float to int16 using an exponent from the input peak, plus the L1 cap for the sparse DFT. It runs the kernel and converts back with `ldexpf`.

- **Scratch memory:** taken from the decoder's existing `fft_buffer`. The aliasing rules are documented in `fxp_bridge.h`.
- **In-place conversion:** the sparse DFT writes `cq31_t` into the float `complex_t` output buffer and converts in place; a static assert checks that the two layouts match in size.
- **Temporary:** the bridge (and its softfloat glue) goes away as neighbouring kernels become fixed point.

## Results

### Quality (host, `results/fxp-fft-suite/{tune,eval}`)

The candidate F is the `codec/` tree, built as `build/host/c2fx`. It is in the `CODECS` table of `scripts/ref-suite.py`.

**Eval set (held out), all speech conditions:**

| metric | SS | FF | UF | UU | target | met |
|---|---:|---:|---:|---:|---|---|
| PESQ-NB | 2.676 | 2.676 | 2.673 | 2.677 | FF within 0.05 of SS | yes (+0.0001) |
| STOI | 0.826 | 0.826 | 0.826 | 0.826 | — | — |
| LSD (dB) | 9.011 | 9.001 | 9.012 | 9.022 | ≤ SS + 0.2 | yes (−0.010) |
| level error (dB) | 0.646 | 0.646 | 0.648 | 0.650 | within ±0.5 of SS | yes (−0.0006) |

- **Per condition:** FF−SS PESQ is within ±0.006 on nominal, quiet, loud and noisy.
  - Transitions is a single item, and its PESQ swings by ±0.03 on tune and ±0.05 on eval for UF−UU as well.
  - LSD is equal or better on every condition except eval transitions (+0.004 dB).
- **Non-speech output levels (eval):**
  - silence: 18.0 (FF) vs 18.1 (SS)
  - −60 dBFS noise: 29.0 vs 29.0
  - −40 dBFS noise: 49.2 vs 49.2
- **Encoder:** unchanged. FvsS frames are 100% identical, and the emulator `packets.bit` is identical to `results/codec-wip-O2`.
- **Waveform agreement with the S decoder:** median 31–45 dB on the same bitstream, up from 16 dB with the FFT-based LPC spectrum. The remaining gap is the postfilter random-phase effect described above.
- **Listening:** not done yet. Listening samples are required before any quality claim is final (doc 03, target 3).

### Cost (emulator, `results/fxp-fft-O2`, same 240-frame corpus as doc 02)

| | float split (doc 02) | this stage | note |
|---|---:|---:|---|
| decode mean / frame | 13,057,631 | 6,893,556 | −47% |
| decode max / frame | — | 9,978,500 | |
| `fxp_rdft512_sparse` (4 calls) | — | 85,212 | 21.3k per call |
| `fxp_irfft512` | — | ≈ 56.6k per call (kbench) | 2 calls per frame |
| bridge glue (`ldexpf`, `lrintf`, `frexpf` softfloat) | — | ≈ 1.0M | disappears when the neighbouring kernels are ported |
| stack high water through decode | 5,552 B | 5,384 B | |
| encode mean / frame | 10,925,126 | 10,925,126 | unchanged |

**Kernel microbenchmarks** (`harness/emu/kbench.c`, `results/kbench-fft-O2`, in instructions):

| kernel | instructions |
|---|---:|
| `fxp_rfft512` | 50.9k–53.8k |
| `fxp_irfft512` | 56.6k |
| `fxp_rdft512_sparse` | 21.3k |

- **Kernel stack:** 240 B.
- **Compiler options:** O3 is about 1% faster; Os is slower.
- **Code size:** `fxp_cfft256` is 1056 B, `fxp_rdft512_sparse` 244 B, `fxp_irfft512` 150 B.
- **Tables:**
  - `fxp_w512`: 1536 B
  - `fxp_digrev256`: 256 B
  - `fxp_cos512_q14`: 1024 B

**FFT optimisation status:** the radix-4 general butterfly costs about 172 instructions, of which about 39 are register spills, because CK803S `-mno-high-registers` leaves only r0–r13 usable. A hand-scheduled or assembly butterfly is a later optimisation, once profiles say the FFT matters relative to the rest; today the bridge's softfloat glue costs more than the kernels.

## Reproduce

```bash
scripts/test-fxp.sh
```
Unit tests with `FXP_CHECK`, UBSan and ASan, no recovery.

```bash
scripts/host-build.sh
```
Builds `build/host/c2fx`.

```bash
.venv/bin/python scripts/ref-suite.py eval results/fxp-fft-suite/eval
```

```bash
.venv/bin/python scripts/suite-summary.py results/fxp-fft-suite/eval
```

Run the next two inside the gcc container:

```bash
scripts/emu-build.sh codec results/fxp-fft-O2
```

```bash
scripts/emu-build.sh codec results/kbench-fft-O2 O2 harness/emu/kbench.c
```

```bash
scripts/emu-run.sh results/fxp-fft-O2 profile
```

```bash
.venv/bin/python scripts/emu-analyze.py results/fxp-fft-O2
```

```bash
.venv/bin/python scripts/profile-report.py results/fxp-fft-O2
```

**Sanitizers on the codec:** the `-DFXP_CHECK -fsanitize=undefined,address` build of the `codec/` tree decoded all 492 corpus items (tune and eval) without a report.
