# 14. Encoder LPC and decoder LPC-spectrum speed-ups (bit-identical)

## Status

Two hot spots are rearranged for speed. The arithmetic is unchanged, so every output is bit-identical to the previous commit:

- `speech_to_lsps_q16` (`codec/src/lpc.c`, encoder)
- `fxp_aks_to_mag2` (`codec/src/fxp_lpcspec.c`, decoder)

The evidence:

- Encode and decode of all 492 corpus items (tune and eval) are byte-identical to a host build of the previous commit.
- The emulator `packets.bit` is identical (md5 `dadc57004a7ed241fac1c2a649a4265e`), and so is the decoder PCM checksum.

On the emulator, encode drops from 366,118 to 355,728 instructions per frame and decode from 449,740 to 438,135. Those are instruction counts, not cycles. There is still no real-time claim; that needs hardware timing.

## What changed

### 1. Autocorrelation: products summed in pairs in int32 (encoder)

The windowed frame x holds int16 values, |x| ≤ 32767 (the window output saturates at ±32767, never −32768). The autocorrelation R[j] = Σ x[i] x[i + j] was accumulated one product at a time into an int64. Each 64-bit add costs a carry chain on the CK803S.

Each product is at most 32767² < 2^30, so the sum of two products is at most 2 · 32767² = 2,147,352,578 < 2^31 − 1. Adding two products in int32 therefore cannot overflow, and the pair goes into the int64 with one 64-bit add. Integer addition is exact, so the total is the same integer. A `_Static_assert` in `lpc.c` checks the bound. An odd-length tail adds its last product alone.

The inner loop is now two `mulsh`, one 32-bit add and one 64-bit add per two samples.

### 2. Inlined log2 and exp2 kernels (decoder)

`fxp_log2_q16` and `fxp_exp2_neg_q30` were called once per spectrum bin. Their bodies now live in `fxp_math.h` as `fxp_log2_q16_inl` and `fxp_exp2_neg_q30_inl`, and the out-of-line functions return those. The code is identical, so the results are too. `fxp_lpcspec.c` uses the inline forms in its per-bin sums and band sums.

### 3. One call site for the per-bin LSP spectrum kernel (decoder)

`lsp_spec_bin_q29` had two call sites in the bin loop: one that also returned a harmonic's amplitude, and one that did not. With two call sites GCC kept it out of line. Now there is one call that passes an output pointer only when a harmonic falls in this bin. It and `log2_plus_floor` are `always_inline`.

The harmonic bookkeeping is unchanged:

- the same bins are hit;
- the same amplitude is stored for every harmonic that rounds to that bin.

## Cost (emulator)

QEMU instruction counts per frame, `-O2`, comparing `results/fxp-hot3-O2` with `results/fxp-fft2-O2`. These are not cycles.

| | before | after |
|---|---:|---:|
| encode mean | 366,118 | 355,728 |
| encode p95 | 388,464 | 378,074 |
| encode max | 390,391 | 380,001 |
| decode mean | 449,740 | 438,135 |
| decode p95 | 484,229 | 474,078 |
| decode max | 489,136 | 478,980 |
| `speech_to_lsps_q16`, inclusive | 97.1k | 86.7k |
| `fxp_aks_to_mag2`, inclusive | 341.4k | 329.8k |
| stack through encode | 672 B | 680 B |
| stack through decode (cumulative) | 1,372 B | 1,444 B |
| `.text` | 18,592 B | 19,276 B |
| `.rodata` (image) | 88,785 B | 88,785 B |
| `.bss` | 6,748 B | 6,748 B |

**Stack and code size.** Inlining raises the decoder's worst-case stack by 72 B, the encoder's by 8 B, and `.text` by 684 B. Those are the costs of this change.

**What remains.**

- `fxp_aks_to_mag2` is still 75% of decode. Almost all of it is the per-bin LSP product: about 20 instructions per factor, 10 factors per bin, 256 bins, two calls per frame for the two interpolated frames.
- `cheb` (about 37k per frame, inside `speech_to_lsps_q16`) can't use a precomputed grid, because the root search restarts off-grid after each root.
- Further large savings would need a different spectrum method, which would no longer be bit-identical.

## Verification

- `scripts/test-fxp.sh` passes, with the LSP and energy agreement figures unchanged.
- Corpus encode and decode are byte-identical on all 492 items against a host build of the previous commit.
- With `-DFXP_CHECK`, UBSan and ASan and no recovery, all 492 items encode and decode with no findings.
- The emulator `packets.bit` and PCM checksum are identical.
- `scripts/emu-build.sh` finds no softfloat, libm or `__muldi3` references.

## Reproduce

```bash
sh scripts/test-fxp.sh
```

```bash
docker run --rm -v "$PWD:/w" -w /w p32-csky-gcc14:14.2.0-arm64 sh scripts/emu-build.sh codec results/fxp-hot3-O2
```

```bash
scripts/emu-run.sh results/fxp-hot3-O2 profile
```

```bash
.venv/bin/python scripts/emu-analyze.py results/fxp-hot3-O2
```

```bash
.venv/bin/python scripts/profile-report.py results/fxp-hot3-O2
```
