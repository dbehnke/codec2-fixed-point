# 13. Faster fixed-point real FFT (bit-identical)

## Status

The 512-point real FFT in `codec/src/fxp_fft.c` (doc 08 and earlier) has been restructured for speed. The encoder uses it in `fxp_rfft512` (analysis) and the decoder in `fxp_irfft512` (synthesis). The arithmetic is unchanged, so every output is bit-identical to the previous commit:

- Encode and decode of all 492 corpus items (tune and eval) are byte-identical to the previous commit.
- 300,000 random and edge-case inputs each go through `fxp_cfft256`, `fxp_rfft512` and `fxp_irfft512` against the previous commit's code: identical outputs and exponents, with UBSan, ASan and `-DFXP_CHECK`. The inputs cover full scale, all −32768, sparse, tiny amplitudes and speech slices.
- The emulator `packets.bit` is identical (md5 `dadc57004a7ed241fac1c2a649a4265e`), and so is the decoder PCM checksum.

On the emulator, encode drops from 389,758 to 366,118 instructions per frame and decode from 479,617 to 449,740. Those are instruction counts, not cycles. There is still no real-time claim; that needs hardware timing.

## What changed

Every change only reorders work or skips passes whose results are already known. The same integer expressions run on the same values.

1. **Group-major butterfly order with a per-stage twiddle table.** The old loop ran j outer and the group inner. It computed three twiddle indices per j and kept six twiddle values live, which spilled them to the stack on the 16-register CK803S.
   - The new loop walks each group of 4q and reads w1, w2, w3 in order from a new generated table, `fxp_tw256[243]`. The table holds copies of the `fxp_w512` entries for q = 4, 16 and 64 (972 B of `.rodata`; `test_fxp_fft` checks it entry by entry).
   - The butterflies of one stage touch disjoint elements, so the order does not change any value.
   - Each stage is inlined with q as a compile-time constant, so p[q], p[2q] and p[3q] become immediate offsets from one pointer. The twiddled butterfly now runs without stack spills.
2. **The last stage's bound is reused.** `fxp_rfft512` re-measured the OR bound of the complex FFT output before the real split. The last radix-4 stage already computes that bound over the same values, so it is passed out instead.
3. **Normalisation is fused with the digit reversal (forward).** `block_normalise` scaled every element in one pass, then `fxp_cfft256` permuted them in another. Now one pass permutes and scales, still once per element, with the same `shl` and `shr_r`.
   - After an exact left normalising shift, the first stage's bound is known to lie in [2^13, 2^14): it is (b << sh) | (2^sh − 1), with b << sh in that range and sh ≤ 13. Only `r4_shift` reads it, and it gives s = 2 there. So the bound pass is skipped, and an `FXP_ASSERT` checks this in test builds.
   - A rounded right shift (input bound ≥ 2^14) still measures the bound.
4. **The inverse uses the same bound rule before the split.** `split()` only compares the bound with 2^13 and 2^14.
5. **The conjugations are fused (inverse).** IFFT(Z) = conj(FFT(conj Z)) used to cost two extra passes over 256 elements. The first conjugation now happens inside the digit-reversal pass, which also measures the bound over the stored values. The second happens in the last radix-4 stage's stores, and that stage skips its unused bound. Both negations are exact, since |i| < 21477.

## Cost (emulator)

QEMU instruction counts per frame, `-O2`, comparing `results/fxp-fft2-O2` with `results/fxp-nlp2-O2`. These are not cycles.

| | before | after |
|---|---:|---:|
| encode mean | 389,758 | 366,118 |
| encode max | 414,032 | 390,391 |
| decode mean | 479,617 | 449,740 |
| decode max | 519,017 | 489,136 |
| `fxp_rfft512`, inclusive (2 calls) | 107.6k | 83.9k |
| `fxp_irfft512`, inclusive (2 calls) | 109.2k | 79.4k |
| stack through decode (cumulative) | 1,432 B | 1,372 B |
| stack through encode | 672 B | 672 B |
| `.text` | 16,352 B | 18,592 B |
| `.rodata` (image) | 87,813 B | 88,785 B |
| `.bss` | 6,748 B | 6,748 B |

The kernel microbenchmark (`harness/emu/kbench.c`, four amplitude classes, identical checksum) gives per call:

| | before | after |
|---|---:|---:|
| `fxp_rfft512` | 50.9k–53.8k | 37.8k–41.9k |
| `fxp_irfft512` | 56.6k | 41.7k |

**Code size.** Inlining the stages with constant spans, plus the separate forward and conjugating last stages, adds 2,240 B of `.text`. That is the price of the spill-free inner loop.

**What remains.**
- The twiddled butterfly is about 125 instructions. Of these, 12 are multiplies, 24 compute the magnitude bound that the next stage's shift needs, and the rest are the adds, shifts, loads and stores of the radix-4 butterfly.
- Removing the bound or fusing stages would change the block scaling, so it would no longer be bit-identical.

## Verification

- `scripts/test-fxp.sh` passes. The FFT SNRs against double and kiss are unchanged, and a new check confirms `fxp_tw256` against `fxp_w512`.
- Corpus encode and decode are byte-identical on all 492 items against a host build of the previous commit.
- The 300,000-iteration exactness harness, comparing against the previous `fxp_fft.c`, found no mismatches.
- With `-DFXP_CHECK`, UBSan and ASan and no recovery, all 492 items encode and decode with no findings.
- The emulator `packets.bit` and PCM checksum are identical.
- `scripts/emu-build.sh` finds no softfloat, libm or `__muldi3` references.

## Reproduce

```bash
sh scripts/test-fxp.sh
```

```bash
docker run --rm -v "$PWD:/w" -w /w p32-csky-gcc14:14.2.0-arm64 sh scripts/emu-build.sh codec results/fxp-fft2-O2
```

```bash
scripts/emu-run.sh results/fxp-fft2-O2 profile
```

```bash
.venv/bin/python scripts/emu-analyze.py results/fxp-fft2-O2
```

```bash
.venv/bin/python scripts/profile-report.py results/fxp-fft2-O2
```
