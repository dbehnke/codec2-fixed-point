# 12. Faster integer NLP (bit-identical)

## Status

The NLP pitch estimator (`codec/src/nlp.c`, doc 10) is rearranged for speed. Its arithmetic is unchanged: every output is bit-identical to the previous commit.

- Encoder packets for all 492 corpus items (tune and eval) are byte-identical to the previous commit.
- The emulator `packets.bit` is identical (md5 `dadc57004a7ed241fac1c2a649a4265e`).
- Encode drops from 498,792 to 389,758 instructions per frame on the emulator (−21.9%). That is an instruction count, not a cycle count, and not DM32 timing.

No real-time claim follows from this. Hardware timing is still needed.

## What changed

There are two changes, both exact rearrangements of integer sums.

### 1. Direct DFT: bins k and 128 − k share one pass

The NLP DFT is an exact int32 direct DFT. It takes 64 int16 samples on a 512-point grid and computes bins 15..128:

re[k] = Σ x[i] C(ik), im[k] = Σ x[i] S(ik), where C and S are the Q14 table cos and sin of 2πm/512, with m = ik mod 512.

The doc 10 preconditions bound every partial sum by Σ|x| · 2^14 < 2^31. So no int32 sum can overflow in any order of addition, and any regrouping gives the same bits.

**Pairing.** i(128 − k) = 128i − ik, and 128i is a multiple of a quarter turn. The table has these exact symmetries:

- C(−m) = C(m)
- C(m + 256) = −C(m)
- S(m) = C(m − 128)

`tests/test_nlp.c` (`test_tables`) checks all of them, with no violations. From them, the twiddles of bin 128 − k come from those of bin k, chosen by i mod 4:

| i mod 4 | cos(i(128−k)) | sin(i(128−k)) |
|---|---|---|
| 0 | C | −S |
| 1 | S | C |
| 2 | −C | S |
| 3 | −S | −C |

Here C = C(ik) and S = S(ik).

So one pass over x computes x[i]·C and x[i]·S once and adds them, with signs, into four accumulators: re and im for bin k, and re and im for bin 128 − k. The loop is unrolled by 4 so the i mod 4 pattern is fixed in the code.

- Bins 15..63 pair with 113..65.
- Bins 64 and 114..128 have no partner in range and keep the single-bin loop.

**Packed table.** A new generated table `fxp_nlp_cs512[512]` holds cos (Q14) in the low 16 bits and sin in the high 16 bits. One 32-bit load then gives both factors. It adds 2,048 B of `.rodata`. `test_tables` checks it against `fxp_cos512_q14`, with sin(m) = C((m + 384) & 511). It comes from `tools/gen/gen-nlp-tables.py`, which leaves the existing tables unchanged.

### 2. Decimating FIR: symmetric taps added first

The 48-tap lowpass is symmetric: h[k] = h[47 − k] (checked in `test_tables`). The FIR keeps its exact split accumulator: a high part (v >> 15) and a low part (v & 0x7fff), summed separately. Each tap pair now adds its two samples' high parts, and its two low parts, before the one multiply by h[k]:

```c
ah += ((p >> 15) + (q >> 15)) * h;
al += ((p & 0x7fff) + (q & 0x7fff)) * h;
```

This halves the multiplies. The sums are the same integers as before, so the result is unchanged.

The doc 10 bound still covers every intermediate. A high part is at most 40961 in magnitude, so a paired term is at most 40961 · 2|h[k]|. The low parts are below 2^15, so a paired low term is below 2^15 · 2|h[k]|. Summed over the 24 pairs, 2 Σ|h[k]| over half the taps is Σ|h| = 45096. So every partial sum, and every single product, is within 40961 · 45096 < 2^31. That is the existing `_Static_assert` in `nlp.c`.

## Cost (emulator)

QEMU instruction counts per frame, `-O2`, `results/fxp-nlp2-O2` against `results/fxp-lpc-O2`. These are not cycles.

| | before | after |
|---|---:|---:|
| encode mean | 498,792 | 389,758 |
| encode p95 | — | 412,106 |
| encode max | — | 414,032 |
| `nlp`, inclusive (2 calls) | 242,557 | 133,523 |
| decode mean | 479,617 | 479,617 |
| stack through encode | 672 B | 672 B |
| `.text` | 15,944 B | 16,352 B |
| `.rodata` (image) | 85,765 B | 87,813 B |
| `.bss` | 6,748 B | 6,748 B |

The image's `.rodata` includes the harness's 76,800 B test corpus. The codec's own `.rodata` grows by the 2,048 B packed table.

The paired DFT inner loop is about 46 instructions per 4 samples × 2 bins. The old loop was 14 per sample per bin.

`nlp` is now about 67k instructions per call, and is still the largest encoder function (34% of encode). The profiler works per function, and `nlp` inlines its helpers, so its cost can't be split by measurement. By instruction count, the DFT is about 44k of the 67k (an estimate). The rest is:

- the DC notch, with two divides by 20 per sample;
- the FIR;
- the window and normalisation;
- the power;
- the searches;
- the memmoves.

The next targets are `fxp_cfft256`, 79k per frame (done in doc 13), and `speech_to_lsps_q16`, 97k inclusive.

## Verification

- `scripts/test-fxp.sh` passes, including the new `test_tables`. FIR exactness and the NLP agreement figures are unchanged from doc 10.
- Corpus encode is byte-identical on all 492 items against a host build of the previous commit.
- With `-DFXP_CHECK`, UBSan and ASan and no recovery, all 492 items encode and decode with no findings.
- The emulator `packets.bit` is identical, and so is decode (479,617).
- `scripts/emu-build.sh` finds no softfloat, libm or `__muldi3` references.

## Reproduce

```bash
sh scripts/test-fxp.sh
```

```bash
docker run --rm -v "$PWD:/w" -w /w p32-csky-gcc14:14.2.0-arm64 sh scripts/emu-build.sh codec results/fxp-nlp2-O2
```

```bash
scripts/emu-run.sh results/fxp-nlp2-O2 profile
```

```bash
.venv/bin/python scripts/emu-analyze.py results/fxp-nlp2-O2
```

```bash
.venv/bin/python scripts/profile-report.py results/fxp-nlp2-O2
```
