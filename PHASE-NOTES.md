# The two sieving phases across register widths

Measurements made on 2026-09-07 with the `-DRP_PHASE_TIMING` instrumentation
added to `sift.c` on this branch.  The question was the one on the TODO list:
does the second phase of the sieve benefit from wide registers the way the
first one does, or would a narrower phase 2 be better?

**Short answer: it depends entirely on how much survives phase 1, and the two
regimes point in opposite directions.**  There is no single best width.

## How it was measured

`-DRP_PHASE_TIMING` brackets the two phases with `rdtsc`; `-DRP_PHASE_COUNTS`
additionally counts the units surviving phase 1.  The exact `gmp` check runs
*inside* phase 2 and is timed separately, because it is the same work for every
build and for point-rich curves it is large.

`rdtsc` counts reference cycles at a constant rate, so it measures time, not
core cycles.  `perf` was not available (`kernel.perf_event_paranoid` is back to
3 after the reboot), so all numbers below are minima over 4-5 rounds with the
variants rotated within each round and pinned to CPU 0.  Repeated rows agree to
about 1% for the shorter runs and to 10-13% for the longest ones, so treat
differences below about 5% as noise unless the spread column was small.  Every
build was checked to produce identical output.

Two curves, deliberately at opposite extremes:

| | `1 0 126 0 441`, height 400000 | Drew Sutherland's record curve, height 200000 |
|---|---|---|
| surviving phase 1, per 64-bit word | 1.07% | 36.1% |
| ... as a fraction of 256-bit arrays | 4.17% | 82.1% |
| exact `gmp` check | 0.25% of the run | 24% of the run |

The second curve has a record number of rational points, so `f` is a square
modulo the small primes very often, the phase-1 primes remove almost nothing,
and 82% of all 256-bit arrays reach phase 2.  These are not two points on a
scale; they are two different problems.

## The sparse regime (a typical curve)

`1 0 126 0 441` at height 400000, TSC units in units of 10^9, minimum of five
rotated rounds:

| build | total | phase 1 | phase 2 | split |
|---|---|---|---|---|
| 64, `CHUNK=1` | 20.34 | 16.59 | 3.75 | 82 / 18 |
| 64, `CHUNK=8` | 12.78 | 9.02 | 3.77 | 70 / 30 |
| 64, `CHUNK=16` | 11.71 | 7.56 | 4.14 | 65 / 35 |
| 128 SSE | 11.10 | 5.80 | 5.30 | 52 / 48 |
| 128 SSE `+long` | 9.58 | 5.68 | 3.90 | 59 / 41 |
| 128 AVX128 | 10.42 | 5.77 | 4.64 | 55 / 45 |
| 128 AVX128 `+long` | 9.70 | 5.71 | 3.99 | 59 / 41 |
| **256 AVX2** | **8.60** | 4.72 | 3.89 | 55 / 45 |
| 256 AVX2 `+long` | 8.50 | 4.65 | 3.85 | 55 / 45 |
| 512 emulated | 20.02 | 14.83 | 5.19 | 74 / 26 |
| 512 emulated `+long` | 18.76 | 14.73 | 4.03 | 79 / 21 |

`+long` is `-DUSE_LONG_IN_PHASE_2`; "512 emulated" is `-DUSE_AVX512` *without*
`-mavx512f`, so gcc lowers the 64-byte vectors -- it says nothing about real
AVX-512 hardware and is included only as a structural check.

**Phase 2 is essentially width-independent here.**  Across every build it costs
between 3.75 and 4.03, a 7% band, while phase 1 varies by a factor of 3.6 (4.72
to 16.59).  The lowest phase-2 numbers of all are the 64-bit ones, but only by
about 3%, which is at the edge of the noise.

The reason is in the counts, and it was predicted before the measurement: with a
per-bit survival probability `q` of order `2^-sp1`, the number of *units* that
reach the `sp2-sp1` loop is `(1-(1-q)^W)/W ~ q` per numerator -- independent of
the width.  Measured surviving arrays: 26.70 M at 64 bits, 26.52 M at 128,
26.16 M at 256, 25.48 M at 512.  A survivor is almost never accompanied by a
second one in the same array, at any width.  So widening the registers divides
the *scanning* work by `W` but leaves the per-survivor work untouched, and the
per-survivor work is most of phase 2.

**Phase 1 does not scale with the width either.**  64 (`CHUNK=16`) -> 128 -> 256
gives 7.56 -> 5.71 -> 4.72, i.e. factors of 1.32 and 1.21 for each doubling.
The structural reason is in `init.c`: each sieve table is stored with
`RBA_PACK` copies of the pattern, so phase 1 reads exactly `sp1` bits of table
per numerator *whatever the width*.  The bytes moved are width-independent; only
the instruction count falls.  That is the memory-bandwidth limit the
documentation warns about, made concrete.

## The dense regime (a curve with very many points)

Drew's curve at height 200000, gmp check excluded, minimum of four rounds:

| build | sieve total | phase 1 | phase 2 | vs 256 |
|---|---|---|---|---|
| 64 `CHUNK=16` | 11.48 | 1.22 | 10.25 | +30.3% |
| 128 SSE | 10.97 | 0.91 | 10.06 | +24.6% |
| 128 SSE `+long` | 10.91 | 0.82 | 10.10 | +23.9% |
| **256 AVX2** | **8.81** | 0.77 | 8.04 | -- |
| 256 AVX2 `+long` | 11.13 | 0.79 | 10.34 | +26.4% |
| 512 emulated | 9.60 | 2.16 | 7.44 | +9.0% |
| 512 emulated `+long` | 12.40 | 2.12 | 10.29 | +40.8% |

Everything is reversed.  Phase 2 is now 90% of the sieve and it *does* scale
with the width: 64 bits is 28% worse than 256, and the emulated 512-bit build
-- which pays two AVX2 operations for every 512-bit one -- still beats 256 in
phase 2 by 7%.  On real AVX-512 hardware phase 2 should gain substantially here.

The reason is the mirror image of the sparse case: when most arrays survive,
phase 2 is not scanning, it is doing the `sp2-sp1` `AND`s (and then the
bit-extraction loop) on nearly every array, and one wide `AND` covers `W`
numerators.  Phase 1, meanwhile, has become irrelevant at 8% of the time.

This is also the regime Drew measured in, which is worth remembering when
reading his 13-14%: on that curve the sieve is about 76% of the run and phase 2
is 90% of the sieve.  A *sparse* curve should show a smaller gain from 512-bit
registers -- the model above predicts roughly 10% -- and it would be a good test
of all this to try one.

## What follows for the build

* **Do not enable `USE_LONG_IN_PHASE_2` at 256 bits.**  It is a wash on the
  sparse curve (-1.2%, inside the noise) and a clear loss on the dense one
  (+26%).  The Makefile's "this is usually slower" is right here.
* **Do enable it at 128 bits.**  There it is -13.7% on the sparse curve and a
  wash on the dense one, and the same holds with `USE_AVX128` (-6.9% sparse).
  The Makefile comment is wrong for this width.
* **`CCFLAGS128` should probably be `-DUSE_AVX128`, not `-DUSE_SSE`.**  The
  `USE_SSE` branch of `rp-private.h` still uses the generic fall-back
  `TEST(a) = (EXT0(a) || EXT(a,1))`, which extracts both halves into general
  registers -- the same weakness that was fixed for AVX-512 in 2.2.3, where
  `USE_AVX128` already has the `_mm_movemask_epi8(_mm_cmpeq_epi8(...))` form.
  It is worth 12% of phase 2 (5.30 -> 4.64) whenever phase 2 runs on
  bit-arrays.  `USE_SSE` also uses `__builtin_ia32_andps`, a floating-point
  `AND` on integer data.
* **`-DRATPOINTS_CHUNK` matters enormously for the 64-bit build**, which
  defaults to `CHUNK=1`: going to 16 takes phase 1 from 16.59 to 7.56, and the
  total from 20.34 to 11.71 (-42%).  If the plain-`unsigned long` build is ever
  used seriously, its default should not be 1.  (Note this also means the
  64-bit column of any width comparison is meaningless unless the chunking is
  matched.)

## Open questions

* Where does the time in phase 2 actually go in the dense regime -- the
  `sp2-sp1` `AND`s, or the bit-extraction loop with its `relprime` call per
  surviving bit?  A counter for surviving *bits* would settle it, and it
  decides whether the width helps for a structural reason or an incidental one.
* The picture above has two points; the survival rate is the parameter, and it
  ought to be swept.  This is the same variable that
  `PERFORMANCE-NOTES.md` (branch `interface-phase1-phase2`) identifies as the
  right one for choosing the phase-2 group size and `sp1`.  Choosing the
  *phase-2 representation* from it belongs in the same fit.
* All of this wants re-measuring with `perf stat -e cpu_core/cycles/` on an
  idle machine; `sudo sysctl kernel.perf_event_paranoid=2` is needed first.
