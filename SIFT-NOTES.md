# The quarter of a long run that was in no phase: TODO item 16

Branch `sift-fill`, off `v2.3` at 60927bf, 2026-09-11.

The item said that adding up every instrumented region of `make testhigh`
gave 74% of the run, that the other 26% was in no timed region, and that the
largest nameable piece was the pass that fills every bit array with the
2-adic pattern before the first phase ANDs anything into it.

Two findings, and only one of them is the one the item expected.

## The instrumentation was most of the missing quarter

`RP_PHASE_COUNTS` counts the bits set on entry to phase 1 by popcounting
every bit array, and that block sits *before* `RP_TIC(_rp_t1)`. So it is in
no timed region, and there are 1.7e10 bit arrays in `make testhigh`.

| `make testhigh`, instrumented build | cycles |
|---|---|
| `-DRP_PHASE_TIMING -DRP_PHASE_COUNTS` | 239.6e9 |
| `-DRP_PHASE_TIMING` alone | 196.5e9 |

`cyc1` and `cyc2` barely move between the two (118.6 and 50.9 against 116.5
and 50.4), so the 43e9 difference is the counting and nothing else: **18% of
the instrumented run, all of it outside every timed region.** That is most of
the 26%, and it says that every accounting of where the time goes has to come
from a timing-only build.

## Where a long run actually goes

Two new counters in the `[phasedata]` line: `cycfill`/`fillarrays` for the
fill, the boundary masking and the padding, and `cycsift`/`siftcalls` for the
whole of `sift()`, so that what is inside `sift()` but in none of the regions
can be seen. Timing-only build, on the merged `v2.3`:

| | test1 | testhigh | testhighmany |
|---|---|---|---|
| phase 1 | 37.3% | 59.3% | 73.9% |
| phase 2, including the checks | 14.6% | 25.7% | 19.8% |
| the fill | 3.9% | 4.1% | 2.1% |
| set-up, including the tables | 7.4% | 1.2% | 0.6% |
| what is left of `sift()` | 10.4% | 3.9% | 2.1% |
| `bp_list`, outside `sift()` | 6.3% | 1.6% | 0.4% |
| everything else | 20.1% | 4.3% | 1.1% |

So at a large height bound there is no missing quarter: the fill is 4%, the
rest of `sift()` is 4%, and everything outside it is 6%. At the small height
bound the picture is quite different --- a fifth of `make test1` is outside
`sift()` altogether, which is `sieving_info` examining thirty primes and
sorting them once per curve, over a thousand curves.

**One trap worth recording.** The first version of `cycsift` came out
*smaller* than the sum of the regions inside it. The cause was an early
return in `sift()` --- `if(b*inter.low > height)`, the one that gives up when
the remaining numerator intervals are empty --- which skipped the closing
`TOC`. If a containing timer reads less than its contents, look for a return
path before doubting the arithmetic.

## The change

The first phase's first prime writes its registers instead of reading them
back:

    ratpoints_bit_array *siv0 = sieves[0].start;
    ratpoints_bit_array reg0 = bits16 & *siv0++;

so the pass that filled the array is gone, and with it one store *and* one
load per bit array. The boundary masking and the padding up to a multiple of
`RATPOINTS_CHUNK` move into `sift0`, after the first phase rather than before
it, which is the same thing because AND is commutative and is now two bit
arrays per call instead of all of them. `_ratpoints_sift0` gains `bits16`,
`mask_low`, `mask_high` and `n_pad`.

Two cases keep a fill. `sp1 == 0`, which `-n 0` asks for, has no first prime
to fold into; and the unchunked arm (`RATPOINTS_CHUNK` outside 2..16) would
have to duplicate three loops to do the same, so it pays for the pass, which
is exactly what it did before.

`RP_PHASE_COUNTS` no longer popcounts the array either: every bit array
starts from the same pattern, so the count is a multiplication. That removes
the 18% above as a side effect, which is why the accounting in this file can
be trusted in a way the item's could not.

## What it is worth

Whole builds against `v2.3`, since the change is structural and cannot be a
runtime switch.  `scratchpad/bench2.sh`, which reports the median of the
per-round ratios and their range as well as the ratio of the medians, because
on this laptop those two disagree.  Every baseline binary was checked against
its reference output before being timed.

| suite | of medians | of ratios | per-round range |
|---|---|---|---|
| test1, random degree 6 at 16383 | 0.9708 | **0.9739** | 0.953 to 0.989 |
| test1many, point-rich at 16383 | 0.9923 | **0.9934** | 0.979 to 1.004 |
| the degree suite at 200000 | 0.9467 | **0.9467** | 0.939 to 0.957 |
| testhighmany, point-rich at 200000 | 0.9645 | **0.9656** | 0.959 to 0.979 |
| testhigh, random at 200000 | 0.9421 | **0.9438** | 0.939 to 0.946 |

The two large-height suites are the tightest measurements in the set --- on
`testhigh` every one of the seven rounds falls between 0.939 and 0.946 ---
and they are also the largest gain.  That is what a change to the innermost
loop looks like when it is real: no overlap with the baseline at all.  The
only suite that barely moves is the point-rich one at the small height bound,
where the run is dominated by the per-curve set-up rather than by sieving.

The accounting afterwards, in the same timing-only build:

| | test1 | testhigh | testhighmany |
|---|---|---|---|
| phase 1 | 38.2% | 59.9% | 73.9% |
| phase 2, including the checks | 15.3% | 27.4% | 21.2% |
| the boundary words and the padding | **1.9%** | **0.9%** | **0.4%** |
| set-up, including the tables | 7.3% | 1.2% | 0.6% |
| what is left of `sift()` | 11.0% | 4.3% | 2.3% |
| `bp_list`, outside `sift()` | 6.3% | 1.6% | 0.5% |
| everything else | 20.1% | 4.6% | 1.1% |

What was 4.1% of `make testhigh` is 0.9%, and the instrumented run itself
went from 196.5e9 cycles to 188.9e9.  What is left of that row is two bit
arrays per call to `sift0` rather than all of them, and there is no obvious
way to remove it: the ends of the numerator interval have to be cleared
somewhere.

## Correctness

Identical output to `v2.3` on `make test`, `make testhigh` and
`make testhighmany`, and separately on `make test1` and `make testdegrees`
built at 64, 128, 256 and 512-bit registers, with `RATPOINTS_CHUNK=1` at 64
and 256 bits, and under `-n 0`, `-n 1`, `-n 2 -N 2`, `-P 0` and `-A 0`.  The
first of those matters most: `-n 0` is the one configuration with no first
prime to fold the pattern into, and it takes the other arm of the branch.
