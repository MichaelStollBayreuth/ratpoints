# A cost model for choosing `sp1` and `sp2`

TODO item 3 asked whether the phase-2 behaviour depends on the *survivor rate*
rather than on `sp1` alone, and whether that could be turned into a way of
choosing the parameters from the input.  It can.  This file records the model,
the measurements it was fitted to, and how well it does.

**Result in one line:** a rule using only quantities `sieving_info()` already
computes chooses `(sp1, sp2)` within 4% of the best available setting on
median, where the shipped default `11/19` is 79% off on median and up to 178%
off on the curves measured here.

## The measurements

Nine curves from `examples/`, spanning a factor of 630 in survivor rate, at
three height bounds, over `sp1` in 8..26 and `sp2` in 14..30: 918 runs of the
256-bit build, pinned to CPU 0, cycles from `perf stat -e cpu_core/cycles/`
apportioned to the stages by the `rdtsc` brackets in `sift.c`.  The raw data is
in `model/sweep-low.csv` and `model/sweep-high.csv`; `model/model2.py` does the
fitting, `model/curves.py` parses the curve lists and computes the densities.

The three stages are the `sp1` primes (stage 1), the scan plus the `sp2-sp1`
primes (stage 2) and the exact `gmp` check (stage 3, which runs inside stage 2
and is subtracted from it).

## What can be known before sieving starts

`sieving_info()` in `find_points.c` already computes, for every prime `p`, the
density of admissible residues

    r = np/p,   np = #{ a mod p : f(a) is a square mod p }

(with a correction when there are points at infinity mod `p`), and sorts the
primes by increasing `r`.  So the survivor rate after the best `k` primes is

    R(k) = product of the k smallest r,

available for free before any sieving happens.  Everything the model needs
follows from `R`, the number of bit-arrays to sweep, and the width `W`:

| quantity | prediction | predicted/measured (median) |
|---|---|---|
| bits surviving stage 1 | `bits_in * R(sp1)` | 0.93 |
| units entering stage 2 | `arrays * (1 - (1-R(sp1))^W)` | 1.06 |
| `AND` steps in stage 2 | `arrays * sum_j (1 - (1-R(sp1+j))^W)` | 0.97 |

The middle line is the one that matters, and it is why the survivor *rate* and
not `sp1` is the right variable: what stage 2 pays for is the number of
bit-arrays that are not empty, and that depends on `R(sp1)` and on `W` together.

## The floor

The number of surviving bits does **not** go to zero as `sp2` grows.  On the
reference curve at height 20000 with `sp2 = 19`, exactly 14348 bits survive
however the primes are split between the phases, and only 229 of them are in
lowest terms.  The other 14119 are the non-reduced representations `(ka, kb)`
of the 44 genuine rational points: they satisfy every congruence condition,
because they *are* points, so no amount of sieving removes them.

So the true survivor count is `bits_in * R(sp2) + F` with a floor `F` set by the
points themselves, and adding phase-2 primes stops paying once `bits_in * R(sp2)`
falls below `F`.  `F` is not known in advance -- it is part of the answer being
computed -- but it only ever makes the model over-estimate the benefit of a
large `sp2`, and it is the reason the measured optimum for `sp2` saturates.

## The fitted costs

Fitted over all 918 runs, in core cycles (256-bit build, i7-1355U):

| what | cycles |
|---|---|
| per `sift0` call | 25 |
| **per bit-array per phase-1 prime** | **1.36** |
| per bit-array scanned in stage 2 | 7.06 |
| **per non-empty unit entering stage 2** | **69.5** |
| per `AND` step in stage 2 | 5.16 |
| per candidate surviving stage 2 | 792 |
| per `check_point` call | 837 |

Stage 1 fits to R² = 0.996, stage 2 to 0.989, stage 3 to 0.972, and the total
to R² = 0.990 with a median relative error of 6.8%.

The 1.36 against the 69.5 is the whole story.  A phase-1 prime costs one load
and one `AND` per bit-array, unconditionally.  A phase-2 prime costs nothing on
empty arrays but 70 cycles the first time an array is not empty -- a
dependent load into a large table, off the streaming path.  So phase 1 is cheap
per prime and phase 2 is cheap per *array*, and which one wins depends entirely
on how many arrays are still alive.

## The rule

Both boundaries follow from marginal cost.  Write `u(k) = 1 - (1-R(k))^W` for
the fraction of bit-arrays still non-empty after `k` primes.  Then

* move a prime from phase 2 into phase 1 while

      1.36  <  69.5 * (u(sp1) - u(sp1+1))  +  5.16 * u(sp1)

  -- i.e. while enough bit-arrays are still alive that sieving them one at a
  time in phase 2 costs more than streaming over all of them in phase 1;

* add a phase-2 prime while it removes more work than it costs,

      5.16 * u(sp2) * arrays  <  792 * bits_in * (R(sp2) - R(sp2+1)) .

Applied to the measurements, and with the cost coefficients fitted on the
*other eight* curves each time (leave-one-curve-out), this gives

| | median | mean | worst |
|---|---|---|---|
| minimising the full model | +3.8% | +4.4% | +16.8% |
| the two marginal rules above | +3.8% | +4.3% | +14.6% |
| the shipped default `11/19` | +79.0% | +98.1% | +278.0% |

all measured against the best setting in the grid.  The two rules do as well as
minimising the model, which is what matters: they need only `R`, `arrays` and
`W`, so they can run inside `sieving_info()` in microseconds.

The optima found range from `8/14` on the sparsest curve to `26/26` and `26/30`
on the densest -- the whole width of the grid.  No fixed pair can cover that,
which is the point.

## Caveats

* **The curve sample is biased.** Both lists in `examples/` collect curves with
  many rational points, so the population here is much denser than a random
  curve.  On the sparsest curve tested (`x^6 + 1`) the shipped default is only
  5-12% off, which is presumably how it was chosen.  The 79% median is a
  statement about point-rich curves, not about all curves.  A sample of random
  curves is needed before changing any default.
* If a fixed default is wanted anyway, `17/26` is the best single choice over
  this sample (1.17x the best on average, against 1.98x for `11/19`).
* The coefficients are for the 256-bit build on this machine and would have to
  be refitted for another width or CPU -- cheaply, since the fit needs only a
  few dozen runs.  The *structure* should carry over; the ratio 69.5/1.36 is
  what sets the phase boundary, and it is a property of the memory system.
* `sp2` is the weaker of the two: the measured optimum is often flat over a
  wide range, and the model tends to overshoot it because it cannot know `F`.
  Every case where the rule loses more than 10% is an `sp2` overshoot at the
  smallest height.
* The grid stops at 30 primes (`PRIME_SIZE=7`); the densest curves want all of
  them in phase 1, so their true optimum may be past the edge.  Note that the
  `test2` target in the `Makefile` already runs the record curve with
  `-n 30 -N 30`, which is exactly what the model predicts.
