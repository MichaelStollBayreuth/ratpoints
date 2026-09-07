# A cost model for choosing `sp1` and `sp2`

TODO item 3 asked whether the phase-2 behaviour depends on the *survivor rate*
rather than on `sp1` alone, and whether that could be turned into a way of
choosing the parameters from the input.  It can.  This file records the model,
the measurements it was fitted to, and how well it does.

**Result in one line:** `sp1` can be chosen from quantities `sieving_info()`
already computes, and doing so is worth nothing at all on random curves --
where the shipped `sp1 = 11` turns out to be right -- but worth a factor of
about two on curves with many rational points, where the default is 79% off on
median and up to 278% off.  `sp2` cannot be chosen the same way; see below.

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
| bits surviving stage 1 | `bits_in * R(sp1)` | 0.93 / 0.99 |
| units entering stage 2 | `arrays * (1 - (1-R(sp1))^m)` | 1.02 / 1.00 |
| `AND` steps in stage 2 | `arrays * sum_j (1 - (1-R(sp1+j))^m)` | 0.99 / 1.00 |

(point-rich sample / random sample)

where `m = bits_in/arrays` is the mean number of bits actually *set* per
bit-array on entry.  `m` is **not** the width: for `which_bits != num_all` only
every other bit is set, and the range masking removes some more.  Measured, `m`
is about 229 on the point-rich curves but only 94 on random ones, so using `W`
in the exponent over-predicts the units entering phase 2 by a factor 2.4 on
random curves.  `ratpoints` knows `m` before sieving.

The middle line is the one that matters, and it is why the survivor *rate* and
not `sp1` is the right variable: what stage 2 pays for is the number of
bit-arrays that are not empty, and that depends on `R(sp1)` and `m` together.

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

## Validation against random curves

Everything above was fitted on the curve lists in `examples/`, which collect
curves with *many* rational points.  `testdata.h` holds 1000 random genus 2
curves with coefficients up to 10 in absolute value, and they are a different
population altogether:

| predicted survivor rate after 11 primes | 5% | 25% | median | 75% | 95% |
|---|---|---|---|---|---|
| random curves (`testdata.h`) | 3.5e-05 | 8.4e-05 | 1.3e-04 | 1.8e-04 | 3.0e-04 |

against 1.4e-05 to 8.5e-03 with median 1.3e-03 for the point-rich sample.  So a
random curve is about an order of magnitude sparser than a typical point-rich
one, and the random population is far more tightly clustered -- a factor of 9
between the 5th and 95th percentiles, against a factor of 600.

Fourteen of them were swept the same way (heights 100000 and 300000, `sp1` in
5..17, `sp2` in 11..26, 616 runs, `model/sweep-random.csv`).  The predictors
carry over unchanged: predicted/measured is 0.99 for all three, using the *same*
formulas and the cost coefficients fitted on the point-rich curves -- a genuine
out-of-population test, not just a held-out curve.

**The shipped `sp1 = 11` is right for random curves.**  It is the optimum in 27
of 28 cases, and the default `11/19` costs +0.0% on median and +1.1% on mean
against the best setting in the grid, worst case +8.3%.  Whatever it was tuned
on, it was tuned well.  The model agrees: it also picks `sp1 = 11` in 27 of 28
cases.  So the model buys nothing here -- and left to itself it *loses* 2.0% on
median, because of `sp2`.

**`sp2` is the weak coordinate.**  Every mistake the model makes on random
curves is an `sp2` overshoot: it picks the largest `sp2` in the grid every time,
where the optimum is 18 in 16 cases and 22 in 10.  The cause is the floor `F`:
the model cannot know how many candidates are non-reduced representations of
genuine points, so it keeps believing that another phase-2 prime will remove
something.  The cost of the overshoot is small because the surface is flat in
`sp2`, but it is real.

**So use the model for `sp1` and a fixed rule for `sp2`.**  Taking `sp1` from
the marginal rule and setting `sp2 = max(19, sp1 + 8)`:

| policy | random curves | point-rich curves |
|---|---|---|
| the shipped default `11/19` | **+0.0% / +1.1%** | +79.0% / +98.1% |
| minimising the full model | +2.0% / +2.4% | +3.8% / +4.4% |
| model `sp1`, `sp2 = max(19, sp1+8)` | **+0.0% / +1.0%** | **+4.4% / +8.4%** |
| default unless the model predicts >10% | **+0.0% / +1.1%** | +4.6% / +5.7% |

(median / mean excess over the best setting in the grid).  The third row is the
recommendation: it is never worse than the default where the default is good,
and it removes almost all of the penalty where the default is bad.  The fourth
row -- keep the default unless the model predicts a large gain -- is an
equally good and even more conservative alternative.

## Caveats

* The 79% median penalty of the default is a statement about **point-rich
  curves only**.  On random curves the default is optimal, as the section above
  shows.  Anyone reading only the first half of this file would draw the wrong
  conclusion.
* If a fixed default is wanted for point-rich work, `17/26` is the best single
  choice over that sample (1.17x the best on average, against 1.98x for
  `11/19`); but it would be a poor default for random curves.
* The coefficients are for the 256-bit build on this machine and would have to
  be refitted for another width or CPU -- cheaply, since the fit needs only a
  few dozen runs.  The *structure* should carry over; the ratio 84/1.36 is what
  sets the phase boundary, and it is a property of the memory system.  The
  coefficients fitted independently on the two samples agree to within about
  15% for `sp1arr`, `arrays` and `units_1`, which is the evidence that they are
  machine constants and not curve-population artefacts.
* The grid stops at 30 primes (`PRIME_SIZE=7`); the densest curves want all of
  them in phase 1, so their true optimum may be past the edge.  Note that the
  `test2` target in the `Makefile` already runs the record curve with
  `-n 30 -N 30`, which is exactly what the model predicts.
* The random-curve runs are smaller than the point-rich ones (many of these
  curves have a negative leading coefficient, so the numerator range is bounded
  by the real locus rather than by the height).  Only about three quarters of
  such a run is inside `sift0`, so their timings are noisier; the sweep was
  repeated at larger heights for this reason and the conclusion did not change.
