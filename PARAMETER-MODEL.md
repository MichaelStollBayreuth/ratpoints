# A cost model for choosing `sp1` and `sp2`

TODO item 3 asked whether the phase-2 behaviour depends on the *survivor rate*
rather than on `sp1` alone, and whether that could be turned into a way of
choosing the parameters from the input.  It can.  This file records the model,
the measurements it was fitted to, and how well it does.

**Result in one line:** at the optimal `sp1` the fraction of numerators
surviving phase 1 is essentially a constant -- about 1 in 7500, within a factor
of 2.7 over curves spanning a factor of 630 in density and optimal `sp1` from 8
to 30 -- so the rule is simply *add phase-1 primes until about 2% of the
bit-arrays are still non-empty*.  That single constant chooses `sp1` within
about 1% of the best setting on both populations tested, beating both the
shipped default (79% off on point-rich curves) and a seven-coefficient fitted
cost model.  `sp2` cannot be chosen the same way; use `max(19, sp1+8)`.

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

## The survivor rate at the optimum is essentially constant

The coarse grids above cannot answer this, because `sp1 = 11` covers a whole
plateau.  A fine scan settles it: `sp1` from 6 to 30 in steps of 1 along
`sp2 = min(30, sp1+8)`, for six point-rich and six random curves, minimum of
three `perf` runs each, with the optimum located by a parabola through the
three cheapest points (`model/sweep-fine-sp1.csv`).

| | m | best `sp1` | interpolated | `R(sp1*)` | `u(sp1*)` |
|---|---|---|---|---|---|
| random `td:696` | 75 | 9 | 9.4 | 1.8e-04 | 0.013 |
| random `td:313` | 94 | 10 | 10.4 | 1.6e-04 | 0.015 |
| random `td:894` | 106 | 11 | 11.1 | 1.2e-04 | 0.013 |
| random `td:564` | 83 | 12 | 11.6 | 9.0e-05 | 0.007 |
| random `td:38` | 94 | 12 | 12.0 | 1.9e-04 | 0.017 |
| random `td:560` | 241 | 13 | 12.6 | 7.0e-05 | 0.017 |
| rich `ex:10` | 247 | 8 | 8.5 | 1.1e-04 | 0.025 |
| rich `ex:118` | 173 | 14 | 14.1 | 1.5e-04 | 0.025 |
| rich `ex:490` | 235 | 17 | 16.5 | 1.1e-04 | 0.025 |
| rich `ex:676` | 239 | 18 | 17.8 | 1.6e-04 | 0.038 |
| rich `bc:136` | 251 | 25 | 25.3 | 1.1e-04 | 0.027 |
| rich `bc:1944` | 246 | 30 | 30.0 | 1.5e-04 | 0.036 |

The optimal `sp1` ranges from 8 to 30 and the curves span a factor of 630 in
`R(11)`, but **the fraction of numerators still alive after phase 1 barely
moves**: `R(sp1*)` has median 1.3e-04 and a total spread of a factor 2.7.
About one candidate in 7500 survives phase 1, whatever the curve.

The array-level version `u = 1-(1-R)^m` is tighter still *within* a population
(a factor 1.5 across the point-rich curves, 2.3 across the random ones) but sits
at different levels for the two -- 0.026 against 0.014 -- because `m` differs,
about 240 against about 94.  Overall `u` spreads by a factor 5.1 against 2.7
for `R`.

### This gives a rule with one constant instead of seven

Choose `sp1` as the first `k` with `u(k) <= 0.02` (or, in the bit-level form,
`R(k) <= 1e-04`), and `sp2 = max(19, sp1+8)`:

| rule | point-rich | random |
|---|---|---|
| `u(sp1) <= 0.02` | **+0.4% / +1.1%** | **+0.0% / +1.0%** |
| `R(sp1) <= 1e-04` | +0.5% / +1.4% | +0.5% / +2.8% |
| minimising the full fitted model | +3.8% / +4.4% | +2.0% / +2.4% |
| the shipped default `11/19` | +79.0% / +98.1% | +0.0% / +1.1% |

(median / mean excess over the best setting in the grid.)  The one-constant rule
**beats the seven-coefficient cost model on both populations**, and it is what
should actually be implemented: `sieving_info()` has `R` and `m`, so it is a
loop over at most 30 primes and one comparison.

It is also robust.  Any threshold `u` in 0.013..0.026 stays under 2.5% mean on
both populations, and any `R` in 1e-04..2e-04 does the same; only outside
roughly a factor of two either way does it start to cost.  The constant is a
property of the machine, not of the curves -- it is the point where one more
phase-1 prime (1.36 cycles per bit-array) stops paying for itself against the
84 cycles that the first non-empty array costs on entry to phase 2, so
`u* ~ 1.36 / (84*(1-r) + 2.7) ~ 0.03` for a typical density `r ~ 0.5`.  That
back-of-the-envelope figure is the right size, which is the reassuring part;
the measured 0.02 is what to use.

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
