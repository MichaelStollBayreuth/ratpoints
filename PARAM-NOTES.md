# Choosing the sieving parameters: items 8, 12, 13 and 14 together

Branch `adaptive-parameters`, off `v2.3` at 8640237.

These four TODO items look like four ideas.  They are one idea seen from four
sides: **the rule that picks `sp1`, `sp2` and `sp3` should use a cost model,
and the cost model has terms that the current rule cannot see.**

    marginal prime is worth adding  <=>  what it saves  >  what it costs

The current rule evaluates a degenerate version of that inequality in which

* every prime costs the same (item 8 says the table for `p` costs `O(p^2)`),
* the run is infinitely long, so set-up is free (item 12 says it is not),
* the exact check costs a fixed amount (item 13 says it scales with the
  degree and the coefficient size),
* and the survivor rate is the predicted `R(n)` (item 14 says measure it).

So the four items add four terms, and each of them can be added, measured and
kept or dropped on its own.  That is how this is organised: not one rewrite,
but a sequence of independent changes, every one of which reduces to the
present behaviour when its new term is negligible.

## The model

Per curve, write

* `U` = numerator **words** swept over the whole run,
* `D` = denominators actually sieved,
* `B` = mean bits set per word on entry (`bits_per_word`, already computed),
* `r_i` = density of admissible numerators mod the `i`th prime, ascending,
* `R(n)` = product of the first `n` of them,
* `d`  = cost of one exact check (gmp Horner plus `mpz_sqrtrem`).

| what | cost per numerator word | scales with the run? |
|---|---|---|
| phase 1, one prime | `a` | yes |
| phase 2, one prime | `b * B*R(sp1)` | yes |
| stage 3, one prime | `c * B*R(sp2)` | yes |
| exact check | `d * B*R(sp3)` | yes |
| **one sieve table** | `k*p*min(D,p)/U` | **no** |
| **one `bp_list` step** | `l*D/U` | **no** |

The last two rows are why `make tune` and `make tunehigh` disagree, and the
`p` in the second-to-last is why the rule mis-ranks primes.

## Plan

**Step 0 -- the two quantities the rule is missing: `U` and `D`.**
Both have closed forms from `args->domain`, `height`, `b_low`, `b_high`, all
of which `sieving_info` has before it chooses anything.  Compute them, print
them under `-DRP_PRIME_STATS` beside the counts the instrumentation already
takes (`_rp_arrays_swept`, `_rp_bp_dens`), and check the prediction against
the count on all four suites.  Nothing depends on a prediction that has not
been checked.

**Step 1 -- item 12: the phase-2 offset from `U`.**

        sp2_extra = round(e_inf / (1 + U0/U))

Two machine constants in place of one.  The measured optimum is 3-5 at height
16383 and 9 at 200000, and the shipped flat 5 costs 8.0% of the pair at
200000 and 17.6% of its point-rich half, so this is the largest single number
on the table.  Success = that 8% recovered without losing anything at 16383.

**Step 2 -- item 8: rank the primes by information per unit cost.**
Sort by `-log(r) / (a + k*p*min(D,p)/U + l*D/U)` instead of by `r`.  Reduces
to the present sort as `U -> infinity`.  The test is the one the item states:
beat the current rule on both populations, paired, and close some of the gap
between the incremental prime rule and a fixed `-p 40`.

**Step 3 -- item 14: measure the survivor rate instead of predicting it.**
The scan already visits every bit array; counting the non-empty ones is an
increment on a path taken half a per cent of the time.  Feed the count back
into `sp2` and `sp3` after a window of denominators, and re-check
periodically.  Safe by construction: sieving only ever removes numerators
that cannot be points, so changing the prime count mid-run changes the time
and nothing else.  This replaces the predicted `S` in the `sp3` rule, which
is currently good to about 40%, and should make
`RATPOINTS_SP3_PER_SURVIVOR`/`_PER_DENOM`/`_COPRIME` matter much less.

**Step 4 -- item 13: the degree axis.**
(a) Fix the regression `PRIME_SIZE` 7->8 left behind: degree 8 lost the
division-free Horner in `sieving_info`.  The fix item 13 prefers is the one
item 11 already built -- multiply-high with a generated reciprocal table --
which helps every degree, not just the one that regressed.
(b) Let `d` enter the rules as an estimate from the degree, the coefficient
sizes and the height, rather than as an assumed constant.
(c) Test data at degrees 3, 4, 7 and 8, checked against Magma, since the
repository has nothing but degree 6 and `rptest.c` hardcodes it.

**Step 5 -- the numbers.**  All four suites against `v2.3`, and all four
regimes against the released 2.2.3, including the prime-starved row that was
invalidated last time.

## Rules of engagement, learned the hard way

* Put both variants in one binary behind an environment variable when the
  change does not move any data structure; compare whole builds at three
  alignments when it does.
* Whole-program comparisons at height 200000 have a floor of about 2 points
  on this machine.  Anything smaller needs a controlled measurement and an
  operation count.
* Confirm a binary by its output before timing it.
* Both populations, always, and reported separately.

---

# What was done

## Step 0: `run_shape()` -- the two quantities the rule was missing

`sieving_info` now computes, before it chooses anything,

* `args->run_words` = `U`, the 64-bit words of numerators the run will sweep;
* `args->run_denoms` = `D`, the denominators it will actually sift.

Both come out of the domain, the height bound, the 2-adic masks `den_bits`
and `num_bits`, the forbidden divisors and the Jacobi test, with a 64-point
midpoint sample over the range of `b` for the numerator measure, which is
piecewise linear there.  Three cases: plain, `USE_SQUARES` (`b = k^2`) and
`USE_SQUARES1` (`b = d*k^2`).

Checked against `_rp_arrays_swept` and `_rp_bp_dens` curve by curve, with a
`[runshape]` line printed per curve when both `-DRP_PRIME_STATS` and
`-DRP_PHASE_TIMING` are on:

| suite | U pred/act | D pred/act |
|---|---|---|
| test1 (random, 16383) | 0.87 | 0.93 |
| test1many (rich, 16383) | 1.10 | 0.98 |
| testhigh (random, 200000) | 1.04 | 0.96 |
| testhighmany (rich, 200000) | 1.26 | 0.96 |

Per curve the middle 80% of `U pred/act` on random curves is 0.6 to 1.1.
That is the accuracy the use needs: `U` only ever appears as `U0/U`,
comparing a fixed cost with a per-word one.

**One trap, walked into.**  `den_bits` bit `j` is the denominator `b = j`
mod 64, not `(b-1)` mod 64, even though the set-up shifts by `(b_low-1) & 63`
-- because the loop shifts *before* it tests.  Reading it the other way
inverted the parity for a third of the curves, gave them `keep = 0`, and made
the prediction 2.4 times too small.  It looked like a modelling error and was
an off-by-one.

## Step 1: item 12, the phase-2 offset from `U`

    sp2_extra = RATPOINTS_SP2_EXTRA / (1 + RATPOINTS_SP2_U0/U)

with `RATPOINTS_SP2_EXTRA` raised from 5 to 9 and `RATPOINTS_SP2_U0` =
1.2e6 words, fitted to the two measured optima (offset 3 at 16383, 9 at
200000) using the *predicted* `U`, so that the fit and the use agree.
`-U 0`, or `sp2_u0 = 0`, restores the flat offset, so one binary measures
both arms and there is no code-alignment confound.

## Step 2: item 13(a), the Horner loop in `examine_prime`

`sieving_info` evaluates `f` at every residue for every prime, which is
`O(degree*p)` per prime per curve and ends in a division by `p`.  Two changes:

* the division becomes a **Barrett reduction** -- `m = floor(2^64/p)`, one
  multiply-high, one multiply, one conditional subtraction.  Unlike the
  Lemire reduction that item 11 put into `sift.c`, this one is exact for a
  full-width accumulator, which is what the Horner needs.
* the accumulator is reduced on a **fixed schedule** rather than when a test
  says it must.  `RP_HORNER_STEPS = LONG_LENGTH/PRIME_SIZE - 1` steps fit
  without a reduction, so at the default `PRIME_SIZE` every degree up to 7
  still reduces only at the end, and degree 8 and above reduce every seven
  steps instead of taking a data-dependent branch after every one.

That is the regression TODO item 13 records: raising `PRIME_SIZE` from 7 to 8
in item 6 moved the division-free boundary from degree 8 down to 7, so genus
3 with an even model started dividing inside every Horner step.  It is fixed,
and the fix helps every degree rather than only the one that regressed.

## Step 3: item 8, ranking the primes by cost as well as by information

`sieving_info` sorted the primes by `r` alone.  A prime's sieve table has `p`
rows and is rebuilt for each denominator class that occurs, at most `p` of
them, so over a run of `U` words it costs `COST_TABLE*p*min(D,p)/U` per word
-- quadratic in the prime and inversely proportional to the length of the
run.  Nothing in the old rule could see that.

The key is now `cost/(-log r)`: what a prime costs divided by what it says.
`-log r` and not `1-r`, because the survival rate is a product, so that is
the quantity a greedy choice should maximise per unit cost.  The primes are
ranked twice, because the cost differs by stage:

* for the first phase, `cost = 1 + table + bp`;
* for the second, `cost = COST_PHASE2*rate + table + bp`, with `rate` the
  measured-in-advance survival rate after the first phase.  Since `rate` is
  about 0.006, the table term is two or three times the per-word term at a
  height bound of 16383, so the size of a prime matters far more in the
  second phase than in the first;
* the third stage builds no table at all, so its primes are ranked by `r`,
  which is what it already did.

`-C 0` (or `cost_table = 0`) drops the table term, and the key is then
monotone in `r`: the old order exactly.

## Step 4: item 14, measuring the survivor rate instead of predicting it

The counts are four: the words swept, the bit arrays that survive the first
phase, the survivors of the second, and the ones that pass the test for
common factors.  The first is free (one addition per call to `sift0`), the
second is an increment on a path taken half a per cent of the time, and the
last two on paths taken a few times in a million.

What they are for is the one thing the prediction can never get right.  If
`a/b` is a rational point then so is `ka/kb` -- the same rational number --
and it gives the same value of `f`, so it passes every prime test there is.
A floor of survivors therefore outlives any amount of sieving, and only the
test for common factors removes it.  Writing

        S(n) = floor + chance*R(n)

two measured rates pin both terms, and `adapt_primes()` then applies the
marginal rule to the curve in hand: a second-phase prime is worth adding
while what it removes, at `COST_SURVIVOR` apiece, beats what it costs, which
is `COST_PHASE2` on every array still in play plus the fixed costs it has to
earn back over the run.  The third stage's rule is the one from item 10, with
the measured survivors per denominator in place of the estimate -- which also
retires `RATPOINTS_SP3_COPRIME`, since the count is taken after the
coprimality test rather than estimated through it.

The first correction comes after a million numerator words and the next at
every doubling, so the estimate improves and the corrections thin out.  At a
height bound of 16383 a whole curve sweeps under a million words, so the
correction never fires there; at 200000 it fires seven to ten times a curve.

**Three things this needed that are not about the rule at all.**  Changing
the number of primes mid-run is safe for correctness -- sieving only removes
numerators that cannot be points -- but not for bookkeeping:

* `sieve_list` must hold every informative prime, not just the ones the three
  stages started with;
* `bp_list` must be sized for all of them, with `sp3_valid` saying how many
  entries are current, so that a prime just brought into play is seeded from
  the denominator rather than stepped from an entry it never had;
* the buffer the sieve tables are built in must be sized for every prime
  *looked at*.  This one was a segfault: the third stage may look past the
  primes the first two phases were given, and those primes never build a
  table -- until the correction promotes one into the second phase.

## The cost constants, measured rather than guessed

Built with `-DRP_PHASE_TIMING -DRP_PHASE_COUNTS`, dividing each part's cycles
by the number of times it ran.  Everything is per numerator word and in units
of what one first-phase prime costs there (`cyc1/(and1*RBA_PACK)`, which is
0.16 to 0.26 rdtsc cycles depending on the suite).

| what | test1 | test1many | testhigh | testhighmany |
|---|---|---|---|---|
| one row of a sieve table | 18.3 | 18.7 | 27.8 | 28.7 |
| one step of `bp_list` | 22.0 | 10.6 | 35.7 | 14.4 |
| filling one `sieve_spec` | 29.7 | 8.8 | 28.5 | 11.6 |

The three parts of the second phase cannot be timed apart in one build, so
they were separated with `-DRP_STOP_AFTER=2` (stop after the scan) and `=3`
(stop after the AND loop) and differencing:

| what | test1 | test1many | testhigh | testhighmany |
|---|---|---|---|---|
| the scan, per bit array (cycles) | 1.68 | 2.20 | 1.59 | 1.87 |
| one phase-2 AND (cycles) | 18.2 | 12.7 | 28.0 | 22.5 |
| one survivor of phase 2 (cycles) | 81 | 113 | 317 | 297 |
| `COST_PHASE2` | 78 | 44 | 170 | 105 |
| `COST_SURVIVOR` | 344 | 398 | 1920 | 1388 |

The last two are what the marginal rule for `sp2` turns on, and their
**ratio** is what decides how many primes it wants: 4.4 and 9.0 at the small
height bound, 11.3 and 13.2 at the large one.  So a single pair of constants
cannot serve both, and the pair to compile in is the large-height one, since
that is the only place the correction ever fires.  The shipped values are
110 and 1400.

Two things worth keeping from that.  **Building the sieve tables is 7.2% of
`make test1` and 0.13% of `make testhigh`** -- not the 22% recorded earlier,
which predates item 4 making `sieve_init` 3.7 times cheaper.  And **filling
`sieve_spec` once per denominator and prime is 5.8% of `make test1`**, which
is comparable to the table building and was not on anyone's list; it is a
fixed cost per denominator, so item 12's model already charges for it, but it
is worth knowing it is there.

---

# Results

Each switch measured on its own inside one binary, against
`-U 0 -R 5 -C 0 -A 0`, which is `v2.3` exactly.  Core cycles under
`perf stat`, pinned with `taskset -c 0`, five rounds at the small height
bound and three at the large one, median of the ratio.

| suite | `-U` (item 12) | `-C` (item 8) | `-A` (item 14) | all three |
|---|---|---|---|---|
| test1 (random, 16383) | -0.77% | **-5.00%** | +0.17% | -4.96% |
| test1many (rich, 16383) | -0.84% | -1.16% | -0.69% | -0.47% |
| testhighmany (rich, 200000) | -8.79% | +0.49% | -4.60% | **-9.06%** |
| testhigh (random, 200000) | +0.61% | +0.08% | -0.04% | +0.24% |

(The `-A` column is the correction as it then stood, taking over both `sp2`
and `sp3`; it is now `-A 2`, and the default takes over `sp3` only.)

Three things stand out.

**The cost-aware ranking is worth 5% of `make test1`**, which is far more
than the offset it was expected to play second fiddle to.  That is where the
model says it should be: a short run has a small `U`, so the `p*min(D,p)/U`
term is large and the difference between a prime of 31 and one of 251 is
most of what either costs.

**The two big wins are in different places and do not get in each other's
way.**  The ranking is worth 5% where the run is short, the offset 9% where
it is long and the curve is point-rich, and together they are 4.96% and
9.06% -- each keeps what it had.  The random curves at the large height bound
are flat within the noise of the method, which is about two points there.

**On the point-rich curves at the small height bound the three together are
worse than any of them alone** (-0.47% against -1.16%).  They are not
independent there: the ranking changes *which* primes the second phase gets,
and it gets cheaper ones, so the number worth having is not the number that
was fitted against the old order.  That is what the refit below is for.
Separately, the fitted offset and the marginal rule are two answers to the
same question about `sp2`, and having both switched on meant the later one
silently overrode the earlier; the marginal rule is now opt-in (`-A 2`).

## A note on what item 8 asked for and what it got

Item 8 says the ranking that matters "is not one order but an assignment of
primes to stages, and a prime that is too expensive to tabulate can still be
worth testing".  The two-pass ranking does that, though not by name: a prime
whose table is dear falls in the *second* ranking, which charges for the
table against a much smaller per-word cost, and so drops past `sp2` into the
pool the third stage draws on -- where it is ranked by `r` alone, because
there it builds nothing.  So an expensive prime is not discarded, it is moved
to the stage that does not pay for its table.

What is *not* done is the reverse: nothing asks whether a prime the third
stage is using would be better in the second.  That would need the third
stage's own marginal rule and the second's to be compared directly, which is
what a single cost model over all three stages would give.

## Refitting with the ranking on

The cost of a sieve-table row is the one constant the model could plausibly
have got badly wrong, so it was swept rather than trusted.  On `make test1`,
against ranking by density alone:

| `-C` | 0 | 10 | 20 | 38 | 70 | 140 |
|---|---|---|---|---|---|---|
| ratio | 1.0000 | 0.9603 | 0.9554 | **0.9515** | 0.9540 | 0.9619 |

The compiled-in 38 is the best of them, and the basin is flat from 20 to 70,
which is where the direct measurement puts it (18 to 30).  So the ranking is
worth 4.85% here with the measured costs in place, and the constant does not
need to be fitted -- measuring it lands inside the flat part.

On the point-rich curves at the same height bound the same sweep is flat ---
1.0000, 1.0027, 1.0000, 1.0070, 1.0011, 0.9876 for the same ladder --- which
is the regime where the curve is starved of primes and there is little to
choose among.  So the ranking is worth five per cent where there is a choice
and nothing where there is not, and it never costs anything.

## How far the correction should go, and the answer is: not as far as it can

With the costs measured and the ranking on, on `make testhighmany`:

| | ratio |
|---|---|
| `-A 0`, no correction | 1.0000 |
| `-A 1`, the third stage only | 1.0007 |
| `-A 2`, the second phase as well | 1.0028 |
| `-A 2` with a flat offset (`-U 0 -R 5`) | 1.0500 |

The last row is the informative one.  **The fitted offset of item 12 beats
the measured marginal rule of item 14**: switching the offset off and letting
the correction do the work instead costs five per cent.  And with the offset
on, the correction adds nothing --- the two land in the same place, and the
fit gets there before the first denominator rather than after the first
million numerator words.

That is a real negative result and worth stating plainly.  The measurement is
not *wrong*: the survivor model it fits is right, the floor is real, and the
counts are cheap.  It is that a two-constant fit, tuned on the same machine,
is already as good as a marginal rule with four measured constants, and it
has no warm-up.  What the counting is still good for is the third stage,
where it replaces an estimate that nothing else supplies and retires a fitted
constant (`RATPOINTS_SP3_COPRIME`) -- at no measurable cost, but at no
measurable gain either on any suite in the package.

## Where the offset wants to be, once the ranking has changed the primes

The ranking gives the second phase cheaper primes, so more of them are worth
having, and the offset fitted against the old order is no longer the right
one.  On `make testhighmany`, ranking on, correction off:

| `-R` | 5 | 7 | 9 | 11 | 14 |
|---|---|---|---|---|---|
| ratio | 1.0000 | 0.9370 | 0.9001 | **0.8816** | **0.8813** |

The optimum has moved from 9 to somewhere between 11 and 14, and the two are
indistinguishable, so the curve is flat there.  That is another 2% on top of
what the scaled offset already recovers, and it is the refit the interaction
in the results table was asking for.

The `-C` sweep at the same height bound is flat (1.0000, 1.0007, 1.0007,
1.0013 for 0, 20, 38, 140), which is again what the model says: `U` is a
hundred times larger there, so the table term is a hundredth of what it is at
16383 and the ranking has nothing to bite on.

At the small height bound, ranking on, correction off:

| `-R` | 2 | 3 | 4 | 5 | 7 | 9 | 11 | 14 |
|---|---|---|---|---|---|---|---|---|
| test1 (random) | 1.0000 | **0.9772** | 0.9975 | 0.9991 | 1.0201 | 1.0447 | | |
| test1many (rich) | | 1.0000 | | **0.9961** | 1.0013 | 1.0061 | 1.0180 | 1.0532 |

So the three points to fit are: offset about 3 at `U = 5.9e5`, about 5 at
`U = 5.1e6`, and 12 at `U = 7.1e8`.  The first and third give
`RATPOINTS_SP2_EXTRA = 13` with `RATPOINTS_SP2_U0 = 2.2e6`, which lands on 3,
9 and 13 for the three.

The middle one is the one the model misses: the point-rich curves at the
small height bound want 5 and are given 9, which costs about 0.6%.  That is
not a failure of the arithmetic but of the shape: their survivor rate is
twice the random curves', so a second-phase prime costs them twice as much
per word, and a formula in `U` alone cannot see that.  The marginal rule can,
which is the argument for it that the timings above do not make.
