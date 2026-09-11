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
