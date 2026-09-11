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
