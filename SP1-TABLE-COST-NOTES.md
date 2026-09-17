# The sp1-table-cost branch: TODO item 29

Step 2 of the tuning session (2026-09-17 evening), off `v2.3` at a102ecb.
The rule that ends the first phase weighs what a modulus costs and what its
survivors cost downstream, neither of which it saw before.

## The rule before

`take_entries` took moduli from the ranked pool while
`bits_per_word * rate > target`, target = `RATPOINTS_SURVIVORS_PER_WORD`
(0.0075, fitted by `make tune` at 16383).  The target stands for the ratio
of two costs: what one more AND per word costs (1, plus fixed costs spread
over the run) and what a surviving word costs from there on (the second
phase's ANDs while it survives them, then extraction, coprimality, third
stage, exact check).  Both were frozen at the values of the tuning height.

## What the counters and a sweep said

Counter build of a102ecb, pinned (`wt-tidy-pt.pinned.*`; TODO 29 came out
of item 28's set-up figures):

| height | sp1 (mean) | tables | set-up/denom | phases 1+2 | rows/curve |
|-------:|-----------:|-------:|-------------:|-----------:|-----------:|
| 200    | 13.37      | 19.9%  | 31.5%        |  3.1%      | 7914       |
| 1000   | 13.34      | 17.2%  | 29.3%        |  9.7%      | 10316      |
| 4000   | 11.97      |  8.4%  | 16.1%        | 28.9%      | 9772       |
| 16383  |  9.97      |  4.3%  |  8.5%        | 57.4%      | 20411      |

Sweep of the forced first-phase count (`rptest -h H -n N`, base binary,
pinned cycles, medians of 3 rounds, relative to the automatic choice;
`sweep-sp1.raw`): the best fixed N is 7 at height 200 (0.746), 11 at 1000
(0.862), 12 at 4000 (1.028; the per-curve choice beats every fixed N there)
and 11 at 16383 (1.065, likewise); 20 and 18 on the point-rich curves at
1000 (0.900) and 16383 (0.871).  So the base over-buys by six primes at 200
and two at 1000, and is right from 4000 on.

## Two attempts

**Cost only** (`new`): take while `bpw*rate > target*cost`, cost =
`prime_cost(p, 1, ...)` = 1 + (tables + sieve_spec + bp_list entries)/U,
the per-word cost the ranking already uses.  Fixes 200 (0.744) but
over-corrects where the run is a few thousand words: 4000 1.065, 16383
1.024-1.033, because the target already contains the fixed costs of 16383
(a factor 1.3-1.5) and because the other side of the ratio moves too.
Lowering the target (-r 0.004) brings 4000 and 16383 back to 1.01 but
would move 200000, where the target is right as it is (`make tunehigh`
kept it).

**Cost and downstream factor** (`new2`, the branch): the survivors are
weighted by `downstream_factor` = (what a surviving word costs with the
`extra` second-phase moduli the run will have) / (what it costs in a long
run), from the two constants that exist: COST_PHASE2 (110) per
second-phase AND while it survives them, COST_SURVIVOR (1400) once it
reaches the extraction; the densities of the moduli the second phase would
take are those next in the pool, the long-run value is an endless second
phase of their mean density.  With extra = 11 the factor is 1.00; with
extra = 0 it is COST_SURVIVOR*(1-r)/COST_PHASE2, about 6.  No constant
changed, no constant added.

Pinned cycles, medians of 3 rounds, against the base (`sweep-r2.raw`):

| suite               | cost only | cost only, -r 0.004 | cost + downstream | ... with COST_SURVIVOR 700 |
|---------------------|----------:|--------------------:|------------------:|---------------------------:|
| rptest -h 200       | 0.759     | 0.749               | **0.766**         | 0.762                      |
| rptest -h 1000      | 0.909     | 0.870               | **0.859**         | 0.870                      |
| rptest -h 4000      | 1.049     | 0.994               | **0.990**         | 0.977                      |
| rptest (16383)      | 1.024     | 1.013               | **0.999**         | 1.008                      |
| rptest-many -h 1000 | 0.872     | 0.843               | **0.865**         | 0.833                      |
| rptest-many (16383) | 0.998     | 1.029               | **0.999**         | 0.992                      |

The variant with the survivor cost halved (the constant's comment says
340 to 1900 were measured, the small figures at 16383) is not better
enough to move a constant for.

## The counts the rule now chooses

Counter builds, pinned (`phasedata-pinned.txt`; means over the 879 curves
that print a [primestats] line):

| suite                    | sp1 before | sp1 after | sp2 before | sp2 after | tables | set-up/denom | rows/curve |
|--------------------------|-----------:|----------:|-----------:|----------:|-------:|-------------:|-----------:|
| rptest -h 200            | 13.37      |  6.82     | 13.37      |  6.82     | 19.9% -> 5.0%  | 31.5% -> 11.1% | 7914 -> 1284 |
| rptest -h 1000           | 13.34      |  9.88     | 13.34      |  9.88     | 17.2% -> 8.9%  | 29.3% -> 17.7% | 10316 -> 4109 |
| rptest -h 4000           | 11.97      | 11.49     | 12.00      | 11.52     |  8.4% -> 7.8%  | 16.1% -> 15.2% | 9772 -> 8579 |
| rptest (16383)           |  9.97      | 10.40     | 11.85      | 12.28     |  4.3% -> 4.3%  |  8.5% -> 8.5%  | 20411 -> 20930 |
| rptest-many (16383)      | 16.97      | 16.98     | 23.94      | 23.95     |  2.5% -> 2.7%  |  5.3% -> 5.3%  | 145726 -> 147124 |
| rptest-many -h 1000      | 22.04      | 19.07     | 22.04      | 19.07     | 31.5% -> 25.8% | 43.2% -> 35.5% | 79395 -> 53701 |
| rptest -h 200000         |  9.01      |  8.99     | 18.82      | 18.81     |  0.3% -> 0.3%  |  1.3% -> 1.4%  | 94706 -> 94679 |
| rptest-high-many 200000  | 19.13      | 19.13     | 30.13      | 30.13     |  0.1% -> 0.1%  |  0.5% -> 0.5%  | 260092 -> 260899 |

The sweep's optima were 7 / 11 / 12 / 11 (random, 200 / 1000 / 4000 /
16383) and 20 / 18 (point-rich, 1000 / 16383).  At 16383 the rule takes
0.4 primes more on the random curves than before: the downstream factor is
1.6 there (extra = 1.9) against a cost factor of 1.3-1.5, so the bar moves
a little down.  At 200000 (extra = 11) nothing moves.

## The paired measurement

pair.sh, new against the base a102ecb, 5 rounds, pinned cycles (`m-base-new2`):

| suite        | cycles B/A (median, min, max)  | instructions | branch misses |
|--------------|--------------------------------|-------------:|--------------:|
| test100      | 0.793 (0.788, 0.821)           | 0.822        | 0.919         |
| test200      | 0.754 (0.747, 0.783)           | 0.783        | 0.908         |
| test1000     | 0.854 (0.846, 0.893)           | 0.852        | 1.022         |
| testmany1000 | 0.868 (0.856, 0.876)           | 0.875        | 1.202         |
| test4000     | 0.987 (0.942, 0.994)           | 0.985        | 0.976         |
| test1        | 1.007 (0.993, 1.025)           | 1.005        | 0.940         |
| test1many    | 1.003 (0.925, 1.038)           | 0.999        | 1.019         |
| testhigh     | 1.004 (0.993, 1.017)           | 0.999        | 1.014         |
| testhighmany | 1.003 (0.999, 1.069)           | 1.000        | 1.007         |

The four standard suites are within the noise; test1's +0.5% of
instructions is the 0.4 primes more, and is what the refit of the target
(step 5) will weigh.

## What is left at height 200

Callgrind of `rptest -h 200` (`cg-*-h200.excl`): 424M instructions before,
332M after.  Now `examine_prime` is 22.6% (75M; the Horner evaluation of f
over the residues of thirty primes, 85k instructions per curve), gmp for the
exact checks and the coefficients about 20%, `sift0` 6.5%, `get_2adic_info`
4.3%, `find_points_work_1` 4%, `sift` 3.8%, the qsort of the candidates
2.2%, `run_shape` 1.5%; the `sieve_init_*` functions, 2.5% each for 43, 41
and 37 before, have left the top of the profile.  Looking at fewer primes at
small heights would be the next step there, if one is wanted (not in this
session).

## What was not done

* The factor (1 - r) of the modulus in hand: still absorbed in the target
  (step 4 of the session).
* The rule looks at the pool in key order for the second phase's moduli;
  the second phase re-ranks them.  Good enough for a factor.
