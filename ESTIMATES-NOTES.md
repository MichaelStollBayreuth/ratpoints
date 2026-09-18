# Step 4 of the tuning session: the estimate corrections -- notes

Branch `estimates` off `v2.3` (c9dea1b), 2026-09-18.  Michael: "Merge
d423aa9 and start step 4."  The corrections are taken as one group and
accepted only with the refit of step 5 (`make tune`, `make tunehigh`):
each alone can look worse, because the constants were fitted with the
old estimates and absorbed their errors.  Every correction is measured
alone for the record, on the quantity it corrects (the counter build's
[runshape] and [primestats] lines) and on cycles (pair.sh against the
base, pinned).

## The list

From the Left-overs of TODO.md and the items' accounts:

1. run_shape's Jacobi factor: the symbol is modelled as letting half the
   denominators through; the exact count (item 24's review) is 0.53.
2. U at small height bounds: each denominator's numerators are rounded
   up to whole bit arrays (and to whole 64-bit words per class) while the
   estimate counts words: Uact/Upred 1.6-2.0 at height 1000.
3. bits_per_word on the square-denominator paths (item 26's skeptic): the
   mean over all 64 classes, where a monic odd-degree curve visits only
   the twelve square classes (mean 22% higher).
4. The density convention (item 21): a prime counts the class it divides
   as 1, which is 1/p^2 too much; a prime power counts exactly.
5. The factor (1 - r) of the modulus in hand in the rule that ends the
   first phase (item 29): what the modulus removes is rate*(1 - r), the
   rule compares rate; absorbed in the target.
6. From item 30: the per-call cost of a modulus (15-20 cycles per sift0
   call and modulus: the row pointers' reductions and the tail legs),
   COST_SETUP's kind at half its size, per call rather than per
   denominator.
7. From item 30: the per-denominator fetch of a modulus's row from beyond
   L1, min(p, arrays per denominator) bit arrays, 1-5 cycles per cache
   line depending on where the tables sit (point-rich curves at 16383:
   1.85 cycles per AND against 1.40).

## What was done

Eight commits on the branch, one correction each, the first seven in the
order of the list and the eighth a consistency fix of the second (below);
the default build's output is unchanged by construction (the sieve's
parameters change, the points do not), `make test` passes after each.

1. 5efc63c -- `keep *= 0.53` in run_shape.  The recorded runs of item 30
   confirm the count: Dact/Dpred is 1.07 on nine tenths of the random
   curves at 16383 and 200000 (the tenth without the Jacobi test sits at
   1.00), which is 0.535 for the factor.
2. 6394d77 -- numerators_for counts the non-empty intervals; run_shape
   adds one bit array (RBA_LENGTH bits) per interval and denominator
   after the packing, since each interval is rounded outwards to whole
   bit arrays at every packing; u_words is what the sieve sweeps, and a
   new out-parameter u_pad hands the padding back so that the third
   stage's survivor estimate S keeps counting numerators (the masks clear
   the padding).  The recorded runs: Uact/Upred 1.78 at 1000, 1.29 at
   4000, 1.09 at 16383 and 1.07 at 200000 on the random curves, of which
   the denominator count's error (Dact/Dpred 1.27 / 1.16 / 1.07 / 1.07,
   the Jacobi factor and, at the small heights, more) is a common factor,
   so the rounding itself is 1.4 / 1.1 / 1.02 / 1.00 -- what one bit array
   per interval gives at 1-3 / 10 / 56 / 200 bit arrays per interval.
3. 4ba9b88 -- mean_bits_per_word: the mean over the classes the run
   visits (the enumeration of run_shape: k^2 mod 64 for k < 32 with
   squares as denominators, d k^2 with squares times divisors, every
   class with a pattern otherwise), each as often as it comes up and
   weighted by the words it sweeps (1/2^k of its numerators).  The stride
   table rp_inv_stride moved to file scope for both users.
4. d41932c -- examine_prime: r = (np+1)(p-1)/p^2 when denominators
   divisible by p occur (the class's row admits the p-1 numerators not
   divisible by p), np/p otherwise; was (np(p-1) + p)/p^2.
5. 1067d10 -- phase_1_wants multiplies the survivors by (1 - r_n): the
   modulus removes rate*(1 - r), the rule compared rate.  Halves the
   left side for a typical modulus, so the target has to roughly double
   at the refit.
6. 1a392e2 -- RATPOINTS_COST_CALL (16 at first; 50 since the sweep review,
   see "Reviews") per call of the sieve and first-phase modulus.  run_shape counts the calls (a fourth out-parameter:
   one per interval and denominator while the interval fits array_size
   bit arrays, as many as it needs above), sieving_info turns them into a
   cost per word, phase_1_key adds it to the first-phase cost (a new
   parameter, passed through add_moduli as well).
7. a84a9f3 -- RATPOINTS_COST_LINE (2.5 at first; 4.5 since the sweep
   review) per 64-byte line of the row a denominator walks: phase_1_key
   adds COST_LINE*min(p, A)*(bytes per bit array/64)*D/U with A =
   U/(D*RBA_PACK) the bit arrays per denominator.  A modulus above A costs
   (at 256 bits) COST_LINE/8 of an AND more per word, whatever the height;
   below A the term is proportional to p.

Not done: the second phase's key keeps its shape (its per-AND cost showed
no footprint effect in item 30's runs, and its fixed costs are the same
COST_SETUP and table terms as before); it does see the survivors per swept
word of the eighth commit.

## Measurements (review-2026-09-12/measurements-2026-09-18/estimates/)

Counter builds of the base (c9dea1b) and of the group (a84a9f3), automatic
choice, per-curve lines; medians over the curves of Uact/Upred and
Dact/Dpred, means of the counts:

| suite                  | Uact/Upred base | group | Dact/Dpred base | group | sp1 base | group | sp2 base | group | bpw base | group |
|------------------------|------:|------:|------:|------:|------:|------:|------:|------:|------:|------:|
| rptest 1000            | 1.783 | 1.083 | 1.271 | 1.199 |  9.92 |  9.24 |  9.92 |  9.24 | 34.96 | 33.07 |
| rptest 4000            | 1.290 | 1.058 | 1.162 | 1.096 | 11.55 |  9.91 | 11.59 |  9.95 | 34.96 | 33.07 |
| rptest 16383           | 1.092 | 0.989 | 1.074 | 1.013 | 10.40 |  8.84 | 12.28 | 10.86 | 34.96 | 33.07 |
| rptest 200000          | 1.069 | 1.004 | 1.073 | 1.012 |  8.99 |  7.89 | 18.81 | 17.74 | 34.96 | 33.07 |
| rptest-many 1000       | 1.214 | 0.906 | 1.000 | 1.000 | 19.02 | 16.87 | 19.02 | 16.87 | 56.31 | 55.24 |
| rptest-many 16383      | 1.011 | 0.990 | 1.000 | 1.000 | 16.69 | 14.28 | 23.66 | 21.34 | 56.31 | 55.25 |
| rptest-high-many 2e5   | 1.118 | 1.053 | 1.127 | 1.063 | 19.13 | 16.90 | 30.13 | 27.90 | 61.30 | 60.78 |

The estimates are where they were off: the words at 1000 and 4000, the
denominators everywhere the Jacobi test runs.  What is left: Dact/Dpred
1.2 at 1000 and 1.1 at 4000 (not the Jacobi factor; the periodic counts
over short ranges, or the valuation test), the point-rich curves at 1000
over-predicted by a tenth (their intervals are many and short, and one
bit array per interval is too much when the interval is a fraction of
one), the point-rich curves at 200000 under-predicted by 5%.  bpw fell 5%
on the random curves: the classes with the larger strides pack their
bits more densely and sweep fewer words, and the weighting by words gives
them less.  The counts fell 0.7-2.5 moduli everywhere before the refit,
which is the (1 - r) factor (a half of the left side) and the two new
costs.

Pinned cycles, pair.sh, 3 rounds, each cumulative commit against the base
(median of B/A):

| tree (cumulative)     | test1 | test1many | testhigh | testhighmany | test1000 | testmany1000 | test4000 | test200 |
|-----------------------|------:|------:|------:|------:|------:|------:|------:|------:|
| 1 Jacobi 0.53         | 0.997 | 1.003 | 1.001 | 1.005 | 1.001 | 0.995 | 1.001 | 0.991 |
| 2 + padding           | 1.006 | 1.001 | 1.001 | 0.999 | 1.019 | 1.011 | 1.005 | 1.046 |
| 3 + bits_per_word     | 1.007 | 1.000 | 1.003 | 1.009 | 1.016 | 1.016 | 1.010 | 1.042 |
| 4 + density exact     | 0.999 | 0.998 | 1.008 | 0.997 | 0.977 | 1.001 | 1.004 | 1.039 |
| 5 + (1 - r)           | 1.052 | 1.022 | 1.084 | 1.046 | 1.001 | 0.992 | 1.010 | 1.014 |
| 6 + COST_CALL         | 1.034 | 1.047 | 1.089 | 1.063 | 1.004 | 1.010 | 1.008 | 1.014 |
| 7 + COST_LINE (group) | 1.064 | 1.054 | 1.081 | 1.060 | 1.004 | 0.999 | 1.011 | 1.025 |

As expected before the refit: the (1 - r) factor halves the survivor
side and the two costs raise the other, the phase stops 1-2 moduli early
and the second phase pays (branch misses x1.5-2.4).  The padding
(commit 2) costs 4.5% at height 200 and 2% at 1000 on its own.

The first-phase threshold swept on the group's tree (-r; medians of 3
rounds at the small heights, one run at 200000), relative to the base at
its default 0.0075:

| suite                  | 0.002 | 0.003 | 0.004 | 0.005 | 0.006 | 0.0075 | 0.01 |
|------------------------|------:|------:|------:|------:|------:|------:|------:|
| rptest 200             | 1.064 | 1.053 | 1.043 | 1.037 | 1.027 | 1.022 | 1.017 |
| rptest 1000            | 1.028 | 1.010 | 1.007 | 1.012 | 0.999 | 0.997 | 1.006 |
| rptest-many 1000       | 1.031 | 1.015 | 0.998 | 1.001 | 1.007 | 0.996 | 1.008 |
| rptest 4000            | 1.018 | 1.005 | 1.006 | 1.004 | 1.004 | 1.014 | 1.032 |
| rptest 16383           | 0.985 | 0.988 | 1.004 | 1.017 | 1.026 | 1.054 | 1.091 |
| rptest-many 16383      | 0.994 | 0.998 | 1.007 | 1.045 | 1.036 | 1.069 | 1.123 |
| rptest 200000          |   -   | 0.995 | 1.013 | 1.034 |   -   | 1.091 |   -   |
| rptest-high-many 2e5   |   -   | 0.992 | 0.988 | 1.003 |   -   | 1.055 |   -   |

The group wants a threshold near 0.003 (the (1 - r) factor took a half
out of the left side; the costs the rest), and there it beats the base by
0.5-1.5% at 16383 and 200000 -- but loses 2-6% at 200 whatever the
threshold, and 1% at 1000.  The cumulative table puts the loss at 200 on
the padding: with the padded words in U the fixed costs of a modulus
spread thinner, so the rule takes more moduli at the small heights, where
the fixed-count sweep of item 29 had put the optimum at what the base
takes.  The flaw is in the correction, not in the tuning: the fixed
costs are rightly spread over every word swept, and so is the AND, but
the survivors per word the rule compares them with are per word of
numerators -- the padding is masked and carries none.  Consistent is the
survivors per *swept* word, bits_per_word*(U - pad)/U, in the two rules
of the first phase and in the second phase's key (the third stage's S
counts per denominator and is right already).  Done as an eighth commit
(56ffa47); the sweep repeated on it:

| suite                  | 0.0015 | 0.002 | 0.003 | 0.004 | 0.005 | 0.0075 |
|------------------------|------:|------:|------:|------:|------:|------:|
| rptest 200             | 1.014 | 1.010 | 1.002 | 1.006 | 0.996 | 0.999 |
| rptest 1000            | 1.015 | 1.008 | 1.000 | 0.999 | 1.000 | 1.009 |
| rptest-many 1000       | 1.030 | 1.019 | 1.000 | 1.002 | 0.999 | 1.009 |
| rptest 4000            | 1.032 | 1.015 | 1.015 | 1.006 | 1.005 | 1.029 |
| rptest 16383           | 1.012 | 0.969 | 0.998 | 0.999 | 1.026 | 1.049 |
| rptest-many 16383      | 0.998 | 0.990 | 0.987 | 1.008 | 1.040 | 1.064 |
| rptest 200000          |   -   | 0.988 | 0.992 | 1.012 |   -   |   -   |
| rptest-high-many 2e5   |   -   | 1.000 | 1.002 | 1.001 |   -   |   -   |

The small heights are back to neutral at 0.003 (200: 1.002, 1000: 1.000).
The rounds at 16383 spread by 3% either way (2.84-3.02 G cycles for the
base itself), so the basin 0.002-0.004 is flat there within the noise,
at 0.97-1.00 of the base; the refit will land near 0.003, where the group
is neutral at 200 and 1000, 1.5% worse at 4000, and 1-2% better on
rptest-many at 16383 and rptest at 200000.  What the group buys is not speed but estimates that are right
where they were off by up to 80%, so that the constants fitted at 16383
carry to the other heights on their merits.

## Step 5: the fits

The group is accepted with the refit, so the fits are on this branch too.

* The first-phase threshold: the sweep above puts the corrected model's
  basin at 0.002-0.004 at 16383 (flat within the round noise of 3%) with
  0.003 neutral at 200-4000; the compiled default is set to 0.003 and
  tune.sh's ladder re-bracketed around it (0.0015 0.002 0.0045 0.0075).
* tune.sh gets a fifth stage: RATPOINTS_SP3_PER_DENOM (-Q), ladder 0.003
  0.006 0.025 0.05 around the 0.013 that was fitted when the third stage
  was set up for every denominator (item 24's P11 made it on demand); make
  tunehigh sweeps it by a factor of two either way like the others.
* `make tune` and `make tunehigh` are run on the group with these
  defaults; what they move goes into ratpoints.h (tuning.mk is the
  per-machine file, the compiled defaults are what 2.3 ships).
* RATPOINTS_COST_PHASE2: a -D pair at 60 and 100 against the default 110
  to show the indifference, as item 24 did for COST_BP.
* RATPOINTS_SP3_PER_SURVIVOR stays (the per-survivor test's cost did not
  change); the Left-overs section of TODO.md is emptied into the items'
  accounts when the branch merges.

After the reviews (9709fb4: the constants in the block's unit, the exact
padding), the estimates against the actuals once more, medians over the
curves:

| suite              | Uact/Upred base | group | (U/D)act/(U/D)pred group | Dact/Dpred group | sp1 base | group (at 0.003) |
|--------------------|------:|------:|------:|------:|------:|------:|
| rptest 200         |   -   | 1.411 | 1.015 | 1.390 |   -   |  6.61 |
| rptest 1000        | 1.783 | 1.191 | 0.993 | 1.199 |  9.92 |  9.36 |
| rptest 4000        | 1.290 | 1.084 | 0.989 | 1.096 | 11.55 | 10.37 |
| rptest 16383       | 1.092 | 1.000 | 0.987 | 1.013 | 10.40 |  9.62 |
| rptest-many 1000   | 1.214 | 1.002 | 1.002 | 1.000 | 19.02 | 16.90 |
| rptest-many 16383  | 1.011 | 1.000 | 1.000 | 1.000 | 16.69 | 15.67 |

The words per denominator are right to 1.5% at every height.  What is
left is the denominator count at the small heights (1.39 at 200, 1.20 at
1000, 1.10 at 4000), which is not the Jacobi factor (1.01 at 16383): the
small denominators pass the tests more often than the asymptotic fractions
say.  It enters the model only through the table term, p*min(D,p)/U; the
per-denominator terms use D/U, which is right.

The threshold: the sweep on 9709fb4 read 2-4% worse than the base at the
small heights for every threshold, but its rounds pair one base run with
six settings and the small suites run 50 ms; the interleaved pairs (5
rounds) are what to read.  At 0.003 (the compiled default) against the
base: test200 1.014, test1000 1.009, testmany1000 1.027, test4000 1.000,
test1 0.995, test1many 0.985, and from the sweep testhigh 0.993,
testhighmany 0.981.  At 0.0015: 1.031 / 1.029 / 1.048 / 1.015 / 0.987 /
0.987 -- the 16383 suites like it a little better, the small heights not
at all.  0.003 stands; `make tune` decides between it and its neighbours.

The one suite the group loses on is the point-rich curves at 1000 (+2.7%):
the rule stops at 16.9 moduli where the base took 19.0 and the fixed-count
sweep of item 30 has its flat optimum at 18-20 (16: +7%, 20: +2.7%).  The
cost side is right there -- 3.1 units per word measured (item 30's cycles
per AND at 4 bit arrays per call) against 3.4 modelled with the per-call
and row-fetch terms -- so it is the survivors' side that is low on those
curves at that height, COST_SURVIVOR (measured 340-1900, set at the
200000 value) standing for a survivor that at 1000 meets five third-stage
primes and often the exact check.  Not tuned here: it is a measured
constant and the only place it reaches the first phase is where the run
has no second phase.  Two more probes of that suite (5 rounds each): the
third stage's cost (-Q 0.003 / 0.006 / 0.025, -P 6) moves it between
1.014 and 1.030 against 1.022 at the defaults, and the threshold at
0.0015 makes it worse (1.048), so neither the first phase's count nor the
third stage's is the lever; the 2% stays unexplained and small.

**With every parameter fixed by hand** (-n 10 -N 10 -P 3, 5 rounds) the
group's code costs nothing: rptest at 1000 0.997, at 16383 (-n 10 -N 12
-P 2) 0.992.  On the point-rich curves at 1000 the same fixed count runs
at **0.62** of the base: the ranking is the difference.  The base's table
term, spread over its under-counted U, still drags nearly silent small
primes into the first phase there (7 at density 0.75, 11 at 0.75, 29 at
0.70), the group takes 43, 59, 67 at 0.45-0.5 instead, and the survivors
reaching the extraction fall from 3.4 M to 1.9 M, the exact checks from
532 K to 274 K, on a suite where the extraction and the checks are three
fifths of the time.  At the automatic count the group takes 16.9 moduli
there where the base took 19.0, and the two effects net to the +2.7%.

### make tune (16383, 3 rounds, 2026-09-18 16:35)

Against the current 0.003 / 11 / 1.6e6 / 38 / 0.013 (the sum of the times
of rptest and rptest-many, each candidate timed back to back with the
current settings, medians of 3):

| stage | candidates (ratio to current)                                              | taken |
|-------|-----------------------------------------------------------------------------|-------|
| 1 -r  | current +0.3%, 0.0015 -1.0%, 0.002 +0.2%, 0.0045 +2.7%, 0.0075 +8.6%          | 0.0015 |
| 2 -R  | current -0.8%, 4 -1.1%, 6 -0.6%, 9 -1.7%, 11 -1.0%, 13 +0.2%, 18 +4.0%        | 9 |
| 3 -U  | current +0.1%, 3e5 +0.7%, 6e5 -1.3%, 1.6e6 -0.3%, 2.4e6 -2.4%, 5e6 -0.2%      | 2.4e6 |
| 4 -C  | current -0.0%, 10 -0.3%, 20 +0.1%, 38 -0.7%, 70 -3.3%, 140 -1.9%              | 70 |
| 5 -Q  | current +0.5%, 0.003 -1.1%, 0.006 -1.0%, 0.013 -1.4%, 0.025 -1.8%, 0.05 -0.6% | 0.025 |

Verdict: nothing beat the current settings by the required 3% (the chain
ends at -1.8%, the current settings measured themselves between -0.8% and
+0.5%), so they are kept.  The threshold's stage is the one with a shape:
0.0075, the old default, is 8.6% worse on the corrected model, 0.0045
2.7%, and 0.0015-0.003 are level -- the basin the sweeps showed.  The
other four are flat within the noise at 16383, the third stage's cost
included (0.003 to 0.05 within 1.2% of each other): on these suites
RATPOINTS_SP3_PER_DENOM decides nothing, which is what a cost that the
on-demand set-up removed should look like; it stays at 0.013 as the
value that is not wrong anywhere.

### make tunehigh (200000, ROUNDS=1, 2026-09-18 16:45-17:40)

The neighbourhood run (a factor of two either way, two primes either way
in the offset), one round, each candidate against the current settings:

| stage | candidates (ratio to current)                                  |
|-------|----------------------------------------------------------------|
| 1 -r  | current -0.0%, 0.0015 +0.4%, 0.006 +5.3%                       |
| 2 -R  | current +0.5%, 9 -0.0%, 11 -3.5%, 13 +5.1%                     |
| 3 -U  | current -0.5%, 8e5 -0.7%, 1.6e6 +0.7%, 3.2e6 -0.0%             |
| 4 -C  | current -3.3%, 19 +0.0%, 38 -0.4%, 76 -3.9%                    |
| 5 -Q  | current +10.7%, 0.0065 +0.1%, 0.013 -3.0%, 0.026 -0.7%         |

The run declared itself noisy at the end (the current settings measured
10.7% from themselves in stage 5, past the 2% bound; "e=11" and "c=38",
identical to the current settings, had measured -3.5% and -0.4% against
it earlier), so nothing was written -- one round at 200000 is a single
40-second timing per test, and the laptop drifted over the hour.  What it
says all the same: the threshold's stage is clear and agrees with the
sweep on the same tree (sweep-r9: 0.002 0.981/0.991, 0.003 0.993/0.981,
0.004 1.011/0.997 against the base) -- 0.003 is where 200000 wants it too,
and 0.006 costs 5%; the offset does not want 13; the run length and the
table cost are flat.  A three-round tunehigh takes two and a half hours
of an idle machine and was not run; the constants stay at 0.003 / 11 /
1.6e6 / 38 / 0.013, which "make tune" confirmed at 16383 and the sweeps
at 200-4000 and 200000.

## Reviews

Two Opus reviewers (briefs review-2026-09-12/rev-*-estimates.txt) on
56ffa47.  The sweep review found nothing that changes what the program
finds and one thing that changes the model: **the two new constants were
in the wrong unit**.  The block's unit is the cost of one first-phase AND
per 64-bit word -- a quarter of one AND on a 256-bit array, 0.26 core
cycles (COST_CHECK's comment pins it: 306 rdtsc cycles over 0.16-0.21) --
and COST_CALL and COST_LINE had been written in cycles per array-AND, a
factor RBA_PACK = 4 too small, with comments saying "one AND is one
cycle".  Re-derived: the per-call excess is 20-22 cycles per call and
modulus (the excess over 1.03 cycles per AND times the arrays per call, at
56 and at 10 arrays per call), of which the fetch of the small primes'
rows accounts for 7-8, leaving 12-13 cycles = 50 units (the joint fit's
1-7 cycles for the call term is the collinear split of the same total and
not to be read alone); the fetch is 1.2 cycles per line at 16383 by the
fit and by the point-rich floor (1.85 against 1.40 at 0.8 of a row per
denominator), 2-4 at 200000 where the term is a few per cent of an AND
regardless, so 4.5 units.  A modulus above the arrays per denominator now
costs 4.5/8 = 0.56 of an AND more per word (was 0.31), the call cost a
tenth of an AND per word at 16383 (was 0.03).  The other findings were
comments and prose: examine_power still described the old density
convention; take_entries' header and phase_1_wants' list of fixed costs
lacked the (1 - r) factor and the two new costs; "from here to
RATPOINTS_COST_CHECK" swept the composite bound in with the costs; the
manual's changelog paragraph lacked the caveat README's had; the notes
said seven commits, "the second phase's key is unchanged", and derived the
rounding from the Jacobi factor alone.  The refit the three documents
describe as done is step 5 of this branch, so they describe the merged
state.

The skeptic review (behaviour identical on every suite and on odd inputs
at heights 1-20000, valgrind clean, the new out-parameters, guards and
floors checked, n_calls exact below 30000 -- 2,829,895 calls for 2,829,895
intervals at 16383) found two things that matter and one it shared with
the sweep (the unit).  First, **the padding constant**: one bit array per
interval assumes the interval sits at random against the bit-array grid,
but most intervals are cut by the height bound at one or both ends and
those ends sit at the same bit for every denominator of a class; measured
125-205 bits per interval (288 at 200, 192 at 1000, 157 at 16383) against
the 256 charged, which left the words per denominator over-predicted by
10% at 1000 on both suites and 3% at 4000-16383 -- and the notes' cause
for the point-rich over-prediction (intervals shorter than a bit array)
was wrong: the packing empties 0-2 intervals per curve.  Fixed in 9709fb4:
clipped_padding gives the exact padding of a cut end per class (at a
bound of 2^n - 1 it is nothing; at 1000 some 24 bits an end), an uncut
end pads half a bit array, numerators_for counts the cut ends.  Second,
**bits_per_word's justification** ran the wrong way: on the square paths
the visited classes' unweighted mean is 16-61% above the old 64-class
mean (not "a fifth"), and the weighting by words takes more than that
away, so the delivered value is 3-11% *below* the old one there (and 0-11%
below on the plain path, where only the weighting acts); the formula is
right -- total bits set over total words swept -- and composes with the
padding factor exactly.  The comment now says so.  Minor: the floors of
run_shape fired after the new out-parameters were derived (now they reset
the padding and keep one call); the counter line printed the unswept
value only (now both, and the padding and calls); the USE_SQUARES1 path
weights every divisor's classes alike, as run_shape does, consistent and
biased towards the larger divisors; min(p, mean A) in place of the mean of
min(p, A_b).

### The COST_PHASE2 pair (item 23's Left-over)

60, 100 and 220 against the default 110, on the group (9709fb4), pinned
cycles, 3 rounds, median of B/A:

| COST_PHASE2 | test1 | test1many | testhigh | testhighmany | test1000 | testmany1000 | test4000 | test200 |
|------|------:|------:|------:|------:|------:|------:|------:|------:|
| 60   | 0.985 | 1.007 | 1.000 | 1.000 | 1.014 | 1.015 | 1.002 | 1.013 |
| 100  | 1.000 | 0.995 | 1.000 | 0.997 | 0.998 | 0.995 | 0.996 | 1.002 |
| 220  | 1.024 | 1.002 | 1.000 | 1.001 | 1.010 | 0.995 | 1.021 | 0.998 |

Indifferent within a factor of about two, like COST_BP in item 24: a
factor two up costs 2% at 16383 and 4000, a factor two down 1.5% at
16383 and 1.3% at the small heights, nothing at 200000.  110 stays.

## What the group is, in the end

Eight corrections and two constants (COST_CALL 50, COST_LINE 4.5), the
threshold from 0.0075 to 0.003, tune.sh with a fifth stage; the other
four fitted constants unchanged and confirmed by `make tune`.  Against
the base at the branch's defaults (pinned cycles; the four suites and
4000 from 3-round pairs and the sweep, 200-1000 from 5-round pairs):

| test1 | test1many | testhigh | testhighmany | test200 | test1000 | testmany1000 | test4000 |
|------:|------:|------:|------:|------:|------:|------:|------:|
| 0.995 | 0.985 | 0.993 | 0.981 | 1.014 | 1.009 | 1.027 | 1.000 |

Performance is where it was, within 2% either way -- which is what the
group could achieve at best, the fixed-count sweep of item 30 having
shown the base at the optimum count on every suite already.  What the
group changes is the model: the run shape is right to a few per cent at
every height (words per denominator within 1.5% from 200 to 200000,
denominators within 1% where the Jacobi test runs), the rule that ends the
first phase weighs what a modulus removes against everything it costs,
and one threshold fits 200 to 200000 -- before the corrections, with the
constants in their unit and the padding exact, the group had wanted a
different threshold at every height, and the base's constants had been
right at all of them only because their errors cancelled.  The one
visible consequence is under the hood: at a fixed count the corrected
ranking runs the point-rich curves at 1000 in 0.62 of the time.

Left as found: the denominator count at the small heights (1.2-1.4 low at
200-1000, not the Jacobi factor); the point-rich curves at 1000 at +2.7%
with no single lever; COST_SURVIVOR at its 200000 value where a survivor
at 1000 may cost more; `make tunehigh` at three rounds not run.
