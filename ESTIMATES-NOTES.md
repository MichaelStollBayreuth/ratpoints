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

(results to come)

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
