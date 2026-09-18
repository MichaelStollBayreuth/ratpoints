# TODO item 30: a footprint term in the cost model -- notes

Branch `footprint` off `v2.3` (6efda13), 2026-09-18.  Step 3 of the tuning
session; taken "with the gate": a bounded measurement first, the term only
if the measurement gives a clean answer, else the numbers go into the item
and the step is left for 2.4.

## The question

`prime_cost` charges every first-phase modulus one AND per word whatever
its size, and the second phase's ANDs one COST_PHASE2 each.  Item 21's
uncapped variant measured 1.6 cycles per first-phase AND with rows of
small primes and 2.9 with rows of 220 arrays -- two whole-run numbers.
Before a term is built: is the extra cost a property of the modulus (a
per-AND term) or of the whole set of first-phase rows (a row budget), and
is it large enough to move the choice at 16383 and 200000?

## What touches what

Phase 1 (`_ratpoints_sift0`, the chunked arm) is chunk-outer,
modulus-inner: a chunk of 16 bit arrays sits in registers and each of the
sp1 rows is ANDed into it at its current position, the position advancing
cyclically through the row (p bit arrays of 32 bytes at 256 bits).  One
call handles at most RATPOINTS_ARRAY_SIZE = 256 bit arrays, the calls of
one denominator continue through the same rows.  So

* the **row set** a denominator streams is R1 = sum over the first-phase
  moduli of p * 32 bytes (min(p, A) * 32 for a call of A < p arrays), read
  once per pass over 16 arrays; against L1 (48 KB on core 0, a P-core);
* the **table set** the run cycles through is T1 = sum of p * min(D, p) *
  32 bytes, the row of the next denominator coming from it; against L2
  (1.25 MB) and L3 (12 MB, shared).

Phase 2 reads one bit array per surviving array and modulus from the
tables of the sp2 - sp1 second-phase moduli, at a position the word number
picks: a random access into T2 = sum p * min(D, p) * 32 bytes, so the
table set is what matters there, not a row set.

## The measurement (sweep-fp.sh, analyze-fp.py)

Three counter builds of d423aa9 (`-DRP_PHASE_TIMING -DRP_PHASE_COUNTS
-DRP_PRIME_STATS`, RATPOINTS_COMPOSITE_MAX 1 / 64 / 255), pinned to core
0 under `perf stat`; the random and the point-rich suites at heights 1000,
4000, 16383 and 200000 with `-n` 2..28 and automatic.  The per-curve line
added in d423aa9 gives every curve's cycles of the two phases, ANDs and
moduli; the run's perf cycles over its rdtsc total scale the rdtsc counts
to core cycles.  Reports in review-2026-09-12/measurements-2026-09-18/footprint/.

## Results of the sweep (255 runs, 08:58-10:17)

Core cycles per first-phase AND (c1), AND-weighted over the curves of a
run; A = bit arrays per sift0 call, R1 = the row set of a denominator.

**The AND itself costs one cycle.**  On the random curves at 200000 (A =
199) c1 is 1.03-1.07 for every set with R1 below 25 KB, whatever the
moduli (primes only, or with the composites up to 64).

**Above a row set of some 25 KB the cost rises linearly**, about 0.008
cycles per AND and KB: random curves at 200000, 1.12 at 33 KB, 1.24 at
44, 1.38 at 61, 1.50 at 80; point-rich curves at 200000 from a floor of
1.20-1.22 (R1 up to 35 KB) to 1.31 at 49, 1.44 at 68, 1.63 at 89, 1.78 at
113.  The automatic choice sits at 21 KB (random, c1 1.03) and at 64 KB
(point-rich, c1 1.43) at 200000; at 16383 at 15 KB (random, 1.42) and 47
KB (point-rich, 1.85).  The uncapped composites at the same count move
R1 up by 5-20 KB and c1 with it (200000 random: 1.22 at R1 27 KB for the
six-modulus set, against 1.07 with primes).

**Two fixed costs the model has no term for** show as a floor above one
cycle that depends on the height, not on the row set:
* per call and modulus, some 15-20 cycles (the row pointers' reductions
  and the tail legs): at height 1000 a call sweeps 1-3 arrays and c1 is
  5-40; at 4000 (A about 10) 3-4 with small sets; at 16383 (A 56) 1.40;
  at 200000 (A 199) 1.03.  In the model's units it is a per-denominator
  cost like COST_SETUP's, of about half its size, and calls not
  denominators.
* per denominator and modulus, the fetch of the row (min(p, arrays per
  denominator) bit arrays) from the table set beyond L1: the point-rich
  curves at 16383 (primes of 60-127 with 128 arrays per denominator,
  so no reuse of a row array within a denominator) sit at 1.82-1.88 for
  every count from 4 to 20 (R1 11-53 KB) -- the row set does not matter
  when nothing is reused -- and the random ones at 1.40.  The table set
  T1 (0.5 MB to 16 MB) has no separate effect at either height.

**The second phase shows no footprint effect.**  A least-squares split
of its cycles into a scan cost per array and a cost per AND gives the
same per-AND cost for the three trees at every height (200000 random:
52.5 / 53.3 / 54.2 for table sets of 1.9 / 2.0 / 2.5 MB; point-rich 19.2
/ 19.2 / 19.1 for 1.1 / 1.1 / 3.0 MB).  What the per-AND cost does
depend on is the height (141-151 at 4000, 60 at 16383, 53 at 200000) and
the suite (19 point-rich against 53 random at 200000): the density of
the survivors, not the tables.

**The count rule is already at the fixed-n optimum.**  Total cycles of
the same runs by -n, relative to the automatic run of the same tree: the
automatic choice beats every fixed n in 17 of the 18 (suite, height,
tree) combinations and ties the 18th (rptest-many at 1000, primes only:
fixed 18 at 0.999).  The best fixed n costs 1.6-2.2% more at 200000
(random 12, point-rich 20), 1-6% more at 16383 (12 and 18), 3-5% more
at 1000 and 4000.  So a footprint term cannot gain through the count;
what it could change is which moduli are taken.

**The trees at the automatic choice**: composites up to 64 gain 6% over
primes only at 16383 (random) and 9% at 200000; uncapped gains the same
9.5% at 200000 random (0.5% better than 64, within noise), loses 5.7% on
the point-rich curves at 200000 and 3% at 16383.  This is what item 21
measured, and it is small change for a term.

**Where a size-aware ranking could bite** (rank-estimate.py, from the
recorded pools): at 200000 the table term is 1e-3 of an AND, so the key
ranks by information alone and on the point-rich curves takes primes of
100-250 (rows of 3-8 KB) into a row set of 64 KB.  Ranking by
(1 + beta*p)/information for the same information, with the ramp above
as the cost: point-rich 200000 saves up to 4% of the first phase at
beta 0.002 and loses from 0.01 on; random 200000 loses 9-11% of it at any
beta (more moduli for the same information, rows that fit L1 anyway);
16383 loses 0-6%.  A crude model -- measured directly below.

## The size penalty measured (chain-beta.sh, pair-wt-*)

The simplest per-AND term: a first-phase modulus costs 1 + beta*p per
word instead of 1 (phase_1_key, compile-time RP_BETA, the composite cap
at 64 unless said).  Plain builds of d423aa9 paired against the base,
pinned cycles, 3 rounds, median of B/A:

| tree                | test1 | test1many | testhigh | testhighmany | test1000 | testmany1000 | test4000 | test200 |
|---------------------|-------|-----------|----------|--------------|----------|--------------|----------|---------|
| beta 0.001          | 1.000 | 0.996     | 1.003    | 0.998        | 1.001    | 0.995        | 0.988    | 1.001   |
| beta 0.002          | 1.016 | 0.996     | 1.002    | 0.994        | 1.000    | 1.000        | 1.017    | 1.001   |
| beta 0.005          | 0.979 | 0.992     | 1.015    | 0.998        | 1.002    | 1.009        | 0.993    | 0.999   |
| beta 0.002, cap 255 | 1.004 | 1.014     | 0.963    | 1.045        | 1.003    | 1.000        | 0.999    | 1.007   |

With the cap the penalty moves nothing by more than 2% in either
direction, and not consistently (beta 0.005: test1 -2%, testhigh +1.5%).
Without the cap it finds the 3.7% on the random curves at 200000 that
the uncapped products have to give (item 21 measured the same tree
without a penalty at 0.5% better than the cap, so the penalty is what
lets the model take the right products there), and loses 4.5% on the
point-rich curves as before: there the uncapped tree puts two products
above 64 per curve into the *second* phase (tables of 3 MB against 0.6,
the phase's cycles up 7% on 10% fewer ANDs) and 0.4 into the first (the
row set 66 KB against 64, c1 1.48 against 1.43).

## The gate's verdict

Not passed, by the criterion set beforehand: there is no one-constant
term, and the simplest one is worth nothing with the cap.

* The count rule needs no term: the automatic choice beats every fixed
  count on every suite, height and tree.
* The one place a term would pay is lifting the composite cap at 200000
  on random curves, 3.7% of testhigh.  To take it without the 4.5% loss
  on the point-rich curves the model needs (a) a set-aware first-phase
  cost -- the L1 knee (some 25 KB of the 48 KB) and the slope beyond it,
  with the marginal cost of a modulus including what its row adds to the
  cost of every other modulus's ANDs -- and (b) a table-set cost for the
  second phase (the products' tables against L2).  Three or four machine
  constants for one suite's 3.7%, against a one-line cap that already
  takes 9% of the 12.7% the products have there.  The cap stays as the
  cheap proxy, and the manual's "a per-AND cost that rises with the
  modulus would let the model use the larger products where they pay"
  gets the measured answer.
* What the measurement adds to the model's account, for the estimate
  corrections of step 4: the per-call cost of a modulus (15-20 cycles per
  sift0 call, the row pointers' reductions and the tail legs), which the
  model has no term for and which is COST_SETUP's kind of cost at half
  its size, counted per call (ceil(arrays per denominator / 256)) rather
  than per denominator; and the per-denominator row fetch for moduli
  with p above the arrays per denominator (16383: 1.85 against 1.40
  cycles per AND on the point-rich curves), which is what made the cap
  at 64 right at that height (item 21's "three density conventions"
  found the fitted model absorbing it).

The per-curve counters of d423aa9 stay: they are what made this
measurement possible in a morning.  The RP_BETA hook is removed again.
