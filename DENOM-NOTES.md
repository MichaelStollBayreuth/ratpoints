# The denominator side: the Jacobi test by Legendre symbols, a word walk, and the per-denominator set-up

Branch `denominators` off `v2.3` at 4bfc676, 2026-09-15.  The third group
of `review-2026-09-12/REVIEW-REPORT.md`: P1 (the Jacobi symbol test), P12
(the denominator walk), P5 (`bp_list`), the `relprime` half of P10, P13
(forbidden-divisor arrays beyond the compiled primes) and P11 (the stage-3
set-up); TODO item 24.  They all touch the loop over the denominators in
`find_points.c` and what `sift()` does once per denominator, and none of
them touches the sieve itself.  Everything here was measured on the
i7-1355U, pinned to P-core 0, in core cycles, paired and interleaved with
the step before, medians of the per-round ratios, the outputs of every
binary compared with the references before it was timed (`pair.sh` in
`review-2026-09-12/`).

## What was done

### P1: the Jacobi symbol test without the Jacobi symbol

For even degree and a non-square leading coefficient `lcf`, a denominator
`b` is dropped unless `(lcf/b*) = 1`, where `b*` is `b` with the prime
factors of `2*lcf` taken out; `jacobi1()` (or `jacobi()` when `lcf` does
not fit a long) computed that symbol by a binary gcd-like loop, about 220
instructions and six mispredicted branches, for the 27 per cent of the
denominators that get past the tests on `b mod 64` and the forbidden
divisors: 4400 calls per curve at 16383, 7.8 per cent of the instructions
of `make test1`, 12 per cent of its cycles and a third of its branch misses
(the review's figures).

The review proposed a table of the symbol over all denominators, built once
per curve from its periodicity: on the `b` coprime to `lcf` the symbol is a
character modulo `8*rad_odd(lcf)`, and the others reduce to `b/q`.  Both
skeptics built it, with byte-identical output.  It is not what was done,
for two reasons that the data made plain.  First, the period is the odd
radical of `lcf`, and while that is at most 33 on the thousand random
curves of `testdata.h`, it is beyond 500 on 16 of the 21 point-rich curves
that have a Jacobi test at all (`testdata-many.h`), beyond 10^4 on 13 of
them and up to 10^10, where no pattern of that length can be built and the
prototype fell back on the symbol.  Second, the table costs memory in the height bound (100 KB at
200000), wants a cap and a fallback beyond it, and is filled for every
denominator while only a quarter of them are ever looked up.

Instead the symbol is evaluated, but not by a loop.  Write `lcf = +-2^v *
prod q_i^(e_i)`.  Then

    (lcf/b*) = (-1/b*)^neg * (2/b*)^v * prod_i (q_i/b*)^(e_i)

and by quadratic reciprocity `(q_i/b*) = (b* mod q_i / q_i) * (-1)^(...)`
with a sign that depends on `b* mod 4` alone; the first two factors depend
on `b* mod 8`.  So once per curve (`jacobi_setup`) the odd part of `lcf` is
factored by trial division against `prime[]`, a table of the non-squares
modulo each `q_i` with an odd exponent is built (`q_i` bytes, the squares
marked by stepping `x^2` to `(x+1)^2`), and a table of eight signs is
filled; and per denominator (`jacobi_test`) the primes of `lcf` are taken
out of the odd part of `b` and the parity of the product is one table
look-up and one exclusive or per prime, the residues by the multiply-high
reduction of item 11 (`RP_MULDIV`, the quotient form, new in
`rp-private.h`, where `RP_MULMOD`, its limit and `RP_CTZL` moved from
`sift.c` so that both files can use them).  Typically `lcf` has one or two
odd primes, and the test is some fifteen instructions with one data-
dependent branch (whether `q_i` divides `b`).

**When it applies.**  Every odd prime factor of `lcf` must be in `prime[]`
(below 1024), and the denominators below 2^32; otherwise `jacobi_setup`
says no and the loop calls `jacobi1`/`jacobi` as before.  On the suites
that covers every one of the 915 Jacobi curves of `testdata.h` and the 48
of `testdata-degrees.h`, and 8 of the 21 in `testdata-many.h` (the rest
have a prime of five or six digits in `lcf`, where the old symbol was never
a noticeable cost: those curves spend their time in the sieve).  A leading
coefficient that fits a long always fits the 8 KB of tables: distinct odd
primes below 1024 whose product stays below 2^63 sum to at most 6057, for
the six largest and a 7 (the sweep review corrected an earlier count here),
and only those with an odd exponent need a table; one beyond a long can
have more, and then the mpz symbol is used as before.

**Order of the tests.**  The Jacobi test used to come after the valuation
test on the primes dividing `lcf`, which divides; now that it is the cheap
one it comes first, since it rejects half of what reaches it.  Nothing
else in the loop changed in this step.

**Checked** by a differential test (`jcheck.c`, scratchpad) of
`jacobi_test` against both `jacobi1` and `jacobi` for every `b` up to
100000 and 5751 leading coefficients (every non-zero `|lcf| <= 3000`, 400
random products of primes below 1024 with random signs and powers of two,
two beyond a long): 575 million comparisons, no mismatch; the set-up
refuses a prime beyond the table and a height bound beyond 2^32 as
designed.  All suites byte-identical.

### P12: the denominators a word at a time

The checked loop visited every `b` from `b_low` to `b_high`: load
`num_bits[b & 0xf]`, shift the word of forbidden-divisor bits by one, test
both, branch -- a dozen instructions for each of the 71 per cent of the
denominators that those two tests reject (the review's count), and one
reload of the word per 64.  Two of the three tests depend on `b mod 64`
alone: bit `b mod 64` of `den_bits`, and whether `num_bits[b mod 16]` has
any bit set at all (the value of that array only matters once `b` is
sifted).  So they are folded into one word per curve (`keep_bits`), the
forbidden-divisor arrays are ANDed into it once per word of 64 denominators
as before, the first and the last word are masked at `b_low` and `b_high`,
and the loop walks the set bits with `RP_CTZL`, as the extraction in
`sift.c` does; a rejected denominator costs nothing, and the per-denominator
work starts at the Jacobi test.  The convention is the one PARAM-NOTES
recorded: `b` is bit `b mod 64` of word `b div 64` (the old loop shifted
before it tested, which is what put the two in step), and `den_bits` is
laid out the same way.

The order of the denominators is unchanged, so the points come out in the
same order.  Checked, besides the suites, by running the old and the new
program over 16 denominator ranges that start and end inside a word (`-dl 1
-du 1`, `63 65`, `64 64`, `127 128`, `3999 4000`, ...) on seven curves of
degrees 5 to 8, with and without the Jacobi test and the forbidden divisors
(`-j`, `-F 0`, `-F 1`): 560 runs, identical output.

### P5: `bp_list` computed, not stepped

`bp_list[n]` is the denominator modulo the n-th sieving prime, which the
per-denominator set-up in `sift()` needs for every prime of all three
stages.  All four denominator loops stepped it from the previous
denominator, `bp += d` and then `while(bp >= p) bp -= p` (the two loops
over squares through `mod()`), which is one to three data-dependent
branches per prime and denominator, taken about a quarter of the time and
mispredicted accordingly: an eighth of all the branch misses of `make
test1` by the review's count, 6.6 per cent of its cycles by the program's
own `RP_BP` region.  Now `fill_bp_list()` computes every entry afresh from
the denominator by the multiply-high reduction, `RP_MULMOD(b, p, magics[n])`
with the reciprocals the sieve already keeps per curve, and divides when
the denominator is beyond 2^32.  No branch, no `last_b`, no `d`; and the
bookkeeping of which entries were still valid after `adapt_primes()` had
brought another prime into play (`sp3_valid`, a field of `ratpoints_args`
and a dozen lines in each loop) goes with it, since nothing is stepped any
more.  The four copies of the fill are one function, which also makes the
call to `adapt_primes` that precedes it.  In the loop over squares times
divisors of the leading coefficient the fill now follows the valuation test
instead of preceding it, so that a denominator that test rejects does not
get one.

### P10, the `relprime` half: no branches in the gcd

`relprime(a, b)` in `sift.c` decides whether a surviving numerator is in
lowest terms, once per survivor of the second phase.  Its binary gcd
replaced numbers by their odd parts with `while(!(x & 1)) x >>= 1` -- one
data-dependent branch per bit -- and chose which of the two to subtract with
another; the review measured the two at a tenth of all the branch misses of
`make test1` (27 per cent together with the same idiom in `jacobi1`, which
P1 has made rare).  Now the odd part is one `RP_CTZL` and a shift, and the
subtraction step is branchless: `d = m - n`, its sign mask, `n = min(m, n)`
and `m = |d|` by mask arithmetic, then the odd part of `m`.  Checked against
a plain Euclidean gcd on fifty million random pairs at three heights, every
seventh pair given a common factor, and on every pair with `|a|, b <= 300`:
no mismatch.  The `jacobi1` half of the review's item is not taken:
`jacobi1` is now called only for a leading coefficient with a prime beyond
1024, where the loop is not what costs.

### P13: forbidden-divisor arrays up to the square root of the height bound

For even degree, a prime `p` with `(lcf/p) = -1` may not divide the
denominator, and the loop tests for the primes of that kind with the word
patterns of `sieves0` -- which exist only for the compiled sieving primes,
up to 251 with `PRIME_SIZE` 8.  The Jacobi symbol supplies the product
form of the same condition, and the two agree exactly when every bad prime
up to `sqrt(b_high)` is in the arrays: what remains of a denominator after
those is at most one bad prime, which the symbol sees.  At 16383 that
holds already (127 < 251), so nothing changes there; at 200000 the arrays
stopped at 251 < 447, and the review counted 1.2 per cent of the sifted
denominators as `q1*q2` or `2*q1*q2` with both primes beyond the table.

So the search for bad primes in `sieving_info` goes on past the compiled
table, up to `sqrt(b_high)` or to the end of `prime[]` (1021, so up to a
height bound of a million), and builds the patterns for the primes it
takes there: `p` words for the prime `p`, word `r` for the word numbers
congruent to `r`, bit `j` clear iff `p | 64r + j`, exactly what
`gen_find_points_h.c` puts into `sieves0`; 64 stores per prime.  They live
in a buffer that stays with `args` and grows when a curve needs more (45 KB
for the 16 primes between 251 and 447).  The arrays and the `forbidden`
list are sized for `prime[]` now, and the default of `max_forbidden` goes
from 30 to 64: with the word walk a prime in the arrays costs four
instructions per 64 denominators, and 30 was what the compiled primes alone
already reached at 200000 (the review measured `-F 53` to change nothing
there, because the cap was not what limited the arrays; the table was).
Nothing changes below a height bound of 66049 = 257^2, the first prime
past the compiled table squared (17161 with `PRIME_SIZE` 7; never with 10,
where the table is all of `prime[]`), and no reference output pins the
list of excluded denominators.

Not done, and not to be done: the review's remark that with the arrays
complete the Jacobi factor in `run_shape` would be 0.65 rather than 0.5.
The sweep review counted exactly, over all odd `b <= 200000` free of a bad
prime below the array bound, for five leading coefficients: the symbol keeps
0.51 to 0.53 of them with the arrays to 251 and 0.53 to 0.55 with the arrays
to 443.  So 0.5 stands, and so does the code's comment that the test
rejects half of what reaches it.

### P11, first half: the third stage's set-up on demand

`sift()` filled `check_spec` -- the prime, its square table, the inverse of
the denominator modulo it, the reciprocal and the bias -- for every prime
of the third stage on every denominator, though only a denominator with a
coprime survivor of the second phase ever reaches that stage: one in seven
at 16383 on a random curve, one in forty at 200000 (the review's counts).
Now `accepted()` in `sift.c` calls `fill_checks()` on the first coprime
survivor of a denominator, and `sift()` only clears the flag
(`stage3_filled` in `ratpoints_args`).  The residue of `b` modulo each of
those primes is recomputed there by the multiply-high reduction, which lets
`fill_bp_list()` stop at `sp2`; the list still has `sp3_max` entries,
because `adapt_primes` can raise `sp2` that far.  `accepted`, `relprime`
and `stage3` carry `always_inline` and `fill_checks` is `noinline`: the
review found that any call gcc might inline into `accepted` stopped it from
inlining `accepted` into the five extraction sites of `sift0`, at a cost of
a per cent, and that the attribute alone was worth 0.16 per cent.

Two model constants were overstated by this: `RATPOINTS_COST_BP`, which is
remeasured and set to 8 below, and `RATPOINTS_SP3_PER_DENOM`, which stood
for a fill that no longer happens per denominator and, being a tuned
constant, is left to a tuning session.

## What it is worth

Each step against the one before, five rounds of the four suites, base =
4bfc676 built from clean, every binary's output compared with the
references before it was timed.  Cycles: median (min-max) of the
per-round ratios; instructions and branch misses: medians.  Columns:
test1 | test1many | testhigh | testhighmany.

| step | test1 | test1many | testhigh | testhighmany |
|---|---|---|---|---|
| P1 Jacobi by Legendre symbols, cycles | **0.860** (0.835-0.879) | **0.993** (0.990-1.012) | **0.967** (0.964-0.972) | **0.995** (0.986-1.001) |
| instructions / misses | 0.926 / 0.617 | 0.998 / 0.974 | 0.982 / 0.814 | 0.999 / 0.987 |
| P12 word walk, cycles | **0.979** (0.951-1.002) | **1.000** (0.995-1.028) | **0.997** (0.996-1.000) | **1.001** (0.992-1.006) |
| instructions / misses | 0.987 / 0.918 | 1.000 / 0.980 | 0.997 / 0.969 | 1.000 / 0.997 |
| P5 `bp_list` computed, cycles | **0.931** (0.915-0.942) | **0.990** (0.982-1.028) | **0.985** (0.984-0.987) | **0.996** (0.994-1.002) |
| instructions / misses | 0.983 / 0.738 | 0.991 / 0.932 | 0.997 / 0.926 | 0.998 / 0.980 |
| P10 branchless `relprime`, cycles | **0.978** (0.963-0.985) | **0.991** (0.980-1.018) | **1.010** (1.005-1.014) | **1.000** (0.999-1.003) |
| instructions / misses | 0.995 / 0.808 | 0.995 / 0.861 | 0.994 / 1.002 | 0.995 / 0.795 |
| P13 arrays to sqrt(H), cycles | **1.006** (0.995-1.022) | **1.000** (0.990-1.007) | **0.994** (0.992-0.997) | **1.000** (0.996-1.004) |
| instructions / misses | 1.000 / 1.000 | 1.000 / 0.999 | 0.990 / 0.991 | 0.998 / 0.994 |
| P11 stage-3 set-up on demand, cycles | **0.980** (0.955-0.994) | **0.992** (0.990-1.003) | **0.986** (0.978-0.988) | **0.988** (0.983-0.991) |
| instructions / misses | 0.976 / 0.994 | 0.991 / 0.975 | 0.997 / 0.984 | 0.996 / 0.996 |
| **the six multiplied out** | **0.756** | **0.967** | **0.941** | **0.980** |

The review predicted P1 at 8-15 / 0 / 1.5-3 / 0 per cent, P12 at 1-2 /
0 / 0.1-0.6 / 0, P5 at 2-5 / 1-2.5 / 0.7-1.4 / 0.3-0.5, the `relprime`
half of P10 at 0.5-1.3 / 1.5-3 / 0.3-2 / 1.2-3, P13 at 0 / 0 / 1.0-1.1 /
0.25 and P11 at 1.5-2.5 / 1-1.5 / 0.2-0.5 / 0.  What came in:

* **P1** at the top of its range: 14 per cent of test1, 3.3 of testhigh,
  with the branch misses down by 38 and 19 per cent.  The instruction
  count falls by 7.4 per cent on test1, as the review measured for its
  table, so the Legendre form costs nothing the table would have saved.
  The point-rich suites move by half a per cent, within their noise: 13
  of their 21 Jacobi curves have a prime beyond 1021 in the leading
  coefficient and keep the old symbol, and the other eight spend their
  time in the sieve.
* **P12** as predicted, 2 per cent of test1, nothing at 200000, where
  the denominator loop is a small part of the run.
* **P5** above its range: 7 per cent of test1 and 1.5 of testhigh for
  1.7 and 0.3 per cent of the instructions -- the branch misses fall by
  26 and 7 per cent.  The review's own counter had put the stepping loop
  at 6.6 per cent of test1's cycles, so this is that loop's whole cost.
* **P10** at 2.2 per cent of test1 and 1 per cent of test1many, above the
  0.5-1.3 the review gave the `relprime` half alone; its misses fall by a
  fifth on test1 and testhighmany.  On testhigh it measures 1.0 per cent
  *slower* in every round, with 0.6 per cent fewer instructions and the
  misses unchanged.  The gcd runs once per survivor of the second phase,
  and testhigh has few of those per denominator, so nothing in what the
  step does can cost that; it is the code placement of `sift0` moving
  with the inlined gcd -- the effect the three-alignment measurement of
  the whole below exists for.  (See there for the verdict.)
* **P13** as designed: 0.6 per cent of testhigh for 1.0 per cent of its
  instructions (the review's 1.08), nothing below a height bound of
  63001.  test1's +0.5 per cent has unchanged instructions and misses and
  is inside that suite's noise band of about 1.5 points either way.
* **P11** at 2 per cent of test1 and, above its range, 1.4 and 1.2 per
  cent of the two large-height suites: `bp_list` stopping at `sp2` is a
  part of that, since the third stage used four to five primes there.

Multiplied out, the six steps are 24 per cent of test1, 6 per cent of
testhigh, 3 per cent of test1many and 2 per cent of testhighmany.  The
review's central estimate for the group was 15-20 per cent of test1 and
3-5 of testhigh.

### The whole against the base, at three code alignments

Whole builds of different source differ in where gcc places the hot loops,
which is worth up to 10 per cent at large heights (TODO item 2), so the
final tree was paired with the base directly, three rounds each, at the
default placement and rebuilt at `-falign-loops=32` and `=64`.  Cycles B/A,
medians; test1 / test1many / testhigh / testhighmany.

| alignment | the six together | instructions | branch misses |
|---|---|---|---|
| default | **0.761** / 0.961 / **0.937** / 0.979 | 0.873 / 0.975 / 0.957 / 0.987 | 0.336 / 0.747 / 0.713 / 0.761 |
| 32 | **0.771** / 0.972 / **0.942** / 0.979 | 0.871 / 0.973 / 0.954 / 0.985 | 0.333 / 0.767 / 0.696 / 0.764 |
| 64 | **0.760** / 0.966 / **0.940** / 0.984 | 0.875 / 0.973 / 0.959 / 0.988 | 0.331 / 0.751 / 0.713 / 0.753 |

The three placements agree to within a point in every suite (the stepwise
product at the default placement, 0.756 / 0.967 / 0.941 / 0.980, agrees
with the direct measurement to half a point), so unlike the trio's this
group has no placement story: **24 per cent of test1, 6 of testhigh, 3 to 4
of test1many and 2 of testhighmany**, against the review's 15-20 / 3-5 for
the group.  The branch misses of test1 fall to a third, those of the other
three suites by a quarter.  The instruction counts fall by 12.7 per cent on
test1 and 4.3 on testhigh, so cycles fall about twice as fast as
instructions: what this group removed was mispredicted branches, as the
review said it would.

### P10 alone, at three alignments

The one step that measured a loss anywhere was the branchless gcd, 1.0 per
cent slower on testhigh at the default placement in every one of five
rounds, with 0.6 per cent fewer instructions and the misses unchanged.  The
gcd runs once per survivor of the second phase, and at 200000 a denominator
has half a survivor, so the step's own work cannot cost a per cent there.
The pair `p5 -> p10` was therefore rebuilt at `-falign-loops=32` and `=64`
and measured again, three rounds of testhigh and test1:

| alignment | test1 | testhigh |
|---|---|---|
| default (5 rounds) | 0.978 | 1.010 |
| 32 | 0.970 | 1.002 |
| 64 | 0.979 | 0.997 |

So the testhigh figure was the code placement of `sift0` moving with the
inlined gcd, and the step is worth 2 to 3 per cent of test1, 1 of
test1many and nothing at 200000 -- above the review's 0.5-1.3 for this
half of the item, because the misses it removes (a fifth of those left on
test1) were mispredicted more often than the review's counter suggested.
The whole-branch measurement above, which agrees at all three placements,
already contains it.

### P11, second half: the hoist -- measured and dropped

The prime and the row offset of every `sieve_spec` are the same for every
denominator, and `sift()` stored them afresh each time.  The review put
hoisting them (and the constants of `check_spec`, which the lazy fill has
since taken out of the per-denominator path) at -2.1 to -2.4 per cent of
test1's instructions, "cycles about half", and a third of that once the
lazy fill was in.  It was built (two per-curve arrays in `args`, one for
each parity of the numerators, filled by `sieving_info` for every prime in
`sieve_list`; `sift()` picks the copy and sets only the table pointer and
its end -- 15e1856, kept as
`review-2026-09-12/prototypes/p11-hoist-sieve-spec.patch`) and paired
against 6f71bac, five rounds:

| | test1 | test1many | testhigh | testhighmany |
|---|---|---|---|---|
| cycles | 0.999 (0.997-1.005) | 0.999 (0.995-1.008) | 0.997 (0.990-0.999) | 0.999 (0.996-1.001) |
| instructions / misses | 0.976 / 1.008 | 0.979 / 1.028 | 0.993 / 1.008 | 0.996 / 1.006 |

It removes 2.4 per cent of test1's instructions and gains a tenth of a per
cent in cycles: the set-up loop is bound by the chain of dependent loads
`sieve_list[n] -> se -> se->sieve[bp]`, and the stores it no longer makes
were free.  Dropped, for a per-curve array and a split of the fill between
two functions that would have bought nothing; the instruction count was the
wrong guide once more (P8 in TRIO-NOTES, the group-of-four scan in
PERFORMANCE-NOTES).  The instrumented build agrees: the per-denominator
set-up region, table building excluded, is 74 cycles per denominator for
16 primes on test1 -- 6 per cent of the run -- and 14 cycles for 22 primes
on testhigh, which no count of stores explains and a load chain that is
cold on a short curve and warm on a long one does.

## The cost constants, remeasured

As in PARAM-NOTES and TRIO-NOTES: builds with `-DRP_PHASE_TIMING
-DRP_PHASE_COUNTS` at `RP_STOP_AFTER` 0, 2 and 3, old code (4bfc676) and
new (6f71bac) back to back on the same morning, pinned, the phase-2 parts
by differencing (`costs.sh`, `costs.py` in `review-2026-09-12/`).  The unit
is one first-phase AND per word.  Old -> new; test1 / test1many / testhigh
/ testhighmany.

| | old | new |
|---|---|---|
| the unit, cycles | 0.207 / 0.309 / 0.156 / 0.223 | 0.205 / 0.278 / 0.152 / 0.218 |
| one phase-2 AND, cycles | 17.5 / 9.0 / 25.1 / 14.3 | 17.3 / 8.8 / 24.3 / 13.3 |
| one survivor of phase 2, cycles (checks excluded) | 132 / 109 / 489 / 243 | 88 / 72 / 308 / 156 |
| `COST_PHASE2` (110 compiled) | 84 / 29 / 160 / 64 | 84 / 32 / 160 / 61 |
| `COST_SURVIVOR` | 640 / 354 / 3120 / 1090 | 432 / 260 / 2030 / 717 |
| `COST_TABLE` (38 compiled) | 19.7 / 19.6 / 30.3 / 29.1 | 19.8 / 20.4 / 29.7 / 28.3 |
| `COST_BP`, per prime and denominator (24 compiled) | 30 / 9.9 / 38 / 14.3 | 9.5 / 5.5 / 11.1 / 7.3 |
| the set-up of `sift()`, per denominator, table building included | 571 / 646 / 637 / 689 | 516 / 610 / 541 / 539 |

What moved is what the group touched: the cost of a survivor of the second
phase (extraction, the gcd, the third stage) fell by a third everywhere --
P10 and P11 --, the `bp_list` entry fell to a third on the random curves
and to a half on the point-rich ones -- P5 --, and the per-denominator
set-up lost the third stage's part -- P11.  The sieve's own constants did
not move: `COST_PHASE2` and `COST_TABLE` agree with TRIO-NOTES to the
noise, which also says the two trees did the same work.

**`COST_BP`** is the one compiled constant now outside every regime: 24
against 5.5-11.1 measured.  It enters `prime_cost` as a fixed cost per
prime and denominator, so a value too high makes a short run stop adding
primes too early.  Whether the runs care was asked with a `-D` pair
(`-DRATPOINTS_COST_BP=8` against the compiled 24, same code, so no
placement confound; three rounds): cycles 1.000 / 1.000 / 1.003 / 1.003,
instructions identical at 200000 and 0.1 per cent apart at 16383 -- the
choice of primes barely moves, and the run not at all, which is the flat
basin TODO item 8 found for every constant of this model.  The constant is
set to **8** all the same, so that it means what its comment says it
means; nothing needs retuning for it, and `make tune` does not touch it.

`COST_SETUP` (30 per prime) is 17-32 per prime on the four suites and
stays.  `COST_SURVIVOR` (1400 compiled) is read by the adaptive correction
of `sp2` only (`-A 2`), which is off by default; measured at 260 to 2030 it
straddles the compiled value as it did before, and stays.  The third stage's
`RATPOINTS_SP3_PER_DENOM`, which stood for a per-denominator cost the stage
no longer has, is a tuned constant and is left to a tuning session.

## How much of it is this machine (for item 22)

Nothing in this group depends on a cache size, a register width or a core
type; `ARCHITECTURE-TRIAGE.md` has the rows.

* **P1**: general.  The machine features it uses are a 64x64->128
  multiply (every 64-bit ISA) and `__builtin_ctzl` (every compiler the
  program builds with); the tables are at most 8 KB per curve.  The size of
  the gain is the misprediction penalty times the six mispredicts a symbol
  cost, 15-20 cycles each on every current core, so a core with a weaker
  predictor gains more.
* **P12**: general -- fewer instructions per rejected denominator, and one
  `ctz` per accepted one.
* **P5**: general -- a mispredicted branch chain replaced by a multiply;
  the same trade as the trio's P4, with the same exactness bound.
* **P10**: general.  `ctz` and an arithmetic right shift of a negative
  `long`, which the C standard leaves implementation-defined and gcc,
  clang and every ISA the program runs on define as the sign extension.
* **P13**: general (mathematics); 45 KB of patterns at 200000, streamed one
  word per 64 denominators, so nothing to cache.
* **P11**: general -- work not done.  `always_inline` and `noinline` are
  gcc/clang attributes, which the file needs anyway for its vector types.

## Measurement

`pair.sh A B ROUNDS OUT suites` from `review-2026-09-12/`: for each round
and suite it runs A then B or B then A, alternating, pinned to core 0
under `perf stat`, and reports the median, minimum and maximum over rounds
of the per-round ratio of cycles, plus the medians of the instruction and
branch-miss ratios; it refuses to time a binary whose output differs from
the reference.  Every step was built from its commit in a worktree of its
own under the scratchpad, the chain of six pairs, the whole at three
alignments and the cost constants ran unattended (`chain.sh`) on the
otherwise idle machine from 09:14 to 11:33, and the P10 and hoist checks
followed in the same way.  The per-round spread is about 1.5 points either
way on the two three-second suites and half a point on the two long ones;
a step measured at 0.5 per cent on a short suite is inside the noise, and
the notes say so where it applies.

## The reviews

Two Opus agents in their own worktrees at 1caef25, one sweeping for
completeness and consistency, one told to break it.

**The sweep** built fifteen configurations from clean (every one
reproduces the references with no warning beyond the two known AVX-512 ABI
notes; the dividing fallback is clean now), ran `make test`, the debug
build under valgrind, `--leak-check=full` on a curve at 70000 and on
forty-eight curves reusing one `args` (the buffer of P13 grows and is
freed), compared the run-time patterns of P13 with the compiled ones by
building at `PRIME_SIZE` 10 (identical excluded lists and points), checked
the word walk on ninety-six sub-ranges against the full runs, recomputed
every count and percentage in these notes, the README and the change log
(all reproduce, two excepted: see below), and compiled the manual against
the one of 4bfc676 (the same four overfull boxes, none new).  It found no
fault in the code.  What it found was text the code had left behind: the
`bp_list` entry still called "stepped" in five comments, two of them
claiming a per-denominator cost for the third stage's primes; the two
phase-timing region labels these notes quote describing what the regions
no longer contain; the `COST_SETUP` comment calling the `bp_list` entry
"about the same" two lines under a constant three times smaller; the
caller lists of `RP_MULMOD` in `rp-private.h` and the Makefile without the
reduction that calls it thirty million times a suite; the `check_spec`
and third-stage comments describing `b` as carried along; the manual's
fallback sentence naming one condition of three and using an undefined
symbol; two unwrapped lines; the `always_inline`/`noinline` attributes
written out where the file's convention (`RP_CTZL`) keeps the plain build
free of gcc-isms -- now `RP_ALWAYS_INLINE` and `RP_NOINLINE`; and two
wrong statements in these notes: the "seven primes" arithmetic of
`jacobi_setup`'s comment (fourteen small primes fit below 2^63; the bound
that matters, 6057 bytes, holds), and the count "10^5 to 10^10 on 16 of
21" (13 of 21, the other three between 500 and 10^4).  Its last finding is
the one worth the most: the review's "0.65" for `run_shape`'s Jacobi
factor is 0.53 by exact count, so that parked remark is withdrawn above.
All applied (1c3775e).

**The skeptic** set out to prove the branch wrong and reported that it
could not.  It ran its own differential test of `jacobi_test` against
both old symbol routines (578.8 million comparisons over 5788 leading
coefficients and every `b` to 10^5, plus 3.6 million just below 2^32; no
mismatch; the refusals exactly the documented ones), checked the exactness
of `RP_MULDIV`/`RP_MULMOD` on 34 million values including the ends of the
range, matched the branchless gcd against Euclid on 48 million pairs,
transcribed the old and the new denominator loop side by side and ran them
on 200000 randomized configurations -- 117 million denominators visited in
the identical sequence, every boundary the brief names included --,
compared the run-time forbidden-divisor patterns bit for bit with the
definition and word for word with `sieves0`, ran about 2000 paired
old-against-new invocations (degrees 4 to 8, heights 1 to 5e9, `-F` 0 to
1000, `-A` 0 to 2, degenerate and reversed denominator ranges, windows
beyond 2^32 that take every fallback), showed with the `RP_PRIME_STATS`
build that the two trees do the same work counter for counter (852 runs,
also where `-A 2` demonstrably moves `sp2`), ran everything under ASan and
UBSan including a harness that drives 48 curves of alternating height
bounds through one `args` so that the P13 buffer grows, is reused smaller
and grows again, and confirmed with `objdump` that `accepted`, `relprime`
and `stage3` are inlined and `fill_checks` is not.  Its findings: a failed
`malloc` for the P13 patterns recorded a capacity the buffer did not have
(fixed: the length is set on success only, and without memory the primes
beyond the table are dropped for that curve); the `run_shape` comment
still explaining the bit convention by a shift the loop no longer does
(fixed); the P13 threshold in these notes, 63001, which is 66049 = 257^2
at the default `PRIME_SIZE` and moves with it (fixed above); the manual's
fallback sentence, now naming all four conditions; and the same comment
findings as the sweep, already applied.  Nothing else.
