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
curves of `testdata.h`, it is 10^5 to 10^10 on the point-rich curves that
have a Jacobi test at all (16 of the 21 in `testdata-many.h`), where no
pattern of that length can be built and the prototype fell back on the
symbol.  Second, the table costs memory in the height bound (100 KB at
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
coefficient that fits a long always fits the 8 KB of tables (at most seven
odd-exponent primes below 1024 multiply to less than 2^63, and no seven
take more than 6100 bytes); one beyond a long can have more, and then the
mpz symbol is used as before.

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
Nothing changes below a height bound of 63001, and no reference output
pins the list of excluded denominators.

Not done: the review's remark that with the arrays complete the Jacobi
factor in `run_shape` is 0.65 rather than 0.5.  That is a model constant,
to be tried as a `-D` pair in a tuning session, not in a step whose output
must not change.

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

Two model constants are now overstatements: `RATPOINTS_COST_BP` is charged
for every prime including those of the third stage, whose per-denominator
step is gone, and `RATPOINTS_SP3_PER_DENOM` stood for a fill that no longer
happens per denominator.  Both are compiled-in constants of the parameter
model, for a tuning session.
