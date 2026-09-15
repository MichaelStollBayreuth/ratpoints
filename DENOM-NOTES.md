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
