# A third sieving stage, indexed from the x-coordinate

Notes on TODO item 10, made while implementing it on the branch
`third-sieve-stage`.  What follows is the reasoning and the measurements;
the code and its comments say what the program does.

## What the stage is

After the second phase a denominator has a handful of surviving numerators
left, and each of them goes to `_ratpoints_check_point`, which homogenises
`f`, evaluates it in multi-precision arithmetic and takes an integer square
root.  The third stage tests those survivors against further primes first,
one numerator at a time, by computing the index into the table of admissible
residues instead of reading a bit out of a sieve table:

    f(a/b) is a square mod p   <==>   is_f_square[(a * b^-1) mod p]

which is what the sieve tables encode as well.  So it rejects exactly what a
table for that prime would have rejected, and never a genuine point, because
the table of squares counts zero as a square.

The point of doing it this way is that it needs no table.  A prime used in
phase 1 or 2 costs `O(p)` per denominator to build one, which is what
`sieve_init` does and what item 4 was about; a prime used here costs one
subtraction per denominator to carry `b` along, and one multiplication and
one reduction per survivor.  So the stage can use primes that would never be
worth building a table for, and it can reach up the prime list to where the
tables would be largest.

## Where the numbers said it would pay, and where not

The four test suites differ in one statistic that decides everything: how
many survivors a denominator still has when the stage would begin.

| suite | survivors per denominator | exact check, share of the run |
|---|---|---|
| `test1`         | 0.035 | 1.0% |
| `test1many`     | 1.4   | 4.6% |
| `testhigh`      | 0.42  | 1.8% |
| `testhighmany`  | 9.2   | 8.1% |

The cost of one more prime splits the same way.  Per survivor it is a
multiplication and a reduction; per denominator it is one step of the loop
that keeps `b` modulo each sieving prime, which was measured at about 3.5
reference cycles per prime and denominator, against 250 to 280 for one exact
check beyond the first for its denominator, and another 190 for the part of
the check that is done once per denominator.

So a prime is worth adding while

    S * (1 - q_s - r) > q_d

where `S` is the expected number of survivors per denominator still in play,
`r` the prime's density, and `q_s`, `q_d` the per-survivor and per-denominator
costs as fractions of one exact check.  `S` falls by a factor of `r` with
every prime added, so this stops of its own accord, and where survivors are
thin it stops at once -- which is right, because then the per-denominator
cost is all there is.  On `test1` it stops at none or one, and forcing four
primes there costs 2%.

`S` is not measured but predicted, from the number of numerators a
denominator considers (sampled over the range of denominators from the
positivity intervals and the height bound), the sixteen-fold pre-sieve, and
the densities of the primes of the first two phases.  Checked against the
counters, that estimate is good to about 40% across all four suites, which is
enough: an error of a factor of two moves the chosen number of primes by one.

## Two things that were not obvious

**The coprimality test has to come first.**  It looked at first as though it
should be the other way round, the cheap arithmetic before the gcd.  But the
survivors it removes are precisely the ones the new stage cannot touch: if
`gcd(a, b) = g > 1` then `a/b` is the same rational number as `(a/g)/(b/g)`,
which a smaller denominator has already dealt with, so `f` takes the same
value there and **every** prime accepts it.  Such pairs are shadows of the
curve's actual rational points, and there are about `H/b'` of them for a
point with denominator `b'`, so their number grows like the height bound
while the survivors that die by chance grow like its square.  That shows up
plainly in the counters: on the same thirty curves, the extraction loop runs
13.4 times per exact check at height 16383 but only 2.16 times at 200000.

**The prime limit had to move.**  With the 30 primes the program looks at by
default, a random curve has about 14 informative ones left over past `sp2`,
at density 0.53, which is plenty.  A curve with many rational points has a
median of one, at density 0.93, and 12 of the 30 curves in `testhighmany`
have none at all -- these are exactly the curves that ran out of primes in
item 6, because `f` is a square modulo every residue for the small primes.
Looking at all 53 compiled-in primes instead, those same curves have 22 to 30
informative primes past `sp2` at densities 0.60 to 0.65.  So the suite where
the exact check costs the most is the one with nothing to spend, unless the
stage is allowed to look further -- which costs nothing but the table of
admissible residues, since it builds no sieve table.  The selection therefore
examines further primes when the ones in hand run out, in the same way and
for the same reason as item 6 does for `sp2`.

## The arithmetic

Reducing modulo `p` is the whole per-survivor cost, and a division would
throw away most of what the stage saves.  It is done by multiplying instead:
with `m` the reciprocal `2^64/p` rounded up, the remainder of `u` modulo `p`
is the top half of `(m*u mod 2^64) * p`.  That is exact for every `u` below
`2^32`, which was checked against `%` for all 53 primes over the whole range.
The numerator is shifted by `p*H` first, so that one reduction does the whole
job of `(a * b^-1) mod p` and no sign correction is needed; that keeps `u`
under `2*p*H`, which fits below `2^32` up to a height of several million.
Beyond that the stage divides after all.  `m` costs one division per prime
and curve, in `examine_prime`.

Build with `-DRP_STAGE3_DIVIDE` to put the divisions back and measure the
difference.

## What it is worth

Comparing `-P 0`, which switches the stage off, against the rule:

| suite | height bound | stage off | with the stage | change |
|---|---|---|---|---|
| `test1`         | 16383  | 6245096961   | 6243606879   | -0.02% |
| `test1many`     | 16383  | 6552959821   | 6426129077   | -1.94% |
| `testhigh`      | 200000 | 331276727141 | 326068774491 | -1.57% |
| `testhighmany`  | 200000 | 206640368749 | 192815689505 | -6.69% |

Core cycles from `perf stat -e cpu_core/cycles/u`, one core, median of 15
rounds on the two fast suites and 5 on the two slow ones, the two settings
interleaved within each round.  Only the ratios mean anything: the same suite
built from source that differs only in comments moves by several per cent from
build to build, which is the code-alignment effect noted under item 2.

How many primes the rule picks, over the curves of each suite:

| suite | primes in the third stage |
|---|---|
| `test1`         | 0 on 44% of curves, 1 on 41%, 2 or 3 on the rest |
| `test1many`     | 0 to 8, most of them 2 to 7 |
| `testhigh`      | 1 to 8, most of them 4 or 5 |
| `testhighmany`  | 10 to 12, except two curves that get none |

Forcing a fixed number instead, with `-P`, does no better than the rule on any
suite: on `test1many` the best fixed value matches it to within 0.15%, and on
`testhighmany` to within 0.3%.

## Left undone

* The per-denominator setup of the stage writes four fields per prime, of
  which only the inverse depends on the denominator.  The other three are
  constants of the curve and could be hoisted out, which would take a little
  off the per-denominator cost -- but only in the case the rule already
  declines, so it has not been done.
* When `p` divides the denominator the prime is skipped.  It could instead
  reject every numerator at once when there are no points at infinity mod `p`,
  which is what `is_f_square[p]` says.  Rare, and it would cost a branch.
* The choice is made per curve, from an average over denominators.  It could
  be made per denominator, where the number of numerators is known exactly;
  that is item 14's business.
* `tune.sh` does not sweep `-P`.
