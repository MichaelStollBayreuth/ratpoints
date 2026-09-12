# The `p | lcf` test on the denominators: TODO item 15

Branch `forbidden-divisors`, off `v2.3` at d2c83ba, 2026-09-12.

The item found forty lines in `sieving_info` that had never run: the arm of
the forbidden-divisor set-up for a prime `p` dividing the leading
coefficient.  This branch makes it run, makes it right, and finds it worth
having.

## What the arm was for

For even degree `d`, a denominator `b` with `p | b` (so `p` does not divide
`a`) has `F(a,b) = c[d] a^d mod p`.  If `c[d]` is a non-square mod `p`, so is
`F(a,b)`, and `b` is excluded: that is the working arm, and the bit arrays
`forb_ba` test for it.  If `p | c[d]`, then `F(a,b)` is `0` mod `p`, which
is a square, and the congruence says nothing.  The dead arm was meant to use
the valuation instead.

## Why it was dead, and why it was wrong

Three defects, of which the first hid the other two.

1. **The guard.**  The arm sat inside `if(!is_f_square[p])`, and
   `is_f_square[p] = squares[pn][c[d] mod p]`, where `squares[pn][0] = 1`
   because zero is a square.  So `p | c[d]` forced the guard false.  It had
   never run on any curve.

2. **The minimum over the wrong terms.**  It computed the smallest `n` such
   that for `v_p(b) >= n` the leading term of `F(a,b)` has strictly smaller
   valuation than every other term, by the same formula `setup_us1` uses for
   the odd-degree case -- but only over the coefficients `c[k]` with
   `p | c[k]`.  A coefficient that is a unit at `p` contributes `v_p(c[k]) =
   0`, the *strongest* constraint, and skipping it treats its valuation as
   infinite.  With `p` not dividing `c[d-1]` the right `n` is `v_p(c[d]) + 1
   >= 2`; the arm got `n = 1` and excluded `p` itself.  Enabled as it stood,
   `make test1` lost 23 points on 12 curves; the smallest case is curve 12,
   `-5 5 8 7 -2 2 -3`, where the two top terms at `(a, b) = (2, 3)` are
   `+192` and `-192`, both of valuation 1, and cancel to `2601 = 51^2`.

3. **The sign.**  `valuation()` reduces `|c|`, so the residue it returns is
   that of `|c[d]|/p^v`, and the arm tested it for being a square mod `p`.
   For `c[d] < 0` and `p = 3 mod 4` that inverts the answer: `c[d] = -18`
   has unit part `-2 = 1 mod 3`, a square, and the test would have said
   non-square and excluded `9 | b`.

And a fourth in waiting: `p^n` was formed by `n` multiplications with no
overflow guard, and `n` can be `v_p(c[d]) + 1`, which for a 400-bit
coefficient at `p = 3` is 250.

## The derivation

Let `v_p(b) = m >= 1` and `w_j = v_p(c[d-j])`, the valuation of the
coefficient of `b^j` in `F(a,b) = sum_j c[d-j] a^(d-j) b^j`.  The term with
`b^j` has valuation `w_j + j*m`.  If one term has strictly smaller valuation
than all the others, it sets `v_p(F)` and `F/p^v mod p`, and then

* `v` odd: `F` is not a square, for any `a` and any `b' = b/p^m`;
* `v` even and `j` even: `F/p^v = r_j a^(d-j) b'^j mod p` with `r_j` the unit
  part of `c[d-j]`; `d-j` and `j` are even, so `F` is a square only if `r_j`
  is one mod `p`;
* `v` even and `j` odd: `b'^j` runs through squares and non-squares as `b'`
  varies (the image of `x -> x^j` has odd index in `F_p^*`, so it is not
  inside the squares), and nothing follows.

If two terms tie for the minimum they may cancel -- that is what happens at
curve 12's point -- and nothing follows.

For `m` large the leading term `j = 0` wins alone, so when `w_0` is odd or
`r_0` is a non-square, every large enough `m` is excluded: "no denominator
is divisible by `p^n`".  The commonest cases: `v_p(c[d]) = 1` and `p |
c[d-1]` exclude `p` itself (Michael's estimate: probability `(p-1)/p^3` for
random coefficients, saving `1/p` of the denominators); `v_p(c[d]) = 1` and
`p` not dividing `c[d-1]` tie at `m = 1` and exclude `p^2` (probability
`(p-1)^2/p^3`, saving `1/p^2`).  Below the tail there can be isolated
valuations: `c[d] = 9u` with `p = 3` not dividing `c[d-1]` excludes `v_3(b)
= 1` (the `b^1` term wins with valuation 1) whatever `u` is -- a square
leading coefficient included -- and then `27 | b` if `u` is a non-square.

`forbidden_valuations()` in `find_points.c` does exactly this: for each `m`
with `p^m <= b_high` it finds the minimising term, and sets bit `m` of a
mask when that term is unique and either its valuation is odd or `j` is even
with a non-square unit.  No `n` is computed; the tail is whatever the mask
has.  If every valuation a denominator can have is excluded, `p` goes to the
bit arrays as before; otherwise a `(p, mask)` entry goes to the array that
the denominator loop tests by division, which now asks whether bit
`v_p(b)` of the mask is set (`valuation1(b, p)` costs one division for the
`b` not divisible by `p`, as the old `b % p^n` did).

**Checked against an independent model.**  `scratchpad/lcf_scan2.py`
derives the excluded set from the Newton polygon in Python, mimicking the
program's reversal of the polynomial; `validate_arm.py` runs `ratpoints -v`
on every curve of the four test suites and compares the `denominators
excluded:` line with it.  Over the 1103 even-degree curves of the four
suites (the eight witnesses included), 119 of which stop before sieving --
no real points, or no points modulo some prime -- **587 reach the arm and
there are 0 disagreements.**

## What it touches beyond the arm

* `RATPOINTS_USE_JACOBI`, a new internal flag: the Jacobi test applies iff
  the degree is even, the leading coefficient is not a square and `-j` was
  not given.  Until now a square leading coefficient switched the whole
  checked denominator loop off; the valuation test does not care whether
  `c[d]` is a square, so the loop now stays available and `sieving_info`
  switches it off only when nothing came of it.  That is what reaches the
  point-rich suites: their curves have square leading coefficients by
  construction (points at infinity), and 16 of the 98 in `test1many` and
  7 of the 30 in `testhighmany` get a `v_3(b) = 1`-type exclusion.
* `forbidden_val { p, mask }` replaces the array of moduli `p^n`;
  `run_shape` takes the fraction `sum (p-1)/p^(m+1)` over the set bits off
  the denominator count.
* `valuation1()` used `abs()` on a `long`.
* `-v` lists what is excluded: `5|b 3^2|b v_7(b)=1 (lcf/b) = -1`.

## What the test suites reach

From the model, restricted to the 30 primes the program looks at by default:

| suite | curves | rule applies | of them square lcf | denominators saved |
|---|---|---|---|---|
| test1 (1000) | random | 563 | 3 | 3.6% via `p`, 2.7% via `p^n`, 0.6% via a single valuation |
| test1many (98) | point-rich | 34 | 16 | 0.5% via `p^n`, 1.5% via a single valuation |
| testhighmany (30) | point-rich | 22 | 7 | 1.3% via `p^n`, 0.5% via a single valuation |
| testdegrees (100) | degrees 3,4,7,8 | 30 | 0 | 1.1% + 1.4% + 0.2% |

"Rule applies" is the model's count over every even-degree curve; what the
program does is a little less, because some curves are settled before any
prime is examined (no real points, or none modulo some small prime).  On
the 1000 random curves the program excludes something on **499**, and the
model agrees on exactly those; on test1many 34, testhighmany 22, and 24 of
the degree suite.  The random suite gets there on half its curves because
its coefficients lie in `[-10, 10]`, where 3, 5 and 7 divide the leading
one often.  For large random coefficients Michael's estimate of 3.4% for the
`p` case stands, and the `p^n` case adds about as much again.

## Measured

Paired cycle counts (`perf stat -e cpu_core/cycles/u`, `taskset -c 0`),
three configurations back to back in a rotating order, one warm-up round
discarded, median of the per-round ratios.  "old" is the new binary with the
arm switched off by the temporary `RP_NO_PADIC` environment switch, "v23"
the `v2.3` binary built in a scratch worktree; `v23/old` measures how alike
the two builds are.

| suite | rounds | new/old | min | max | v23/old |
|---|---|---|---|---|---|
| test1 | 15 | **0.9169** | 0.9033 | 0.9382 | 0.9940 |
| test1many | 15 | **0.9763** | 0.9681 | 0.9851 | 0.9997 |
| testdegrees | 15 | **0.9475** | 0.9298 | 0.9748 | 0.9970 |
| testhigh | 7 | **0.9116** | 0.9088 | 0.9198 | 0.9984 |
| testhighmany | 7 | **0.9790** | 0.9546 | 0.9816 | 0.9962 |

The saving is larger than the share of denominators excluded (8.3% against
6.9% on test1, 5.3% against 2.7% on the degree suite).  That is not noise:
**a denominator divisible by `p` is dearer than average when `p` divides the
leading coefficient**, because `p` then passes every numerator for it --
`is_f_square[p]` is 1 for the class `b = 0 mod p` -- so it reaches the later
phases with several times the usual survivors.  Excluding exactly those
denominators saves more than their count.

## Correctness

* `make test1`, `test1many`, `testdegrees`, `test2`, `testhigh` and
  `testhighmany` are byte-identical to their references.
* The debug build under valgrind on three curves that reach the arm: clean.
* **Eight witness curves appended to `testdata.h`** (rows 1001-1008), two
  for each case that admits a point the rule must let through, found by
  random search under the case's congruence conditions and checked three
  ways: by a brute-force pass over every coprime `(a, b)` up to height
  16383 in 128-bit arithmetic (`scratchpad/brute.c`, independent of the
  sieve), against the new binary, and against the `v2.3` binary, which also
  confirms that none of them has a point of height between 16383 and
  200000, so `testbase` serves `make testhigh` as before.

  | rows | leading coefficient | what is excluded at 3 | the witness |
  |---|---|---|---|
  | 1001-1002 | `-72`, `-45` (`9u`, `-u = 1 mod 3` a square) | `v_3(b) = 1` only | `b = 54`, `b = 27`: a dropped sign would exclude `27 \| b` |
  | 1003-1004 | `45`, `99` (`9u`, `u` a non-square) | `v_3(b) = 1` and `27 \| b`; tie at `v = 2` | `b = 9` |
  | 1005-1006 | `9`, `36`, square, with a square constant term | `v_3(b) = 1`; the checked loop on a curve with points at infinity | `b = 18`, `b = 9` |
  | 1007-1008 | `-27`, `27` with `9 \| c[5]`, `c[4] = 1 mod 3` | `9 \| b`; at `v = 1` the `b^2` term wins with an even valuation and a square unit | `b = 3` |

  The fifth case, `c[4] = 2 mod 3` in the last row, excludes `v_3(b) = 1`
  through that unit and every larger valuation through the leading term, so
  it forbids 3 outright and has no witness by construction.

## The review

Five Opus agents read commit b77b342 through different lenses -- the
mathematics as implemented, C defects, control flow and flags,
documentation against code, tests and build -- and every medium or high
finding was then handed to two or three skeptics told to refute it.  Eleven
findings, none of them a defect in the code.  One was confirmed: the
documents said the test "fires on 563" of the random curves, which is the
model's count of curves the rule *applies* to, whereas the program excludes
something on 499 of them (the rest are settled before any prime is looked
at); one was refuted (the wording of the sign case in `testdata.h`, which
was consistent but read two ways, and is now spelled out); the rest were
low-severity polish, all taken: the README stated the rule without the
even-power condition, the `-F` paragraph read as scoped to a non-square
leading coefficient, the loop that enumerates the possible valuations in
`sieving_info` lacked the shift bound its twin has, two new `-Wextra`
sign-compare warnings, the search for further forbidden primes past
`num_primes` running to no purpose on a square leading coefficient, the
`-v` list not wrapped, and the Makefile's target summary omitting
`testdegrees`.  The one thing not taken: `make dist` and `make doc` still
ship and build only the 2.2 documentation, which belongs to the version
pass (TODO item 17).

## Not done, and why

* **The primes past `num_primes`.**  The loop that looks for more forbidden
  primes up to `sqrt(b_high)` tests `kronecker(c[d], p) == -1` only; for
  `p | c[d]` it could run the valuation test too.  Those primes are above
  127, so the saving would be under 1% with probability under `1/p^2`.
  (The review noticed that this loop now runs, and can never succeed, on a
  square leading coefficient; it is skipped there.)
* **Ties at `m = 1`.**  When `v_p(c[d]) = 1` and `p` does not divide
  `c[d-1]`, the two top terms tie at `v_p(b) = 1` and `F/p = a^(d-1)(r_0 a +
  r_1 b') mod p`.  That is a unit, and `F` not a square, unless `a = -r_1
  b'/r_0 mod p`: so for those denominators only one numerator class mod `p`
  can give a point.  The sieve cannot use it as it stands, because its
  tables are indexed by `b mod p` and this depends on `b mod p^2`: it would
  take `p-1` further tables for such a prime, indexed by `b/p mod p`, and a
  hook where the class-0 table is selected.  **Estimated 2026-09-12** (after
  the merge question came up): for those denominators, `(p-1)/p^2` of all,
  the prime passes `(p-1)/p` of the numerators now and would pass `1/p`,
  against `(p+1)/2p` for a typical class, so on a curve affected at `p` the
  refinement removes 11% (`p = 3`), 15% (`p = 5`), 14% (`p = 7`), 12%
  (`p = 11`) ... of the survivor-dependent work; a random curve is affected
  at `p` with probability `(p-1)^2/p^3`.  Summed over the phase-1 primes
  that is about 8% of the survivor-dependent work on random curves, which
  is **1.1% of a run at height 16383 and 2.1% at 200000** (phase 2 with the
  checks being 14.6% and 25.7% of those), and nothing on curves with a
  square leading coefficient (the analogous tie at `v_p(b) = 2` for
  `v_p(c[d]) = 2` is worth 0.3% of that work).  Run-time cost negligible;
  the cost is the code.  Not done on this branch.
* **`p = 2`.**  Handled by `get_2adic_info` modulo 16, which covers the
  analogous cases up to `v_2(b) = 3` exactly and treats `v_2(b) >= 4` as
  one class through `c[d] mod 16`; what it misses (`v_2(c[d]) >= 5` odd)
  has probability 1/64 and saves 1/64.
* **Odd degree** has its own valuation argument in `setup_us1`, which loops
  over all coefficients and is right.  Untouched.
