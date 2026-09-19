# The test suite that exercises every branch (TODO item 20)

Branch `testsuite` off `v2.3` (a134a08), 2026-09-19.

## What was built

- **`test4.sh` / `testbase4` / `make test4`.**  452 invocations of
  `./ratpoints`, organised by what they exercise, each announced by a line
  `# <arg> <arg> ...` so that a failure locates itself in the diff and so
  that the checker knows what was run.  Part 1 (points and messages, the
  same in 2.2.4 and 2.3), part 2 (the `-v` reports, version-specific
  texts), part 3 (the options 2.2.4 does not have: `-r -R -U -C -A -P -Q
  -W`).  1.7 s on this laptop at `-O2`.  The `-v` reports are filtered by
  an awk function `v()` that drops what depends on the register width or
  the tuning (the counts of moduli chosen, the lists of moduli, the cost of
  an exact check); runs without `-q` go through `m()`, which drops the
  banner's first line and the report on the moduli used; runs that look at
  one line of a report go through `f()`, which announces the filter with
  the arguments, so that the checker sees them too (the first version
  piped `t` itself into grep, which ate the announce line of 18
  invocations and hid them from the checker -- the reviews caught it).
  The `-x` runs, which print the sieve's survivors, pin their primes with
  `-n 15 -N 30 -p 30`: what the sieve leaves depends on the moduli chosen
  and so on the tuning; nothing else in the reference does.
- **`verify-test4.py`.**  Reads `testbase4` and, for every invocation whose
  output is a list of points, searches all coprime `(a, b)` in the range by
  brute force in Python's arbitrary precision -- no sieve -- emulating the
  program's conversion of the search intervals to numerator ranges (double
  arithmetic, `ceil(b*low)` and `floor(b*up)`) so that the two agree at
  the ends.  `-x` and `-y` are checked as sets of `(a, b)`; `-1` as one
  point of the curve; `-z` by the count; a run refused as not squarefree
  is checked to be one (polynomial gcd over Q, plus the genus drop when
  leading zeros are stripped from an odd degree); the error tests are
  checked to have been refused.  Searches too large for brute force (the
  two curves at heights 30000-100000 in part 3, 19 invocations) go to
  PARI/GP's `hyperellratpoints`, another implementation of the search.
  An announce glued to the end of a format test's output (no final
  newline) is split off, so every invocation is seen.  Result on the
  final reference: 366 checked, 0 failed, 86 skipped (error messages,
  output formats, filtered reports), 6 s.  The skeptic review ran its own
  brute force on invocations of its choosing and 800 randomised runs
  against the program, and mutated the reference ten ways: every mutation
  was caught.
- **`rpapi.c` / `testbase-api` / `make testapi`.**  The library interface
  the program cannot reach: `cof == NULL`, `height <= 0`, `domain ==
  NULL`, degree 0 and -1 and the genus drop (the error codes -2, -3, -1), the
  fields the library normalises (`num_inter < 0`, `b_low`/`b_high` out of
  range, `array_size`, `sturm`, `num_primes`, `sp1 > sp2`,
  `max_forbidden`), the input fields coming back as they went in, the
  flags through the API, a callback that declines or stops, two intervals,
  and the wrapper `find_points`.
- **`test4-configs.sh` / `make test4configs`.**  Builds the library in
  `build-test4-<name>/` from symbolic links to the sources, for the
  register widths 64, 128 (`USE_AVX128` and `USE_SSE`), 512 (`USE_AVX512`
  without `-mavx512f`, emulated), `RATPOINTS_CHUNK=1`,
  `USE_LONG_IN_PHASE_2`, `PRIME_SIZE=7` and `=9`, `RATPOINTS_COMPOSITE_MAX=1`
  and `=1023`, and runs `test4.sh` against each; all ten compare equal to the
  one `testbase4` -- the `PRIME_SIZE=7` build without part 2, since two of
  the `-v` reports say what the primes beyond 127 say ("no points mod p =
  131", the third stage and the forbidden divisors of the square-everywhere
  curve), which a table of 30 primes cannot.  About 40 s.  The matrix
  caught one thing the first version of the `-v` filter had missed: the
  second line of the "bits set per word" message ("and N primes more in
  the third stage") depends on the width too, through the run shape.
- **`coverage.sh` / `make coverage`.**  Builds instrumented (`--coverage
  -O0`) in `build-coverage/`, runs test1, test1once (all three forms),
  test1many, testdegrees, test3, test4 and testapi through it, compares
  each with its reference, prints gcov's summary per source file and
  leaves the annotated sources there.  test2 and timing are left out: 60 s
  and 7 s at `-O0`, and they add nothing to the coverage.

## What the suite found

Three bugs, fixed on this branch (commit "what the suite found"), and two
more edge cases plus a missing check found by the two reviews of the branch
(commit "what the reviews turned up"):

1. **`sift.c`, `fill_checks`:** the third stage's bias was decided by
   `(double)p*(double)(2*height) < RP_STAGE3_LIMIT`; at a height bound of
   `2^62` the product `2*height` overflows a `long`, the comparison comes
   out true, `p*height` overflows too, and the stage rejects every
   survivor.  Found by the planted Pythagorean points at height `2^62`
   (denominators `2^32 - 1`, `2^32` and `2^62 - 2^31`), which vanished
   whenever the third stage had primes; heights up to `2^62 - 1` were
   fine.  Now `(double)p*2.0*(double)height`.  2.2.4 has no bias and no
   bug here.
2. **`find_points.c`, `find_points_work_1`:** `lcfsq =
   mpz_perfect_square_p(c[degree])` stood in the declarations, before the
   test `if(c == NULL) return(RATPOINTS_BAD_ARGS)`, so a `NULL` coefficient
   pointer crashed instead of being refused (found by `rpapi`).  It is
   now set after the checks and after the leading zeros are stripped
   (where it also used the given degree rather than the stripped one,
   harmlessly: a drop that keeps the genus goes from even to odd degree,
   and odd degrees do not read it).  2.2.4 has the same order.
3. **`main.c`, `read_input`:** `if(degree == 0) error(5)` after
   `degree--`, so an empty coefficient string (`ratpoints '' 10`) gave
   degree -1, reached the library and printed "Bug no. 1 - please
   report!" (exit 9) instead of "The polynomial must have degree at least
   1." (exit 5).  Now `degree <= 0`.  2.2.4 has the same bug; part 1 of
   the suite agrees between the versions on everything else.
4. **`find_points.c`, `sift()`:** the bit interval was made half-open by
   `high++` and its bit arrays counted as `(high + RBA_LENGTH - 1) >>
   RBA_SHIFT`; both overflow a `long` when the last bit is within
   `RBA_LENGTH` of `LONG_MAX`, which a height bound within 256 (64 with
   64-bit registers) of `LONG_MAX` reaches with a search interval that
   runs to the top -- the interval was then dropped without a word (the
   skeptic review bisected the first bad height to exactly `LONG_MAX -
   RBA_LENGTH + 1`).  Now the interval stays closed, `w_high = (high >>
   RBA_SHIFT) + 1`, and the mask of the last bit array is `RBA_LENGTH - 1 -
   (high & (RBA_LENGTH - 1))`.  The same in 2.2.4.
5. **`find_points.c`, the denominator loops:** `for(b = b_low; b <=
   b_high; b++)` never ends when `b_high` is `LONG_MAX` (the increment
   wraps), and `bb = b*b` and `bb = d*b*b` in the loops over the squares
   overflow for denominators within `2*10^9` of it (the sweep review).
   The plain loop (and in 2.2.4 the checked one, which walks denominators
   one at a time too) breaks at `LONG_MAX`; the square loops test `b <=
   b_high/b` and `b <= (b_high/d)/b`, which is the same inequality without
   the product.  The first two have tests at `LONG_MAX`; the square loops
   cannot be tested there in a fast suite, since they start at `b = 1`
   whatever `b_low` is (a pre-existing limitation: with `-dl` near
   `LONG_MAX` the loop over the squares runs `3*10^9` empty iterations),
   so those two conditions rest on the reasoning alone.
6. **`find_points.c`, `find_points_work`:** `long c_long[degree+1]` is
   sized before any check, so a negative degree from a library caller was
   undefined behaviour (it did not crash here, but it is the shape of bug
   2).  `find_points_work` now refuses a negative degree with
   `RATPOINTS_BAD_ARGS` before calling `find_points_work_1`; `rpapi`
   tests it.  Not done in 2.2.4, whose `find_points_work` has no wrapper
   to put the check in.

Also learnt while looking for inputs, not bugs:

- The numerator strides above 8 cannot occur.  For the odd denominators
  the set of admissible `k` is what `f(k)` makes a square mod 64; the odd
  squares mod 64 are exactly the residues 1 mod 8, so if any admissible
  `k` has `f(k)` odd, its whole class mod 8 is admissible and the stride is
  at most 8; if every admissible `k` has `f(k)` even, `f(k)` is in {0, 4,
  16, 36} and the Taylor expansion `f(k + t) = f(k) + t f'(k) + t^2 f_2(k)
  mod 64` shows that `k + 16` or `k + 32` is admissible too (the skeptic
  review supplied this half of the argument and confirmed the conclusion
  by exhaustion over every polynomial mod 64 of degree up to 3 and 70
  million samples of higher degree).  For an even denominator `b` with
  `v_2(b) = j` the map `a -> b a^-1` has fibres that are classes mod
  `2^(6-j)`, so `j >= 3` gives stride at most 8 at once, and for `j = 1, 2`
  expanding `g(t + 16), g(t + 32), g(t + 48)` shows that `gsq` never meets
  a class of `t` in one point only.  `RP_NUM_STRIDES` is 7 (strides up to
  64); the entries for 16, 32 and 64 are unreachable but harmless.
- A prime power can only be useless (`examine_power` returning 0) when
  no denominator divisible by its prime occurs: with `pinf` the classes
  `p | b` count at what their rows admit, which is less than everything.
  `-7x^4 + 18x^3 - 20x^2 - 6x + 19` is a square modulo 9 at every `x`
  with a leading coefficient that is not a square modulo 3, and does it.
- Two moduli sharing a prime are never both taken (`take_entries`), so
  the reuse of a prime power's sieve entry by a second modulus
  (`power_se[idx] != NULL` in `make_modulus`) cannot happen -- as long
  as the mask of the primes a modulus involves has a bit for every prime,
  which it has up to `PRIME_SIZE` 8 (53 primes in an `unsigned long`).
  With `PRIME_SIZE` 9 (96 primes) the mask is 0 for the primes past the
  64th, the exclusion silently stops working for them, and that arm
  becomes reachable; the points are unaffected (a redundant modulus
  sieves correctly), the matrix runs that build and it agrees with the
  reference, but it is a thing to know if the table is ever enlarged.
- `adapt_primes` in mode 2 returns before measuring when `sp2 == sp1`
  (there is no second phase to measure, and `r1 - r2` would be zero), so
  `-A 2 -R 0` never adds a prime; `-R 1` with `-r 0.5 -C 0` does.
- In `jacobi1` the reduction `f %= b` after dividing the gcd out is
  commented out (as in 2.2.4), which makes the `if(f == 0) return(1)`
  after it dead code; the multi-precision `jacobi` has the reduction and
  both outcomes are reached (`-F 0` on the 17-prime leading coefficient,
  denominator 45).

## Coverage

`make coverage` after the suite, gcov 14, `-O0`; before it the existing
tests (test1, test1once, test1many, testdegrees, test2, test3, timing)
stood at 92.9 / 90.9 / 100 / 100 / 57.3 per cent of the lines.

| file            | lines executed   | branches taken at least once |
|-----------------|------------------|------------------------------|
| find_points.c   | 99.94% of 1566   | 95.2% of 1274                |
| sift.c          | 95.72% of 187    | 90.9% of 186                 |
| init.c          | 100% of 60       | 98.1% of 52                  |
| sturm.c         | 100% of 144      | 99.3% of 142                 |
| main.c          | 98.73% of 473    | 97.7% of 300                 |

What is left, and why it is left:

- **Guards against states no caller can produce** (dead by construction,
  kept as guards): `ensure_ba_buffer`'s second condition; the zero
  arguments of `valuation`/`valuation1` (a zero coefficient is skipped
  before); `mpz_fits_slong_p` in `setup_us1` (a remaining factor below
  `1021^2` fits); `compare_entries` on equal keys; `prime_cost` with
  `tabled == 0`; `info <= 0` in the keys; `root < 1` in `check_cost`;
  the floors and zero checks of `run_shape`, `mean_bits_per_word`,
  `phase_2_offset`, `forbidden_fraction`, `modinv`; `check_rel <= 0`;
  `want < 0` for the second phase; `sp3_max < sp3`; the malloc-failure
  arm for the patterns of forbidden divisors beyond the table; the
  `m < LONG_LENGTH - 1` bounds of the valuation loops (`p^m` exceeds the
  height first); `den_bits == 0` with a point at infinity (a square
  leading coefficient always leaves the class `b = 0 mod 64`
  admissible); in `main.c` the return of `read_input` (always 0) and its
  `message(1)`, and `error(9)` with its text "Bug no. 1 - please report!",
  which fix 3 made unreachable from the program (the library's error codes
  other than "not squarefree" no longer have a path from `main`); the
  `a < 0` arm of `mod_mul` in `sift.c` (the row offsets carry
  `RP_ROW_BIAS`, so its argument is never negative); the `q + i < words`
  bound of `lay_out_pattern` in `init.c`.
- **Other compile-time settings**: `pn >= LONG_LENGTH` in the prime masks
  (`PRIME_SIZE >= 9`); the end of the composite loop at
  `RATPOINTS_MAX_PRIME_EVEN`, three prime factors and the candidate
  bound (`RATPOINTS_COMPOSITE_MAX >= 105`, run by `test4configs` as
  `comp1023`); `sift.c`'s subtraction ladder in `mod()` (values below
  `16p`, which the `!small` path that calls it never has).
- **Impossible by arithmetic**: strides above 8 (above); the reuse of a
  power's entry (above); `f == 0` after the division in `jacobi1`
  (above); `relprime` with both arguments even (even denominators take
  odd numerators); `sl == 0` at an interval end in the Sturm bisection
  (the signs passed down are never zero).
- **Not reached, reachable in principle**: `adapt_primes` returning for
  too few survivors (`n_bits < 200` with `n_arrays >= 1000`, a curve with
  very few survivors of the second phase at a million words); a composite
  modulus at the top of the second phase when mode 2 walks it down (the
  second-phase ranking puts the primes there).

## How to keep it up to date

A change to the sieve that is right leaves `testbase4` alone: nothing in
it depends on the register width, the tuning or the constants, except
that the `-x` runs print the sieve's survivors and so pin their primes
(the configuration matrix is the check of the first claim; a retuning
does not touch the reference).  A change to what `-v`
or the messages print changes part 2 and the `m()` runs; regenerate with
`./test4.sh > testbase4` and look at the diff.  A new branch in the code
wants an invocation that reaches it: `make coverage`, look at the
`#####` lines and the `taken 0%` branches in `build-coverage/*.c.gcov`,
add the invocation to the section it belongs to, regenerate, and run
`verify-test4.py testbase4` (a few seconds; the two large searches need
`gp`).

## The two reviews

Two Opus agents, briefs in `review-2026-09-12/rev-common/sweep-testsuite.txt`
and `skeptic-testsuite.txt`.  The sweep reproduced everything, confirmed
the three fixes and the mathematics of every comment and of the checker
(its `squarefree()` against sympy on 4000 polynomials, its homogenisation
and interval emulation against the C line by line, five mutations of the
reference caught), and listed 20 items: a wrong forbidden-valuation
comment (`3x^2 + 2x + 1` excludes `3^2 | b`, not every valuation; the
bit-array arm is reached by `3x^2 + 4`), a negative-definite curve filed
under the 2-adic classes (replaced by `-53x^2 + 28x + 4`, one numerator in
four for odd denominators), the 18 piped invocations invisible to the
checker, the announces glued to format output, the VLA before the checks
(fix 6), the latent overflows near `LONG_MAX` in the loops (fix 5), and
wording (degrees "1 to 9 and up to 100", "odd primes below 1024", 19
invocations checked by gp, the coverage figure).  The skeptic held every
claim it attacked, with caveats that are now fixed or written down: it
found the `w_high` overflow (fix 4) and bisected its threshold, showed
that the `-x` runs depended on the tuning (pinned), that the strides
argument was half an argument (completed above), and that `PRIME_SIZE` 9
is a cheap tenth configuration (added).  Not taken up: making `make
test4` fail the build (the house convention is `echo "Test failed!"`),
and a test of `RATPOINTS_VERBOSE` through the API (its report depends on
the register width, and `testbase-api` must not).

## The 2.2.4 adaptation

Branch `testsuite-224` off `main`: `test4.sh` without part 3, without
`-P` in part 2, and with the same filters (2.2.4's lists "use N primes for
first/second stage:" match the same pattern), `testbase4` from 2.2.4
(part 1 identical to 2.3's line for line, part 2 its own), `rpapi.c`
without the 2.3 fields, the `sp*_used` print and the negative-degree
check, `testbase-api`, `verify-test4.py`, `test4-configs.sh` with the
switches 2.2.4 has (widths, chunk, phase-2 mode, `PRIME_SIZE` 8),
`coverage.sh` with its tests, the four targets, and the fixes 2, 3, 4
and 5.  The files are derived from the 2.3 ones by a script
(`adapt224.py`, session-only).  Michael tags and pushes `main` as v2.2.4
himself.
