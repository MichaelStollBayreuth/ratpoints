# The test suite that exercises every branch (TODO item 20)

Branch `testsuite` off `v2.3` (a134a08), 2026-09-19.

## What was built

- **`test4.sh` / `testbase4` / `make test4`.**  About 430 invocations of
  `./ratpoints`, organised by what they exercise, each announced by a line
  `# <arg> <arg> ...` so that a failure locates itself in the diff and so
  that the checker knows what was run.  Part 1 (points and messages, the
  same in 2.2.4 and 2.3), part 2 (the `-v` reports, version-specific
  texts), part 3 (the options 2.2.4 does not have: `-r -R -U -C -A -P -Q
  -W`).  1.1 s on this laptop at `-O2`.  The `-v` reports are filtered by
  an awk function `v()` that drops what depends on the register width or
  the tuning (the counts of moduli chosen, the lists of moduli, the cost of
  an exact check); runs without `-q` go through `m()`, which drops the
  banner's first line and the report on the moduli used.
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
  Result on the final reference: 356 checked, 0 failed, 74 skipped (error
  messages, output formats, `-v` reports without points), 4 s.
- **`rpapi.c` / `testbase-api` / `make testapi`.**  The library interface
  the program cannot reach: `cof == NULL`, `height <= 0`, `domain ==
  NULL`, degree 0 and the genus drop (the error codes -2, -3, -1), the
  fields the library normalises (`num_inter < 0`, `b_low`/`b_high` out of
  range, `array_size`, `sturm`, `num_primes`, `sp1 > sp2`,
  `max_forbidden`), the input fields coming back as they went in, the
  flags through the API, a callback that declines or stops, two intervals,
  and the wrapper `find_points`.
- **`test4-configs.sh` / `make test4configs`.**  Builds the library in
  `build-test4-<name>/` from symbolic links to the sources, for the
  register widths 64, 128 (`USE_AVX128` and `USE_SSE`), 512 (`USE_AVX512`
  without `-mavx512f`, emulated), `RATPOINTS_CHUNK=1`,
  `USE_LONG_IN_PHASE_2`, `PRIME_SIZE=7`, `RATPOINTS_COMPOSITE_MAX=1` and
  `=1023`, and runs `test4.sh` against each; all nine compare equal to the
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

Three bugs, fixed on this branch (commit "what the suite found"):

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

Also learnt while looking for inputs, not bugs:

- The numerator strides above 8 cannot occur: the odd squares modulo 64
  are exactly the residues 1 mod 8, so whether `f(k)` is an odd square
  depends on `k` mod 8 only, and the even squares reduce to squares mod
  16 and mod 4 the same way.  `RP_NUM_STRIDES` is 7 (strides up to 64);
  the entries for 16, 32 and 64 are unreachable but harmless.
- A prime power can only be useless (`examine_power` returning 0) when
  no denominator divisible by its prime occurs: with `pinf` the classes
  `p | b` count at what their rows admit, which is less than everything.
  `-7x^4 + 18x^3 - 20x^2 - 6x + 19` is a square modulo 9 at every `x`
  with a leading coefficient that is not a square modulo 3, and does it.
- Two moduli sharing a prime are never both taken (`take_entries`), so
  the reuse of a prime power's sieve entry by a second modulus
  (`power_se[idx] != NULL` in `make_modulus`) cannot happen.
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
| find_points.c   | 99.94% of 1563   | 93.9% of 1270                |
| sift.c          | 95.72% of 187    | 90.9% of 186                 |
| init.c          | 100% of 60       | 98.1% of 52                  |
| sturm.c         | 100% of 144      | 99.3% of 142                 |
| main.c          | 98.73% of 473    | 97.3% of 300                 |

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
  `message(1)`.
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
it depends on the register width, the tuning or the constants (the
configuration matrix is the check of that claim).  A change to what `-v`
or the messages print changes part 2 and the `m()` runs; regenerate with
`./test4.sh > testbase4` and look at the diff.  A new branch in the code
wants an invocation that reaches it: `make coverage`, look at the
`#####` lines and the `taken 0%` branches in `build-coverage/*.c.gcov`,
add the invocation to the section it belongs to, regenerate, and run
`verify-test4.py testbase4` (a few seconds; the two large searches need
`gp`).

## The 2.2.4 adaptation

Branch `testsuite-224` off `main`: `test4.sh` without part 3 and with
`v()` dropping 2.2.4's lists ("use N primes for first/second stage:" and
the line after), `testbase4` from 2.2.4 (part 1 identical to 2.3's but
for the empty coefficient string, part 2 its own), `rpapi.c` without the
2.3 fields and `sp*_used`, `testbase-api`, the `test4`/`testapi`
targets, and the `read_input` fix.  Michael tags and pushes `main` as
v2.2.4 himself.
