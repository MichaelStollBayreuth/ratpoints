# The tidy branch: step 1 of the tuning session

Three Left-overs of TODO.md that touch no estimate and no constant, taken
together as the first step of the tuning session (2026-09-17 evening), off
`v2.3` at 7f2239e.  `make test` only; nothing here changes what the default
build computes or how fast.

## 1. The unchunked arm's `&survivors[range - p]` (Left-over of 23)

`sift.c`, the `#else` arm for `RATPOINTS_CHUNK` outside 2..16.  The packet
loop ended on a pointer `surv_end = &survivors[range - p]`, which points
before the array whenever the range is shorter than the prime -- never
dereferenced, but undefined as pointer arithmetic, and a wart to read.  It
now counts: `left = range - r` bit arrays remain after the first `r`;
packets of `p` while `left >= p`; then the `left` that are left.  In the
`else` branch `range >= r` holds (the `if` took `range < r`), so `left >= 0`.

Checked with `CCFLAGS1='${CCFLAGS256} -DRATPOINTS_CHUNK=1'` (alt.sh's
`chunk1` variant): test1, test1many, testdegrees against their references;
`rptest -h 100/200/1000/4000` and the `timing` curve byte-identical to the
chunked default build.  Nothing builds this arm by default.

## 2. The Makefile's `SHELL` (Left-over of 18)

`SHELL = /bin/bash` was file-scope in v2.3; 2.2.4 on `main` scopes it to the
targets whose recipes use the shell's `time` (a built-in dash lacks).  Now
`test1 test1many testhigh testhighmany testdegrees test2 timing: SHELL =
/bin/bash`; every other recipe runs under whatever `/bin/sh` is.  The seven
are all the recipes that use `time`; no other recipe uses a bash-only
construct (grep for `[[`, brace expansion, `$((`, `<(`, `&>`, `set -o`).
`tune.sh` and `test3.sh` have their own `#!`.  Verified by `make test`:
bash's `real/user/sys` lines appear for each timed target.

## 3. The manual's "a fifth of the run is table construction" (Left-over of 18)

Measured on the merged tree 7f2239e with `-DRP_PHASE_TIMING -DRP_PHASE_COUNTS
-DRP_PRIME_STATS` (counters-tidy.txt; cycles by rdtsc, shares of `cyctot`):

| suite                      | tables | set-up/denominator | fill | phase 1 | phase 2 | phase 3 | phases 1+2 | outside sift() |
|----------------------------|-------:|-------------------:|-----:|--------:|--------:|--------:|-----------:|---------------:|
| test1 (16383, random)      |  4.3%  |  8.5%              | 2.5% | 38.8%   | 18.5%   | 1.3%    | 57.2%      | 16.2%          |
| test1many                  |  2.6%  |  5.4%              | 0.9% | 58.5%   | 21.1%   | 1.1%    | 79.6%      |  5.6%          |
| rptest -h 1000             | 18.0%  | 29.7%              | 1.8% |  6.3%   |  3.3%   | 1.2%    |  9.6%      | 48.2%          |
| testhigh (200000, random)  |  0.3%  |  1.3%              | 0.9% | 61.3%   | 28.6%   | 0.1%    | 89.9%      |  2.8%          |
| testhighmany               |  0.1%  |  0.5%              | 0.4% | 78.0%   | 17.4%   | 0.2%    | 95.4%      |  1.4%          |

"tables" is `cyctab` (the `sieve_init_*` calls), "set-up/denominator" is
`cycsetup`, which contains it (the rows for the phase-2 primes, the table
calls); both lie inside `cycsift`.  Callgrind on the plain build of the same
tree at 16383 (cg-tidy-rptest.*, instructions): `_ratpoints_sift0` 60% self,
`sift` 15% self, `find_points_work_1` 6%, `fill_bp_list` 5.5%, the
`sieve_init_*` family 4.3%, gmp (the exact checks) 2.5%, Sturm 2.2%,
`examine_prime` 1.2%.

So the tables are 4% of test1 now, not a fifth (they were 22% when the
high suites were added, before items 23 and 24 reworked the set-up), and
0.3% of testhigh, not 0.5%; the two sieving phases are 57% and 90% of the
run, not 48% and 87%.  What is still outside the two sieving phases at
16383 is two fifths of the run: the prime choice, the set-up per
denominator, the `b mod p` reductions, the exact checks, the third stage.
Five passages said the old numbers; all now say the new ones with the old
as history where the sentence is historical:

* `ratpoints-doc-2.3.tex`, the test-suite section (`22%`/`0.5%`/`48%`/`87%`);
* the `make tune` section ("a fifth of the time ... building sieve tables");
* the 2.3 narrative on the high suites ("a fifth of the run to be table
  construction"; kept as what was true when they were added, with the
  present figure);
* `Makefile`, the tunehigh comment and the testhigh comment.

The README has no such figure (its "fifth" sentences are gains of items 25
and 27, not shares).  The manual compiles with 0 errors and the same four
overfull boxes as before.

## What was not done

No measurement of the default build: none of the three changes is compiled
into it (1), affects a recipe's work (2) or the code (3).
