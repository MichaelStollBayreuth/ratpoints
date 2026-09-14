# The sieve-side trio: sentinel scan, phase-2 row by the reciprocal, unswitched loops

Branch `sieve-trio` off `v2.3` at 281178e, 2026-09-14.  Items P7, P4 and
P8 of `review-2026-09-12/REVIEW-REPORT.md`, taken up as one group because
they touch the same two loops of `sift.c` and the build; TODO item 23.
Everything here was measured on the i7-1355U, pinned to P-core 0, in core
cycles (`perf stat -e cpu_core/cycles/u`), paired and interleaved, medians of
the per-round ratios; the outputs of every binary were compared with the
references before it was timed.

## What was done

### P7: the scan runs into a sentinel

After the first phase, `_ratpoints_sift0` steps over the bit arrays looking
for the few that are non-zero.  The loop was

    while(i < w_high && !TEST(*surv0)) { surv0++; i++; base++; }

nine instructions and two branches per bit array, with `w_high` reloaded
from the stack.  Now a non-zero bit array is written just past the range
(`*surv_end = ~zero`) and the loop is

    while(TESTZ(*surv0)) { surv0++; }

with the position recovered from the pointer once per survivor.  Under AVX
that is `add`, `vmovdqa`, `vptest`, `je`: four instructions, one branch.
`find_points_work` allocates `array_size + 2` bit arrays instead of `+ 1`:
one spare pays for the alignment, and since the range can be all of
`array_size`, the sentinel needs the other.

`TESTZ` is new in `rp-private.h`: "the bit array is all zero", defined for
AVX (and AVX-512) with the `testz` intrinsic and as `!TEST` everywhere else.
`TEST` itself is unchanged; it tests values that are already in registers,
where the two forms cost the same, and the review found the `vptest` form
to be a regression there when the short-circuit was still in the loop.

Two forms of the AVX test were tried.  `_mm256_testz_si256(ones, a)` should
take the bit array straight from memory (one `vptest` with a memory operand,
which is what the review's F1 verdict described), but gcc 14 does not fold
it: it loads the array, then rebuilds the all-ones constant *inside* the
loop with `vpcmpeqd`, five instructions per array.  `_mm256_testz_si256(a,
a)` gives the four above and is what the review's own prototype measured.

### P4: the second phase computes its row

The table for a prime has `p` bit arrays, and the one for word number `i`
is at index `(i + offset) mod p`.  The second phase used to reach it from a
pointer set up at the head of every call,

    ptr = ssp->start + base; while(ptr >= ssp->end) { ptr -= p; }

which is one to three data-dependent, mostly mispredicted branches per AND
(2.1 billion iterations on `make testhigh`), plus `sp2 - sp1` reductions at
the head of every call to set the `start` fields, whether or not a single
bit array survived.  Now it is

    AND(nums, ssp->ptr[small ? RP_MULMOD(i + ssp->offset, p, magic) : mod(i + ssp->offset, p)]);

with `RP_MULMOD` the multiply-high reduction of item 11 (exact below 2^32)
and `magic` from the per-curve `magics` array that already existed.

**No sign fix, no new array.**  The reduction needs a non-negative value,
and `i` can be negative.  The review's prototype added a per-prime bias from
a second array; here the bias is folded into the offset that the loop reads
anyway: `sift()` sets `ssp[n].offset` to the odd-numerator shift *plus*
`se->bias = p * ceil(2^31 / p)`, a multiple of `p` in `[2^31, 2^31 + p)`
(`RP_ROW_BIAS`, a new field of the sieve entry, computed once per prime and
curve beside `magic`).  Every reduction of `w + offset` is unchanged modulo
`p`, and one test at the head of `sift0`,

    small = (w_low >= -2^31 && w_high <= 2^31 - 2*RATPOINTS_MAX_PRIME)

says that every value the call reduces lies in `[0, 2^32)`.  That reaches a
height bound above 10^11; beyond it `mod()` divides, as the third stage does
beyond its own limit.  The first phase's start loop uses the same offsets
through `mod_mul`, whose sign branch is now never taken.

**What went with it.**  The start fields are set for the first-phase primes
only, in both arms of `sift0`; the unchunked arm's second start loop is gone;
`small` and `magics` moved to function scope (the unchunked arm had its own
copies); `base` is gone from the second phase; and the commented-out
`BASE_REPEAT` experiment of 2009, which was about exactly this loop, is
gone.  The `USE_LONG_IN_PHASE_2` arm, which used to walk `start` one word at
a time with its own subtraction loop (the review's prototype lost points
there because it left that loop reading a `start` nobody set), finds the row
the same way and takes word `k` of it.

**The fallbacks.**  With `-DRP_MULMOD_DIVIDE`, or on a compiler without
`__int128`, `RP_MULMOD` is a division, and the second phase divides per AND.
The review advised keeping the old pointer loop for that case; it is not
kept, because it would be a second copy of the loop for a compiler that
does not exist in practice (every 64-bit gcc, clang and icc has `__int128`),
and `RP_MULMOD_DIVIDE` is a development switch whose purpose is to measure
what the reciprocal is worth.  The `RP_MOD_CHOICE` runtime switch covers
the start loop only, as before.

### P8: `-funswitch-loops` -- measured and dropped

One flag in `CCFLAGS0`, tried as the third step.  It duplicates a loop whose
body contains a loop-invariant test, one copy per outcome: the "which
reduction" test in the start loop and now in the second-phase loop of
`sift.c`, and the `which_bits` tests in the per-denominator loops of
`find_points.c`.  `sift0` grows from 5450 to 6673 bytes.  It removes 1 to 3
per cent of the instructions and gains nothing in cycles at any of three code
alignments (table below), so it is not used, and the Makefile records that.

### A fourth stage in `make tune`

`ARCHITECTURE-TRIAGE.md` found the `COST_*` constants to be the
machine-bound part of the model, and `RATPOINTS_COST_TABLE` the one that
matters at small heights; `-C` existed as a runtime option but nothing
tuned it.  `tune.sh` now sweeps it after the run length (`C_VALUES`,
default `10 20 70 140` around the compiled-in 38; `C_FACTORS='0.5 2'` in
`make tunehigh`), passes the current value to every other stage, and writes
`-DRATPOINTS_COST_TABLE` into `tuning.mk`.  Thirteen settings a round in
`tunehigh` rather than ten, twenty-one in `tune` rather than sixteen.

## What it is worth

Each step against the one before, five rounds of the four suites, base =
281178e built from clean.  Cycles: median (min-max) of the per-round ratios;
instructions and branch misses: medians.

| step | test1 | test1many | testhigh | testhighmany |
|---|---|---|---|---|
| P7 sentinel, cycles | **0.985** (0.977-0.998) | **0.984** (0.974-0.990) | **0.965** (0.963-0.970) | **0.975** (0.971-0.985) |
| instructions / misses | 0.935 / 1.023 | 0.952 / 1.048 | 0.871 / 1.026 | 0.914 / 1.005 |
| P4 row, cycles | **0.985** (0.975-1.008) | **0.983** (0.963-0.984) | **0.959** (0.939-0.961) | **0.947** (0.937-0.965) |
| instructions / misses | 0.974 / 0.966 | 0.942 / 0.903 | 0.959 / 0.815 | 0.955 / 0.716 |
| P8 flag, cycles | 0.998 (0.976-1.015) | 1.001 (0.980-1.003) | 1.014 (1.010-1.019) | 1.020 (1.012-1.027) |
| instructions / misses | 0.975 / 1.037 | 0.969 / 0.996 | 0.989 / 1.044 | 0.990 / 0.998 |

The review predicted P7 at 1.5 / 2 / 4.0 / 2.0 per cent (its own paired
measurement, with the `vptest` form) and P4 at 2-3 / 4-6 / 4-6 / 4-6; both
came in as predicted, P4 a little better at the large bound.  The two
together are 3 per cent at 16383 and 7.5 per cent (random) to 7.7 per cent
(point-rich) at 200000, and the mechanisms are as the report said: the
sentinel removes instructions from the highest-IPC loop of the program (the
cycle gain is a quarter of the instruction gain, and the misses rise by the
one per call that now falls off the sentinel), the reciprocal removes
mispredicted branches (18 to 28 per cent of all misses at 200000, where the
fixup loop ran three times per AND).

P8 removed 1 to 3 per cent of the instructions and gained nothing: a wash
at 16383, and 1.4 and 2.0 per cent *slower* on the two large-height suites,
every round.  Whether that is the flag or the code placement it changes is
what the alignment check below is for.

### The alignment check, and the verdict on P8

Whole builds of different source differ in where gcc places the hot loops,
which is worth up to 10 per cent here (TODO item 2), so the pairs were
rebuilt at `-falign-loops=32` and `=64` and measured again, three rounds
each: base against v3 (everything), and v2 against v3, where v2 is the same
tree built with `-fno-unswitch-loops`, i.e. the flag alone.  Cycles B/A,
medians; test1 / test1many / testhigh / testhighmany.

| alignment | v3 against base | the flag alone (v3 against v2) | P7 + P4 (derived) |
|---|---|---|---|
| default | 0.978 / 0.968 / 0.937 / 0.944 | 0.998 / 1.001 / 1.014 / 1.020 | 0.980 / 0.967 / 0.924 / 0.925 |
| 32 | 0.971 / 0.982 / 0.942 / 0.958 | 1.007 / 1.002 / 1.006 / 1.018 | 0.964 / 0.980 / 0.937 / 0.941 |
| 64 | 0.988 / 0.984 / 0.957 / 0.952 | 1.005 / 0.997 / 1.002 / 1.001 | 0.983 / 0.987 / 0.955 / 0.951 |

(The stepwise measurement at the default alignment gives P7 + P4 = 0.970 /
0.967 / 0.925 / 0.923, in agreement with the derived row.)

**The flag never gains.**  At 16383 it is within a per cent of a wash at
every alignment; at 200000 it costs 0 to 2 per cent, with the instruction
count down 1 per cent and the branch misses up to 4 per cent up.  The review
predicted 1.5-4 and 0.5-1.6 per cent from the instruction counts and said
the cycles were not measured; now they are, and the instruction count was
the wrong guide, as it was for the group-of-four scan in PERFORMANCE-NOTES.
Dropped.  The five-line hand unswitch the review mentioned for non-gcc
builds is therefore not wanted either.

**What P7 and P4 are worth together** is then the last column: about 3 per
cent at 16383, and between 4.5 and 7.5 per cent at 200000 depending on the
code placement, 7.5 at the default one.  The spread at the large bound is
the alignment floor of about 2 points that `ratpoints-measuring` warns of,
not noise between rounds (the rounds agree to half a point).

## Checks

* Outputs byte-identical after each step on `test1`, `test1many` and
  `testdegrees`, compared before any timing (pair.sh refuses to time a
  binary whose output differs).
* The division path (`small == 0`): `ratpoints "a0^2 -2a0 1 0 0 0 1" 6e11
  -dl 1 -du 1 -l 559999999000 -u 560000001000` with `a0 = 5.6e11`, so that
  the word numbers exceed 2^31 and `f(a0) = a0^6` is a square: old and new
  binary both find `(a0 : a0^3 : 1)`, and both find nothing for `-dl 3 -du 3`.
* Every build variant reproduces the references on `test1`, `test1many`,
  `testdegrees` and `test3` (`alt.sh`, worktree at 1f00249): plain 64-bit,
  `USE_AVX128`, `USE_SSE`, `USE_AVX512` emulated (its usual two ABI
  warnings, no others), `USE_AVX` with `-mavx` only (the `vptest` TEST and
  TESTZ), `RATPOINTS_CHUNK=1` (the unchunked arm), `USE_LONG_IN_PHASE_2`,
  `RP_MULMOD_DIVIDE` (the second phase dividing), `RP_PHASE_TIMING` with
  `RP_PHASE_COUNTS` (identical on the three suites; `test3` differs there
  as it must, since the report goes to stderr), and the default.  On the
  default build also `test1once`, `test2`, `testhigh` and `testhighmany`,
  and valgrind clean on the debug build and on the optimised one (a random
  curve at 16383 and the sparse curve at 20000).

## `COST_PHASE2` remeasured

The review asked for it: the constant is what one AND on a surviving bit
array costs in units of one first-phase AND per word, and both ends of that
ratio moved.  Measured as in PARAM-NOTES (builds with `-DRP_PHASE_TIMING
-DRP_PHASE_COUNTS` at `RP_STOP_AFTER` 0, 2 and 3, the phase-2 parts by
differencing, pinned, old code and new back to back on the same evening;
the counters of the two trees agree to the last digit, so they did the same
work).  Old -> new:

| | test1 | test1many | testhigh | testhighmany |
|---|---|---|---|---|
| the unit: one phase-1 AND per word, cycles | 0.227 -> 0.208 | 0.335 -> 0.306 | 0.165 -> 0.158 | 0.241 -> 0.232 |
| one phase-2 AND, cycles | 17.9 -> 18.4 | 11.0 -> 9.3 | 29.5 -> 24.4 | 22.5 -> 14.0 |
| `COST_PHASE2` | 79 -> 88 | 33 -> 30 | 179 -> 155 | 93 -> 61 |
| the scan, cycles per bit array | 1.70 -> 1.74 | 2.34 -> 2.34 | 1.59 -> 1.45 | 1.91 -> 1.84 |
| `COST_TABLE` | 19 -> 21 | 18 -> 20 | 29 -> 30 | 28 -> 29 |
| `COST_BP` | 29 -> 31 | 9.7 -> 9.9 | 38 -> 38 | 14 -> 15 |

The unit fell by 4 to 9 per cent, because the reductions for the phase-2
primes at the head of every call were in the phase-1 timed region and are
gone; that is why every other constant rose by that much.  The AND itself
fell by 5 and 8.5 cycles at the large bound, one to three mispredicted
fixup iterations' worth, and by nothing measurable at the small one, where
the fixup ran 1.5 times per AND and the differencing of two 1.3-second
runs is not precise to a cycle.  (The instrumented build sees the direction
of P7 in the scan row but is not the measurement of it; the paired cycle
counts above are.  The per-survivor row of PARAM-NOTES is not repeated
here: it is the small difference of two large numbers from separate runs
and moved by 30 per cent between two measurements of the *same* code.)

**The constant stays at 110.**  It was compiled in as a value between the
two large-height measurements, 170 and 105 at the time, because that is the
regime in which the `-A 2` correction fires; it still lies between them
(155 and 61), and the basin of every constant in this model is flat (TODO
item 8).  Whether the ranking of the phase-2 primes would like it nearer 60
is a question a `-D` pair can answer without any alignment confound, since
only a constant changes; it is left for a tuning session.

## `make tune` with the fourth stage

Run on the final tree, idle machine, default three rounds (six minutes).
The fourth stage ran and the current settings were kept, which is the
designed outcome of a tuned build:

    stage 1 (threshold):   current +0.7%   0.003 +11.5%   0.005 +2.0%   0.012 -4.2%   0.02 -1.8%
    stage 2 (offset, r=0.012):  current -0.4%   4 +4.0%   6 +8.1%   9 +0.6%   13 +0.3%   18 +4.8%
    stage 3 (run length):  current -1.2%   3e5 +1.6%   6e5 +1.5%   2.4e6 -1.3%   5e6 -0.6%
    stage 4 (table cost):  current -1.6%   10 -0.2%   20 +0.2%   70 +2.0%   140 +2.5%
    Nothing beat the current settings (0.0075, 11, 1.6e6, 38.0) by the required 3%.

Two things came out of reading that report.  The `current` row is the
settings timed against themselves and lands within 1.6 per cent of 1 in
every stage, which is the noise the script tolerates (`NOISE` 0.02) and the
reason it demands 3 per cent before it writes anything: the -4.2 per cent of
the threshold 0.012 in stage 1 did not survive into stage 4, where every
candidate carried it.  And it could not have survived as such, because the
search had a gap that predates this branch: each stage carries the winners
of the stages before it into every candidate, but the constant it sweeps
appears only with its ladder values, never with its current one, so the
combination "earlier winners, this constant unchanged" is never timed and
can never win.  `tune.sh` now adds the current value to the candidates of
stages 2 to 4 when the ladder lacks it (one more setting per stage and
round: twenty-four a round in `tune`, fifteen in `tunehigh`), and was run
again to validate that.  The second run showed the new candidates in every
stage (`e=11`, `u=1.6e6`, `c=38.0`) and then declared the machine too noisy
-- the settings measured 6 per cent away from themselves in stage 1 -- and
wrote nothing.  That verdict was right: a stray test run of mine on another
core fell into that stage, and the package power budget did the rest.  So
the script's own guard was exercised as well.  A clean third run, and the
first ever `make tunehigh`, follow the reviews; their results are below.
TUNE3-PENDING

## How much of it is this machine (for item 22)

* **P7 sentinel**: general.  It removes a bound test and two induction
  variables from a loop; every core executes fewer instructions faster.
  The size of the gain is inversely related to the core's width, since the
  loop was already at nine instructions in 1.5 cycles here.  The `vptest`
  form of `TESTZ` is the AVX-specific half, worth at most 0.6 per cent on
  this P-core and up to 1 per cent on the E-core by the review's
  measurement; other widths get the sentinel with `!TEST`, which is most of
  the gain.
* **P4 row**: general.  It trades one to three mispredicted branches per
  AND for a multiply-high, which every 64-bit ISA has; the size scales with
  the misprediction penalty, 15-20 cycles on every current core.  The
  `magics`/`offset` reads are the same lines the loop touched before.
* **P8 flag**: compiler, not architecture -- and measured to be worth nothing
  in cycles at three alignments, so dropped.
* Nothing here depends on a cache size or a core type.

## Measurement

`pair.sh A B ROUNDS OUT suites` in the session scratchpad (a copy is in
`review-2026-09-12/`'s `bench3.sh` lineage): for each round and suite it
runs A then B or B then A, alternating, pinned to core 0 under `perf stat`,
and reports the median, minimum and maximum over rounds of the per-round
ratio of cycles, plus the medians of the instruction and branch-miss ratios.
The per-round spread on the two three-second suites is about 1.5 points
either way, on the two long ones about half a point; a 1.4 per cent
difference seen in every one of five rounds on `testhigh` is therefore real
for *these two binaries*, and the alignment check is what says whether it is
real for the two *programs*.
