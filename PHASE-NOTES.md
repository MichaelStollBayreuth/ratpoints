# The two sieving phases across register widths

Measurements of 2026-09-07 with the `-DRP_PHASE_TIMING` instrumentation in
`sift.c` on this branch.  The question is the one on the TODO list: does the
second phase of the sieve benefit from wide registers the way the first one
does, or would a narrower phase 2 be better?

**This supersedes the first version of these notes.**  That one was taken with
`sp1 = 11`, `sp2 = 19` -- the fixed defaults of 2.2.3 -- and most of what it
concluded was an artefact of those values rather than a property of the
register width.  Everything below uses the automatic choice, so each curve is
sieved at the survivor rate the parameters aim at.

**Short answer: phase 2 is never faster on narrow registers.**  Where it
depends on the width at all it prefers wide registers, and where it does not,
it is flat to within 4%.  The hunch on the TODO list was wrong, and so was the
first measurement, which appeared to support half of it.

## How it was measured

Core cycles from `perf stat -e cpu_core/cycles/u`, split into the three stages
by the `rdtsc` ratios of the same run, pinned to CPU 0 with `taskset`.  Nine
builds and two curves are rotated within each round and the minimum of five
rounds is kept; the observed spread between the fastest and slowest round is
in the last column, and is 1-3% except for one entry.  Where a difference
mattered it was re-measured as an explicit ratio, the two builds run back to
back, seven times, median kept -- this machine slows by up to a quarter as it
warms up, and unpaired numbers taken minutes apart cannot be compared.

The exact `gmp` check runs *inside* phase 2 and is timed separately, because
it is the same work for every build.  "setup" is the rest of the run: the
Sturm bounds, the 2-adic information and the sieve tables.

The two curves are deliberately at opposite ends:

| | `1 0 126 0 441`, height 400000 | Drew Sutherland's record curve, height 200000 |
|---|---|---|
| chosen automatically | `sp1 = 12`, `sp2 = 17` | `sp1 = 23`, `sp2 = 23` |
| surviving phase 1, per 64-bit word | 0.54% | 1.75% |
| bits set per surviving word | 1.003 | 1.008 |
| exact `gmp` check | 1% of the run | 28% of the run |

Two things are already visible here.  The automatic rule picks **the same
`sp1` and `sp2` at every width** -- by construction, since the threshold is
per 64-bit word and the mean number of bits set per word does not depend on
the register width -- so this comparison is between builds doing identical
arithmetic, which the earlier one was not.  And a surviving word essentially
never carries a second bit (1.003 bits per surviving word), which is what
makes the model below work.

The second curve is prime-starved rather than dense: see the last section.

## The sparse regime (a typical curve)

`1 0 126 0 441` at height 400000, `sp1 = 12`, `sp2 = 17`.  Core cycles in
units of 10^9:

| build | total | phase 1 | phase 2 | gmp | setup | 1/2 split | vs 256 | spread |
|---|---|---|---|---|---|---|---|---|
| 64, `CHUNK=16` | 21.31 | 13.49 | 5.40 | 0.13 | 2.30 | 71/29 | +74.1% | 0.9% |
| 128 SSE | 16.72 | 10.04 | 5.43 | 0.12 | 1.12 | 65/35 | +36.6% | 1.4% |
| 128 SSE `+long` | 17.29 | 10.10 | 5.92 | 0.12 | 1.15 | 63/37 | +41.2% | 0.8% |
| 128 AVX128 | 15.15 | 10.08 | 3.81 | 0.12 | 1.13 | 73/27 | +23.7% | 0.5% |
| 128 AVX128 `+long` | 16.40 | 10.04 | 5.09 | 0.12 | 1.14 | 66/34 | +34.0% | 1.4% |
| **256 AVX2** | **12.24** | 8.33 | 3.06 | 0.12 | 0.73 | 73/27 | -- | 3.4% |
| 256 AVX2 `+long` | 14.11 | 8.42 | 4.85 | 0.13 | 0.72 | 63/37 | +15.3% | 1.7% |
| 512 emulated | 31.97 | 26.54 | 4.32 | 0.13 | 0.98 | 86/14 | +161.2% | 1.4% |
| 512 emulated `+long` | 32.84 | 26.29 | 5.45 | 0.13 | 0.98 | 83/17 | +168.2% | 2.1% |

`+long` is `-DUSE_LONG_IN_PHASE_2`; "512 emulated" is `-DUSE_AVX512` *without*
`-mavx512f`, so gcc lowers the 64-byte vectors -- it says nothing about real
AVX-512 hardware and is a structural check only.

**Phase 2 now scales with the width**, 5.40 -> 3.81 -> 3.06 from 64 to 128
(AVX128) to 256.  The first version of these notes found it width-independent;
that was because the phase-2 loop then re-entered its body for every empty
bit-array, and the tight skip loop added since has made the scan cheap enough
that its cost, which is proportional to `N/W`, is visible.

Two constants describe phase 2 across the whole table.  Writing `A` for the
bit-arrays swept and `S` for the arrays surviving phase 1,

    phase 2  =  1.17 * A  +  173 * S    cycles

fitted on the 128 AVX128 and 256 rows, predicts the 64-bit one to within 2.3%
(5.28 against a measured 5.40).  `A` halves with every doubling of the width;
`S` does not move at all (13.56, 13.51, 13.42, 13.24 million at 64, 128, 256,
512), because a survivor almost never has a companion in the same array.  So
widening the registers divides the scan and leaves the per-survivor work
alone, and at this rate the per-survivor term is 76% of phase 2 at 256 bits
(45% at 64, 61% at 128).  That is the ceiling on what any wider register can buy here.

**Phase 1 still does not scale with the width**: 13.49 -> 10.08 -> 8.33 is a
factor of 1.34 and then 1.21 per doubling.  The reason is unchanged and
structural, in `init.c`: each sieve table is stored with `RBA_PACK` copies of
the pattern, so phase 1 reads exactly `sp1` bits of table per numerator
whatever the width.  The bytes moved are width-independent; only the
instruction count falls.  That is the memory-bandwidth limit the documentation
warns about, made concrete.

## The point-rich regime

Drew's curve at height 200000, `sp1 = sp2 = 23`:

| build | total | phase 1 | phase 2 | gmp | setup | 1/2 split | vs 256 | spread |
|---|---|---|---|---|---|---|---|---|
| 64, `CHUNK=16` | 11.60 | 5.27 | 3.42 | 2.43 | 0.48 | 61/39 | +32.1% | 0.7% |
| 128 SSE | 9.29 | 3.09 | 3.42 | 2.44 | 0.33 | 47/53 | +5.8% | 2.3% |
| 128 SSE `+long` | 9.05 | 2.94 | 3.39 | 2.40 | 0.33 | 46/54 | +3.1% | 1.1% |
| 128 AVX128 | 9.13 | 2.91 | 3.44 | 2.43 | 0.34 | 46/54 | +4.0% | 1.9% |
| 128 AVX128 `+long` | 9.13 | 2.93 | 3.44 | 2.42 | 0.33 | 46/54 | +4.0% | 1.1% |
| **256 AVX2** | **8.78** | 2.60 | 3.42 | 2.48 | 0.29 | 43/57 | -- | 1.6% |
| 256 AVX2 `+long` | 8.54 | 2.59 | 3.26 | 2.41 | 0.28 | 44/56 | -2.7% | 1.5% |
| 512 emulated | 13.29 | 6.98 | 3.56 | 2.42 | 0.33 | 66/34 | +51.4% | 3.1% |
| 512 emulated `+long` | 12.92 | 6.87 | 3.32 | 2.40 | 0.33 | 67/33 | +47.2% | 9.0% |

**Phase 2 is flat**: 3.26 to 3.56 across every build, a 4% band, against a
factor of 2.7 in phase 1.  Here `sp2 = sp1`, so phase 2 does no `AND` steps at
all (the counter reads exactly zero) and consists of the scan -- which is small
because there are few numerators -- plus per-survivor work that does not depend
on the width.  The first version of these notes had phase 2 *scaling* with the
width in this regime, 64 bits being 28% worse than 256; that was entirely the
old `sp1 = 11`, at which 82% of all 256-bit arrays reached phase 2 and phase 2
really was doing wide `AND`s on nearly every array.  At the parameters the
program now chooses, that regime does not occur.

The `gmp` check is 2.4e9 cycles in every build, 28% of the run.  On this curve
it, not the sieve, is what a further optimization would have to attack.

## What follows for the build

* **Do not enable `USE_LONG_IN_PHASE_2` at any width.**  It costs 3.4% at
  128 SSE, 8.3% at 128 AVX128 and 15.3% at 256 on the sparse curve, and is a
  wash to -2.7% on the point-rich one.  The first version of these notes
  recommended enabling it at 128 bits, on a measurement that predates the
  skip loop; that recommendation is withdrawn.  The Makefile's "this is
  usually slower" is right at every width.
* **`CCFLAGS128` should be `-DUSE_AVX128`, not `-DUSE_SSE`.**  It is 9.4% on
  the total and 30% on phase 2 of the sparse curve (5.43 -> 3.81), and a wash
  on the point-rich one.  `USE_SSE` still uses the generic fall-back
  `TEST(a) = (EXT0(a) || EXT(a,1))`, which extracts both halves into general
  registers, and `__builtin_ia32_andps`, a floating-point `AND` on integer
  data; `USE_AVX128` has the `_mm_movemask_epi8(_mm_cmpeq_epi8(...))` form
  that was introduced for AVX-512 in 2.2.3.  Despite the name it needs only
  SSE2 intrinsics and no `-m` flag beyond the x86-64 baseline, so it is as
  portable as the `USE_SSE` build it would replace.
* **`-DRATPOINTS_CHUNK` matters enormously for the 64-bit build**, which
  defaults to `CHUNK=1`.  The 64-bit rows above are all at `CHUNK=16`; at the
  default the build is far slower, and the 64-bit column of any width
  comparison is meaningless unless the chunking is matched.

## Answers to the TODO item

* *"Phase 2 may well be fastest at 64 bits."*  No.  It is fastest at 256 in
  the sparse regime and width-independent in the point-rich one.  There is no
  regime in which narrow registers win phase 2.
* *"A wider `nums` takes more primes to clear, so the early exit fires later."*
  Not measurable, because at the rate the parameters now aim at a surviving
  array holds 1.003 bits on average at every width.  The effect needs several
  bits per array, which was the old `sp1 = 11` and is not the current regime.
* *"The optimal `sp1` should rise with the register width."*  It should not,
  and with the per-word rule it does not: the choice is width-independent by
  construction and nothing here argues against that.
* *"A mixed build -- phase 1 wide, phase 2 on 64-bit words."*  That is exactly
  `USE_LONG_IN_PHASE_2`, and it loses at every width.

## What the measurement turned up instead

**The prime table, not the register width, is what limits point-rich curves.**
Drew's curve gets `sp1 = sp2 = 23` not because 23 primes reach the target rate
but because 23 is every prime below 128 that carries any information at all:
on a curve with that many points `f` is a square modulo every residue for the
smallest primes, and `sieving_info` discards those.  The rule stops at 1.75%
survivors per word against a target of 0.75%.

Raising the table to `PRIME_SIZE=8` (primes to 251) and allowing 40 of them
gives `sp1 = 19`, `sp2 = 24`, 12.7 times fewer calls to the exact check, and

    ratio to the default build: 0.448 0.430 0.441 0.430 0.510 0.437 0.454

seven paired runs -- **a factor of 2.3**.  Larger primes are what a point-rich
curve needs, because `np/p` tends to 1/2 as `p` grows however many points the
curve has, while for small `p` it is close to 1.

It does not follow that the table should simply be enlarged.  Of the 98 curves
in `testdata-many.h`, 19 are starved this way and 79 are not, and the two
groups move in opposite directions (height 100000, median of three paired
runs):

| | starved (6 curves) | not starved (6 curves) |
|---|---|---|
| ratio, `PRIME_SIZE=8 -p 40` vs default | 0.64-0.73 | 1.03-1.49 |

and over the whole suites the two cancel: `test1many` 0.97, `test1` 1.22 --
random curves, which are never starved, only lose.  The reason is that the
rule ranks primes by information alone: given more to choose from it takes
larger ones, whose tables are bigger and colder, and it pays for that whenever
it did not need them.

The height bound is the second half of it.  A prime `p` needs `p` sieve tables
of `p` bits, so the tables grow quadratically in `p` and are paid for once per
denominator whatever the height bound, while the sieving they serve grows with
its square.  Ratios of `PRIME_SIZE=8 -p 40` to the default build, median of
three paired runs:

| curve | h=5000 | h=20000 | h=80000 | h=320000 |
|---|---|---|---|---|
| `1 0 126 0 441` (random) | 1.55 | 1.08 | 1.12 | 1.13 |
| point-rich, starved | 1.17 | 0.91 | 0.67 | 0.56 |
| point-rich, not starved | 1.65 | 1.51 | 1.51 | 1.22 |

The starved curve crosses over at about `h = 15000` and is nearly twice as fast
by `h = 320000`; the other two never cross over, they only become less bad.
The suites are all at small heights, which is the other reason the effect
cancels there.

Both criteria are sharp and available before any sieving starts:
`sieving_info` knows whether the loop over `prec[]` ends without the rate
falling below the target -- visible as `sp1 == sp2` in the verbose output --
and `args->height` is known.  Choosing the range from the two is TODO item 6.

## Does a larger prime table cost anything by itself?

Asked in preparation for TODO item 6: if the rule is to be allowed to reach
for larger primes when it needs them, the tables for those primes have to be
in the build, and the question is what merely having them costs.

It is a clean experiment, because `-p 30` makes every build use exactly the
same thirty primes (3 to 127; `prime[]` in `primes.h` is one fixed list and
`RATPOINTS_NUM_PRIMES` only says how far into it a build may look).  `sift.c`
does not mention the prime size at all, so the sieve is textually identical;
`find_points.c` mentions it twice, both in `sieving_info`.  Verified: the same
`sp1` and `sp2` are chosen at every prime size, and all four builds reproduce
`testbase` and `testbase-many`.  Only the size and stride of the tables differ.

Ratios to `PRIME_SIZE=7`, each candidate run immediately after a
`PRIME_SIZE=7` run, median of seven pairs for the suites and five for the
curves:

| | `PRIME_SIZE=8` | 9 | 10 |
|---|---|---|---|
| `make test1` | 1.006 | 1.014 | 1.004 |
| `make test1many` | 1.020 | 1.017 | 0.995 |
| `1 0 126 0 441`, h=400000 | 1.025 | 1.007 | 0.992 |
| record curve, h=200000 | 0.965 | 0.967 | 0.966 |

**There is no cost.**  On the suites and on the sparse curve everything is
within 2.5% of 1 and not monotone in the prime size, which is what no effect
looks like at this precision.

The point-rich curve is 3.4% *faster* at every larger prime size, which is
consistent across all five pairs and spread evenly over all four components of
the run -- including the exact `gmp` check, which the prime tables cannot
touch.  Two controls: `PRIME_SIZE=7` paired against itself gives 0.997, so it
is not an artefact of always running the baseline first; and rebuilding
`PRIME_SIZE=7` with `-falign-functions=32` and `=64`, which moves the code
without changing it, gives 0.989 and 1.007 on the total and leaves the `gmp`
component at 1.002, so ordinary code-layout luck is about half the size and
does not touch the part that moved most.  The likely explanation is the heap:
`find_points_init` allocates tens of megabytes more, and the `mpz` temporaries
that the check works in are placed around it.  Recorded rather than explained.

What a larger table does cost is address space.  `find_points_init` sizes
`ba_buffer` from **all** `RATPOINTS_NUM_PRIMES`, not from `args->num_primes`,
so `-p 30` does not reduce it:

| `PRIME_SIZE` | primes | largest | `ba_buffer` malloc (256-bit) | peak RSS |
|---|---|---|---|---|
| 7 | 30 | 127 | 5.0 MB | 7.2 MB |
| 8 | 53 | 251 | 33.2 MB | 9.4 MB |
| 9 | 96 | 509 | 240.6 MB | 10.1 MB |
| 10 | 171 | 1021 | 1668.5 MB | 8.8 MB |

The resident set barely moves, because the pages belonging to primes that are
never used are never touched, and the figures double at 512-bit registers.
Still, 1.7 GB of address space per `ratpoints_args` is not something to ship,
and it is per *thread* under TODO item 1.

**Fixed since** (`v2.3`): `find_points_init` reserves for
`RATPOINTS_DEFAULT_NUM_PRIMES` and `find_points_work` enlarges the block when
a call asks for more.  It cannot read `args->num_primes` at init, because the
documented use of the library sets that field between the two.  The last
column of the table above becomes 8, 9, 11 and 18 MB, and `-p 40` at
`PRIME_SIZE=10` costs 26 MB.

## How large should the prime table be?

The previous section asked what a larger table costs when its primes are not
used.  This one asks what it is worth when they are, and therefore what
`PRIME_SIZE` should default to.

**Raising `PRIME_SIZE` alone does nothing at all.**
`RATPOINTS_DEFAULT_NUM_PRIMES` is a fixed 30 whatever the table holds, so
without `-p` the rule never looks past the thirtieth prime: at 7, 8, 9 and 10
the same `sp1` and `sp2` are chosen and the cycle counts agree to 0.5%.  The
two constants are coupled, and only the pair means anything.

**The benefit is captured by about forty primes.**  For the five most starved
curves, the primes the rule asks for as it is allowed to see more (`sp1/sp2`;
`sp1 == sp2` is the signature of running out):

| | `-p 30` | 35 | 40 | 45 | 53 | 70 | 96 | 171 |
|---|---|---|---|---|---|---|---|---|
| record curve | 23/23 | 21/26 | 19/24 | 18/23 | 17/22 | 15/20 | 14/19 | 13/18 |
| four others | 22-26 starved | all unstarved | 17-18/22-23 | | | | | 13/18 |

By `-p 35` -- primes to 151 -- none of them is starved any more.  Everything
beyond that only lets the rule swap in larger primes for their information,
which it does happily, because it ranks by density and never by cost.

**And forty is where the time is best.**  Core cycles at height 200000,
relative to what the program does today (`PRIME_SIZE=7`, `-p 30`), median of
three paired runs; the last two columns need `PRIME_SIZE=9`:

| curve | 8, `-p 30` | 8, 35 | **8, 40** | 8, 45 | 8, 53 | 9, 70 | 9, 96 |
|---|---|---|---|---|---|---|---|
| record curve | 0.967 | 0.471 | **0.431** | 0.447 | 0.467 | 0.530 | 0.569 |
| starved #2 | 0.971 | 0.521 | **0.501** | 0.535 | 0.556 | 0.618 | 0.694 |
| starved #3 | 0.981 | 0.603 | **0.576** | 0.616 | 0.647 | 0.717 | 0.806 |
| starved #4 | 0.971 | 0.623 | **0.572** | 0.610 | 0.664 | 0.726 | 0.816 |
| starved #5 | 0.957 | 0.549 | **0.511** | 0.552 | 0.583 | 0.675 | 0.713 |

**1.7 to 2.3 times faster**, with the optimum at `-p 40` for every one of them
and a shallow minimum -- 35 and 45 are within a few per cent.  Forty primes
reach to 173, so **`PRIME_SIZE = 8` is enough and 9 and 10 buy nothing**:
every configuration that needs them is worse than `8, -p 40` on every curve.

**What it must not be is a blanket default.**  The same `-p 40`, in cycles:

| | ratio |
|---|---|
| random curve, h=400000 | 1.109 |
| random curve, h=50000 | 1.096 |
| point-rich but *not* starved, h=200000 | 1.260 |
| `make test1` | 1.243 |
| `make test1many` | 1.006 |

Ten per cent on random curves at any height, and twenty-six on a point-rich
curve that did not need the extra primes.  Nor does the penalty amortise the
way a pure set-up cost would: at `-p 40` a random curve is given `sp1 = 11`
instead of 12, out of larger primes with larger tables, so phase 1 itself gets
slower.  Going further is worse again: on `make test1`, `-p 53` costs 64% and
the full `PRIME_SIZE=10` table 6.6 times.

So the two constants want to move differently:

* **`PRIME_SIZE` to 8.**  It is free -- within 2.5% on everything measured in
  the previous section, and about 1 MB of address space now that the
  reservation follows `num_primes`.  It is the smallest table that contains
  the whole win.  (One side effect: the fast Horner path in `sieving_info`
  needs `(degree+1)*PRIME_SIZE <= 64`, so degree 8 loses it.  Measured on a
  degree-8 curve at `-p 30`: 1.038 at h=2000 and 0.970 at h=50000, i.e.
  nothing.)
* **`RATPOINTS_DEFAULT_NUM_PRIMES` stays at 30**, and the extra primes are
  looked at only for the curves that run out -- TODO item 6.  The criterion
  costs nothing to evaluate and the two populations are cleanly separated by
  it: nothing in `testdata.h` is ever starved, and the starved fifth of
  `testdata-many.h` is exactly the set that gains.

## Open questions

* The bit-extraction loop `for(a = a0; nums; a += d, nums >>= 1)` walks from
  bit 0 to the highest set bit.  On Drew's curve it runs 185.6 million
  iterations to find 5.8 million set bits -- 32 iterations per bit, the
  expected position of a single random bit in a word.  `__builtin_ctzl` would
  make it one iteration per bit.  This is a few per cent of the run in the
  point-rich regime and negligible in the sparse one.
* Whether real AVX-512 hardware changes the sparse picture.  The model says
  phase 2 would go from 3.06 to about 2.66, since only the scan term halves;
  phase 1 is where the gain would have to come from.
