# Notes on where the time goes, and on some things that did not work

Working notes from a round of profiling in September 2026, kept so that the
measurements do not have to be repeated and the dead ends are not walked into
again.  Nothing here describes how to use `ratpoints`; see the documentation
for that.

## Where the time goes

Two workloads behave very differently, and it is worth being explicit about
which one an optimisation is aimed at.  Profiles from `gprof`, AVX-256 build:

| | `rptest`, 1000 curves at height 16383 | one curve at height 400000 |
|---|---|---|
| `_ratpoints_sift0` | 55.4% | 96.7% |
| `sieve_init_<p>` | 26.5% | ~0% |
| `find_points_work` (self) | 12.5% | 0.2% |
| `sift` (self) | 5.7% | 3.0% |
| `_ratpoints_check_point` | ~0% (14955 calls) | 0.2% |
| `_ratpoints_compute_sturm` | ~0% | ~0% |

Two things are worth noting.  The exact check with `gmp` costs essentially
nothing in either case -- the sieve is doing its job well enough that
verification is free -- so there is nothing to gain there.  And for large
heights everything is in `_ratpoints_sift0`; any optimisation for serious users
has to be inside it, or has to reduce the work it is given.

Inside `_ratpoints_sift0` at height 400000, with the defaults `sp1 = 11`,
`sp2 = 19`:

    phase 1 (and with sp1 primes) : 3.14 s (53%)
    phase 2 (scan + sp2-sp1)      : 2.74 s (47%)
    bit-arrays swept: 627,200,000   surviving phase 1: 26,162,862 (4.17%)

So the second phase spent most of its time discovering that 96% of the
bit-arrays were zero.

**The load budget is the thing to keep in mind.**  The first phase reads `sp1`
sieve bit-arrays for every bit-array of survivors; at `sp1 = 13` that is about
8.2 G of the 11.2 G L1 loads of the whole run, roughly three quarters.  Any
change confined to the second phase is moving a small denominator.  Both of the
mis-estimates recorded below came from reasoning about phase 2 in isolation.
## Skipping runs of empty bit-arrays in phase 2 (adopted)

The second phase now looks at four bit-arrays at a time and skips the group
when their `or` is zero.  This needed an `ORR` macro next to `AND` for each
register width.

**Where the gain really comes from.**  Not from the grouping.  The saving is in
having a *tight inner loop* that walks over empty bit-arrays without re-entering
the large outer loop body -- which otherwise recomputes `ssp = &sieves[sp1]`,
sets `n = sp2-sp1`, and evaluates `TEST` twice, for every one of the ~96% that
are empty.  A group size of one, which does no `or`s at all, already gets
-13.6% at the default `sp1 = 11`, where the group of four gets only -1.3%.

**The group size interacts with `sp1`, and not in one direction.**  Cycles at
height 400000, against the same build with no skip loop:

| group | 11/19 | 12/19 | 13/19 | 14/19 | 15/19 |
|---|---|---|---|---|---|
| 1 | **-13.6%** | **-13.1%** | -12.1% | -12.9% | -- |
| 2 | -7.1% | -11.9% | **-13.7%** | -15.2% | -- |
| 4 | -1.3% | -8.9% | -13.3% | **-15.8%** | **-17.8%** |
| 8 | -- | -- | -11.7% | -15.4% | -17.7% |
| 16 | -- | -- | -9.7% | -12.4% | -16.1% |

This is a ridge, not a slope.  At a low `sp1` the survival rate is high, and a
group that contains a survivor has wasted its whole `or` tree, so small groups
win.  At a high `sp1` the survival rate is low, nearly every group is empty, and
the cost per bit-array is amortised over the group, so larger groups win.  The
crossover is at about `sp1 = 13`.  Sixteen is never best; eight approaches four
from below but does not pass it.

Four is kept because the settings that are actually good are `sp1 = 13..15`,
where it is at or near the top; but note that at the *shipped* default of 11 it
is worth almost nothing, and a group of one would be worth 13%.

**Instructions are a bad guide here.**  At 12/19 the group of four executes
1.8 G *fewer* instructions than the group of one and needs 0.64 G *more*
cycles: four independent load-test-branch iterations are replaced by four loads
feeding a two-level `or` tree that has to resolve before the branch can.  Fewer
operations, longer critical path.  Callgrind instruction counts were actively
misleading for this change and only `perf` cycles settled it.

**Best against best**, over both `sp1` and the group size: about -13%
(15.38 -> 13.38 G cycles at height 400000).  That ratio came out the same in two
sweeps taken hours apart whose absolute numbers differ by 4%, which is the usual
pattern on this machine -- trust ratios measured close together, never absolutes.

**Precision.**  The rows for 13/19 and 14/19 were measured in two separate
sweeps and disagree by up to 1.5 points.  Differences below that are not
meaningful, which puts 2 and 4 in a tie at `sp1 = 13`, and 4 and 8 in a tie at
14 and 15.

The effect is much smaller for many small curves: instruction counts for
`rptest` at height 4000 fall by a flat 6.1-7.1% across all settings tried,
because that workload is dominated by `sieve_init` and the denominator loop.

Incidentally, the baseline's own optimum is not the shipped default either --
at height 400000 it uses about 16.2 G cycles at 11/19 against 15.38 G at 12/19
and 13/19.

**This makes the tuning problem three-dimensional.**  `sp1`, `sp2` and the group
size are not independent: the best `sp1` depends on how expensive the second
phase is, and the best group size depends on the survival rate, which is set by
`sp1`.  The open question is what combination is best for a given distribution
of height bounds and curves.  `sp2` is the least important of the three, being
close to neutral over 19..23 at this height.

## Per-chunk emptiness flags from phase 1 (tried, rejected)

The first phase holds a whole chunk of `RATPOINTS_CHUNK` bit-arrays in
registers just before storing them, so or-ing them together costs about one
operation per bit-array, and the second phase could then skip an empty chunk
with a single byte load instead of `RATPOINTS_CHUNK` vector loads.

It works exactly as designed -- at height 400000 with `sp1 = 13`, 84.0% of the
chunks are empty and all 84.0% are skipped -- and it is still slower:

| `sp1`/`sp2` | with group skip | + chunk flags | |
|---|---|---|---|
| 12/19 | 13,624,503,340 | 14,344,797,429 | +5.3% |
| 13/19 | 12,918,048,031 | 13,361,157,236 | +3.4% |
| 14/19 | 12,968,315,580 | 13,200,400,895 | +1.8% |
| 15/19 | 13,345,331,632 | 13,457,683,905 | +0.8% |

The load counts say why: skipping 84% of the chunks removes only about 3% of
all L1 loads (11.20 -> 10.88 G at 13/19), because of the load budget noted
above.  Against that, the or-accumulation costs some 0.6 G instructions across
all chunks, empty or not, and the extra bookkeeping in phase 2 costs about as
much again.  The group-of-four skip had already made scanning cheap enough that
removing it entirely saves less than the flags cost to produce.

The penalty does shrink as `sp1` grows (+5.3% -> +0.8%), so for a very low
survival rate this might break even; not at the settings that are good.

**If anyone re-tries this**: the obvious implementation drifts out of chunk
alignment, because the group-of-four skip can step across a chunk boundary and
then the aligned test never matches again.  That version skipped only 36% of
the 84% and looked like a clear loss for the wrong reason.  The group skip has
to be bounded by the end of the current chunk.

## The set-up phase (`sieve_init`)

Rebuilding the sieve tables is about 21% of `make test1` and essentially 0% of
a single curve with a large height bound, so it only matters for many curves at
moderate heights.  Within it, `sieving_info` -- the coefficient reduction, the
Horner loop over all residues, the sorting of primes -- is only about 1.5%;
everything is in `sieve_init_<p>` in `init.c`, and 95% of that is in the primes
above `LONG_LENGTH`, which cost roughly 720 ns per table against 92 ns for the
smaller ones.  `bench_init` measures this in isolation and checks the tables
against a reference computed straight from the definition.

Things that were tried and gave nothing (each measured against a noise floor of
about 3%):

* replacing the stride-`p` replication loop by a doubling `memcpy`: -2%
* keeping the rotation state in registers instead of the `help[]` array: 0%
  (gcc had already promoted it)
* writing the table into an L1-hot buffer instead of the 5.3 MB one: -5%,
  so it is not memory-bound
* reformulating the rotation loop for `p < LONG_LENGTH` as independent
  double-word shifts: **+17%, i.e. worse** -- the existing two-operation
  recurrence is already efficient

What did look promising in an isolated test was interleaving several
independent `(p, b)` table constructions, which gave -24% at eight-way: the
loops are serial recurrences and distinct `(p, b)` pairs are independent, so
the chains overlap.  It would need the lazy per-`(p, b)` call in `sift()` to be
restructured into batches, and it is worth at most a fifth of `sieve_init`,
which is itself a fifth of one workload and none of the other.

## Measuring on a laptop

The development machine is an i7-1355U: 2 P-cores with hyperthreading (logical
CPUs 0-3) and 8 E-cores (4-11), in a 15 W package.  Wall-clock timings on it
are not trustworthy without care.  Running the same benchmark pinned to CPU 0
while a spinner runs elsewhere:

| spinner on | slowdown | P-core clock |
|---|---|---|
| nothing | -- | 3.15 GHz |
| CPU 1 (SMT sibling) | +68% | 2.39 GHz |
| CPU 2 (other P-core) | +22% | 2.68 GHz |
| one E-core | +13% | |
| all eight E-cores | +88% | 1.80 GHz |

Note that loading only E-cores, which share no execution resources with the
measured core, still costs 88%: the package power budget is global, so having
spare cores does not help, and using them actively hurts.  For anything below
about 10%, prefer

* `perf stat -e cpu_core/cycles/` -- cycles are immune to frequency drift
  (needs `sysctl kernel.perf_event_paranoid=2`), or
* `valgrind --tool=callgrind` -- instruction counts, bit-for-bit reproducible,
  but blind to cache and latency effects, and about 80 times slower.

Wall-clock comparisons should be paired (both versions measured adjacently) and
taken as a minimum over several runs, and an outlier that disagrees with its
neighbours should be re-run rather than believed.
