# The 2-adic stride on the numerator

Branch `stride` off `v2.3` at cd66310, 2026-09-17.  TODO item 28, the
review's P3 (verdict F2: skeptic 1 uncertain about the size of the gain
only, skeptic 2 confirmed).  Measured on the i7-1355U, pinned to P-core 0,
in core cycles, paired and interleaved with the base, medians of the
per-round ratios, outputs compared with the references before anything
was timed (`pair.sh`).

## What was wrong

For a class of the denominator mod 64 the admissible numerators -- those
with b^D f(a/b) a square mod 64, item 26's pattern -- are a union of
residue classes mod 64, and on a random curve they often lie in a single
class modulo 4 or 8: computed from the pattern over the 1008 test curves
(the thousand random ones and the eight added in 2026), the odd
denominators take their numerators from one class mod 4 on 223 curves and
mod 8 on 155 (one class mod 2 on 54, no restriction on 486, none
admissible on 90; over the thousand random ones alone 220/155/54/481/90);
the even denominators reach 8 as well.  Nothing reaches 16, and the stride
is the same for every odd class and within each class of v_2(b) on every
curve, as the theory says (the pattern of b = u b' with u a unit mod 64 is
u times that of b', and multiplication by a unit keeps a set of residues
within one class mod 2^k).  Against the mod-16 information (491/53/232/147
and 85 with nothing admissible) the modulus 64 raises the stride on 11 of
the 1008 curves -- one from 1 to 2, two from 1 to 4, two from 1 to 8, six
from 4 to 8 -- and excludes five more that mod 16 would have sieved at
stride 4 (the sweep reviewer's transition table); nothing reaches 16
either way.

The program knew only the two-fold packing: `which_bits` said whether the
odd denominators take all, only even or only odd numerators, and an even
denominator always took the odd ones.  So on two fifths of the random
curves the bit arrays held two or four times as many candidates as they
need to, and the first phase swept, and the scan walked, two or four
times as many bit arrays for those classes.

## What was done

**Per class of b mod 64 a packing.**  `get_2adic_info` (find_points.c)
builds the unpacked pattern of every class as before (the two words fsq
and gsq of item 26, scattered), then reads off each class the largest k
with all admissible numerators in one class a0 mod 2^k, and packs: bit t
of a bit array stands for the numerator a0 + 2^k t, the packed pattern has
period 64/2^k in t and is one word repeated as before.  The class record
(`rp_num_class` in rp-private.h: `bits`, `k`, `a0`, `offset`) replaces
`num_bits[64]` and the `bit_selection` enum, which is gone from every
interface.

**The table row.**  Bit t wants the pattern for the residue
(a0 + 2^k t) b^-1 = (t + a0 2^-k) (b 2^-k)^-1 mod p, so the row for the
denominator is the one for b 2^-k mod p, read a0 2^-k bits further on,
that is a0 (2^k RBA_LENGTH)^-1 mod p bit arrays further on -- what the
compiled `offsets[] = (2 RBA_LENGTH)^-1` did for k = 1, a0 = 1 (it is gone
from gen_find_points_h.c).  `fill_bp_list` computes b 2^-k mod p directly,
by one multiply with 2^-k mod p (kept per prime in the sieve entry,
`dinv[k]`) inside the multiply-high reduction; the row shifts are one
table per curve, `class_offsets()`, one row of sp3_max entries per
distinct (k, a0) with the RP_ROW_BIAS multiple built in, on
find_points_work's stack, so sift()'s per-prime loop is a load instead of
the halving of bp and the offset select.  `class_offsets` also fills
`dinv[k]` and the inverses of 2^k RBA_LENGTH it needs, by halving from
1 and from RBA_LENGTH^-1, for the strides actually in use.  The p | b row
(`sieves0`, "bit index not divisible by p") needs no special case, as the
verdict said.

**The interval and the extraction.**  sift() maps the numerator interval
[low, high] to bits ceil((low - a0)/2^k) .. floor((high - a0)/2^k) by
shifts; sift0's two extraction sites set d = 2^k, a0 + word * 2^k
RBA_LENGTH, da = 2^k LONG_LENGTH from the class.  sift.c changes in six
lines plus the signature.

**The model.**  run_shape's factor on the numerators is the mean over the
classes kept of 1/2^k (item 27's half-width factor was the k <= 1 case);
`bits_per_word` rises by itself, so sp1 rises through the existing rule
(10.91 -> 11.68 on test1, 10.07 -> 10.65 at 200000).

**Verbose output** names the stride for the seven kinds of denominator
(b odd, 2 mod 4, ..., 32 mod 64, 0 mod 64); the `[runshape]`
instrumentation line has `kodd=` (the odd classes' k) instead of `wb=`.

## Correctness

`make test` passes; test1, test1many and testdegrees byte-identical with
the references at 64 and 128 bits and with RATPOINTS_CHUNK=1 as well as
at 256; alt.sh: all ten build variants (64bit, avx128, sse, avx512emu,
avx-only, chunk1, long2, mulmod-divide, phase-timing, default) pass the
three suites and test3 (the phase-timing build's test3 differs by its
counter lines on stderr, as it did for items 26 and 27), test1once,
test2, testhigh and testhighmany identical, valgrind clean on the debug
and the optimised binary.

**Reviewed** by two Opus agents (briefs rev-common/sweep/skeptic-stride.txt
in review-2026-09-12/), no correctness fault.  The sweep derived the row
shift from scratch and checked the halvings over every odd prime below
1024, every register width and every stride; a brute force over all
coprime (a, b) up to height 250 and 2000 on 86 curves covering 86 distinct
stride patterns agreed with the program, as did eight option sets on a
subset; eleven build variants (the ten of alt.sh and PRIME_SIZE 7 and 9)
warning-free and byte-identical; the manual's overfull boxes are the base's
four; every number in the notes recomputed from the raw reports.  It found
two statements in the notes that their own numbers contradicted (the
mod-16 comparison, the curve count), a cost figure that mixed two builds,
stale `offsets[]` and `which_bits` in gen_find_points_h.c, the Makefile
and ratpoints.h, and a dozen nits (c614d86: the rows sized by the number
of packings, `kodd=-1` for an empty odd class, the verbose report's dash
explained, test3 pins the report).  The skeptic packed 219 random curves
of degrees 1-10 with coefficients up to 2^70 independently of fsq/gsq
and compared with the DEBUG build's per-class dump (0 mismatches),
checked 2148 (denominator, prime) rows of the DEBUG trace against the
row formula, ran 900 randomized -l/-u interval cases and 960 sieve-free
cases (-x -n 0 -N 0 -P 0 -j -F 0) against brute force, 3320 differential
cases against the base over five builds and both reduction paths with
identical points (only -x survivor lists differ, legitimately), ASan and
UBSan silent, the U drift attributed to the array rounding by stride
class (zero at stride 1, +0.7 at stride 8 at height 1000), no curve with
its exact checks doubled.  It found the guard of fill_bp_list ignoring k:
above 2^24 a stride-1 class took the division where the base multiplied
(5.5% more instructions on such a curve at such heights); fixed with a
middle tier of two reductions up to 2^32.

## What it is worth

Counters (`-DRP_PHASE_TIMING -DRP_PHASE_COUNTS -DRP_PRIME_STATS`), new
over base:

| suite | arrays swept | phase-1 ANDs | sift0 calls | survivors of phase 1 | exact checks | mean sp1 base -> new |
|---|---|---|---|---|---|---|
| test1 | 0.850 | 0.881 | 1.000 | 0.849 | 0.939 | 10.91 -> 11.68 |
| test1many | 0.978 | 0.983 | 1.000 | 0.983 | 0.989 | 16.95 -> 17.03 |
| testdegrees | 0.782 | 0.808 | 1.000 | 0.804 | 0.956 | 12.01 -> 13.12 |
| test1 at 1000 | 0.897 | 0.928 | 0.999 | 0.963 | 0.876 | 13.01 -> 13.61 |
| test1many at 1000 | 0.984 | 0.990 | 1.000 | 0.993 | 0.997 | |
| testhigh | 0.850 | 0.878 | 0.887 | 0.807 | 0.891 | 10.07 -> 10.65 |
| testhighmany | 0.989 | 0.991 | 0.990 | 0.970 | 1.004 | |

The verdict predicted arrays x 0.86 at both heights before the chunk
padding (which item 25 removed since) and 0.99 point-rich: reproduced
(0.98-0.99 point-rich).  Instructions (perf stat, whole run, final code;
setup-instructions.txt in the reports): test1 0.921, test1many 0.972,
testdegrees 0.907 -- the point-rich suite gains more than its
arrays because the per-denominator loop got simpler.

Cycles, pinned, new over base (median of the per-round ratios; 5 rounds
at the default placement, 3 at each of -falign-loops=32 and 64):

| suite | default placement (5 rounds) | -falign-loops=32 (3) | -falign-loops=64 (3) |
|---|---|---|---|
| test1 | 0.913 (0.898-0.915) | 0.926 (0.906-0.942) | 0.920 (0.911-0.950) |
| test1many | 0.991 (0.983-0.994) | 0.988 (0.961-0.989) | 0.989 (0.981-0.994) |
| testhigh | 0.875 (0.870-0.879) | 0.881 (0.879-0.882) | 0.887 (0.884-0.887) |
| testhighmany | 0.982 (0.980-0.985) | 0.984 (0.982-0.986) | 0.990 (0.987-0.990) |
| test200 | 1.048 (1.039-1.062) | | |
| test1000 | 1.035 (1.028-1.046) | 1.039 (1.035-1.046) | 1.040 (1.031-1.106) |
| test4000 | 0.970 (0.963-0.976) | | |
| testmany1000 | 0.998 (0.989-1.002) | | |

Instructions at the default placement: test1 0.921, test1many 0.972,
testhigh 0.871, testhighmany 0.985, test200 1.043, test1000 1.024,
test4000 0.978, testmany1000 0.988.  The first form (6716603, before the
set-up cuts) measured test1 0.902, test200 1.081, test1000 1.044 in the
same protocol (first-chain-6716603/ in the reports): the cuts took three
points off the small heights, and test1's difference is the noise between
two chains.

The verdict said 3-5% of test1 and 9-12% of testhigh.  It discounted for
the chunk padding, which item 25 removed since, and it did not count the
simpler per-denominator loop (2.2k instructions per curve at height 100,
a few per cent of test1) or the scan, which walks the arrays too.

The run-length estimate: log(Uact/Upred) on test1 mean 0.071 sd 0.077
before, 0.088 sd 0.100 after; at height 1000 the geometric mean of
Uact/Upred goes 1.72 -> 2.04.  With fewer
bits per denominator the rounding of each denominator's range up to whole
bit arrays weighs more (the Left-over "U at small height bounds" of item
27); at 200000 unchanged (1.051 -> 1.053).

## The set-up, and the small height bounds

The first form of the change (6716603) cost 26.3k instructions per curve
more than the base at height 100 (+7.6%; perf stat over the 1008 curves
of `rptest -h 100`, setup-instructions.txt in the reports), 28.9k at 200
and 21.1k at 1000, against the 11k the sieve saves there.  Callgrind on
`rptest -h 100` split it (taken with the quadratic search for a shared
row already replaced by a key look-up, hence a total of +21.2k there;
setup-callgrind.txt): the seven halvings per prime for `dinv` and the row
inverses in examine_prime 3.7k, the packing and the offset rows in
find_points_work 4.8k, a malloc and free per curve for the rows 4.2k, a
floating-point division per class in run_shape 0.6k, the tables of the
extra first-phase prime 6-7k; sift saves 2.2k (its per-denominator loop
lost the halving and the offset select).  The second commit (00673e9)
takes what is cheap to take: the stride search runs upwards and stops at
the first failure (most classes have stride 1 or 2), a class with stride
1 is not repacked, the rows are deduplicated by the key 2^k + a0 and live
on the stack, the inverses are computed once per prime and only up to the
largest stride in use, run_shape reads 1/2^k from a table.  Now +14.5k per
curve at height 100 (+4.2%), +17.3k at 200 (+4.3%), +15.2k at 1000
(+2.4%) by perf stat; callgrind's split of the +14.7k at 100:
find_points_work +7.3k (the packing of 64 classes, about a dozen rows of
30 shifts), run_shape +0.6k, sift0 +0.4k, sift -2.2k, and the tables
+8.9k.

The tables are not this item's doing but the sp1 rule's: `bits_per_word`
rises with the packing, so `primes_for_phase_1` takes 0.6 primes more at
height 200 (12.8 -> 13.4), and at that height a first-phase prime costs
about thirty tables and saves nothing measurable -- the `sieve_init_*`
functions are a quarter of the instructions there, sift0 six per cent.
The base already over-buys primes at small heights for the same reason;
the rule balances an AND per word against the second phase and knows
nothing of the tables.  That is TODO item 29, for the tuning session.
What remains of the set-up, about 2% of a run at height 100, is the price
of deciding 64 classes and building their rows; nothing lazy would help
at that height, where every class is visited.

## Not done, and why

* Rows of sp2 instead of sp3_max entries (adapt_primes can raise sp2
  during the run, so the rows would have to grow with it): 1.7k
  instructions per curve at height 200, not worth the plumbing.
* A bit-compress instead of the loop over set bits for the repacking
  (cheaper for the dense stride-2 patterns, dearer for the sparse ones):
  a kilo-instruction per curve either way.
* The per-class generality (a stride per class of b mod 64, up to 64)
  costs nothing over the verdict's "per class of v_2(b)"; the theory is
  used only by the verbose report.
