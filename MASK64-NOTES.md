# The 2-adic mask modulo 64

Branch `mask64` off `v2.3` at ae47ef0, 2026-09-16.  The review's P6
(`review-2026-09-12/REVIEW-REPORT.md`; verdict F16 in `verdicts-final.md`,
both skeptics confirmed it and each built a prototype), the second item of
the sieve group; TODO item 26.  Measured on the i7-1355U, pinned to P-core
0, in core cycles, paired and interleaved with the base, medians of the
per-round ratios, outputs compared with the references before anything was
timed (`pair.sh` in `review-2026-09-12/`).

## What was there

`get_2adic_info` (find_points.c) decided modulo 16 which numerators a
denominator can have: `num_bits[b mod 16]` was a word of period 16 (8 in
the num_odd and num_even packings) that the first prime of the first phase
ANDs into every bit array it sieves, and `den_bits` said which classes of
the denominator mod 16 have any admissible numerator at all.  The odd
denominators were decided from a table "f(a) is a square mod 16" (a = 0..15)
through the inverse of b mod 16; the even ones from hand-derived
congruences for f(odd/2), f(odd/4), f(odd/8) and f(odd/2^n), n >= 4, with
shortcuts like "k^d = k^2 mod 16" that hold only up to a harmless
discrepancy.  The manual said "higher powers of 2 would be possible, but
not very likely to give a significant improvement" -- a hunch, never
measured.

## The mathematics

With D the degree of f rounded up to an even number, F(a,b) = b^D f(a/b)
is a form of even degree D with integer coefficients, so F(a,b) mod 64
depends on a and b mod 64 only, and a point needs F(a,b) to be a square.
The squares mod 64 are the twelve residues 0, 1, 4, 9, 16, 17, 25, 33, 36,
41, 49, 57; mod 16 they are four of sixteen.  The finer modulus changes
nothing for odd F (an odd number is a square mod 64 exactly when it is one
mod 8) and halves what is accepted where F = 0 mod 4: mod 16 the classes 0
and 4 pass wholesale, mod 64 only 0, 4, 16 and 36 of the eight residues mod 64
they contain.  Modulo 256 the squares are 44 of 256, which would remove one
admissible class in twelve more -- a quarter of what the step from 16 to 64
removes -- for sixteen times the set-up; 64 is where it stops paying.

Two words decide every class of the denominator exactly, with no lifting
and no hand-derived congruences:

* **odd b**: b^D is the square of a unit, so F(a,b) is a square mod 64
  exactly when f(a b^-1 mod 64) is.  Bit k of `fsq` says that f(k) is a
  square mod 64 (64 evaluations of f by Horner's rule).
* **even b**: the numerator is odd, a^D is the square of a unit, and
  F(a,b) = a^D frev(b/a) with frev(t) = t^D f(1/t) = c_0 t^D + c_1 t^(D-1)
  + ... + c_D (c_D = 0 when the degree is odd).  Bit t of `gsq` says that
  frev(t) is a square mod 64 (32 evaluations, t even).

The verdict's two skeptics disagreed on the even denominators.  Skeptic 1
built exactly this form; skeptic 2 held the uniform rule to be "vacuous for
even b" because "2^6 | b^d makes F = 0 mod 64" -- which treats f(a/b) as an
integer.  It is not: F(a,b) = sum c_j a^j b^(D-j) has the term c_D a^D with
no factor of b at all (F(1,2) = 39 mod 64 on the first random curve of the
check below).  Skeptic 2's prototype, which lifts k = a (b/2^v)^-1 from
mod 2^(6-v) to mod 64 and accepts if any lift gives a square, is exact
too, but only because all lifts give the same F(a,b) mod 64 -- the lifting
is spurious.  The old mod-16 arms were exact mod 16 as well (checked arm by
arm: every shortcut differs from the exact value by 8 on odd residues,
which swaps 1 and 9 and 3 and 11 and so on and preserves squareness mod
16), so the new rule refines the old one class by class.

Checked by brute force (`check64.py`, in the measurements directory):
300 random polynomials of degrees 1 to 10, every pair (a mod 64, b mod 64)
not both even -- 921600 pairs -- F(a,b) mod 64 against the fsq and gsq
rules and against the scatter construction below: no mismatch.

## What was done (8fda20f, 0c638c7)

`get_2adic_info` is about half its former length.  The coefficients are
reduced mod 64 (`mpz_get_si` keeps the low bits and the sign of a
coefficient too large for a long, which is all that is needed); the two
words are built by Horner's rule in unsigned long arithmetic, which is
exact mod 2^64 so that the residue is taken once at the end; the packing
(num_all/odd/even/none) is chosen from the parities present in `fsq`
(bit a of an odd denominator's pattern is bit a b^-1 of fsq, so the
parities are the same for every odd b) -- on the finer information, a free
refinement; the 64 patterns are built by *scattering*: for every set bit k
of fsq and every odd b the numerator a = k b gets its bit in b's pattern,
and for every set bit t of gsq and every odd a the numerator a gets its bit
in the pattern of b = t a (t = 0 is the class b = 0 mod 64, which admits
every odd numerator when frev(0) = c_D is a square -- always for odd
degree).  That is 32 iterations per set bit of the two words, no inverses,
some 2000 iterations on a random sextic.  `den_bits` gets bit b exactly when
class b has a pattern: that is what it meant before too (the old bit for
the odd classes was "some numerator works for some odd b", which is the
same thing since a b^-1 runs through all residues), so the word walk's
second mask over `num_bits` and run_shape's AND of the two are redundant
and went; run_shape counts the bits of `den_bits`.

Elsewhere: `num_bits[16]` -> `[64]`; of the eight index sites `b & 0xf`
(the review counted seven; the denominator walk of item 24 added one) six
became `b & 0x3f` and two -- run_shape's AND and the walk's second mask --
went, run_shape's square-denominator loops to period 32 (k^2 mod 64
has period 32 in k), `bits_per_word` averaged over the 64 classes, the
messages and comments, the DEBUG dump of the patterns.  The squares-mod-16
table is gone.  sift.c does not change: `sift.o` is byte-identical to the
base (the parameter `bits16` was renamed `bits64` in 0c638c7, which changes
no code).  test3.sh's line 43 grepped the verbose output for the wording
"mod 16"; it and the two reference lines say "mod 64" now -- the only
change to a reference file, and a wording one.

## The set-up cost

The verdict's prototype (uniform Horner over all 64x64 pairs first, then a
reversed form) cost 16.7k instructions per curve more than the mod-16
set-up and lost 4.5% at height 100, 1.3% at 400, with break-even near
600-1000; the checklist said to make it lazy or height-conditional if that
held.  It does not hold: the scatter construction costs **2260
instructions per curve** more than the mod-16 code (perf stat, `rptest -h
100`, 1008 curves): +0.66% of a run at height 100, and at height 200 the
run is already 0.3% shorter, 1.2% at 400, 1.7% at 1000, 2.4% at 4000
(instructions, default build, one reading each).  Nothing is lazy or
conditional.

## What it is worth

Chain `chain-mask64.sh` (16:04-16:45): the counter builds, `alt.sh`, then
`pair.sh` -- 0c638c7 against ae47ef0, 5 rounds of eight suites at the
default code placement, 3 rounds of five at `-falign-loops=32` and `=64`
(the worktrees built from the same sources with that one flag added).
Cycles new/base, medians of the per-round ratios, [min, max] in the raw
reports (`m-base-new*` in `review-2026-09-12/measurements-2026-09-16/
mask64/`):

| suite | default (5 rounds) | -falign-loops=32 (3) | -falign-loops=64 (3) | instructions |
|---|---|---|---|---|
| test1 | 0.970 | 0.956 | 0.988 | 0.972 |
| test1many | 0.984 | 0.985 | 0.998 | 0.982 |
| testhigh | 0.968 | 0.964 | 0.971 | 0.967 |
| testhighmany | 0.985 | 0.993 | 0.990 | 0.987 |
| height 200 | 0.980 | | | 0.997 |
| height 1000 | 0.970 | 0.972 | 0.980 | 0.982 |
| height 4000 | 0.974 | | | 0.976 |
| point-rich at 1000 | 1.001 | | | 0.990 |

**Worth 3% of `make test1` (1.2 to 4.5 at the three placements) and of
`make testhigh` (2.9 to 3.6), 1.5% of `make test1many`, about 1% of
`make testhighmany`, 2 to 3% at height bounds from 200 to 4000, nothing
measurable for the point-rich set at height 1000.**  The verdict had
predicted 2.5% of test1, 3% of testhigh and 1.3 to 1.9% for the point-rich
suites.  Branch misses are unchanged within noise (the change removes
work, not mispredictions).

The counters (`RP_PHASE_TIMING -DRP_PHASE_COUNTS -DRP_PRIME_STATS` builds,
`prelim-mask64.txt`), base -> new:

| suite | phase-1 ANDs | bits entering | arrays swept | table rows | exact checks | denominators |
|---|---|---|---|---|---|---|
| test1 | -4.5% | -21.7% | -0.76% | -6.4% | 49659 -> 48819 | -0.9% |
| test1many | -2.0% | -8.4% | -1.1% | -1.2% | 104901 -> 95953 | -1.6% |
| testdegrees | -4.9% | -22.1% | -0.15% | -11.2% | 2558 -> 2428 | -0.2% |
| height 1000 | -3.8% | -21.9% | -0.74% | -9.4% | 3766 -> 3576 | -0.9% |
| testhigh | -4.5% | -21.5% | -0.77% | -2.5% | 318914 -> 299898 | -0.9% |
| testhighmany | -1.7% | -3.8% | -1.4% | -0.9% | 250903 -> 496629 | -2.4% |

The mechanism is the one the verdict described: the bits entering the
first phase fall by a fifth on the random curves (a quarter of the mod-16
admissible classes are gone), so the phase needs about half a prime less
to reach the same survivor density -- `sp1` 11.34 -> 10.87 on test1,
10.48 -> 10.06 at height 200000, 12.54 -> 11.94 on testdegrees (52 of 87
curves lose one prime), 17.01 -> 16.87 on test1many (15 of 98) -- and the
survivors of the phase are essentially unchanged.  Three curves more of
test1 (879 against 882 of 1008 reach the sieve) return before sieving
because no class of the denominator admits a numerator mod 64; they had
no points before either, the output being byte-identical.
The one thing to watch is the last column's last row: on the point-rich
height-200000 suite the exact checks double (the verdict's skeptic saw the
same, 252252 -> 497970).  Four of the thirty curves get fewer primes for
their third stage: `bits_per_word` feeds the `may_extend` rule of
`sieving_info`, which stops looking for further primes sooner when the mask
is sharper; three of them lose one prime and the first curve of the suite
(`456976 -448032 -255200 208380 61033 -12834 81`) has 29 instead of 37.  The suite still
comes out 1 to 1.5% faster because exact checks are cheap relative to
sieving there, but it is the existing extension rule responding to a
smaller `bits_per_word`, not a fault of the pattern -- a matter for the
tuning of the rule, noted in TODO.md.

`alt.sh` (`alt-mask64.log`): the plain 64-bit build, AVX-128, SSE, the
emulated AVX-512, AVX without AVX2, `RATPOINTS_CHUNK=1`, `USE_LONG_IN_PHASE_2`,
`RP_MULMOD_DIVIDE`, the phase-timing build and the default all reproduce
the references (the timing build's test3 differs by its report, as
always), test2, testhigh and testhighmany reproduce theirs, valgrind is
silent on the debug build and on two optimised runs.

`make tune` runs on the merged tree, as it did for the trio and the tail
leg; `sp1` falls by half a prime here, so the constants deserve the check.
The result is recorded in TODO.md item 26.
