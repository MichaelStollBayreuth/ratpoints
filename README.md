# ratpoints

## The ratpoints library and command-line program

This is a program that uses an optimized quadratic sieve algorithm in order
to find rational points on hyperelliptic curves.

The program is distributed under the GNU GPL, version 2 (or later).

Read the [documentation](https://www.mathe2.uni-bayreuth.de/stoll/programs/ratpoints-doc-2.2.pdf).

The current version is ratpoints-2.2.3 from September 6, 2026. This version can use 256-bit AVX
registers and has been optimized further, so that it now runs considerably faster than ratpoints-2.1.3.

The number of primes used in each of the two sieving stages is no longer a fixed default: unless
`-n` and `-N` are given, both are chosen from the curve, from the densities of admissible numerators
modulo the small primes that the program computes anyway. That is worth about 14% on random curves
and about 36% on curves with many rational points, which sieve very differently. A curve with very
many rational points can also leave too few of the small primes saying anything at all for that
choice to have the primes it wants; the program then looks at a few more of its own accord, which
is worth up to a factor of two on such curves and is never asked for otherwise.

There is now a third sieving stage as well. The two sieving stages work on bit arrays, which is the
right shape while a whole array's worth of numerators is still in play; by the time they are done
there is about one candidate left per denominator, and each of those went straight to the exact
test, which works in multi-precision arithmetic and takes an integer square root. The third stage
tests those candidates against further primes first, one at a time, working out where to look in
the table of admissible residues rather than reading a bit out of a sieving table. So it builds no
table, which is what makes it worth doing at that point, and it can use primes far up the list for
which building one would be out of the question. How many it uses is decided per curve; on a random
curve at a small height bound it uses none, and on curves with many rational points at a large
height bound it is worth about 7%.

How many primes each stage uses, and which ones, now depends on how long the run is going to be
and on what each prime costs, both of which the program works out before it sieves anything. A
prime's sieving table has `p` rows and is rebuilt for each denominator class that turns up, so it
costs `O(p^2)` spread over the whole run: at equal quality a smaller prime is strictly better, and
at a small height bound it is much better. Ranking the primes by what they cost rather than by what
they say alone is worth 5% of `make test1`, and scaling the number of second-stage primes by the
length of the run is worth 12% of `make testhighmany`. While it runs, the program also counts what
the sieve is actually finding and corrects the third stage from the count — which gets at the one
thing no prediction can, namely that the non-reduced forms `(ka, kb)` of a rational point pass
every prime test, so a floor of survivors outlives any amount of sieving.

The first sieving phase no longer reads the bit arrays it is about to overwrite. Every bit array
starts from the same pattern of admissible numerators modulo 64, and that pattern used to be
written into all of them on a pass of its own; the first prime now ANDs it in as it sieves, so the
pass is gone and with it one store and one load per bit array. That is worth 5.6% of
`make testhigh`, 5.3% of the degree suite at a height bound of 200000, 3.4% of
`make testhighmany` and 2.6% of `make test1`.

Two loops of the second sieving phase have been rewritten. The scan that steps over the empty bit
arrays after the first phase used to carry a bound test and two counters, nine instructions per bit
array; it now runs into a non-zero sentinel written just past the range and is a load, a test, a
branch and the pointer step. And the phase used to reach the table row of a surviving bit array
from a pointer set up at the start of every call, subtracting the prime until the pointer was back
inside the table -- one to three data-dependent, mostly mispredicted branches for every AND; it now
computes the row from the word number with a precomputed reciprocal, with a multiple of the prime
folded into the offset so that nothing is ever negative. Together they are worth 3% of
`make test1` and `make test1many` and 7.5% of `make testhigh` and `make testhighmany` at the
default code placement, between 4.5% and 7.5% at others.

The loop over the denominators, and what is done once per denominator, have been reworked. The
Jacobi symbol test on the denominators is evaluated as a product of Legendre symbols read from
tables built once per curve, instead of by a binary algorithm run for each denominator, which was
12% of `make test1` and a third of its branch mispredictions; the tests on `b mod 64` are applied
to a word of 64 denominators at once and only the denominators that pass them are visited; the
residues of `b` modulo the sieving primes are computed with a precomputed reciprocal instead of
being stepped from the previous denominator by mispredicted conditional subtractions; the gcd that
tests a candidate for lowest terms is branchless; the forbidden-divisor bit arrays cover every
excluded prime up to the square root of the height bound, not only the compiled sieving primes
(the default of `-F` is 64 now); and the third stage's per-denominator set-up happens only when a
candidate reaches it. Together: 24% of `make test1`, 4% of `make test1many`, 6% of
`make testhigh` and 2% of `make testhighmany` at the default code placement, and within a point of
that at two others, with two thirds of the branch mispredictions of `make test1` gone.

The first sieving phase no longer sieves bit arrays that are then thrown away. It works on 16 bit
arrays at a time, one per vector register, and the range of every denominator and interval used to
be padded up to a multiple of that; the padding was sieved with every prime of the phase, walked by
the scan of the second phase and zeroed -- a seventh of all the bit arrays swept in `make test1`,
four fifths of them at a height bound of 1000. The chunk loop now stops at the last whole chunk, and
the bit arrays left over are sieved in legs of 8, 4, 2 and 1 registers, one leg for each set bit of
their number, which neither wrap the table pointers around nor store them back. Worth 10% of
`make test1` (9.5% to 11.4% at three code placements), 6.5% of `make test1many`, 2% of
`make testhigh`, between nothing and 2.5% of `make testhighmany`, and a fifth of the run at height
bounds of 1000 and 4000.

The information from the polynomial modulo a power of 2 is taken modulo 64 instead of 16. The
pattern of admissible numerators that the first prime ANDs into every bit array is one repeated
64-bit word whatever its period, so the finer modulus costs the sieve nothing, and the squares are
12 of the 64 residues modulo 64 against 4 of 16: where `b^D f(a/b)` is 0 mod 4 the coarser modulus
let twice as much through. The even denominators are decided exactly as well: with `a` odd,
`F(a,b) = a^D frev(b/a)`, where `frev(t) = t^D f(1/t)` is the reversed polynomial, an integer
polynomial in `t`, and `a^D` is a unit square, so one table of `frev` modulo 64 replaces the
hand-derived congruences. On random curves the first phase needs half a prime less and does 4.5%
fewer ANDs; the set-up costs 2260 instructions more per curve, so a run at height 100 is 0.7%
longer and the gain begins at height 200. Worth 3% of `make test1` and of `make testhigh`, 1.5% of
`make test1many`, about 1% of `make testhighmany`, and 2 to 3% at height bounds from 200 to 4000.

The choice of primes for the third sieving stage no longer starves. The stage used to take primes
from those already looked at while each paid for itself and stop at the first that did not, asking
for more only when none were left; on a curve with very many rational points the small primes say
nothing and the informative ones come late, so one such curve of `make testhighmany` was left with
no third stage at all and 250000 exact checks instead of 7900. Now the stage looks at a further
prime whenever none in hand pays but a prime of the density the curve has been offering would (the
same test now also decides whether an empty pool is refilled at all; a fixed number of primes is
taken as before). And
the estimate of the run's length counts the numerators of even denominators at half width, which was
what the sieve then swept for them: it was too large by up to a third on such curves, and its spread
over the random curves halves. Worth 1% of `make test1many`, half a per cent of `make testhighmany`
and nothing measurable elsewhere; the exact checks of `make testhighmany` fall 4.6 times.

The bit arrays hold only the numerators a denominator can have modulo a power of two. The
paragraph before last says which numerators mod 64 each class of the denominator mod 64 admits;
on two fifths of the random test curves those of the odd denominators lie in a single class
modulo 4 or 8, and those of the even denominators often do as well. The program used to notice
only whether they were all odd or all even and pack every second numerator; now every class gets
the largest stride `2^k` with all its admissible numerators congruent modulo `2^k`, and a bit of
its arrays stands for every `2^k`-th integer. The sieve tables serve unchanged: the denominator's
residue is multiplied by the inverse of `2^k` and the row is read at a shift that depends on the
class, which is what the two-fold packing already did, and the sieving loops do not change. On
the random curves the first phase sweeps 15% fewer bit arrays and does 12% fewer ANDs (the rule
then wants three quarters of a prime more in the first phase), and the work done once per
denominator got simpler. Worth 9% of `make test1` (7.4 to 8.8% at three code placements), 12.5% of
`make testhigh`, 1% of `make test1many`, 1.5% of `make testhighmany` and 3% at height bound
4000. A run at height 1000 is 3.5% longer and one at 200 5% longer: deciding the 64 classes costs
about 7000 instructions per curve, and the extra first-phase prime brings tables that a run of
that length cannot use -- the rule that counted the first-phase primes knew nothing of their
tables, which predated this change and is fixed two paragraphs below.

The sieve can use composite moduli. A table row only has to be periodic in the bit index with the
modulus as its period and be selected by the denominator's residue, and nothing in the sieving
loops asks for a prime; the program used primes only, and the small ones say little for what they
cost: modulo 3 two residues in three are admissible on a random curve, so that AND removes a third
of the candidates where a large prime's removes half. Modulo 9 the admissible residues are those
with `f(x)` a square mod 9, and the squares are 4 of the 9 residues: over the thousand random test
curves that carries twice what 3 does, for a table of nine rows, and 25, 27 and 49 do the like for
5, 3 and 7. The row of a product of coprime moduli is the AND of its factors' rows, so one AND per
word carries them all: 45 = 9*5 carries what two large primes do. Every odd composite up to 64 is
now a candidate in the ranking that chooses the primes, with the product of its factors' densities
and its own table cost, and moduli sharing a prime exclude one another; the third stage keeps to
primes, and below a height bound of a few hundred no table pays and nothing is offered (at 1000
the modulus 9 enters on a third of the random curves, and the run does not change measurably). The
first phase of a random curve at height 16383 uses 10 moduli instead of 11.7 primes and does 14%
fewer ANDs. Worth 7% of `make test1` (6 to 8% in two measurements), 11% of `make testhigh` and 3 to
4% of `make testhighmany`, and nothing measurable on `make test1many`, where the small moduli are nearly
silent. Larger moduli save more instructions
and lose time: a row of 220 bit arrays does not stay in the first-level cache, and the cost model
prices an AND the same whatever the modulus, which is on the list.

The rule that ends the first phase weighs what a modulus costs and what its survivors cost. It
used to take moduli while the expected survivors per word exceeded a fitted constant, which stands
for the ratio of two costs -- one more AND per word against what a surviving word costs from there
on -- both frozen at the values of the height bound the constant was fitted at. Now the modulus's
side carries its tables and its per-denominator entries spread over the words of the run, the
same per-word cost the ranking already charges, and the survivor's side is weighted by what one
costs with the second phase the run will have: in a long run eleven cheap ANDs kill it, in a run of
a few thousand words there is no second phase and every survivor reaches the extraction, which
costs several times more. Both factors are near one where the constant is fitted, so nothing was
retuned or added. At height 200 the phase takes 7 moduli instead of 13 (a sweep of fixed counts
puts the optimum at 7), the tables fall from a fifth of the run to a twentieth, and the run is
25% shorter; 21% at height 100, 15% at 1000 (13% on the point-rich curves), 1% at 4000; at
16383 and 200000 the counts and the four suites do not move (instructions within 0.5%). What is
left at height 200 is looking at the thirty primes, a quarter of the instructions, and the exact
checks.

A test on the denominators that had never run now does. When a prime `p` divides the leading
coefficient, the congruence modulo `p` says nothing about a denominator divisible by `p`, but the
valuation does: for each valuation the denominator can have at `p`, the program asks whether a single
term of `F(a,b)` sets the valuation of the whole, and rules the valuation out if that is odd, or
even with an even power of `b` and a non-square unit. Typically that forbids `p` itself, or `p^2`, or one particular valuation. The
code for this case existed but sat behind a guard it could never pass, and it was wrong besides --
it would have lost 23 points on 12 of the thousand test curves. The corrected test excludes something on 499 of
those curves and is worth 8.3% of `make test1`, 5.3% of `make testdegrees` and 2.4% of
`make test1many`, 8.8% of `make testhigh` and 2.1% of `make testhighmany`, more than the share of denominators it excludes, because those
are the denominators for which that prime could not sieve anything.

How many primes the third stage uses is settled by weighing what a prime removes against what it
costs, and what it removes is exact tests — so it depends on what one exact test costs, which is a
property of the curve and not only of the machine. The test evaluates a binary form of the given
degree and takes an integer square root, so it costs about 15 cycles per degree plus a step for
each further limb the square root needs, and over the range from degree 3 with small coefficients
to degree 14 with 400-bit ones it varies by a factor of ten. The program estimates it from the
degree, the largest coefficient and the height bound before it sieves anything, and `-W` overrides
the estimate. Measured against the old assumption that every curve costs the same, this changes no
running time by as much as half a per cent on any test set, because the rule it feeds sits at its
own optimum; what it buys is that the constants mean what they are documented to mean for curves
outside the ones they were tuned on.

Four machine-dependent constants govern the choice. They can be set with the `-r`, `-R`, `-U` and
`-C` options; reasonable values are compiled in, and `make tune` will measure better ones for your
machine if you want it to: it minimises the sum of the times of `make test1` (random curves) and
`make test1many` (curves with many points), and writes what it finds to `tuning.mk`, which the
Makefile picks up. No source file is touched, and deleting that file restores the compiled-in
values.

It is a separate step rather than part of `make all`, because it takes several minutes (three to
eight, depending on the register width), wants an otherwise idle machine, and must not be run under
`make -j`. Timing this reliably is the hard part — the cost surface is flat while a laptop under
load slows by a quarter as it warms up — so each candidate is timed back to back with the current
settings and only the ratio is kept, and nothing is written unless the winner is clearly better and
the machine measured consistently. A run that reports that nothing beat the current settings has
done its job. See the documentation for the details.

There is also a variant that uses 512-bit AVX registers, which needs a CPU with AVX512F capability;
see the documentation for how to enable it. One caveat: it has not been tested completely, since I
have no such CPU available (it has only been checked indirectly, by having the compiler express the
512-bit operations through narrower ones). It was also unclear whether it would be any faster at all,
since the sieving loop is limited by how fast it loads bit arrays from the first-level cache, and by how
much of that cache the sieve tables occupy, rather than by arithmetic; but one data point
kindly provided by [Drew Sutherland](https://github.com/andrewvsutherland), for a curve with many
rational points on a Zen 5 CPU, shows a speedup of 13-14% over the 256-bit version. It is still
advisable to run `make test` and compare the timings on your own machine.

A review of the code in September 2026 found five bugs, fixed here and (the older ones) in 2.2.4:
degree 1 crashed the Sturm sequence code; the points at infinity were lost when the positivity
region of `f` missed the search domain; an array of numerator patterns was read before it was
filled when no denominator class survived mod 16 (as the modulus then was); the run-length estimate collapsed when only
even denominators survived, switching two sieving stages off; and the library wrote into input
fields of `ratpoints_args`, so a program that fills the structure once and loops over curves ran
every curve after the first with the first one's forbidden divisors, primes and search region. The
input fields now come back as they went in, `sp1_used`, `sp2_used` and `sp3_used` report what was
used, and two test targets cover this: `make test3` runs one invocation per bug against
`testbase3`, `make test1once` runs the test curves with the fields set once.

The code now assumes a 64-bit `long` and refuses to compile without one, so the 32-bit
architectures are no longer supported; version 2.2.4 is the last one that accommodates a 32-bit
`long`.

`make test` now includes `make testdegrees`, a hundred curves of degree 3, 4, 7 and 8. Everything
else in the package is degree 6, and every constant was tuned there, so nothing was known about the
rest; one regression had already slipped through because of it. The points on those curves were
checked once against a brute-force search over every coprime pair within a small height bound,
written independently of the sieve.

There is now [ratpoints-gpu](https://github.com/wgxli/ratpoints-gpu) by [Samuel Li](https://github.com/wgxli), which has similar functionality, but does the sieving on a GPU, which makes it much faster. His code is independent from what is in this repository.