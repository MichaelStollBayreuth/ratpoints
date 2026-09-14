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
starts from the same pattern of admissible numerators modulo 16, and that pattern used to be
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
filled when no denominator class survived mod 16; the run-length estimate collapsed when only
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