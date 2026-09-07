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

Two machine-dependent constants govern the choice. They can be set with the `-r` and `-R` options,
reasonable values are compiled in, and `make tune` will measure better ones for your machine if
you want it to: it minimises the sum of the times of `make test1` (random curves) and
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
since the sieving loop is limited by memory bandwidth rather than by arithmetic; but one data point
kindly provided by [Drew Sutherland](https://github.com/andrewvsutherland), for a curve with many
rational points on a Zen 5 CPU, shows a speedup of 13-14% over the 256-bit version. It is still
advisable to run `make test` and compare the timings on your own machine.

There is now [ratpoints-gpu](https://github.com/wgxli/ratpoints-gpu) by [Samuel Li](https://github.com/wgxli), which has similar functionality, but does the sieving on a GPU, which makes it much faster. His code is independent from what is in this repository.