# ratpoints

## The ratpoints library and command-line program

This is a program that uses an optimized quadratic sieve algorithm in order
to find rational points on hyperelliptic curves.

The program is distributed under the GNU GPL, version 2 (or later).

Read the [documentation](https://www.mathe2.uni-bayreuth.de/stoll/programs/ratpoints-doc-2.2.pdf).

The current version is ratpoints-2.2.4 from September 13, 2026, which fixes bugs found in a
review of the code (see the change log in the documentation): a crash for polynomials of degree 1,
lost points at infinity when the positivity region of `f` misses the search domain, an
uninitialised read that could print wrong points under `-j -F 0`, and input fields of
`ratpoints_args` that the library overwrote, so that a program filling the structure once and
looping over curves ran every curve after the first with the first one's forbidden divisors,
primes and search region. This version can use 256-bit AVX
registers and has been optimized further, so that it now runs considerably faster than ratpoints-2.1.3.

The test suite written for the next version has been adapted to this one: `make test4` runs a few
hundred invocations of `ratpoints` chosen so that every branch of the code runs at least once
(every degree and shape of curve, every option, restricted ranges, height bounds from 1 to 2^63-1,
every error message) against a reference that was checked by brute force, independently of the
sieve, and `make testapi` covers the library interface that the program cannot reach; `make
test4configs` runs the suite on the library built with the other compile-time switches, and `make
coverage` measures what the tests execute. Writing it found a few small things here, all fixed: the
library tested whether the leading coefficient is a square before it checked the coefficient
pointer for `NULL`; an empty coefficient string made the program report "Bug no. 1" instead of
refusing the input; and at the very top of the range of a `long` a search interval reaching the
last few hundred numerators below 2^63 was dropped without a word, while the loop over the
denominators never ended when the bound was 2^63-1. Two more followed: the loops over the square
denominators started at 1 whatever the lower bound `-dl` said, which at the top of the range meant
three thousand million useless squares before the first one in the range, and now start at the
first square in the range; and the test on the valuation of such a denominator at a prime dividing
the leading coefficient looked at its low 32 bits only (`abs` where `labs` was meant), so that
from 2^31 on it could exclude a denominator that carries a point.

A test that fails now fails its `make` target, and `make test` runs all of its suites whatever the
earlier ones did and fails at the end if any of them failed, so that a script can tell whether the
build passed. The scripts behind `make test4configs` and `make coverage` exit with status 1 when an
output differs from its reference and with 2 when a build did not succeed.

There is also a variant that uses 512-bit AVX registers, which needs a CPU with AVX512F capability;
see the documentation for how to enable it. One caveat: it has not been tested completely, since I
have no such CPU available (it has only been checked indirectly, by having the compiler express the
512-bit operations through narrower ones). It was also unclear whether it would be any faster at all,
since the sieving loop is limited by how fast it loads bit arrays from the first-level cache, and by how
much of that cache the sieve tables occupy, rather than by arithmetic; but one data point
kindly provided by [Drew Sutherland](https://github.com/andrewvsutherland), for a curve with many
rational points on a Zen 5 CPU, shows a speedup of 13-14% over the 256-bit version. It is still
advisable to run `make test` and compare the timings on your own machine.

There is now [ratpoints-gpu](https://github.com/wgxli/ratpoints-gpu) by [Samuel Li](https://github.com/wgxli), which has similar functionality, but does the sieving on a GPU, which makes it much faster. His code is independent from what is in this repository.