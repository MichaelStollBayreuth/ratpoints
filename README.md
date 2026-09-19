# ratpoints

## The ratpoints library and command-line program

This is a program that uses an optimized quadratic sieve algorithm in order
to find rational points on hyperelliptic curves.

The program is distributed under the GNU GPL, version 2 (or later).

Read the [documentation](https://www.mathe2.uni-bayreuth.de/stoll/programs/ratpoints-doc-2.2.pdf).

The current version is ratpoints-2.2.4 from September 19, 2026, which, compared to version 2.2.3,
fixes a number of bugs found in a review of the code and when running the new test suite (see below),
some of which could lead to incorrect results in certain cases.
See the change log in the documentation pdf for more details.

This version can use 256-bit AVX registers (when available) and has been optimized further,
so that it now runs considerably faster than ratpoints-2.1.3.

On a CPU with AVX512F capability, the program can use 512-bit AVX registers as well;
see the documentation for how to enable it. One caveat: it has not been tested completely, since I
have no such CPU available (it has only been checked indirectly, by having the compiler express the
512-bit operations through narrower ones). It was also unclear whether it would be any faster at all,
since the sieving loop is limited by how fast it loads bit arrays from the first-level cache, and by how
much of that cache the sieve tables occupy, rather than by arithmetic;, so it is still
advisable to run `make test` and compare the timings on your own machine.

A test suite originally written for the next version has been adapted to this one: `make test4` runs a few
hundred invocations of `ratpoints` chosen so that every branch of the code runs at least once
against a reference that was checked by brute force, independently of the
sieve, and `make testapi` covers the library interface that the program cannot reach;
`make test4configs` runs the suite on the library built with the other compile-time switches,
and `make coverage` measures what the tests execute.

A test that fails now fails its `make` target, and `make test` runs all of its suites whatever the
earlier ones did and fails at the end if any of them failed, so that a script can tell whether the
build passed. The scripts behind `make test4configs` and `make coverage` exit with status 1 when an
output differs from its reference and with 2 when a build did not succeed.

### AI Declaration

The bugs this version fixes compared to version 2.2.2 were found and the changes making the
512-bit version work (in the emulation that my machine allows) were made by Claude Code (Anthropic),
in multiple sessions supervised by myself.

### ratpoints-gpu

There is now [ratpoints-gpu](https://github.com/wgxli/ratpoints-gpu) by [Samuel Li](https://github.com/wgxli), which has similar functionality, but does the sieving on a GPU, which makes it much faster. His code is independent from what is in this repository.