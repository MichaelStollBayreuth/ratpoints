# ratpoints

## The ratpoints library and command-line program

This is a program that uses an optimized quadratic sieve algorithm in order
to find rational points on hyperelliptic curves.

The program is distributed under the GNU GPL, version 2 (or later).

Read the [documentation](https://www.mathe2.uni-bayreuth.de/stoll/programs/ratpoints-doc-2.2.pdf).

The current version is ratpoints-2.2.3 from September 6, 2026. This version can use 256-bit AVX
registers and has been optimized further, so that it now runs considerably faster than ratpoints-2.1.3.

There is also a variant that uses 512-bit AVX registers, which needs a CPU with AVX512F capability;
see the documentation for how to enable it. Two caveats: it has not been tested completely, since I
have no such CPU available (it has only been checked indirectly, by having the compiler express the
512-bit operations through narrower ones), and it may well turn out to be *slower* than the 256-bit
version, because the sieving loop is limited by memory bandwidth rather than by arithmetic. Please
run `make test` and compare the timings before using it.

There is now [ratpoints-gnu](https://github.com/wgxli/ratpoints-gpu) by [Samuel Li](https://github.com/wgxli), which has similar functionality, but does the sieving on a GPU, which makes it much faster. His code is independent from what is in this repository.