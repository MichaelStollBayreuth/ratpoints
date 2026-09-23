# ratpoints

## The ratpoints library and command-line program

This is a program that uses an optimized quadratic sieve algorithm in order
to find rational points on hyperelliptic curves.

The program is distributed under the GNU GPL, version 2 (or later).

Read the [documentation](https://www.mathe2.uni-bayreuth.de/stoll/programs/ratpoints-doc-3.1.pdf).

The current version is **ratpoints-3.1.0** from XXXX XX, 2026. New in it:

* **The sieve can run on several threads.** `ratpoints ... -t n` sieves on
  n threads (`-t 0`: as many as the machine has processors), and programs
  using the library ask for them through the field `num_threads`. What is
  printed, and what the library's callback receives, is the same in the
  same order whatever the number of threads, and the callback is called
  from the calling thread only. Threads pay for long runs: on my laptop a
  search with height bound 400000 takes 1.9 s on one thread and 0.4 s on
  ten. `make THREADS=0` builds without threads. Programs using the library
  must be recompiled, since `ratpoints_args` has gained fields.
* **The 512-bit build has been verified** on AVX-512 hardware; see below.
* `make testthreads` checks that the output does not depend on the threads,
  and `make tsan` runs the threaded tests under ThreadSanitizer.

Version 3.0.0 from September 21, 2026, compared to version 2.2.4, works
out for itself what it used to be told, sieves between two and ten times faster, and comes with a test
suite that exercises nearly every line of the code. It needs a 64-bit `long`; **version 2.2.4 is the last one for 32-bit machines**.

The main improvements of version 3.0.0 over version 2.2.x:

* **The sieving parameters are chosen from the curve.**
  How many moduli each sieving stage uses, and which ones, follows from the
  densities of admissible residues modulo the small primes, from how
  long the run is going to be and from what each modulus costs.
  The five machine-dependent constants behind this can be set on the command line
  (`-r`, `-R`, `-U`, `-C`, `-Q`), and `make tune` measures them for your machine
  and writes them to `tuning.mk`.
* **A third sieving stage**
  tests the few candidates that survive the two bit-array sieving stages against
  further primes, one candidate at a time and without sieving tables.
  How many primes it uses is decided per curve, from an estimate of what one
  exact test costs on that curve, and corrected while
  the program runs from what the sieve is actually finding.
* **Composite moduli.**
  The odd composite numbers below 64 sieve alongside the primes;
  this improves the efficiency of the first two sieving stages.
* **Finer 2-adic information.**
  The numerator patterns are taken modulo 64 instead of 16,
  and the bit arrays of each denominator class pack the numerators
  with the largest stride a power of two allows, not merely by parity.
* **More denominators are excluded.**
  A test at the primes dividing the leading coefficient, derived from the
  Newton polygon, now runs; the forbidden divisors cover every excluded prime
  up to the square root of the height bound; the Jacobi symbol test is a product
  of table-driven Legendre symbols, and the tests on `b mod 64` are applied
  to 64 denominators at a time.
* **Leaner sieving loops.**
  Beginning and end of the first sieving stage are sped up;
  the sieving tables are built four times faster, and the memory reserved
  for them follows the number of primes in use.
* **A test suite.**
  `make test` runs random and point-rich curves, curves of other degrees, one
  invocation per bug fixed, a suite of a few hundred invocations chosen so
  that nearly every line of the code runs at least once (its reference checked
  by brute force, independently of the sieve) and the library interface,
  and it fails when a test fails; `make test4configs` repeats the suite on the
  other compile-time configurations, and `make coverage` measures what
  the tests execute.

What it comes to, against the released 2.2.4 on the same curves and the same laptop
(cycles pinned to one core, medians of three interleaved rounds;
2.2.4 at its own defaults; both find the same points on every suite):

| suite                                    | 2.2.4, Gcycles | 3.0.0 | ratio | faster by |
|------------------------------------------|------:|------:|------:|------:|
| `make test1` (random curves, 16383)      |   8.45 |   2.89 | 0.34 |  2.9x |
| `make test1many` (point-rich, 16383)     |  11.27 |   4.80 | 0.43 |  2.3x |
| `make testhigh` (random, 200000)         | 412.7  | 190.8  | 0.46 |  2.2x |
| `make testhighmany` (point-rich, 200000) | 764.4  | 138.3  | 0.18 |  5.5x |
| random curves at height 4000             |   3.06 |   0.54 | 0.18 |  5.6x |
| random curves at height 1000             |   2.06 |   0.21 | 0.10 | 10.1x |
| point-rich curves at height 1000         |   0.44 |   0.15 | 0.34 |  3.0x |
| random curves at height 200              |   0.96 |   0.12 | 0.13 |  8.0x |

Instructions fell to between 0.16 and 0.47 of 2.2.4's and mispredicted branches
to between 0.03 and 0.40.

On a CPU with AVX512F capability, the program can use **512-bit** AVX registers
as well; see the documentation for how to enable it. I have no such CPU and can
test that variant only in the emulation the compiler allows;
[Drew Sutherland](https://github.com/andrewvsutherland) has kindly run the
test suite of version 3.0.0 on an AMD Ryzen 9 9950X3D (Zen 5), where every
test passes, and timed the two variants against each other, both with the
constants that `make tune` finds for that machine: the 512-bit version is 14%
faster on a long run (the curve of `make timing`, height bound 400000), but
8% slower on the point-rich test curves at height 16383, whose sieve tables
are twice as large to build, and equal on the random ones. Whether it pays
thus depends on the machine and on the kind of run; run `make test` and
`make tune` and compare the timings on your own machine. (`make tune` alone
was worth 6-7% on that machine at either register width: the constants
compiled in were measured on my laptop.)

### AI Declaration

The bugs fixed in versions 2.2.3 and 2.2.4 were found
and the changes making the 512-bit version work (in the emulation that
my machine allows) were made by Claude Code (Anthropic),
in multiple sessions supervised by myself.

The same holds for the improvements of version 3.0.0 over 2.2.4 described above,
for the test suite and for the threads of version 3.1.0: they were worked out,
implemented and measured by Claude Code in many such sessions, with the
decisions taken by myself.

### ratpoints-gpu

There is now [ratpoints-gpu](https://github.com/wgxli/ratpoints-gpu)
by [Samuel Li](https://github.com/wgxli), which has similar functionality,
but does the sieving on a GPU, which makes it much faster.
His code is independent from what is in this repository.
