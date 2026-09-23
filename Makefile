#   ratpoints-3.1.0
#    - A program to find rational points on hyperelliptic curves
#   Copyright (C) 2008, 2009, 2022, 2023, 2026  Michael Stoll
#
#   This program is free software: you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation, either version 2 of
#   the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of version 2 of the GNU General
#   Public License along with this program.
#   If not, see <http://www.gnu.org/licenses/>.
#
#
#   Makefile
#
#   Michael Stoll, September 21, 2009; January 7, 2022; September 6-21, 2026
#   with changes by Bill Allombert, December 29, 2021

# The main targets are
#   all         the library, the program and this documentation
#   test        test1, test1once, test1many, testdegrees, test2, test3,
#               test4, testapi and timing (see below)
#   test1once   the curves of test1 with the input fields of args set once
#               before the loop: the library must leave them alone
#   test3       the invocations of test3.sh, regression tests for bugs
#               that were fixed, against testbase3
#   test1       1000 random genus 2 curves and eight chosen ones, checked
#               against testbase
#   test1many   curves with many rational points, against testbase-many;
#               the two cover the two regimes that behave differently, and
#               both should be used when judging a change to the sieve
#   testhigh    the curves of test1 at a height bound of TESTHEIGHT below,
#               where the sieve and not the set-up decides the running time
#   testhighmany  the thirty curves of test1many with the most points, at the
#               same height; testhigh and testhighmany take about a minute each
#   test4       the invocations of test4.sh, the suite that exercises nearly
#               every line of the code -- every degree and shape of curve,
#               every option, restricted ranges, heights from 1 to 2^63-1,
#               every error message -- against testbase4; a few seconds.
#               verify-test4.py checks that reference by brute force
#   testapi     rpapi, the tests of the library interface that the program
#               cannot reach (the argument checks and their error codes,
#               the callback, the flags), against testbase-api
#   test4configs  test4 again on the library built with the other
#               compile-time switches: register widths 64, 128, 512,
#               RATPOINTS_CHUNK=1, USE_LONG_IN_PHASE_2, PRIME_SIZE 7 and
#               composite moduli off and unbounded (test4-configs.sh);
#               builds each in a directory of its own, about a minute
#   coverage    build the sources instrumented for gcov in build-coverage/,
#               run the tests through that build and print how much of each
#               source file they executed (coverage.sh); a minute or so
#   tune        measure the four machine-dependent constants that decide how
#               many primes each sieving stage uses and which, and write them
#               to tuning.mk (see tune.sh); takes several minutes, wants an idle
#               machine, and must not be run under -j
#   bench_init  correctness check and benchmark for the sieve table set-up
#   bench_check benchmark for the exact check, and how to remeasure the four
#               RATPOINTS_CHECK_* constants
#   clean       remove the intermediate files; distclean also the executables
#               and any tuning.mk

# The program sieves with the odd primes below 2^PRIME_SIZE.  The default is 8
# rather than 7 not because 53 primes are normally used -- the choice below
# still starts from RATPOINTS_DEFAULT_NUM_PRIMES of them -- but so that the few
# primes past that are there for the curves that run out, which are the ones
# with very many rational points; on those it is worth up to a factor of two.
# Having the larger table costs nothing when it is not used.
PRIME_SIZE = 8
VERSION = 3.1.0

# The height bound for "make testhigh" and "make testhighmany"; see there.
TESTHEIGHT = 200000

CC = gcc
RM = rm -f
INSTALL = cp

INSTALL_DIR = /usr/local

# -funswitch-loops, the -O3 optimisation that compiles a loop with an
#  invariant test inside it twice, once per outcome, is deliberately not
#  used.  It removes 1 to 3% of the instructions -- the "which reduction" test
#  in the per-call start loop and in the second-stage loop of sift.c, and the
#  tests on the numerator packing in the per-denominator loops of
#  find_points.c -- and gains nothing in cycles: a wash on the two suites at
#  height 16383 and 0 to 2% slower on the two at 200000, measured at three
#  code alignments.  Instruction counts predict a gain; cycles are what count.
# Threads: the library sieves on several threads when asked to (the field
# num_threads of ratpoints_args, the option -t), with POSIX threads.
# "make THREADS=0" builds without them, for a system that has none; the
# field and the option are then ignored.
ifeq (${THREADS},0)
THREADFLAGS = -DRATPOINTS_NO_THREADS
else
THREADFLAGS = -pthread
endif

CCFLAGS0 = -Wall -O2 -fomit-frame-pointer -DRATPOINTS_MAX_BITS_IN_PRIME=${PRIME_SIZE} ${THREADFLAGS}
# For gcc on Apple, may have to add '-fnested-functions' to CCFLAGS0.
# Add "-DUSE_LONG_IN_PHASE_2" to sieve the survivors of the first stage one
#  64-bit word at a time instead of a whole bit-array at a time. The first
#  stage and the scan for survivors stay at the full register width either
#  way. It is a wash to 4% slower at 128 and 256 bits, and the reason is worth
#  knowing: the number of AND steps is the same either way, because the other
#  words of a surviving bit-array were already zero, and a narrow and a wide
#  read of the same table come from one cache line. At 64 bits the two are the
#  same code.
# Add "-DRP_PHASE_TIMING" to have sift.c time the two stages of the sieve
#  separately and write a report to stderr when the program exits; add
#  "-DRP_PHASE_COUNTS" as well to count the bit-arrays surviving stage 1, or
#  "-DRP_STOP_AFTER=<n>" to cut the pipeline short after a chosen stage, so
#  that the cost of a stage can be had as a difference of two runs.
#  These are development aids; see the comment at the top of sift.c.
# Add "-DRP_INIT_BRANCH" or "-DRP_INIT_ONEWAY" to select, one at a time, the
#  two simpler forms of the sieve table set-up: testing is_f_square with a branch
#  rather than shifting the value into place, and stepping the residue one row at
#  a time rather than four.  Between them they cost 13% of "make test1" on the
#  machine this was written on; the flags are there to measure that again on
#  yours.  "-DRP_INIT_NOACC" and "-DRP_INIT_NOREP" leave out a stage of
#  sieve_init altogether, so that its share can be had as a difference of two
#  runs; they produce wrong tables, and bench_init will say so.
# Add "-DRP_MULMOD_DIVIDE" to reduce modulo a sieving prime by dividing rather
#  than by multiplying with the reciprocal, which is what the third sieving
#  stage, the start-of-sieve computation, the row look-up of the second
#  stage, the reduction of the denominator modulo each sieving prime and the
#  Jacobi symbol test cost without that.  "-DRP_MOD_CHOICE" instead builds both forms of
#  the start-of-sieve computation into one binary, selected by the
#  environment variable RP_MOD_MUL, so that they can be timed against each
#  other without the code-alignment difference two builds would bring (the
#  multiplication is worth 2 to 3%).  "-DRP_MOD_COUNTS" reports how often the
#  chain of conditional subtractions actually divides.
# When comparing two builds whose *source* differs, be aware that where gcc
#  happens to place the hot loops is worth about 10% here, reproducibly, so
#  repeating the runs will not reveal it. Rebuild both with, say,
#  "-falign-loops=32" and "-falign-loops=64" and check that the difference
#  survives. Those flags are not set by default because no value is reliably
#  better: measured against the default over three register widths and two
#  curves, -falign-loops=32 lands between 0.992 and 1.010 with no consistent
#  sign, and 64 is worse.
# Add "-DRATPOINTS_CHUNK=<n>" to force the use of 2 <= n <= 16 registers
#  in stage 1 of sieving. For n=1, the loop is left to the compiler.
#  If SSE/AVX registers are used and this is not set, 16 registers will be used.
#  In some cases, using n=8 may be faster
#  (e.g., Intel(R) Xeon(R) CPU E3-1220 V2 with -DUSE_AVX -mavx).

# The following uses 64-bit registers, i.e., plain unsigned longs.
# This works on any machine the library builds on at all (a 64-bit long is
# what the code needs in any case; see rp-private.h).
CCFLAGS64 =
# The following uses 128-bit registers. In spite of its name, USE_AVX128 needs
# only SSE2, which every x86-64 machine has, so this is as portable as the
# USE_SSE variant below and never slower: its test for an empty register
# compares the whole register against zero instead of extracting both halves
# into general registers. That is worth 30% of the second sieving stage, hence
# 9% of a long run, but only about 3% of "make test1", whose height bounds are
# small enough that the set-up is a large share of the time.
CCFLAGS128 = -DUSE_AVX128
# The older variant, using the SSE intrinsics directly. Kept for comparison.
CCFLAGS128s = -DUSE_SSE
# The following uses 256-bit AVX registers.
# Change "-mavx2" to "-mavx" when your processor has AVX, but no AVX2.
# This may be a bit slower compared to using AVX2 instructions.
CCFLAGS256 = -DUSE_AVX -mavx2
# To use 512-bit AVX registers, use the following
# (if your processor has AVX512f capability).
# Whether it is faster than the 256-bit build depends on the machine and on
# the run (on a Zen 5: 14% faster on a long run, 8% slower on point-rich
# curves at a small height bound; see the manual), so time the two against
# each other, and run "make tune" after switching.  Leaving out "-mavx512f"
# makes gcc emulate the 64-byte vectors with narrower ones; that is slower,
# but it runs anywhere and exercises the same code path.  (gcc then warns
# "AVX512F vector argument without AVX512F enabled changes the ABI"; this is
# harmless here, since the whole library is built with the same flags.)
CCFLAGS512 = -DUSE_AVX512 -mavx512f

# This will be the default. Change as appropriate.
# CCFLAGS1 = ${CCFLAGS128}
CCFLAGS1 = ${CCFLAGS256}

# Machine-dependent tuning of the constants that decide how many primes each
# sieving stage uses and which.  "make tune" writes tuning.mk; without it the
# values compiled into ratpoints.h are used.  This needs GNU make for the
# conditionals; if yours is not GNU make, delete the block and either leave
# tuning.mk out or trust it unconditionally.
-include tuning.mk

TUNE_CONFIG = ${CCFLAGS1} / PRIME_SIZE=${PRIME_SIZE}
ifdef TUNED_FOR
ifneq (${TUNED_FOR},${TUNE_CONFIG})
$(warning tuning.mk was measured for "${TUNED_FOR}", but this build is)
$(warning "${TUNE_CONFIG}" -- ignoring it; run "make tune" again)
TUNEFLAGS =
endif
endif

# CCFLAGS_H is the part that decides what the generated headers must contain
# (the register width and PRIME_SIZE); the tuning flags do not affect them.
CCFLAGS_H = ${CCFLAGS0} ${CCFLAGS1}
CCFLAGS_0 = ${CCFLAGS_H} ${TUNEFLAGS}

# Further compiler flags for linking
CCFLAGS2 = -lgmp -lgcc -lc
CCFLAGS3 = -L. -lratpoints -lm
# Further flags that can be added by calling "make CCFLAGS=..."
CCFLAGS =

# Files that make up the distribution
DISTFILES = Makefile ratpoints.h rp-private.h primes.h \
            gen_find_points_h.c gen_init_sieve_h.c \
            sift.c init.c sturm.c find_points.c \
            main.c rptest.c testdata.h testbase ratpoints-doc-3.1.tex \
            README.md gpl-2.0.txt testbase2 testdata-many.h testbase-many \
            testdata-high-many.h testbase-high-many \
            testdata-degrees.h testbase-degrees \
            test3.sh testbase3 \
            test4.sh testbase4 verify-test4.py test4-configs.sh test4-threads.sh tsan.sh \
            rpapi.c testbase-api coverage.sh \
            bench_init.c bench_check.c tune.sh

# Temporary files that are generated during build and test
# and can be removed afterwards
TEMPFILES = sift.o init.o sturm.o find_points.o \
            sift.s sift.i init.s find_points.h init_sieve.h \
            gen_find_points_h gen_init_sieve_h \
            rptest.out rptest-many.out rptest-high.out \
            rptest-high-many.out rptest-degrees.out config.stamp build.stamp \
            sift-debug.o find_points-debug.o main.o test2.out \
            test3.out rptest-once.out rptest-once2.out rptest-once3.out \
            test4.out testapi.out test4-*.out \
            rptest-t2.out rptest-t5.out rptest-many-t3.out \
            rptest-degrees-t4.out test3-t3.out rptest-threads-*.out

# Executables and library produced when building
TARGETFILES = ratpoints libratpoints.a rptest rptest-many rptest-high-many \
              rptest-degrees rpapi ratpoints-debug \
              bench_init bench_check ratpoints-doc-3.1.pdf

FAILED = "Test failed!"
# what a test does when its output differs from the reference: print the
# message and fail the recipe, so that the exit status of make reports it
FAIL = { echo ${FAILED}; false; }

all: ratpoints libratpoints.a doc

doc: ratpoints-doc-3.1.pdf

# The suites "make test" runs, each a target below.  A test whose output
# differs from its reference prints "Test failed!" and fails its target
# (make then exits with status 2 and names the target); "make test" runs
# every suite whatever the earlier ones did and fails at the end if any of
# them failed.
TESTS = test1 test1once test1many testdegrees test2 test3 test4 testapi timing

.PHONY: test
test:
	@status=0; for t in ${TESTS}; do \
	   ${MAKE} --no-print-directory $$t || status=1; done; \
	 exit $$status

# Measure good values for the four machine-dependent constants and write them
# to tuning.mk; see tune.sh.  Deliberately not part of "make all": it takes a
# couple of minutes, wants an otherwise idle machine, and must not run under
# "make -j".  Run "make all" afterwards to rebuild with the result.
.PHONY: tune
tune: rptest rptest-many
	@TUNE_CONFIG='${TUNE_CONFIG}' ./tune.sh

# The same, measured on the large-height suites instead (see testhigh and
# testhighmany below).  Which one to use depends on the runs that matter: at
# the 16383 of "make test" two fifths of the time on the random curves goes
# into work other than the two sieving stages -- choosing the primes, the
# set-up per denominator, the reductions of the denominators, the exact checks
# -- on which the threshold and the offset have little or no effect (the sieve
# tables alone are 4%, and the cost of building a table is the one constant
# that bears on them), so a short-run tuning judges them partly on work they
# do not touch.
#
# One timing here is two minutes against three seconds there, so this does not
# sweep the whole ladder of candidates.  It starts from the settings in force
# -- which is to say from what "make tune" found, since it reads the same
# tuning.mk -- and asks only whether a factor of two in the threshold, in the
# run length, in the table cost or in the third stage's per-denominator cost
# either way, or two more or fewer primes in the second stage, is better.
# That is eighteen settings a round rather than twenty-nine, and it rests on
# the two regimes not wanting wildly different
# values.  If it moves a value, run it again from there.
#
# Running this after "make tune" is what pins RATPOINTS_SP2_U0, the third
# constant: the offset it scales is the number of extra primes an
# arbitrarily long run wants, and only a pair of tunings at very different
# run lengths can separate that from the length at which the scaling bites.
# Expect a couple of hours at the default three rounds;
# "ROUNDS=1 make tunehigh" is the short version.
.PHONY: tunehigh
tunehigh: rptest rptest-high-many
	@TUNE_CONFIG='${TUNE_CONFIG}' \
	 TUNE_TESTS='./rptest:testbase ./rptest-high-many:testbase-high-many' \
	 TUNE_HEIGHT='${TESTHEIGHT}' \
	 R_FACTORS='0.5 2' E_DELTAS='-2 0 2' U_FACTORS='0.5 2' C_FACTORS='0.5 2' \
	 Q_FACTORS='0.5 2' \
	 ./tune.sh

# The timed test targets time their runs with the shell's "time".  That is a
# shell built-in which dash, /bin/sh on Debian and Ubuntu, does not have, so
# those targets run under bash; the build itself runs under whatever /bin/sh
# is.  Override on the command line if bash lives elsewhere on your system.
test1 test1many testhigh testhighmany testdegrees test2 timing: SHELL = /bin/bash

# Run ratpoints on a set of 1008 test cases -- 1000 random genus 2 curves and
# eight chosen to reach the test on the denominators at a prime dividing the
# leading coefficient, see testdata.h -- and check the output
test1: rptest testbase
	time ./rptest > rptest.out
	cmp -s testbase rptest.out || ${FAIL}

# The same, but for curves with many rational points, which sieve very
# differently: the small primes say little about them, so the first stage
# needs about sixteen moduli where a random curve needs ten, and five times
# as many candidates per denominator survive the second.  The two tests take
# about the same time and should be used together whenever the constants
# that choose sp1 and sp2 (see ratpoints.h) are retuned; a change that
# helps one regime can easily hurt the other.
test1many: rptest-many testbase-many
	time ./rptest-many > rptest-many.out
	cmp -s testbase-many rptest-many.out || ${FAIL}

# The same two regimes at a height bound of ${TESTHEIGHT} instead of the 16383
# that test1 and test1many use.  What a suite measures depends a good deal on
# that bound: the sieve tables are built lazily, once for each pair (prime,
# denominator mod that prime), and then reused, so their total cost is bounded
# by the primes and does not grow with the height, while the sifting does.  On
# the random curves the tables are 4% of "make test1" (22% when these suites
# were added) but 0.3% here, the whole set-up per denominator 8.5% against
# 1.3%, and the two sieving stages go from 57% to 90% of the run.  So these
# are the suites to judge a change to the sieving loops by, and the ones to
# point "make tune" at if the runs that matter are long ones.  Each takes
# about a minute.

# The curves of test1 at the larger height bound.  None of them has a
# rational point of height between 16383 and 200000, so the list of
# points is the same one and testbase is the reference for both -- which is
# itself worth checking.  Raise TESTHEIGHT and that stops being true; the new
# output has to be looked at and kept as a reference of its own.
testhigh: rptest testbase
	time ./rptest -h ${TESTHEIGHT} > rptest-high.out
	cmp -s testbase rptest-high.out || ${FAIL}

# The thirty curves of test1many that run out of sieving primes, at the same
# height; see testdata-high-many.h for what they are and why just those.
testhighmany: rptest-high-many testbase-high-many
	time ./rptest-high-many -h ${TESTHEIGHT} > rptest-high-many.out
	cmp -s testbase-high-many rptest-high-many.out || ${FAIL}

# Run ratpoints on the curve with the record number of known
# rational points, with a fairly large height bound,
# time it and compare with the expected output
# Curves of degree 3, 4, 7 and 8.  Everything else in the package is degree
# 6, and every constant in ratpoints.h was tuned there, so this is the suite
# that says whether anything is peculiar to genus 2.  It found one thing
# already: raising PRIME_SIZE from 7 to 8 had moved the division-free path in
# sieving_info from degree 8 down to degree 7.
#
# testbase-degrees was produced by this program, so the suite is a regression
# test and not an independent one.  What makes it trustworthy is that the
# points were checked once against a brute-force search over every coprime
# pair (a, b) within a small height bound, written independently of the
# sieve.
testdegrees: rptest-degrees testbase-degrees
	time ./rptest-degrees > rptest-degrees.out
	cmp -s testbase-degrees rptest-degrees.out || ${FAIL}

# The curves of test1 again, with the input fields of args set once before
# the loop instead of once per curve: the library must leave them alone
# (rptest prints a line whenever one has changed), and the points must be
# the same.  The second run, with a lower bound on the denominator, makes
# the library take its own decision not to reverse the polynomial, which
# must not show in the caller's flags; the third gives values that the
# library normalises, which must not show in the fields.
test1once: rptest testbase
	./rptest -O > rptest-once.out
	cmp -s testbase rptest-once.out || ${FAIL}
	./rptest -O -dl 2 -z > rptest-once2.out
	! grep -q changed rptest-once2.out || ${FAIL}
	./rptest -O -dl 0 -du 1000000 -S 100 > rptest-once3.out
	cmp -s testbase rptest-once3.out || ${FAIL}

# Regression tests for bugs that were fixed: a list of invocations of
# ratpoints in test3.sh, against testbase3.
test3: ratpoints testbase3 test3.sh
	./test3.sh > test3.out 2>&1
	cmp -s testbase3 test3.out || ${FAIL}

# The suite that exercises nearly every line (see test4.sh): a few hundred
# invocations of ratpoints, against testbase4.  The reference was checked
# by brute force, independently of the sieve: verify-test4.py does that,
# and can be run again whenever testbase4 changes.
test4: ratpoints testbase4 test4.sh
	./test4.sh > test4.out 2>&1
	cmp -s testbase4 test4.out || ${FAIL}

# The tests of the library interface (see rpapi.c), against testbase-api.
testapi: rpapi testbase-api
	./rpapi > testapi.out 2>&1
	cmp -s testbase-api testapi.out || ${FAIL}

# The suites with the sieve on several threads (the option -t): the same
# references, since what the program prints must not depend on the number
# of threads; then test4 and test1 on builds with fixed block lengths, so
# that the blocks of the loop over the denominators end elsewhere
# (test4-threads.sh; a couple of minutes).
.PHONY: testthreads
testthreads: rptest rptest-many rptest-degrees ratpoints testbase \
             testbase-many testbase-degrees testbase3 testbase4 \
             test3.sh test4.sh test4-threads.sh
	./rptest -t 2 > rptest-t2.out
	cmp -s testbase rptest-t2.out || ${FAIL}
	./rptest -t 5 > rptest-t5.out
	cmp -s testbase rptest-t5.out || ${FAIL}
	./rptest-many -t 3 > rptest-many-t3.out
	cmp -s testbase-many rptest-many-t3.out || ${FAIL}
	./rptest-degrees -t 4 > rptest-degrees-t4.out
	cmp -s testbase-degrees rptest-degrees-t4.out || ${FAIL}
	RPOPTS='-t 3' ./test3.sh > test3-t3.out 2>&1
	cmp -s testbase3 test3-t3.out || ${FAIL}
	RPOPTS='-t 3' ./test4.sh > test4-t3.out 2>&1
	cmp -s testbase4 test4-t3.out || ${FAIL}
	RPOPTS='-t 7' ./test4.sh > test4-t7.out 2>&1
	cmp -s testbase4 test4-t7.out || ${FAIL}
	./test4-threads.sh

# The threaded runs under ThreadSanitizer (tsan.sh: a build with
# -fsanitize=thread in build-tsan/, the drivers with threads, rpapi, and
# test4 with -t 3; every output must match its reference and the sanitizer
# must report nothing; a few minutes).
.PHONY: tsan
tsan: testbase testbase-many testbase-degrees testbase-api testbase4 tsan.sh
	./tsan.sh

# test4 on the library built with the other compile-time switches, each in
# a build directory of its own (see test4-configs.sh; the script exits with
# 1 when an output differs from the reference, with 2 when a build failed).
.PHONY: test4configs
test4configs: testbase4 test4.sh test4-configs.sh
	./test4-configs.sh

# How much of the code the tests execute, by gcov (see coverage.sh; the
# script exits with 1 when an output differs from its reference, with 2
# when the instrumented build failed).
.PHONY: coverage
coverage: coverage.sh
	./coverage.sh

test2: ratpoints testbase2
	time ./ratpoints '247747600 -985905640 567207969 2396040466 52485681 -470135160 82342800' 1000000 -n 30 -N 30 -p 30 -q > test2.out
	cmp -s testbase2 test2.out || ${FAIL}

# Time a call to ratpoints with a largish height parameter.
# This can be helpful to assess modifications to the sieving process.
timing: ratpoints
	time ./ratpoints '1 0 126 0 441' 400000 -q ${RPOPTS} > /dev/null

install-bin: ratpoints
	${INSTALL} ratpoints ${INSTALL_DIR}/bin/
	chmod 755 ${INSTALL_DIR}/bin/ratpoints

install-lib: ratpoints.h libratpoints.a
	${INSTALL} ratpoints.h ${INSTALL_DIR}/include/
	chmod 644 ${INSTALL_DIR}/include/ratpoints.h
	${INSTALL} libratpoints.a ${INSTALL_DIR}/lib/
	chmod 644 ${INSTALL_DIR}/lib/libratpoints.a

install: install-bin install-lib

# To generate the documentation, run pdflatex twice
# to get the cross-references right.
ratpoints-doc-3.1.pdf: ratpoints-doc-3.1.tex
	pdflatex ratpoints-doc-3.1.tex
	pdflatex ratpoints-doc-3.1.tex

dist: ${DISTFILES}
	mkdir -p ratpoints-${VERSION}
	cp ${DISTFILES} ratpoints-${VERSION}/
	tar --create --file=ratpoints-${VERSION}-`date --rfc-3339=date`.tar.gz \
	    --gzip --dereference ratpoints-${VERSION}
	rm -r ratpoints-${VERSION}

clean:
	${RM} ${TEMPFILES}
	${RM} -r build-coverage build-test4-* build-threads-* build-tsan

distclean: clean
	${RM} ${TARGETFILES} tuning.mk

debug: ratpoints-debug

libratpoints.a: sift.o init.o sturm.o find_points.o
	ar rs libratpoints.a sift.o init.o sturm.o find_points.o

ratpoints: libratpoints.a main.c ratpoints.h build.stamp
	${CC} main.c -o ratpoints ${CCFLAGS_0} ${CCFLAGS2} ${CCFLAGS3} ${CCFLAGS}

main.o: main.c ratpoints.h build.stamp
	${CC} main.c -c -o main.o ${CCFLAGS_0} -O3 ${CCFLAGS}

ratpoints-debug: sift-debug.o init.o sturm.o find_points-debug.o main.o build.stamp
	${CC} sift-debug.o init.o sturm.o find_points-debug.o main.o \
              -o ratpoints-debug ${CCFLAGS_0} ${CCFLAGS2} -lm ${CCFLAGS}

sift.o: sift.c ratpoints.h rp-private.h build.stamp
	${CC} sift.c -c -o sift.o ${CCFLAGS_0} -funroll-loops ${CCFLAGS}

sift-debug.o: sift.c ratpoints.h rp-private.h build.stamp
	${CC} sift.c -c -o sift-debug.o ${CCFLAGS_0} -funroll-loops -DDEBUG ${CCFLAGS}

sift.s: sift.c ratpoints.h rp-private.h build.stamp
	${CC} sift.c -S -o sift.s ${CCFLAGS_0} -funroll-loops ${CCFLAGS}

sift.i: sift.c ratpoints.h rp-private.h build.stamp
	${CC} sift.c -E -o sift.i ${CCFLAGS_0} -funroll-loops ${CCFLAGS}

init.o: init.c ratpoints.h rp-private.h init_sieve.h build.stamp
	${CC} init.c -c -o init.o ${CCFLAGS_0} -funroll-loops -O3 ${CCFLAGS}

init.s: init.c ratpoints.h rp-private.h init_sieve.h build.stamp
	${CC} init.c -S -o init.s ${CCFLAGS_0} -funroll-loops -O3 ${CCFLAGS}

sturm.o: sturm.c ratpoints.h rp-private.h build.stamp
	${CC} sturm.c -c -o sturm.o ${CCFLAGS_0} ${CCFLAGS}

find_points.o: find_points.c ratpoints.h rp-private.h primes.h find_points.h build.stamp
	${CC} find_points.c -c -o find_points.o ${CCFLAGS_0} ${CCFLAGS}

find_points-debug.o: find_points.c ratpoints.h rp-private.h primes.h find_points.h build.stamp
	${CC} find_points.c -c -o find_points-debug.o ${CCFLAGS_0} -DDEBUG ${CCFLAGS}

# Correctness check and benchmark for the sieve table set-up (see bench_init.c).
bench_init: libratpoints.a bench_init.c ratpoints.h rp-private.h primes.h build.stamp
	${CC} bench_init.c -o bench_init ${CCFLAGS_0} -O3 -funroll-loops \
              ${CCFLAGS3} ${CCFLAGS2} ${CCFLAGS}

# Benchmark for the exact check (see bench_check.c).  It links nothing from
# the library -- it reproduces the gmp calls the check makes -- but it reads
# the four RATPOINTS_CHECK_* constants out of ratpoints.h to print beside
# what it measures.
bench_check: bench_check.c ratpoints.h build.stamp
	${CC} bench_check.c -o bench_check ${CCFLAGS_0} ${CCFLAGS2} -lm ${CCFLAGS}


rptest: libratpoints.a rptest.c ratpoints.h testdata.h build.stamp
	${CC} rptest.c -o rptest ${CCFLAGS_0} ${CCFLAGS2} ${CCFLAGS3} ${CCFLAGS}

rptest-many: libratpoints.a rptest.c ratpoints.h testdata-many.h build.stamp
	${CC} rptest.c -o rptest-many -DRATPOINTS_TESTDATA='"testdata-many.h"' \
	      ${CCFLAGS_0} ${CCFLAGS2} ${CCFLAGS3} ${CCFLAGS}

rptest-high-many: libratpoints.a rptest.c ratpoints.h testdata-high-many.h \
                  build.stamp
	${CC} rptest.c -o rptest-high-many \
	      -DRATPOINTS_TESTDATA='"testdata-high-many.h"' \
	      ${CCFLAGS_0} ${CCFLAGS2} ${CCFLAGS3} ${CCFLAGS}

rptest-degrees: libratpoints.a rptest.c ratpoints.h testdata-degrees.h \
                build.stamp
	${CC} rptest.c -o rptest-degrees \
	      -DRATPOINTS_TESTDATA='"testdata-degrees.h"' \
	      ${CCFLAGS_0} ${CCFLAGS2} ${CCFLAGS3} ${CCFLAGS}

rpapi: libratpoints.a rpapi.c ratpoints.h build.stamp
	${CC} rpapi.c -o rpapi ${CCFLAGS_0} ${CCFLAGS2} ${CCFLAGS3} ${CCFLAGS}

# What is compiled depends on flags, which make cannot see by itself: change
# CCFLAGS1 or PRIME_SIZE, or retune, and every file must be built again even
# though no source file has changed.  Two stamp files record the flags in use,
# so that a change to either shows up as an out-of-date prerequisite.
#
#   config.stamp  the flags the generated headers depend on.  Without this a
#                 stale find_points.h would be compiled against a different
#                 RBA_PACK.  (The headers also carry a compile-time check of
#                 their own, for the case that they are used outside this
#                 Makefile.)
#   build.stamp   everything that reaches the compiler, including the tuning
#                 flags.  Every rule that runs ${CC} depends on it: without it
#                 a width change would rebuild the generated headers and the
#                 two objects that include them while leaving the others at
#                 the old register width -- which links and then crashes --
#                 and "make tune" would write tuning.mk without anything
#                 being recompiled to use it.
#
# They are separate so that a retune, which changes only build.stamp, does not
# regenerate the headers: their contents do not depend on the tuning.
.PHONY: FORCE
config.stamp: FORCE
	@echo '${CCFLAGS_H}' > $@.tmp
	@cmp -s $@.tmp $@ || \
	  { mv $@.tmp $@; echo "configuration changed; regenerating headers"; }
	@rm -f $@.tmp

build.stamp: FORCE
	@echo '${CCFLAGS_0} ${CCFLAGS}' > $@.tmp
	@cmp -s $@.tmp $@ || \
	  { mv $@.tmp $@; echo "compilation flags changed; rebuilding"; }
	@rm -f $@.tmp

gen_init_sieve_h: gen_init_sieve_h.c ratpoints.h rp-private.h primes.h config.stamp
	${CC} gen_init_sieve_h.c -o gen_init_sieve_h  ${CCFLAGS_0} ${CCFLAGS2} ${CCFLAGS}

gen_find_points_h: gen_find_points_h.c ratpoints.h rp-private.h primes.h config.stamp
	${CC} gen_find_points_h.c -o gen_find_points_h  ${CCFLAGS_0} ${CCFLAGS2} ${CCFLAGS}

init_sieve.h: gen_init_sieve_h
	./gen_init_sieve_h > init_sieve.h

find_points.h: gen_find_points_h
	./gen_find_points_h > find_points.h
