#   ratpoints-2.2.3
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
#   Michael Stoll, September 21, 2009; January 7, 2022; September 6, 2026
#   with changes by Bill Allombert, December 29, 2021

# The main targets are
#   all         the library, the program and this documentation
#   test        test1, test1many, test2 and timing (see below)
#   test1       1000 random genus 2 curves, checked against testbase
#   test1many   curves with many rational points, against testbase-many;
#               the two cover the two regimes that behave differently, and
#               both should be used when judging a change to the sieve
#   tune        measure the two machine-dependent constants that decide how
#               many primes each sieving stage uses, and write them to
#               tuning.mk (see tune.sh); takes several minutes, wants an idle
#               machine, and must not be run under -j
#   clean       remove the intermediate files; distclean also the executables
#               and any tuning.mk

# The program sieves with the odd primes below 2^PRIME_SIZE.  The default is 8
# rather than 7 not because 53 primes are normally used -- the choice below
# still starts from RATPOINTS_DEFAULT_NUM_PRIMES of them -- but so that the few
# primes past that are there for the curves that run out, which are the ones
# with very many rational points; on those it is worth up to a factor of two.
# Having the larger table costs nothing when it is not used.
PRIME_SIZE = 8
VERSION = 2.2.3

# The test targets time their runs with the shell's "time".  That is a shell
# built-in, and dash -- which is /bin/sh on Debian and Ubuntu -- does not have
# it, so name a shell that does rather than relying on whatever /bin/sh is.
# Override on the command line if bash lives elsewhere on your system.
SHELL = /bin/bash

CC = gcc
RM = rm -f
INSTALL = cp

INSTALL_DIR = /usr/local

CCFLAGS0 = -Wall -O2 -fomit-frame-pointer -DRATPOINTS_MAX_BITS_IN_PRIME=${PRIME_SIZE}
# For gcc on Apple, may have to add '-fnested-functions' to CCFLAGS0.
# Add "-DUSE_LONG_IN_PHASE_2" to sieve the survivors of the first phase one
#  64-bit word at a time instead of a whole bit-array at a time. The first
#  phase and the scan for survivors stay at the full register width either
#  way. It is a wash to 4% slower at 128 and 256 bits, and the reason is worth
#  knowing: the number of AND steps is the same either way, because the other
#  words of a surviving bit-array were already zero, and a narrow and a wide
#  read of the same table come from one cache line. At 64 bits the two are the
#  same code. See PHASE-NOTES.md on the phases-by-register-width branch.
# Add "-DRP_PHASE_TIMING" to have sift.c time the two phases of the sieve
#  separately and write a report to stderr when the program exits; add
#  "-DRP_PHASE_COUNTS" as well to count the bit-arrays surviving phase 1.
#  These are development aids; see the comment at the top of sift.c.
# Add "-DRATPOINTS_CHUNK=<n>" to force the use of 2 <= n <= 16 registers
#  in phase 1 of sieving. For n=1, this reverts to the code used previously.
#  If SSE/AVX registers are used and this is not set, 16 registers will be used.
#  In some cases, using n=8 may be faster
#  (e.g., Intel(R) Xeon(R) CPU E3-1220 V2 with -DUSE_AVX -mavx).

# The following uses word-length registers.
# This should work on essentially every machine.
CCFLAGS64 =
# The following uses 128-bit registers. In spite of its name, USE_AVX128 needs
# only SSE2, which every x86-64 machine has, so this is as portable as the
# USE_SSE variant below and never slower: its test for an empty register
# compares the whole register against zero instead of extracting both halves
# into general registers. That is worth 30% of the second sieving phase, hence
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
# This has not been run on an AVX512F machine yet; see the comment at the
# top of the USE_AVX512 branch in rp-private.h .  Leaving out "-mavx512f"
# makes gcc emulate the 64-byte vectors with narrower ones; that is slower,
# but it runs anywhere and exercises the same code path.  (gcc then warns
# "AVX512F vector argument without AVX512F enabled changes the ABI"; this is
# harmless here, since the whole library is built with the same flags.)
CCFLAGS512 = -DUSE_AVX512 -mavx512f

# This will be the default. Change as appropriate.
# CCFLAGS1 = ${CCFLAGS128}
CCFLAGS1 = ${CCFLAGS256}

# Machine-dependent tuning of the two constants that decide how many primes
# each sieving stage uses.  "make tune" writes tuning.mk; without it the
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
            main.c rptest.c testdata.h testbase ratpoints-doc-2.2.tex \
            gpl-2.0.txt testbase2 testdata-many.h testbase-many tune.sh

# Temporary files that are generated during build and test
# and can be removed afterwards
TEMPFILES = sift.o init.o sturm.o find_points.o \
            sift.s sift.i init.s find_points.h init_sieve.h \
            gen_find_points_h gen_init_sieve_h \
            rptest.out rptest-many.out \
            sift-debug.o find_points-debug.o main.o test2.out

# Executables and library produced when building
TARGETFILES = ratpoints libratpoints.a rptest rptest-many ratpoints-debug \
              ratpoints-doc-2.2.pdf

FAILED = "Test failed!"

all: ratpoints libratpoints.a doc

doc: ratpoints-doc-2.2.pdf

test: test1 test1many test2 timing

# Measure good values for the two machine-dependent constants and write them
# to tuning.mk; see tune.sh.  Deliberately not part of "make all": it takes a
# couple of minutes, wants an otherwise idle machine, and must not run under
# "make -j".  Run "make all" afterwards to rebuild with the result.
.PHONY: tune
tune: rptest rptest-many
	@TUNE_CONFIG='${TUNE_CONFIG}' ./tune.sh

# Run ratpoints on a set of 1000 test cases
# and check the output
test1: rptest testbase
	time ./rptest > rptest.out
	cmp -s testbase rptest.out || echo ${FAILED}

# The same, but for curves with many rational points, which sieve very
# differently: about one numerator in 10^4 survives the first phase on a
# random curve, up to a hundred times more on these.  The two tests take
# about the same time and should be used together whenever the constants
# that choose sp1 and sp2 (see ratpoints.h) are retuned; a change that
# helps one regime can easily hurt the other.
test1many: rptest-many testbase-many
	time ./rptest-many > rptest-many.out
	cmp -s testbase-many rptest-many.out || echo ${FAILED}

# Run ratpoints on the curve with the record number of known
# rational points, with a fairly large height bound,
# time it and compare with the expected output
test2: ratpoints testbase2
	time ./ratpoints '247747600 -985905640 567207969 2396040466 52485681 -470135160 82342800' 1000000 -n 30 -N 30 -p 30 -q > test2.out
	cmp -s testbase2 test2.out || echo ${FAILED}

# Time a call to ratpoints with a largish height parameter.
# This can be helpful to assess modifications to the sieving process.
timing: ratpoints
	time ./ratpoints '1 0 126 0 441' 400000 -q > /dev/null

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
ratpoints-doc-2.2.pdf: ratpoints-doc-2.2.tex
	pdflatex ratpoints-doc-2.2.tex
	pdflatex ratpoints-doc-2.2.tex

dist: ${DISTFILES}
	mkdir -p ratpoints-${VERSION}
	cp ${DISTFILES} ratpoints-${VERSION}/
	tar --create --file=ratpoints-${VERSION}-`date --rfc-3339=date`.tar.gz \
	    --gzip --dereference ratpoints-${VERSION}
	rm -r ratpoints-${VERSION}

clean:
	${RM} ${TEMPFILES}

distclean: clean
	${RM} ${TARGETFILES} tuning.mk

debug: ratpoints-debug

libratpoints.a: sift.o init.o sturm.o find_points.o
	ar rs libratpoints.a sift.o init.o sturm.o find_points.o

ratpoints: libratpoints.a main.c ratpoints.h
	${CC} main.c -o ratpoints ${CCFLAGS_0} ${CCFLAGS2} ${CCFLAGS3} ${CCFLAGS}

main.o: main.c ratpoints.h
	${CC} main.c -c -o main.o ${CCFLAGS_0} -O3 ${CCFLAGS}

ratpoints-debug: sift-debug.o init.o sturm.o find_points-debug.o main.o
	${CC} sift-debug.o init.o sturm.o find_points-debug.o main.o \
              -o ratpoints-debug ${CCFLAGS_0} ${CCFLAGS2} ${CCFLAGS}

sift.o: sift.c ratpoints.h rp-private.h
	${CC} sift.c -c -o sift.o ${CCFLAGS_0} -funroll-loops ${CCFLAGS}

sift-debug.o: sift.c ratpoints.h rp-private.h
	${CC} sift.c -c -o sift-debug.o ${CCFLAGS_0} -funroll-loops -DDEBUG ${CCFLAGS}

sift.s: sift.c ratpoints.h rp-private.h
	${CC} sift.c -S -o sift.s ${CCFLAGS_0} -funroll-loops ${CCFLAGS}

sift.i: sift.c ratpoints.h rp-private.h
	${CC} sift.c -E -o sift.i ${CCFLAGS_0} -funroll-loops ${CCFLAGS}

init.o: init.c ratpoints.h rp-private.h init_sieve.h
	${CC} init.c -c -o init.o ${CCFLAGS_0} -funroll-loops -O3 ${CCFLAGS}

init.s: init.c ratpoints.h rp-private.h init_sieve.h
	${CC} init.c -S -o init.s ${CCFLAGS_0} -funroll-loops -O3 ${CCFLAGS}

sturm.o: sturm.c ratpoints.h rp-private.h
	${CC} sturm.c -c -o sturm.o ${CCFLAGS_0} ${CCFLAGS}

find_points.o: find_points.c ratpoints.h rp-private.h primes.h find_points.h
	${CC} find_points.c -c -o find_points.o ${CCFLAGS_0} ${CCFLAGS}

find_points-debug.o: find_points.c ratpoints.h rp-private.h primes.h find_points.h
	${CC} find_points.c -c -o find_points-debug.o ${CCFLAGS_0} -DDEBUG ${CCFLAGS}

rptest: libratpoints.a rptest.c ratpoints.h testdata.h
	${CC} rptest.c -o rptest ${CCFLAGS_0} ${CCFLAGS2} ${CCFLAGS3} ${CCFLAGS}

rptest-many: libratpoints.a rptest.c ratpoints.h testdata-many.h
	${CC} rptest.c -o rptest-many -DRATPOINTS_TESTDATA='"testdata-many.h"' \
	      ${CCFLAGS_0} ${CCFLAGS2} ${CCFLAGS3} ${CCFLAGS}

gen_init_sieve_h: gen_init_sieve_h.c ratpoints.h rp-private.h primes.h
	${CC} gen_init_sieve_h.c -o gen_init_sieve_h  ${CCFLAGS_0} ${CCFLAGS2} ${CCFLAGS}

gen_find_points_h: gen_find_points_h.c ratpoints.h rp-private.h primes.h
	${CC} gen_find_points_h.c -o gen_find_points_h  ${CCFLAGS_0} ${CCFLAGS2} ${CCFLAGS}

init_sieve.h: gen_init_sieve_h
	./gen_init_sieve_h > init_sieve.h

find_points.h: gen_find_points_h
	./gen_find_points_h > find_points.h
