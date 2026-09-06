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

PRIME_SIZE = 7
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
# Add "-DUSE_LONG_IN_PHASE_2" to work with unsigned long's instead of bit-arrays
#  in phase 2 of the sieving. This is usually slower.
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
# The following uses 128-bit SSE-registers.
CCFLAGS128 = -DUSE_SSE
# A variant of the above.
CCFLAGS128a = -DUSE_AVX128
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

CCFLAGS_0 = ${CCFLAGS0} ${CCFLAGS1}

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
            gpl-2.0.txt testbase2 testdata-many.h testbase-many \
            bench_init.c

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
	${RM} ${TARGETFILES}

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

# Correctness check and benchmark for the sieve table set-up (see bench_init.c).
bench_init: libratpoints.a bench_init.c ratpoints.h rp-private.h primes.h
	${CC} bench_init.c -o bench_init ${CCFLAGS_0} -O3 -funroll-loops \
              ${CCFLAGS3} ${CCFLAGS2} ${CCFLAGS}

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
