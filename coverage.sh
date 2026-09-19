#!/bin/sh
# Measure how much of the library and the program the tests exercise, with
# gcov: builds the sources instrumented and without optimisation in
# build-coverage/ (from symbolic links, so that the build in the working
# directory is left alone), runs the tests through that build -- test1,
# test1once, test3, test4 and testapi; test2 and timing add nothing but
# minutes at -O0 -- and prints the summary gcov gives for each source
# file: the lines executed and the branches taken at least once.  The
# annotated sources, <file>.c.gcov, are left in build-coverage/ for a look
# at what was missed: a line marked ##### was never executed, and "branch
# N taken 0%" after a line says one of its outcomes never happened.  "make
# coverage" runs this; a minute or so.  The suite was written for 2.3 and
# adapted; the notes of the 2.3 sources say what is left uncovered there.
# The exit status is 1 when the output of a test differs from its
# reference and 2 when the instrumented build failed.

dir=build-coverage
rm -rf "$dir"; mkdir "$dir"
for f in Makefile *.c *.h *.sh testbase* ; do ln -s "../$f" "$dir/$f"; done
rm -f "$dir/find_points.h" "$dir/init_sieve.h"
cd "$dir" || exit 2
# --coverage is -fprofile-arcs -ftest-coverage at compile time and the
# profiling library at link time; -O0 last overrides the Makefile's -O2/-O3,
# so that every branch of the source is a branch of the code
if ! make -s ratpoints rptest rpapi CCFLAGS='--coverage -O0' > make.log 2>&1
then echo "build failed, see $dir/make.log"; exit 2; fi
fail=0
run() { # name, command, reference
  eval "$2" > "$1.out" 2>&1
  if [ -n "$3" ] && ! cmp -s "$3" "$1.out"; then echo "$1: differs from $3"; fail=1; fi
}
run test1 './rptest' testbase
run test1once './rptest -O' testbase
run test1once2 './rptest -O -dl 2 -z' ''
grep -q changed test1once2.out && { echo "test1once2: a field changed"; fail=1; }
run test3 './test3.sh' testbase3
run test4 './test4.sh' testbase4
run testapi './rpapi' testbase-api
[ $fail = 0 ] && echo "all tests agree with their references"
echo
# gcov names the data after the object: main.c was compiled straight into
# the executable, so its notes are ratpoints-main.gcno
summary() { # the four lines gcov prints for the file named, less one
  awk -v f="File '$1'" '$0 == f { p = 4 } p > 0 { print; p-- }' | grep -v 'Branches executed'
}
for f in find_points sift sturm; do gcov -b "$f.c" 2>&1 | summary "$f.c"; done
# the code of init.c is in the generated header, which is what gcov reports
gcov -b init.c 2>&1 | summary init_sieve.h
gcov -b ratpoints-main.gcno 2>&1 | summary main.c
echo "(annotated sources in $dir/*.c.gcov)"
exit $fail
