#!/bin/sh
# Measure how much of the library and the program the tests exercise, with
# gcov: builds the sources instrumented and without optimisation in
# build-coverage/ (from symbolic links, so that the build in the working
# directory is left alone), runs the tests through that build -- test1,
# test1many, testdegrees, test1once, test3, test4 and testapi; test2 and
# timing add nothing but minutes at -O0 -- and prints the summary gcov gives
# for each source file: the lines executed and the branches taken at least
# once.  The annotated sources, <file>.c.gcov, are left in build-coverage/
# for a look at what was missed: a line marked ##### was never executed,
# and "branch N taken 0%" after a line says one of its outcomes never
# happened.  "make coverage" runs this; a minute or so.
#
# What the suite of item 20 left uncovered on 2026-09-19 is listed in
# test4.sh's notes: guards against states the callers cannot produce,
# arms reached only with other compile-time settings (test4-configs.sh
# covers those builds, but this script measures one), and two conditions
# the arithmetic makes impossible.

dir=build-coverage
rm -rf "$dir"; mkdir "$dir"
for f in Makefile *.c *.h *.sh testbase* ; do ln -s "../$f" "$dir/$f"; done
rm -f "$dir/find_points.h" "$dir/init_sieve.h"
cd "$dir" || exit 1
# --coverage is -fprofile-arcs -ftest-coverage at compile time and the
# profiling library at link time; -O0 last overrides the Makefile's -O2/-O3,
# so that every branch of the source is a branch of the code
if ! make -s ratpoints rptest rptest-many rptest-degrees rpapi \
          CCFLAGS='--coverage -O0' > make.log 2>&1
then echo "build failed, see $dir/make.log"; exit 1; fi
fail=0
run() { # name, command, reference
  eval "$2" > "$1.out" 2>&1
  if [ -n "$3" ] && ! cmp -s "$3" "$1.out"; then echo "$1: differs from $3"; fail=1; fi
}
run test1 './rptest' testbase
run test1many './rptest-many' testbase-many
run testdegrees './rptest-degrees' testbase-degrees
run test1once './rptest -O' testbase
run test1once2 './rptest -O -dl 2 -z' ''
grep -q changed test1once2.out && { echo "test1once2: a field changed"; fail=1; }
run test1once3 './rptest -O -dl 0 -du 1000000 -S 100' testbase
run test3 './test3.sh' testbase3
run test4 './test4.sh' testbase4
run testapi './rpapi' testbase-api
[ $fail = 0 ] && echo "all tests agree with their references"
echo
# gcov names the data after the object: main.c was compiled straight into
# the executable, so its notes are ratpoints-main.gcno
for f in find_points sift init sturm; do
  gcov -b "$f.c" 2>&1 | grep -A3 "^File '$f.c'" | grep -v 'Branches executed'
done
gcov -b ratpoints-main.gcno 2>&1 | grep -A3 "^File 'main.c'" | grep -v 'Branches executed'
echo "(annotated sources in $dir/*.c.gcov)"
exit $fail
