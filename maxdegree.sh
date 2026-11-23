#!/bin/sh
# The refusals at the bound on the degree, which the default build cannot
# reach: its RATPOINTS_MAX_DEGREE is far beyond any list of coefficients a
# test can give.  The program and rpapi are built with RATPOINTS_MAX_DEGREE
# 100 in build-maxdegree/ (symbolic links to the sources, as
# test4-configs.sh does).  There a list of 101 coefficients is searched, as
# the program of the default build searches it; one of 102 is refused with
# error 3, "Too many coefficients."; and 102 integers followed by a word
# that is not one with error 11, since such a word is refused wherever it
# stands.  rpapi must print its reference, which shows that the degree
# RATPOINTS_MAX_DEGREE + 1, LONG_MAX and the negative ones are refused by
# find_points, find_points_work and find_points_init, and that a search at
# degree 100 runs.  "make testmaxdegree", which "make test" runs, runs
# this.  Exits with 1 on a difference, with 2 when the build failed; a few
# seconds besides the build.  KEEP=1 keeps the build directory with the
# outputs.
dir=build-maxdegree
rm -rf "$dir"; mkdir "$dir"
for f in Makefile *.c *.h; do ln -s "../$f" "$dir/$f"; done
rm -f "$dir/find_points.h" "$dir/init_sieve.h"
if ! (cd "$dir" && make -s ratpoints rpapi CCFLAGS=-DRATPOINTS_MAX_DEGREE=100 \
        > make.log 2>&1)
then echo "build failed, see $dir/make.log"; exit 2; fi
status=0
# name, reference: the output is $dir/name.out
check() {
  if cmp -s "$2" "$dir/$1.out"
  then echo "$1 ok"
  else echo "Test failed! ($1: diff $2 $dir/$1.out)"; status=1
  fi
}
# a run and its exit status
run() { "$@"; echo "exit $?"; }
# y^2 = x^100 + 1, 101 coefficients; with one more, y^2 = x^101 + x^100 + 1
s=1; i=1; while [ $i -lt 100 ]; do s="$s 0"; i=$((i+1)); done; s="$s 1"
run ./ratpoints "$s" 10 -q > "$dir/coeffs101.ref"
run "./$dir/ratpoints" "$s" 10 -q > "$dir/coeffs101.out"
check coeffs101 "$dir/coeffs101.ref"
printf '\nToo many coefficients.\n\nexit 3\n' > "$dir/coeffs102.ref"
run "./$dir/ratpoints" "$s 1" 10 -q > "$dir/coeffs102.out"
check coeffs102 "$dir/coeffs102.ref"
printf '\nIncorrect coefficients (not integers separated by blanks).\n\nexit 11\n' \
  > "$dir/coeffs102x.ref"
run "./$dir/ratpoints" "$s 1 x" 10 -q > "$dir/coeffs102x.out"
check coeffs102x "$dir/coeffs102x.ref"
"./$dir/rpapi" > "$dir/rpapi.out" 2>&1
check rpapi testbase-api
[ $status = 0 ] && [ -z "$KEEP" ] && rm -rf "$dir"
exit $status
