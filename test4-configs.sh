#!/bin/sh
# Run test4.sh against the library built with the other compile-time
# switches: the register widths (64, 128 by the SSE intrinsics and by
# SSE2, 512, 256 being what the tree is built with), one register in the
# first phase (RATPOINTS_CHUNK=1), the second phase one word at a time
# (USE_LONG_IN_PHASE_2), and the larger prime table (PRIME_SIZE 8, the
# default of 2.3).  The points found must not depend on any of these, and
# test4.sh prints nothing that does, so every build is compared with the
# same testbase4.  "make test4configs" runs this.
#
# Each configuration is built in a directory of its own, build-test4-<name>,
# from symbolic links to the sources, so that the build in the working
# directory is left alone; the directories are removed afterwards unless
# KEEP is set.  About a quarter of a minute for all of them.  The exit
# status is 1 when the output of a build differs from the reference and 2
# when a build failed (whatever the other builds did).
#
# The configurations, as name:flags, the flags separated by commas; they
# replace CCFLAGS1 of the Makefile (and PRIME_SIZE, for the one that
# changes it).  A third field "nopart2" compares the output without the
# second part of test4.sh, the reports of -v: two of those look at the
# primes beyond 127 when the table has them.
CONFIGS=${CONFIGS:-'64: 128:-DUSE_SSE 128a:-DUSE_AVX128 512:-DUSE_AVX512
  chunk1:-DUSE_AVX,-mavx2,-DRATPOINTS_CHUNK=1
  long2:-DUSE_AVX,-mavx2,-DUSE_LONG_IN_PHASE_2
  prime8:-DUSE_AVX,-mavx2,PRIME_SIZE=8:nopart2'}

status=0
for cfg in $CONFIGS; do
  name=${cfg%%:*}
  rest=${cfg#*:}
  part2=yes
  case "$rest" in *:nopart2) part2=no; rest=${rest%:nopart2} ;; esac
  flags=$(echo "$rest" | tr ',' ' ')
  psize=7
  case "$flags" in
    *PRIME_SIZE=*) psize=$(echo "$flags" | sed 's/.*PRIME_SIZE=\([0-9]*\).*/\1/')
                   flags=$(echo "$flags" | sed 's/PRIME_SIZE=[0-9]*//') ;;
  esac
  dir=build-test4-$name
  rm -rf "$dir"; mkdir "$dir"
  for f in Makefile *.c *.h; do ln -s "../$f" "$dir/$f"; done
  # the generated headers must not be the links: remove them, so that the
  # build in the directory makes its own
  rm -f "$dir/find_points.h" "$dir/init_sieve.h"
  # tuning.mk, if there is one, was measured for the tree's configuration
  # and would be ignored with a warning; leave it out
  echo "== $name: CCFLAGS1='$flags' PRIME_SIZE=$psize"
  if ! (cd "$dir" && make -s ratpoints CCFLAGS1="$flags" PRIME_SIZE=$psize > make.log 2>&1)
  then echo "build failed, see $dir/make.log"; status=2; continue; fi
  RP="./$dir/ratpoints" ./test4.sh > "test4-$name.out" 2>&1
  ref=testbase4
  if [ $part2 = no ]
  then sed '/^==== part 2/,/^==== part 3/{/^==== part 3/!d;}' testbase4 > "$dir/testbase4-nopart2"
       sed '/^==== part 2/,/^==== part 3/{/^==== part 3/!d;}' "test4-$name.out" > "$dir/out-nopart2"
       mv "$dir/out-nopart2" "test4-$name.out"
       ref="$dir/testbase4-nopart2"
  fi
  if cmp -s "$ref" "test4-$name.out"
  then echo "ok"; [ -n "$KEEP" ] || rm -rf "$dir" "test4-$name.out"
  else echo "Test failed! (diff $ref test4-$name.out)"; [ $status = 2 ] || status=1
  fi
done
exit $status
