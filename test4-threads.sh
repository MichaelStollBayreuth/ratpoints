#!/bin/sh
# test4 and test1 with the sieve on several threads and the blocks of the
# loop over the denominators ending elsewhere: the library is built with
# -DRP_BLOCK_LEN=n for a few n, each in a build directory of its own
# (build-threads-n/, symbolic links to the sources, as test4-configs.sh
# does), and test4 and rptest are run on it with two thread counts, against
# the one reference each.  What the program prints must depend neither on
# the threads nor on where the blocks end.  "make testthreads" runs this
# after the suites on the tree's own build.  Exits with 1 when an output
# differs, with 2 when a build failed.  A couple of minutes.
#
# LENGTHS and COUNTS in the environment choose the block lengths and the
# thread counts; KEEP=1 keeps the build directories and the outputs.
LENGTHS=${LENGTHS:-'1 3 17'}
COUNTS=${COUNTS:-'3 6'}
status=0
for len in $LENGTHS; do
  dir=build-threads-$len
  rm -rf "$dir"; mkdir "$dir"
  for f in Makefile *.c *.h; do ln -s "../$f" "$dir/$f"; done
  # the generated headers must not be the links: remove them, so that the
  # build in the directory makes its own
  rm -f "$dir/find_points.h" "$dir/init_sieve.h"
  echo "== block length $len"
  if ! (cd "$dir" && make -s ratpoints rptest CCFLAGS="-DRP_BLOCK_LEN=$len" > make.log 2>&1)
  then echo "build failed, see $dir/make.log"; status=2; continue; fi
  for t in $COUNTS; do
    RP="./$dir/ratpoints" RPOPTS="-t $t" ./test4.sh > "test4-threads-$len-$t.out" 2>&1
    if cmp -s testbase4 "test4-threads-$len-$t.out"
    then echo "test4 -t $t ok"; [ -n "$KEEP" ] || rm -f "test4-threads-$len-$t.out"
    else echo "Test failed! (diff testbase4 test4-threads-$len-$t.out)"; [ $status = 2 ] || status=1
    fi
    "./$dir/rptest" -t $t > "rptest-threads-$len-$t.out"
    if cmp -s testbase "rptest-threads-$len-$t.out"
    then echo "rptest -t $t ok"; [ -n "$KEEP" ] || rm -f "rptest-threads-$len-$t.out"
    else echo "Test failed! (diff testbase rptest-threads-$len-$t.out)"; [ $status = 2 ] || status=1
    fi
  done
  [ -n "$KEEP" ] || rm -rf "$dir"
done
exit $status
