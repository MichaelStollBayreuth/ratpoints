#!/bin/sh
# The threaded runs under ThreadSanitizer.  The library and the drivers are
# built with -fsanitize=thread in build-tsan/ (symbolic links to the sources,
# as test4-configs.sh does; -O1 and -g, so that a report names the lines);
# then rptest, rptest-many and rptest-degrees run with threads, rpapi runs
# its own threaded cases, and test4 runs with -t 3.  Every output must
# match its reference, and the sanitizer must report nothing: its reports
# go to the .log file next to each output.  "make tsan" runs this.  Exits
# with 1 on a difference or a report, with 2 when the build failed.  A few
# minutes.  KEEP=1 keeps the build directory with the outputs and logs.
dir=build-tsan
rm -rf "$dir"; mkdir "$dir"
for f in Makefile *.c *.h; do ln -s "../$f" "$dir/$f"; done
rm -f "$dir/find_points.h" "$dir/init_sieve.h"
if ! (cd "$dir" && make -s rptest rptest-many rptest-degrees ratpoints rpapi \
                        CCFLAGS='-fsanitize=thread -g -O1' > make.log 2>&1)
then echo "build failed, see $dir/make.log"; exit 2; fi
# every report, not just the first
TSAN_OPTIONS="halt_on_error=0${TSAN_OPTIONS:+ $TSAN_OPTIONS}"
export TSAN_OPTIONS
status=0
# name, reference: the output and the log are $dir/name.out and .log
check() {
  if cmp -s "$2" "$dir/$1.out" && ! grep -q 'WARNING: ThreadSanitizer' "$dir/$1.log"
  then echo "$1 ok"
  else echo "Test failed! ($1: diff $2 $dir/$1.out; reports in $dir/$1.log)"; status=1
  fi
}
"$dir/rptest" -t 2 > "$dir/rptest-t2.out" 2> "$dir/rptest-t2.log"
check rptest-t2 testbase
"$dir/rptest" -t 4 > "$dir/rptest-t4.out" 2> "$dir/rptest-t4.log"
check rptest-t4 testbase
"$dir/rptest-many" -t 4 > "$dir/rptest-many-t4.out" 2> "$dir/rptest-many-t4.log"
check rptest-many-t4 testbase-many
"$dir/rptest-degrees" -t 3 > "$dir/rptest-degrees-t3.out" 2> "$dir/rptest-degrees-t3.log"
check rptest-degrees-t3 testbase-degrees
"$dir/rpapi" > "$dir/rpapi.out" 2> "$dir/rpapi.log"
check rpapi testbase-api
RP="./$dir/ratpoints" RPOPTS='-t 3' ./test4.sh > "$dir/test4-t3.out" 2> "$dir/test4-t3.log"
check test4-t3 testbase4
[ $status = 0 ] && [ -z "$KEEP" ] && rm -rf "$dir"
exit $status
