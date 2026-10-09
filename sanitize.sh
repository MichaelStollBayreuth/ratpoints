#!/bin/sh
# test4 and rpapi under AddressSanitizer and UndefinedBehaviorSanitizer, for
# what the other tests cannot see: a read or write beyond an allocation,
# memory that is never freed (AddressSanitizer looks for leaks when a run
# ends), a signed overflow or a conversion out of range, which can all
# leave every output as it should be.  The program and rpapi are built in
# build-sanitize/ (symbolic links to the sources, as test4-configs.sh does;
# -O1 and -g, so that a report names the lines) with both sanitizers and
# float-cast-overflow, which -fsanitize=undefined leaves out, and with every
# report fatal; then test4 runs on one thread and with -t 3, and rpapi
# runs.  Every output must match its reference, and the sanitizers must
# report nothing: their reports go to the .log file next to each output.
# AddressSanitizer fills every allocation up to 1 MB with non-zero bytes
# rather than the first 4096 bytes only, so that a scan for a terminating
# zero that was never written runs off its array every time.  "make
# sanitize" runs this.  Exits with 1 on a difference or a report, with 2
# when the build failed (a compiler without the sanitizers, say).  A minute
# or so.  KEEP=1 keeps the build directory with the outputs and logs.
dir=build-sanitize
rm -rf "$dir"; mkdir "$dir"
for f in Makefile *.c *.h; do ln -s "../$f" "$dir/$f"; done
rm -f "$dir/find_points.h" "$dir/init_sieve.h"
if ! (cd "$dir" && make -s ratpoints rpapi \
        CCFLAGS='-O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined,float-cast-overflow -fno-sanitize-recover=all' \
        > make.log 2>&1)
then echo "build failed, see $dir/make.log"; exit 2; fi
ASAN_OPTIONS="max_malloc_fill_size=1048576${ASAN_OPTIONS:+:$ASAN_OPTIONS}"
UBSAN_OPTIONS="print_stacktrace=1${UBSAN_OPTIONS:+:$UBSAN_OPTIONS}"
export ASAN_OPTIONS UBSAN_OPTIONS
status=0
# name, reference: the output and the log are $dir/name.out and .log
check() {
  if cmp -s "$2" "$dir/$1.out" && ! grep -q 'Sanitizer\|runtime error' "$dir/$1.log"
  then echo "$1 ok"
  else echo "Test failed! ($1: diff $2 $dir/$1.out; reports in $dir/$1.log)"; status=1
  fi
}
RP="./$dir/ratpoints" ./test4.sh > "$dir/test4.out" 2> "$dir/test4.log"
check test4 testbase4
RP="./$dir/ratpoints" RPOPTS='-t 3' ./test4.sh > "$dir/test4-t3.out" 2> "$dir/test4-t3.log"
check test4-t3 testbase4
"$dir/rpapi" > "$dir/rpapi.out" 2> "$dir/rpapi.log"
check rpapi testbase-api
[ $status = 0 ] && [ -z "$KEEP" ] && rm -rf "$dir"
exit $status
