#!/bin/sh
# What the library and the program do when memory runs out (see rpalloc.c).
# rpalloc and ratpoints-alloc are built in build-alloc/ (symbolic links to
# the sources, as test4-configs.sh does), with the calls of malloc, calloc
# and realloc that the library and main.c make routed through rpalloc.c by
# GNU ld's option --wrap, so that the n-th of them can be made to fail, or
# every one from the n-th on.  rpalloc runs its searches on one thread and
# on three, for n = 1, 2, 3, ... until nothing fails any more, and judges
# each against the same search without failures: it must deliver the same
# points and return the same value, or return RATPOINTS_NO_MEMORY having
# delivered the first so many of those points.  Then the program runs a few
# searches the same way, each n in a run of its own: it must print what it
# prints without failures, or the first lines of that and then "Not enough
# memory." with exit status 10, followed by the hint that -s saves memory
# when it was the library that ran out of memory and the Sturm step was on.
# A crash, or a run that has not ended after a minute, fails the test as
# well.  The memory gmp takes for its integers is not affected (gmp ends
# the program when that runs out).
#
# "make testalloc" runs this.  It needs a linker with --wrap: GNU ld, gold
# and lld have it, the linker of macOS has not; without it, the script says
# so and exits with 2, as it does when the build failed; with 1 when a
# search went wrong.  Half a minute or so.  CC names the compiler (gcc by
# default; make passes its own); KEEP=1 keeps the build directory with the
# outputs.
dir=build-alloc
CC=${CC:-gcc}
rm -rf "$dir"; mkdir "$dir" "$dir/ld"
# whether the linker sends a call of malloc to __wrap_malloc
cat > "$dir/ld/wrap.c" << 'EOF'
#include <stdlib.h>
void *__real_malloc(size_t n);
void *__wrap_malloc(size_t n);
void *__wrap_malloc(size_t n) { return(n == 12345 ? NULL : __real_malloc(n)); }
int main(void) { void *p = malloc(12345); int r = (p != NULL); free(p); return(r); }
EOF
if ! (cd "$dir/ld" && $CC -O0 wrap.c -o wrap -Wl,--wrap=malloc > wrap.log 2>&1 && ./wrap)
then echo "testalloc.sh needs a linker with GNU ld's option --wrap, and that of $CC has not"
     rm -rf "$dir"; exit 2
fi
for f in Makefile *.c *.h; do ln -s "../$f" "$dir/$f"; done
rm -f "$dir/find_points.h" "$dir/init_sieve.h"
if ! (cd "$dir" && make -s rpalloc ratpoints-alloc CC="$CC" > make.log 2>&1)
then echo "build failed, see $dir/make.log"; exit 2; fi
status=0

# the library: rpalloc says what it tried and whether every search was right
for t in 1 3; do
  for mode in one from; do
    out="$dir/rpalloc-$t-$mode.out"
    if "./$dir/rpalloc" $t $mode > "$out" 2>&1
    then tail -n 1 "$out"
    else cat "$out"; echo "Test failed! (rpalloc $t $mode: see $out)"; status=1
    fi
  done
done

# the program: name, then the arguments of ratpoints.  The run without
# failures is the reference; then n = 1, 2, ... for each of the two ways
# of failing, until n has passed the last allocation, which the line
# rpalloc.c writes to standard error when the program exits says.  The
# program makes its own allocations before it calls the library, so that
# the message with the hint on -s (see nomem) comes for the n beyond them
# and the one without it for the n up to them: a run with the hint must not
# come before one without it, and some run must have it when the Sturm
# step is on (no -s among the arguments), none when it is off.
prog() {
  name=$1; shift
  out="$dir/$name"
  sturm=1; for a; do [ "$a" = -s ] && sturm=0; done
  if ! "./$dir/ratpoints-alloc" "$@" > "$out.ref" 2> "$out.log"
  then echo "Test failed! ($name: no success without failures)"; status=1; return; fi
  for var in RPALLOC_FAIL RPALLOC_FROM; do
    n=1; bad=0; hinted=0
    while :; do
      env "$var=$n" "./$dir/ratpoints-alloc" "$@" > "$out.out" 2> "$out.log"; st=$?
      failed=$(sed -n 's/^rpalloc: [0-9]* allocations, \([0-9]*\) failed$/\1/p' "$out.log")
      if [ $st = 0 ] && cmp -s "$out.ref" "$out.out"; then :
      elif [ $st = 10 ] && nomem "$out.out" "$out.ref" && form; then :
      else echo "  $name, $var=$n: exit status $st, output in $out.n$n"
           cp "$out.out" "$out.n$n"; bad=1
      fi
      [ "$failed" = 0 ] && break
      n=$((n + 1))
      if [ $n -gt 10000 ]; then echo "  $name, $var: the allocations do not end"; bad=1; break; fi
    done
    if [ $sturm = 1 ] && [ $hinted = 0 ]
    then echo "  $name, $var: no run had the hint on -s"; bad=1
    fi
    if [ $bad = 0 ]
    then echo "$name, $var=n for n = 1 to $((n - 1)): ok"
    else echo "Test failed! ($name, $var: compare with $out.ref)"; status=1
    fi
  done
}
# whether the form of the message that nomem found is one the run may
# print at this point (see prog); hinted records that a run had the hint
form() {
  if [ $hint = 1 ]; then [ $sturm = 1 ] && hinted=1
  else [ $hinted = 0 ]
  fi
}
# whether the output $1 is the first lines of $2 and then the message of
# the lack of memory, "\nNot enough memory.\n\n" (every line the program
# prints before it ends with a newline), alone (hint=0) or followed by the
# hint "At large degree, the option -s saves memory (no Sturm sequence).\n\n"
# (hint=1), which the program adds when the library ran out of memory with
# the Sturm step on.  When -1 has stopped the search (stopped=1), the point
# it stopped at must not be among those lines: a point that could not be
# kept comes after it, and the search must not fail then.
nomem() {
  l=$(($(wc -l < "$1")))
  hint_text='At large degree, the option -s saves memory (no Sturm sequence).'
  if [ $l -ge 5 ] && [ "$(tail -n 5 "$1" | tr '\n' '|')" = "|Not enough memory.||$hint_text||" ]
  then hint=1; k=5
  elif [ $l -ge 3 ] && [ "$(tail -n 3 "$1" | tr '\n' '|')" = "|Not enough memory.||" ]
  then hint=0; k=3
  else return 1
  fi
  [ $stopped = 0 ] || [ $((l - k)) -lt $(($(wc -l < "$2"))) ] || return 1
  head -n $((l - k)) "$1" > "$1.head"; head -n $((l - k)) "$2" > "$1.refhead"
  cmp -s "$1.head" "$1.refhead"
}
# the program's own allocations (the arrays of the coefficients and of the
# search intervals, the strings of -f, -fs, -fm and -fe, the digits of an
# interval end it reads exactly, and what it hands to the callback), and
# then the library's, on one thread, also without the Sturm step (-s); the
# banner and the report of a run without -q; on three threads, also with a
# stop at the first point
stopped=0
prog strings '1 0 1' 300 -q -f '%x %y %z\n' -fs 'start\n' -fm '' -fe 'end\n'
# the list of the coefficients read from a file (@path) of 10000 bytes,
# for which the buffer is taken and then doubled twice
printf '1 0%9995s1\n' '' > "$dir/coeffs"
prog file "@$dir/coeffs" 300 -q
prog nosturm '1 0 1' 300 -q -s
prog banner '5 3 7' 70000 -dl 66000 -du 66100 -l 0 -u 0.01
prog threads '1 0 1' 1000 -q -t 3
stopped=1
prog stop '1 0 1' 1000 -q -t 3 -1 -i
[ $status = 0 ] && [ -z "$KEEP" ] && rm -rf "$dir"
exit $status
