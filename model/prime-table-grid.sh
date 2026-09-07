#!/bin/bash
# Core cycles for (PRIME_SIZE, num_primes) on both populations, each candidate
# paired with PRIME_SIZE=7 -p 30 and only the ratio kept.
set -u
S=/tmp/claude-1000/-home-mstoll-software-git-ratpoints/83b9b35a-2986-4d7f-a258-6c32978b4583/scratchpad
PAIRS=${PAIRS:-3}
CFG=${CFG:-"8:30 8:40 8:53 9:30 9:53 9:70 9:96 10:30 10:53 10:96 10:171"}

cyc() { perf stat -x, -e cpu_core/cycles/u -- taskset -c 0 "$@" 2>&1 >/dev/null \
        | grep 'cpu_core/cycles/u' | cut -d, -f1; }

for suite in rptest rptestmany; do
  echo "=== $suite, cycles relative to PRIME_SIZE=7 -p 30 ==="
  for cfg in $CFG; do
    ps=${cfg%%:*}; np=${cfg##*:}
    rs=""
    i=1
    while [ $i -le $PAIRS ]; do
      a=`cyc $S/bin4/$suite-ps7 -p 30`
      b=`cyc $S/bin4/$suite-ps$ps -p $np`
      rs="$rs `awk -v a=$a -v b=$b 'BEGIN{printf "%.4f", b/a}'`"
      i=`expr $i + 1`
    done
    med=`echo $rs | tr ' ' '\n' | sort -n | sed -n "$(( (PAIRS+1)/2 ))p"`
    printf "  PRIME_SIZE=%-3s -p %-4s %s   (%s)\n" $ps $np "$med" "$rs"
  done
done
