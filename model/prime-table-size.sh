#!/bin/bash
# Point-richest curves: core cycles at each (PRIME_SIZE, num_primes), each run
# paired with PRIME_SIZE=7 -p 30 (what the program does today).
set -u
S=/tmp/claude-1000/-home-mstoll-software-git-ratpoints/83b9b35a-2986-4d7f-a258-6c32978b4583/scratchpad
H=${H:-200000}; PAIRS=${PAIRS:-3}
CFG="8:30 8:35 8:40 8:45 8:53 9:70 9:96"

cyc() { perf stat -x, -e cpu_core/cycles/u -- taskset -c 0 "$@" 2>&1 >/dev/null \
        | grep 'cpu_core/cycles/u' | cut -d, -f1; }

printf "%-7s" "curve"; for c in $CFG; do printf "%10s" "$c"; done; echo
i=1
while read -r crv; do
  printf "#%-6s" $i
  for cfg in $CFG; do
    ps=${cfg%%:*}; np=${cfg##*:}
    rs=""; k=1
    while [ $k -le $PAIRS ]; do
      a=`cyc $S/bin4/rp-ps7 "$crv" $H -p 30`
      b=`cyc $S/bin4/rp-ps$ps "$crv" $H -p $np`
      rs="$rs `awk -v a=$a -v b=$b 'BEGIN{printf "%.4f", b/a}'`"
      k=`expr $k + 1`
    done
    printf "%10s" "`echo $rs | tr ' ' '\n' | sort -n | sed -n "$(( (PAIRS+1)/2 ))p"`"
  done
  echo
  i=`expr $i + 1`
done < $S/rich.txt
