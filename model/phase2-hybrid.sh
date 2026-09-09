#!/bin/bash
# Hybrid phase 2 (wide phase 1 and scan, 64-bit AND and extraction) against the
# all-wide baseline of the same register width, over a range of sp2 - sp1.
set -u
S=/tmp/claude-1000/-home-mstoll-software-git-ratpoints/83b9b35a-2986-4d7f-a258-6c32978b4583/scratchpad
PAIRS=${PAIRS:-3}
C1="1 0 126 0 441"
C2="247747600 -985905640 567207969 2396040466 52485681 -470135160 82342800"
cyc() { perf stat -x, -e cpu_core/cycles/u -- taskset -c 0 "$@" 2>&1 >/dev/null \
        | grep 'cpu_core/cycles/u' | cut -d, -f1; }

for spec in "sparse:$C1:150000:12" "point-rich:$C2:200000:22"; do
  n=${spec%%:*}; r=${spec#*:}; c=${r%%:*}; r=${r#*:}; h=${r%%:*}; n1=${r##*:}
  echo "=== $n, h=$h, sp1=$n1 -- hybrid / baseline, median of $PAIRS pairs ==="
  printf "  %-6s"; for d in 0 1 2 3 5 8; do printf " %9s" "d=$d"; done; echo
  for w in 128 256 512e; do
    printf "  %-6s" $w
    for d in 0 1 2 3 5 8; do
      rs=""; k=1
      while [ $k -le $PAIRS ]; do
        a=`cyc $S/hyb/rp-$w-base "$c" $h -n $n1 -N $((n1+d))`
        b=`cyc $S/hyb/rp-$w-hyb  "$c" $h -n $n1 -N $((n1+d))`
        rs="$rs `awk -v a=$a -v b=$b 'BEGIN{printf "%.4f", b/a}'`"
        k=`expr $k + 1`
      done
      printf " %9s" `echo $rs | tr ' ' '\n' | sort -n | sed -n "$(( (PAIRS+1)/2 ))p"`
    done
    echo
  done
done
