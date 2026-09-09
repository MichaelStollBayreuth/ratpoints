#!/bin/bash
# Is any -falign setting systematically better?  Each candidate is compared
# with the default build of the same width on the same curve; a setting worth
# adopting has to win across widths and curves, not just somewhere.
set -u
cd /home/mstoll/software/git/ratpoints
S=/tmp/claude-1000/-home-mstoll-software-git-ratpoints/83b9b35a-2986-4d7f-a258-6c32978b4583/scratchpad
C1="1 0 126 0 441"
C2="247747600 -985905640 567207969 2396040466 52485681 -470135160 82342800"
declare -A W=( [64]="" [128]="-DUSE_AVX128" [256]="-DUSE_AVX -mavx2" )
declare -A A=( [default]="" [loops32]="-falign-loops=32" [loops64]="-falign-loops=64" \
               [both32]="-falign-functions=32 -falign-loops=32" )

for w in 64 128 256; do
  for a in default loops32 loops64 both32; do
    make clean >/dev/null 2>&1
    make ratpoints CCFLAGS1="${W[$w]}" CCFLAGS="${A[$a]}" >/dev/null 2>&1 \
      && cp ratpoints $S/al/rp-$w-$a || echo "FAILED $w $a"
  done
done
make clean >/dev/null 2>&1

cyc() { perf stat -x, -e cpu_core/cycles/u -- taskset -c 0 "$@" 2>&1 >/dev/null \
        | grep 'cpu_core/cycles/u' | cut -d, -f1; }

printf "%-8s %-11s %10s %10s %10s\n" width curve loops32 loops64 both32
for w in 64 128 256; do
  for spec in "sparse:$C1:200000" "point-rich:$C2:300000"; do
    n=${spec%%:*}; c=${spec#*:}; c=${c%:*}; h=${spec##*:}
    printf "%-8s %-11s" $w $n
    for a in loops32 loops64 both32; do
      rs=""; k=1
      while [ $k -le 3 ]; do
        x=`cyc $S/al/rp-$w-default "$c" $h`
        y=`cyc $S/al/rp-$w-$a "$c" $h`
        rs="$rs `awk -v a=$x -v b=$y 'BEGIN{printf "%.4f", b/a}'`"
        k=`expr $k + 1`
      done
      printf " %10s" `echo $rs | tr ' ' '\n' | sort -n | sed -n 2p`
    done
    echo
  done
done
