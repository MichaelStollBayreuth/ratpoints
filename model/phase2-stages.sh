#!/bin/bash
# The five levels are run back to back inside each round, so that the
# differences between them are taken on a machine in the same thermal state;
# the median of the per-round differences is reported.  Levels below 0 share
# the set-up and the whole first phase, so both cancel exactly.
set -u
S=/tmp/claude-1000/-home-mstoll-software-git-ratpoints/83b9b35a-2986-4d7f-a258-6c32978b4583/scratchpad
ROUNDS=${ROUNDS:-5}
C1="1 0 126 0 441"
C2="247747600 -985905640 567207969 2396040466 52485681 -470135160 82342800"

cyc() { perf stat -x, -e cpu_core/cycles/u -- taskset -c 0 "$@" 2>&1 >/dev/null \
        | grep 'cpu_core/cycles/u' | cut -d, -f1; }

for spec in "sparse:$C1:200000:12:17" "point-rich:$C2:300000:22:27"; do
  n=${spec%%:*}; r=${spec#*:}; c=${r%%:*}; r=${r#*:}; h=${r%%:*}; r=${r#*:}
  n1=${r%%:*}; n2=${r##*:}
  echo "=== $n, h=$h, sp1=$n1 sp2=$n2 ==="
  printf "  %-6s %9s %9s %9s %9s %9s\n" build "phase 1" scan AND extract check
  for w in 64 64s 128 256 512e; do
    : > $S/lad/d.txt
    k=1
    while [ $k -le $ROUNDS ]; do
      set --
      for lvl in 1 2 3 4 0; do set -- "$@" `cyc $S/lad/rp-$w-L$lvl "$c" $h -n $n1 -N $n2`; done
      awk -v a=$1 -v b=$2 -v d=$3 -v e=$4 -v f=$5 \
        'BEGIN{printf "%d %d %d %d %d\n", a, b-a, d-b, e-d, f-e}' >> $S/lad/d.txt
      k=`expr $k + 1`
    done
    printf "  %-6s" $w
    for col in 1 2 3 4 5; do
      v=`awk -v k=$col '{print $k}' $S/lad/d.txt | sort -n | sed -n "$(( (ROUNDS+1)/2 ))p"`
      sp=`awk -v k=$col '{print $k}' $S/lad/d.txt | sort -n | awk 'NR==1{a=$1} {b=$1} END{printf "%.2f", (b-a)/1e9}'`
      printf " %7s+-%-4s" `awk -v v=$v 'BEGIN{printf "%.3f", v/1e9}'` $sp
    done
    echo
  done
done
