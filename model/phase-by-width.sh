#!/bin/bash
# Phase 1 / phase 2 / gmp-check split by register width, with the parameters
# chosen automatically (i.e. at the survivor rate the tuning aims at).
# Variants are rotated inside each round and the minimum over rounds is kept,
# because this machine slows by up to a quarter as it warms up.
set -u
S=/tmp/claude-1000/-home-mstoll-software-git-ratpoints/83b9b35a-2986-4d7f-a258-6c32978b4583/scratchpad
cd $S
ROUNDS=${ROUNDS:-5}

V=(64c16-lo 128-ba 128-lo 128a-ba 128a-lo 256-ba 256-lo 512e-ba 512e-lo)
C1="1 0 126 0 441"
C2="247747600 -985905640 567207969 2396040466 52485681 -470135160 82342800"

: > $S/phase2.csv
echo "curve,variant,round,perf,cyc1,cyc2,cyc3,cyctot,sp1,sp2" >> $S/phase2.csv

for r in $(seq 1 $ROUNDS); do
  n=${#V[@]}
  for k in $(seq 0 $((n-1))); do
    i=$(( (k + r) % n ))            # rotate the order every round
    v=${V[$i]}
    for cur in sparse dense; do
      if [ $cur = sparse ]; then c="$C1"; h=400000; else c="$C2"; h=200000; fi
      out=$(perf stat -x, -e cpu_core/cycles/u -- \
            taskset -c 0 ./bin2/rp-$v-time "$c" $h 2>&1 >/dev/null)
      pc=$(echo "$out" | grep 'cpu_core/cycles/u' | cut -d, -f1)
      pd=$(echo "$out" | grep phasedata)
      get() { echo "$pd" | tr ' ' '\n' | grep "^$1=" | cut -d= -f2; }
      echo "$cur,$v,$r,$pc,$(get cyc1),$(get cyc2),$(get cyc3),$(get cyctot),$(get sp1),$(get sp2)" \
        >> $S/phase2.csv
    done
  done
  echo "round $r done  ($(date +%H:%M:%S))"
done
