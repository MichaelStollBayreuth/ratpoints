#!/bin/bash
# Per-phase cost of a larger PRIME_SIZE at a fixed -p 30, so that exactly the
# same primes are used and only the memory layout differs.  Each candidate is
# run immediately after a PRIME_SIZE=7 run and only the ratios are kept.
set -u
S=/tmp/claude-1000/-home-mstoll-software-git-ratpoints/83b9b35a-2986-4d7f-a258-6c32978b4583/scratchpad
C1="1 0 126 0 441"; C2="247747600 -985905640 567207969 2396040466 52485681 -470135160 82342800"
PAIRS=${PAIRS:-5}

row() { taskset -c 0 "$1" "$2" "$3" -p 30 2>&1 >/dev/null \
        | tr ' ' '\n' | grep -E "^(cyc1|cyc2|cyc3|cyctot)=" | cut -d= -f2 | tr '\n' ' '; }

for spec in "sparse:$C1:400000" "dense:$C2:200000"; do
  n=${spec%%:*}; rest=${spec#*:}; c=${rest%:*}; h=${rest##*:}
  echo "=== $n (h=$h), ratio to PRIME_SIZE=7, median of $PAIRS pairs ==="
  printf "  %-6s %8s %8s %8s %8s %8s\n" build phase1 phase2 gmp setup total
  for ps in 8 9 10; do
    : > /tmp/claude-1000/pr.txt
    i=1
    while [ $i -le $PAIRS ]; do
      set -- $(row $S/bin3/rp-ps7-time "$c" $h);    a1=$1 a2=$2 a3=$3 at=$4
      set -- $(row $S/bin3/rp-ps$ps-time "$c" $h);  b1=$1 b2=$2 b3=$3 bt=$4
      awk -v a1=$a1 -v a2=$a2 -v a3=$a3 -v at=$at -v b1=$b1 -v b2=$b2 -v b3=$b3 -v bt=$bt \
        'BEGIN{printf "%.4f %.4f %.4f %.4f %.4f\n", b1/a1, b2/a2, b3/a3,
                      (bt-b1-b2-b3)/(at-a1-a2-a3), bt/at}' >> /tmp/claude-1000/pr.txt
      i=`expr $i + 1`
    done
    printf "  ps%-4s" $ps
    for col in 1 2 3 4 5; do
      printf " %8s" "$(awk -v k=$col '{print $k}' /tmp/claude-1000/pr.txt | sort -n | sed -n "$(( (PAIRS+1)/2 ))p")"
    done
    echo
  done
done
