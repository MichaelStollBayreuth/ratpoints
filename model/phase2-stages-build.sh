#!/bin/bash
# The RP_STOP_AFTER ladder at each register width.
set -u
cd /home/mstoll/software/git/ratpoints
S=/tmp/claude-1000/-home-mstoll-software-git-ratpoints/83b9b35a-2986-4d7f-a258-6c32978b4583/scratchpad
declare -A W=( [64]="" [64s]="" [128]="-DUSE_AVX128" [256]="-DUSE_AVX -mavx2" [512e]="-DUSE_AVX512" )
declare -A X=( [64]="-DRATPOINTS_CHUNK=16" [64s]="-DRATPOINTS_CHUNK=16 -DRP_LONG_SKIP" \
               [128]="" [256]="" [512e]="" )
for w in 64 64s 128 256 512e; do
  for lvl in 1 2 3 4 0; do
    make clean >/dev/null 2>&1
    if make ratpoints CCFLAGS1="${W[$w]}" CCFLAGS="${X[$w]} -DRP_STOP_AFTER=$lvl" \
            >$S/lad/build.log 2>&1
    then cp ratpoints $S/lad/rp-$w-L$lvl
    else echo "FAILED $w L$lvl"; tail -4 $S/lad/build.log; fi
  done
  echo "built $w"
done
make clean >/dev/null 2>&1
