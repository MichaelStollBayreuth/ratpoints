#!/bin/bash
# Rebuild every register-width variant from the current HEAD, which has the
# automatic sp1/sp2 choice and the phase-2 skip loop.  The generated headers
# depend on CCFLAGS1, and this branch has the config.stamp guard, so a plain
# "make ratpoints" is enough -- but "make clean" is cheap and removes all doubt.
set -u
cd /home/mstoll/software/git/ratpoints
S=/tmp/claude-1000/-home-mstoll-software-git-ratpoints/83b9b35a-2986-4d7f-a258-6c32978b4583/scratchpad
mkdir -p $S/bin2

# name -> CCFLAGS1 ; the 64-bit build needs chunking to be competitive at all
declare -A W=( [64c16]="" [128]="-DUSE_SSE" [128a]="-DUSE_AVX128" \
               [256]="-DUSE_AVX -mavx2" [512e]="-DUSE_AVX512" )
declare -A X=( [64c16]="-DRATPOINTS_CHUNK=16" [128]="" [128a]="" [256]="" [512e]="" )

for w in 64c16 128 128a 256 512e; do
  for p2 in ba lo; do
    # at 64 bits rp-private.h forces USE_LONG_IN_PHASE_2; skip the duplicate
    [ "$w" = 64c16 ] && [ "$p2" = ba ] && continue
    p2flag=""; [ "$p2" = lo ] && p2flag="-DUSE_LONG_IN_PHASE_2"
    for inst in plain time count; do
      iflag=""
      [ "$inst" = time  ] && iflag="-DRP_PHASE_TIMING"
      [ "$inst" = count ] && iflag="-DRP_PHASE_TIMING -DRP_PHASE_COUNTS"
      make clean >/dev/null 2>&1
      if make ratpoints CCFLAGS1="${W[$w]}" CCFLAGS="${X[$w]} $p2flag $iflag" \
              >$S/build2.log 2>&1
      then cp ratpoints $S/bin2/rp-$w-$p2-$inst; echo "built $w-$p2-$inst"
      else echo "BUILD FAILED $w-$p2-$inst"; tail -5 $S/build2.log; fi
    done
  done
done
make clean >/dev/null 2>&1
