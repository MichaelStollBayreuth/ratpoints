#!/bin/sh
# Regression tests for the bugs found in the review of September 2026: each
# line is one invocation of ./ratpoints whose output is compared with
# testbase3 by "make test3".  All point lists were checked against a
# brute-force search over all coprime (a, b) with |a|, b <= H.
#
#  1-7  degree 1, given as such or reached by stripping a zero leading
#       coefficient or by reversal of y^2 = c2 x^2 + c1 x: used to read past
#       the end of the Sturm arrays and crash
#  8-12 the positivity region misses the search domain: the points at
#       infinity used to be dropped (8-10 in general, 11-12 for even degree
#       with a square leading coefficient when reversal is suppressed)
# 13-15 no denominator class admits a numerator mod 16: num_bits[] used to be
#       read uninitialised, which with -j -F 0 printed points outside the
#       height bound after a 5000-fold blow-up.  13 and 15 print nothing
#       either way; 14 shows the message of the early return, which the
#       unfixed program does not have
# 16-17 no odd denominator admits a numerator, but the even ones do.  On the
#       2.3 line this used to put the run-length estimate on its floor and
#       switch the second and third sieving stages off (the points were
#       right, the run slow): 16 pins the point list, 17 counts the -v lines
#       saying that a stage uses no primes, which must be 0 (2.2.4 has no
#       such estimate and no line 17)
# 18    not squarefree: the Sturm chain reaches a zero remainder, and the
#       loop that found it used to index one below the array
# 19-20 the report the program prints after a run (the primes used, the
#       reversal, the search intervals), which no other test looks at; the
#       numbers of primes are pinned so that retuning does not change it
RP=./ratpoints
$RP '1 2' 20 -q
$RP '1 2 0' 50 -q
$RP '0 19 1' 45 -q
$RP '0 1 -1' 45 -q
$RP '0 3 3' 45 -q
$RP '0 -9' 45 -q
$RP '-4 1' 45 -q
$RP '-1000000000 0 0 1' 45 -q
$RP '0 1 0 0 -1000000000' 45 -q
$RP '-10000000 0 0 1' 100 -q
$RP '-1000000 0 0 0 1' 10 -q -k
$RP '-1000000 0 0 0 1' 10 -q -l -5 -u 5
$RP '2 0 3' 200 -q -j -F 0 -x
$RP '2 0 3' 200 -v -j -F 0 -x | grep 'mod 16'
$RP '2 0 3' 300000 -q -j -F 0
$RP '10 10 5 -7 0 3 -2' 16383 -q
$RP '10 10 5 -7 0 3 -2' 16383 -v 2>&1 | grep -c 'use 0 primes for second stage\|use 0 primes for third stage'
$RP '1 2 1' 20 -q
$RP '1 0 126 0 441' 100 -n 5 -N 8 -P 2 -z | sed -n '/primes used/,$p'
$RP '10 10 5 -7 0 3 -2' 1000 -n 5 -N 8 -P 2 -z | sed -n '/primes used/,$p'
