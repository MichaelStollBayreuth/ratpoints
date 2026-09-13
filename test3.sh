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
# 13-14 no denominator class admits a numerator mod 16: num_bits[] used to be
#       read uninitialised, which with -j -F 0 printed points outside the
#       height bound after a 5000-fold blow-up
# 15    no odd denominator admits a numerator, but the even ones do: U was put
#       on its floor, which switched the second and third stages off (the
#       points were right; the run was slow)
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
$RP '2 0 3' 300000 -q -j -F 0
$RP '10 10 5 -7 0 3 -2' 16383 -q
