#!/bin/sh
# The suite that exercises every branch of the code (TODO item 20).  Each
# line runs ./ratpoints once; "make test4" compares the output of the whole
# script with testbase4.  Unlike test1 and its relatives, which run a
# thousand curves through one setting, this one runs a few hundred settings
# through a few curves: every degree from 1 up, odd and even, square and
# non-square leading and constant coefficients, every reason to reverse the
# polynomial, curves with no points modulo some prime, curves with no
# admissible numerator modulo 64, curves with and without forbidden
# divisors of the denominator, every command-line option, restricted
# denominator ranges and search intervals, height bounds from 1 to 2^63 - 1,
# and every error the program reports.  It runs in a few seconds; the
# heights are small, and the large ones come with a narrow search interval
# and a denominator range of one.
#
# The reference was checked independently: verify-test4.py reads testbase4,
# and for every invocation whose output is a list of points it searches all
# coprime pairs (a, b) in the range by brute force -- arbitrary precision,
# no sieve -- and compares; for the other invocations it checks what can be
# checked (that -1 printed a point of the curve, that the count -z reports
# is the number of points, ...).  A second, independent check is that
# ratpoints 2.2.4 prints the same points for the whole first part.
#
# Before each invocation the script prints its arguments, one per <...>, so
# that a failing test can be found from the diff of test4.out against
# testbase4, and so that verify-test4.py knows what it is checking.
#
# The script takes the program from $RP, ./ratpoints by default, so that it
# can be run against a build with other flags: test4-configs.sh does that
# for the register widths and the other compile-time switches.
#
# This is the copy for 2.2.4.  The first part prints points and messages
# that are the same in 2.2.4 and 2.3 (and the reference agrees with 2.3's
# there, line for line); the second part looks at the reports of -v, whose
# texts differ between the versions, so the reference is 2.2.4's own.  The
# script in the 2.3 sources has a third part for the options 2.3 added
# (-r -R -U -C -A -P -Q -W), and its second part pins the third sieving
# stage with -P, which 2.2.4 does not have; those options are left out
# here.

RP=${RP:-./ratpoints}
# no program to run: exit 2 (the comparison with the reference is make's,
# whose target fails with 1 when they differ)
[ -x "$RP" ] || { echo "$RP: not an executable" >&2; exit 2; }

# run ratpoints, the arguments announced first
t() { printf '#'; for a; do printf ' <%s>' "$a"; done; echo; "$RP" "$@"; }
# the same, with the output run through a filter (the first argument, a
# shell command reading its standard input), which is announced with the
# arguments: for the runs that look at one line of a report, or count them
f() { filt=$1; shift
      printf '#'; for a; do printf ' <%s>' "$a"; done; echo " | $filt"
      "$RP" "$@" | eval "$filt"
}
# the same, and report the exit status: for the tests of the error
# messages, with the usage text taken out (it names the version and lists
# the options, which differ between 2.2.4 and 2.3)
e() { printf '#'; for a; do printf ' <%s>' "$a"; done; echo
      out=$("$RP" "$@"); st=$?
      printf '%s\n' "$out" | grep -v '^This is ratpoints-\|^Usage: \|^ *\['
      echo "exit $st"
}
# the same as t, for runs without -q: the first line of the banner names
# the version, and the report on the moduli used depends on the width and
# the tuning, as below
m() { t "$@" | grep -v '^This is ratpoints-\|used for\|further primes\|primes used'; }
# the same as t, for the tests of the report -v prints, with the lines that
# depend on the register width, the prime table and the machine-dependent
# tuning taken out: how many bits are set per word and how many moduli
# each stage was given (the mean over the classes visited depends on
# nothing else, but the numbers chosen do, and the message continues on a
# second line), the lists of the moduli each stage uses (the ranking
# depends on the width of a bit array), and what the exact check is put at
# (the compiled-in constants).  Two of the reports look at primes beyond
# 127 when the table has them, so they differ in a build with PRIME_SIZE 8
# (the table ends at 127 with the 7 of 2.2.4); test4-configs.sh compares
# that build without this part
v() { t "$@" | awk '
  /bits set per word/ { skip = 1; next }
  /^  use [0-9]+ (moduli|primes) for (first|second|third) stage:/ { skip = 1; next }
  /one exact check is put at/ { skip = 1; next }
  skip { skip = 0; next }
  { print }'
}

echo '==== part 1: points and messages ===='

echo '---- degree 1: y^2 = c1 x + c0 ----'
# monic (the denominators are squares); a rational root x = -1/2 gives a
# point with y = 0, printed once
t '1 2' 20 -q
t '1 2' 20 -q -k
# lcf not +-1: squares times divisors of the leading coefficient
t '3 5' 100 -q
t '-3 7' 60 -q
t '2 -4' 60 -q
t '1 12' 50 -q
# zero constant term: y^2 = x has the point (0 : 0 : 1); y^2 = 6x
t '0 1' 30 -q
t '0 6' 30 -q
t '0 -1' 30 -q
# a leading zero coefficient is dropped; from an odd degree that lowers the
# genus (the form of even degree gets a double root at infinity), which is
# refused
t '1 2 0' 50 -q
t '1 2 3 4 0' 50 -q
e '0 19 1 0' 45 -q

echo '---- degree 2 ----'
# leading coefficient a square: points at infinity, no Jacobi test, no
# forbidden divisors, so every denominator is tried
t '1 0 1' 50 -q
t '-1 0 1' 50 -q
t '1 0 -1' 50 -q
# negative everywhere: no real points
t '-1 0 -1' 50 -q
t '-1 0 -1' 50 -q -s
# constant term zero: reversed to degree 1
t '0 1 2' 40 -q
t '0 1 2' 40 -q -k
# lcf a square, constant term not: reversed
t '2 3 1' 50 -q
t '2 3 1' 50 -q -k
t '1 3 2' 50 -q
# a non-square leading coefficient: the Jacobi symbol test and the primes
# with (lcf/p) = -1 as forbidden divisors
t '5 3 7' 100 -q
t '5 3 7' 100 -q -j
t '5 3 7' 100 -q -F 0
t '5 3 7' 100 -q -F 0 -j
t '5 3 7' 100 -q -F 1
t '4 0 3' 50 -q
t '4 0 3' 50 -q -F 0
t '9 0 -2' 50 -q
# the prime 3 divides the leading coefficient: which valuations of the
# denominator at 3 are excluded depends on the Newton polygon (see
# forbidden_valuations): v_3(b) = 1 only; v_3(b) = 1 and 3^4 | b; 3^3 | b;
# v_3(b) = 1 and 3^3 | b; 3^2 | b (for 3x^2 + 2x + 1 the two top terms tie
# at v_3(b) = 1).  When every valuation is excluded, as for 3x^2 + 4 above,
# 3 | b is tested by a bit array instead
t '1 1 9' 100 -q
t '1 1 27' 100 -q
t '1 3 27' 100 -q
t '1 1 18' 100 -q
t '1 2 3' 100 -q
t '1 2 3' 100 -q -F 0
# no denominator admits a numerator mod 64, and no point at infinity
t '2 0 3' 200 -q
# not squarefree
e '0 0 1' 20 -q
e '1 2 1' 20 -q

echo '---- degree 3 ----'
t '1 0 0 1' 100 -q
t '-1 0 0 1' 100 -q
t '1 0 0 -1' 100 -q
t '1 0 0 -1' 100 -q -k
# a rational root
t '1 1 1 1' 100 -q
# squares times divisors of the leading coefficient
t '3 2 -1 12' 100 -q
t '3 2 -1 12' 100 -q -k
t '2 0 0 5' 100 -q
# constant term zero, coefficient of x +-1, lcf not +-1: reversed
t '0 1 0 2' 50 -q
t '0 -1 0 2' 50 -q
t '0 1 0 2' 50 -q -k
# the leading coefficient 1031*1033 has no prime factor in the table of
# odd primes below 1024, so it is not factored and the denominators are
# not restricted
t '5 0 0 1065023' 100 -q
# x^3 is not squarefree
e '0 0 0 1' 20 -q

echo '---- degree 4 ----'
t '1 0 126 0 441' 200 -q
t '1 0 126 0 441' 200 -q -k
t '1 0 126 0 441' 200 -q -j
t '1 0 126 0 441' 200 -q -j -J
t '1 0 126 0 441' 200 -q -k -K
t '1 0 126 0 441' 200 -q -x -X -y -Y -i -I -z -Z -k -K -j -J -s -S
t '1 0 0 0 -1' 50 -q
t '-1 0 0 0 1' 100 -q
t '1 1 1 1 1' 100 -q
# constant term zero: reversed to degree 3, which is monic, so the
# reversed curve restricts its denominators to squares
t '0 1 0 0 2' 50 -q
t '0 1 0 0 2' 50 -q -k
# lcf a square and constant term not: reversed
t '2 0 0 0 1' 50 -q
t '-4 0 0 0 9' 50 -q
t '-4 0 0 0 9' 50 -q -k
# no points modulo 3 (and none modulo 64 either); no points modulo 3 but
# numerators modulo 64, so the primes are looked at
t '2 0 1 0 2' 100 -q
t '2 63 1 0 2' 100 -q
# no points modulo 131, which is beyond the 30 primes looked at: the run
# sieves and finds none (part 2 makes the third stage look that far)
t '2 131 4 0 2' 100 -q
# a leading coefficient beyond 2^63: 10^20 + 1 = 73 * 137 * 1676321 *
# 5964848081, so the Jacobi symbol is computed on gmp integers
t '1 0 0 0 100000000000000000001' 100 -q
t '1 0 0 0 100000000000000000001' 100 -q -j
t '100000000000000000001 0 0 0 1' 100 -q
# divisible by x^2
e '0 0 1 0 2' 20 -q

echo '---- degrees 5 to 8 ----'
t '1 0 0 0 0 1' 60 -q
# -x^5 - 38x^4 - 33x^3 + 35x^2 - 38x + 22 admits no numerator for an odd
# denominator, so the squares of odd numbers are skipped
t '22 -38 35 -33 -38 -1' 100 -q
t '-3 5 -7 11 -13 17' 60 -q
t '1 1 1 1 1 1 1' 100 -q
t '-3 5 -7 11 -13 17 -19' 100 -q
t '1 0 0 0 0 0 0 4' 40 -q
t '1 0 0 0 0 0 0 0 1' 40 -q
t '2 0 0 0 0 0 0 0 3' 40 -q

echo '---- higher degrees ----'
t '1 0 0 0 0 0 0 0 0 1' 30 -q
t '4 0 0 0 0 0 0 0 0 0 0 0 1' 20 -q
t '1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1' 15 -q
# degree 100, the largest allowed
t '1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1' 10 -q
# degree 101
e '1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1' 10 -q

echo '---- curves with many rational points ----'
# these run out of informative primes, so the program looks at more of them
t '10 10 5 -7 0 3 -2' 1000 -q
t '0 -223813512 -29715903 47242068 5804812 -1242416 3136' 300 -q
t '0 -87865180 -234179851 125769766 107149081 -18406280 547600' 300 -q

echo '---- the numerator classes mod 64 ----'
# what the bit arrays hold, by class of the denominator mod 64 (the report
# of -v says it in part 2): odd denominators only, with 24 of 64 numerators;
# odd denominators one numerator in two, denominators 2 mod 4 one in four,
# 8 mod 16 every odd one; odd denominators one in four; one in eight; a
# cubic with odd denominators one in eight and six kinds of even ones; and
# no class at all (like 3x^2 + 2 above, without a point at infinity)
t '-36 -54 -5' 100 -q
t '44 -2 -12' 100 -q
t '4 28 -53' 100 -q
t '-52 -28 43' 100 -q
t '32 48 55 -44' 100 -q
t '-40 56 45' 100 -q

echo '---- the sieving primes ----'
# no first phase at all; one modulus; second phase empty; third stage
# empty; every prime in the first phase; more primes than the table has
t '1 0 126 0 441' 200 -q -n 0
t '1 0 126 0 441' 200 -q -n 0 -N 0
t '1 0 126 0 441' 200 -q -n 1
t '1 0 126 0 441' 200 -q -n 2 -N 2
t '1 0 126 0 441' 200 -q -n 3 -N 8
t '1 0 126 0 441' 200 -q -n 30 -N 30 -p 30
t '1 0 126 0 441' 200 -q -p 5
t '1 0 126 0 441' 200 -q -p 100
t '1 0 126 0 441' 200 -q -N 40 -p 10
t '1 0 126 0 441' 200 -q -n 8 -N 4
t '3 2 -1 12' 200 -q -n 0
t '3 2 -1 12' 200 -q -n 1 -N 1
t '5 3 7' 200 -q -n 2 -N 6
t '10 10 5 -7 0 3 -2' 300 -q -n 4 -N 6
t '10 10 5 -7 0 3 -2' 300 -q -p 8
# few primes to look at, with the choice free to look further
t '1 0 126 0 441' 200 -q -p 3
t '-3 5 -7 11 -13 17 -19' 200 -q -p 3
t '3 2 -1 12' 200 -q -p 2
# f = (x^3 + x + 1)^2 + 3*5*...*127 is a square modulo every one of the
# first thirty odd primes, so none says anything and the choice has to
# look beyond them; the constant term does not fit a long
t '2007238469666518094547220599513022568322942623866 2 1 2 2 0 1' 100 -q
t '2007238469666518094547220599513022568322942623866 2 1 2 2 0 1' 100 -q -p 40
# no exact check: the survivors of the sieve are printed as they are; with
# every prime of the table sieving, they are the points.  Every -x run
# pins its primes like this: what the sieve leaves depends on the moduli
# chosen, and so on the tuning, which nothing else in this file does
t '1 0 126 0 441' 200 -q -x -n 15 -N 30 -p 30
t '3 2 -1 12' 200 -q -x -n 15 -N 30 -p 30
t '1 0 1' 100 -q -x -n 15 -N 30 -p 30
t '1 0 1' 100 -q -x -X

echo '---- forbidden divisors ----'
t '1 1 9' 100 -q -F 0
t '1 1 9' 100 -q -F 1
t '1 1 27' 100 -q -F 1000
t '5 3 7' 100 -q -F 2
t '5 3 7' 100 -q -F 2 -j
# a denominator range beyond 251^2, where primes past the table of sieving
# primes become forbidden divisors and their patterns are built on the spot
t '5 3 7' 70000 -q -dl 66000 -du 66100 -l 0 -u 0.01
t '5 3 7' 70000 -q -dl 66000 -du 66100 -l 0 -u 0.01 -F 0

echo '---- the isolation of the positivity region ----'
# without it; with 0, 1, 2 and 100 iterations
t '-1 0 0 0 1' 100 -q -s
t '-1 0 0 0 1' 100 -q -S 0
t '-1 0 0 0 1' 100 -q -S 1
t '-1 0 0 0 1' 100 -q -S 2
t '-1 0 0 0 1' 100 -q -S 100
t '-1 0 0 0 1' 100 -q -S
# three intervals: f = -(x^2-1)(x^2-4)(x^2-9) is positive between the
# roots, alternately
t '36 0 -49 0 14 0 -1' 100 -q
t '36 0 -49 0 14 0 -1' 100 -q -s
t '36 0 -49 0 14 0 -1' 100 -q -S 3
# positive only beyond the height bound (|x| > 1000), so no affine point,
# but the points at infinity are there
t '-1000000000000 0 0 0 1' 100 -q
t '-1000000000000 0 0 0 1' 100 -q -k
t '-1000000000000 0 0 0 1' 100 -q -i
t '-1000000000000 0 0 0 1' 100 -q -l -5 -u 5
# ... and not there
t '-1000000 0 0 0 3' 100 -q
t '-1000000000 0 0 0 3' 100 -q
t '-1000000000 0 0 0 3' 100 -q -s
# f = x^4 - 55*3*5*...*127 is negative up to |x| = 10^12, so there are no
# affine points below the bound, but its constant term is 1 mod 64 and 0
# modulo every sieving prime: without the isolation the numerator 0
# survives every prime and the exact check sees a negative value
t '-110398115831658495200097132973216241257761844312575 0 0 0 1' 10 -q
t '-110398115831658495200097132973216241257761844312575 0 0 0 1' 10 -q -s
t '-110398115831658495200097132973216241257761844312575 0 0 0 1' 10 -q -s -n 0 -N 0

echo '---- the search intervals ----'
t '1 0 126 0 441' 200 -q -l 0 -u 1
t '1 0 126 0 441' 200 -q -l -1
t '1 0 126 0 441' 200 -q -u 0
t '1 0 126 0 441' 200 -q -l -2 -u -1 -l 0 -u 0.1
t '1 0 126 0 441' 200 -q -l 0.5 -u 0.5
t '1 0 126 0 441' 200 -q -l 0.05 -u 0.1
t '1 0 126 0 441' 200 -q -l -0.1 -u 0.1 -l 1 -u 2 -l 5
t '36 0 -49 0 14 0 -1' 100 -q -l 0 -u 100
t '36 0 -49 0 14 0 -1' 100 -q -l -2.5 -u 2.5
t '36 0 -49 0 14 0 -1' 100 -q -l -2.5 -u -1.5 -l 1.5 -u 2.5
t '36 0 -49 0 14 0 -1' 100 -q -l 0 -u 0.5 -l 2 -u 2.5
t '36 0 -49 0 14 0 -1' 100 -q -l 3 -u 50
t '3 2 -1 12' 100 -q -l 0 -u 1
t '3 2 -1 12' 100 -q -l -1 -u 0
t '3 2 -1 12' 100 -q -l -0.6 -u -0.4
# intervals outside the height bound
t '1 0 126 0 441' 200 -q -l 300 -u 400
t '1 0 126 0 441' 200 -q -l -400 -u -300
t '1 0 126 0 441' 200 -q -l -400 -u 400

echo '---- the denominators ----'
t '1 0 126 0 441' 200 -q -dl 2
t '1 0 126 0 441' 200 -q -du 10
t '1 0 126 0 441' 200 -q -dl 21 -du 21
t '1 0 126 0 441' 200 -q -dl 22 -du 22
t '1 0 126 0 441' 200 -q -dl 0 -du 0
t '1 0 126 0 441' 200 -q -dl 5 -du 3
t '1 0 126 0 441' 200 -q -dl 200 -du 200
t '1 0 126 0 441' 200 -q -dl 201 -du 300
t '1 0 126 0 441' 200 -q -dl -5 -du -5
# no square in the range; the range starts above the first squares
t '1 0 0 1' 100 -q -dl 5 -du 8
t '1 0 0 1' 100 -q -dl 5
# the loop over the squares starts at the first square in the range,
# ceil(sqrt(b_low))^2, and the loops over the squares times a divisor d of
# the leading coefficient at d*ceil(sqrt(b_low/d))^2; y^2 = x^3 + 17 has
# points with the denominators 1, 4, 9, 25, 64 and 81 below height 300,
# 12x^3 - x^2 + 2x + 3 with 1, 2, 3, 4, 6, 12, 16, 25, 27, 48, 49, 50, ...
t '17 0 0 1' 300 -q -dl 4 -du 9
t '17 0 0 1' 300 -q -dl 5 -du 9
t '17 0 0 1' 300 -q -dl 2 -du 4
t '17 0 0 1' 300 -q -dl 9 -du 9
t '17 0 0 1' 300 -q -dl 10 -du 80
t '17 0 0 1' 300 -q -dl 26 -du 63
t '17 0 0 1' 300 -q -dl 65
t '3 2 -1 12' 300 -q -dl 2 -du 2
t '3 2 -1 12' 300 -q -dl 3 -du 6
t '3 2 -1 12' 300 -q -dl 7 -du 12
t '3 2 -1 12' 300 -q -dl 13 -du 26
t '3 2 -1 12' 300 -q -dl 49 -du 49
# a divisor of the leading coefficient beyond the denominator bound
t '3 2 -1 12' 100 -q -du 5
# a denominator range whose squares times divisors are 50 = 2 * 5^2 and
# 54 = 6 * 3^2
t '3 2 -1 12' 100 -q -dl 50 -du 60
# the prime dividing the leading coefficient is beyond the denominator bound
t '1 2 3' 2 -q
# a denominator range beyond 1021^2: the search for forbidden divisors runs
# to the end of the table of odd primes below 1024
t '5 3 7' 1050000 -q -dl 1042441 -du 1042500 -l 0 -u 1e-5
t '5 3 7' 1050000 -q -dl 1042441 -du 1042500 -l 0 -u 1e-5 -F 1000
t '3 2 -1 12' 100 -q -dl 4 -du 50
t '3 2 -1 12' 100 -q -dl 5 -du 24
t '1 0 1' 100 -q -dl 64 -du 64
t '1 0 1' 100 -q -dl 63 -du 65
t '5 3 7' 100 -q -dl 64 -du 64
t '5 3 7' 100 -q -dl 63 -du 65
t '5 3 7' 100 -q -dl 1 -du 1
t '5 3 7' 100 -q -dl 100 -du 100
# a denominator range and a search interval together
t '5 3 7' 1000 -q -dl 500 -du 520 -l 0.5 -u 0.7

echo '---- small height bounds ----'
for h in 1 2 3 4 7 8 9 63 64 65 127 128 129 255 256 257 511 512 513; do
  t '1 0 126 0 441' $h -q
done
for h in 1 2 3 64 256; do
  t '3 2 -1 12' $h -q
  t '1 0 1' $h -q
  t '5 3 7' $h -q
done
t '1 2' 1 -q
t '0 1' 1 -q
t '2 0 3' 1 -q

echo '---- large height bounds, a denominator at a time ----'
# y^2 = x^2 + 1 has the point (a : y : b) for every primitive Pythagorean
# triple; one is planted inside each window of about 10^5 numerators
# (see verify-test4.py for the triples).  The denominators are 2^32 - 1
# and 2^32, where the reduction modulo the sieving primes changes its
# method, and 2^62 - 2^31, where every word number is far beyond 2^32
t '1 0 1' 4611686018427387904 -q -dl 4294967295 -du 4294967295 -l 9544371.76666664 -u 9544371.766689923
t '1 0 1' 4611686018427387904 -q -dl 4294967296 -du 4294967296 -l 1073741823.9999999 -u 1073741824.0000234
t '1 0 1' 4611686018427387904 -q -dl 4611686016279904256 -du 4611686016279904256 -l 0.7499999994179233 -u 0.749999999417945
t '1 0 1' 4611686018427387904 -q -dl 4611686016279904256 -du 4611686016279904256 -l 0.7499999994179233 -u 0.749999999417945 -n 3 -N 6
t '1 0 1' 4611686018427387904 -q -dl 4611686016279904256 -du 4611686016279904256 -l -0.749999999417945 -u -0.7499999994179233
# a non-square leading coefficient, so the Jacobi symbol test runs on
# denominators beyond 2^32 (by jacobi1, the table-driven form not being
# available there) and beyond 2^59, and once with the denominator dividing
# the leading coefficient
t '1 1 6' 4611686018427387904 -q -dl 8589934593 -du 8589934600 -l 0 -u 1e-8
t '1 1 17179869186' 4611686018427387904 -q -dl 8589934593 -du 8589934600 -l 0 -u 1e-8
t '1 1 2305843009213693954' 4611686018427387904 -q -dl 1152921504606846977 -du 1152921504606846980 -l 0 -u 1e-17
t '1 1 6' 4611686018427387904 -q -dl 1152921504606846977 -du 1152921504606846980 -l 0.5 -u 0.50000000000001
# and with a leading coefficient beyond 2^63
t '1 1 100000000000000000001' 4611686018427387904 -q -dl 8589934593 -du 8589934600 -l 0 -u 1e-8
# a denominator above 2^31 in the loop over the squares times the divisors
# of the leading coefficient, 3 * 59827^2, which is not divisible by 9 but
# whose low 32 bits are: the test on the valuation at 3 (an even positive
# v_3(b) is excluded here) must look at the whole denominator.  The point
# is x = 1/(3k^2) on y^2 = 3x^3 + 9k^6 + 2 again
t '412691823703593288194398668803 0 0 3' 1099511627776 -q -dl 10737809787 -du 10737809787 -l 0 -u 1.4e-10
# the leading coefficient 3*5*...*61 has 17 odd prime factors, too many
# for the table-driven Jacobi symbol test; 971*...*1021 has nine whose
# non-square tables would not fit
t '1 1 58644190679703485491635' 300 -q
t '1 1 58644190679703485491635' 300 -q -F 0
t '1 1 979798255664800613767220501' 300 -q
t '1 1 8209' 300 -q
# 8209 * 15: the table-driven test is not available, so jacobi1 computes
# the symbol for every denominator, including the ones sharing a factor
# with the leading coefficient (which -F 0 lets through), and with a
# negative leading coefficient
t '1 1 123135' 300 -q -F 0
t '1 1 -123135' 300 -q -F 0
t '1 1 123135' 300 -q
# the largest height bound a long holds, 2^63 - 1: a window at the very top
# of the range (the last bit array of the last word; a monic linear
# polynomial, so that the denominators are the squares up to 1), the
# denominator 2^63 - 1 itself in the loop over every denominator (which
# has to stop without incrementing past it), and the last eight
# denominators in the loop that tests them
t '-9223372036854769900 1' 9223372036854775807 -q -dl 1 -du 1 -l 9223372036854767616 -u 1e30
t '1 0 1' 9223372036854775807 -q -dl 9223372036854775807 -du 9223372036854775807 -l 0.5 -u 0.500000000000001
t '5 3 7' 9223372036854775807 -q -dl 9223372036854775800 -du 9223372036854775807 -l 0 -u 1e-18
# the loops over the squares at 2^63 - 1: the range starts at the largest
# square below it, k^2 for k = 3037000499, which is the only square in the
# range and carries a planted point, x = 1/k^2 on y^2 = x^3 + k^6 + 2
# (y = (k^6 + 1)/k^3); then with the leading coefficient 3 and the range
# starting at 3k^2 for k = 1753413056, the largest with 3k^2 <= 2^63 - 1,
# where x = 1/(3k^2) lies on y^2 = 3x^3 + 9k^6 + 2 and the loop over the
# plain squares finds none in the range and must stop at once.  Before
# the loops started at the first square in the range they ran through
# three thousand million useless squares here, six seconds
t '784637715410305245771861851903983139945652963971781747003 0 0 1' 9223372036854775807 -q -dl 9223372030926249001 -du 9223372036854775807 -l 0 -u 1e-18
t '261545905470885580590835526458680687693017248327589167106 0 0 3' 9223372036854775807 -q -dl 9223372034853777408 -du 9223372036854775807 -l 0 -u 1e-18
# a height bound just below and above 2^31, with few denominators
t '1 0 1' 2147483647 -q -dl 1 -du 3 -l 0 -u 1e-8
t '1 0 1' 2147483648 -q -dl 1 -du 3 -l 0 -u 1e-8
# a leading coefficient of eleven primes: 2^11 divisors, more than the
# table of squarefree divisors holds once the denominator bound admits
# them, so the restriction is not used
t '1 0 0 3710369067405' 100000 -q -l 0 -u 0.00005
t '1 0 0 3710369067405' 1000 -q

echo '---- output options ----'
t '1 0 126 0 441' 200 -q -i
t '1 0 126 0 441' 200 -q -i -I
t '1 0 126 0 441' 200 -q -1
t '1 0 126 0 441' 200 -q -1 -i
t '3 2 -1 12' 100 -q -1
t '3 2 -1 12' 100 -q -1 -i
t '5 3 7' 100 -q -1
t '1 0 0 1' 100 -q -1
t '1 0 0 1' 100 -q -1 -i
t '1 0 1' 100 -q -1 -i
t '5 0 0 1065023' 100 -q -1 -i
t '1 0 126 0 441' 200 -q -1 -x -n 15 -N 30 -p 30
t '2 3 1' 50 -q -1
t '-1 0 -1' 50 -q -1
t '1 0 126 0 441' 200 -q -y
t '1 0 126 0 441' 200 -q -y -Y
t '1 0 126 0 441' 200 -q -x -n 15 -N 30 -p 30
t '1 0 126 0 441' 200 -q -x -y -n 15 -N 30 -p 30
f 'grep found' '1 0 126 0 441' 200 -z
f 'grep found' '1 0 126 0 441' 200 -z -i
t '1 0 126 0 441' 200 -z -Z -q
f 'grep found' '1 0 126 0 441' 200 -z -y
f 'grep found' '1 0 126 0 441' 200 -z -x -n 15 -N 30 -p 30
f 'grep found' '1 0 126 0 441' 200 -z -1
f 'grep found' '-1 0 -1' 50 -z
f 'grep found' '2 0 3' 200 -z
t '1 0 126 0 441' 200 -q -f '%x/%z\n'
t '1 0 126 0 441' 200 -q -f '[%x, %y, %z] '
t '1 0 126 0 441' 200 -q -f '%x\t%y\t%z\n'
t '1 0 126 0 441' 200 -q -f '%x %z%' -y
t '1 0 126 0 441' 200 -q -f '%x %z\'
t '1 0 126 0 441' 200 -q -f '%q %x \% \\ \q %z\n'
t '1 0 126 0 441' 200 -q -f '%y\n' -y
t '1 0 126 0 441' 200 -q -fs '[' -fm ', ' -fe ']' -f '(%x:%y:%z)'
t '1 0 126 0 441' 200 -q -fs 'points:\n' -fe 'end\n' -fm ''
t '1 0 126 0 441' 200 -q -fs 'start\t\\\%\n' -fe '\'
t '1 0 126 0 441' 200 -q -fs 'a\qb\n' -fm '\' -fe '\n'
t '+1 +0 +126 +0 +441' 200 -q -f '%x %z\n'
t '-1 0 -1' 50 -q -fs '[' -fe ']\n'
t '1 0 126 0 441' 200 -q -f '%x %y %z\n' -fs '{' -fe '}\n' -1
# without -q: the banner and the messages after the run
m '1 0 126 0 441' 200
m '1 0 126 0 441' 200 -l 0 -u 1
m '36 0 -49 0 14 0 -1' 100 -l -2.5 -u -1.5 -l 1.5 -u 2.5
m '36 0 -49 0 14 0 -1' 100 -l -2.5 -u 2.5 -z
m '0 1 2' 40
m '1 0 0 1' 50
m '3 2 -1 12' 50
m '5 3 7' 50 -j
m '5 3 7' 50 -x -n 15 -N 30 -p 30
m '5 3 7' 50 -y
m '5 3 7' 50 -z -y
m '-1 0 -1' 50
m '2 0 3' 20
# the curve equation as the program prints it
f "grep 'Curve equation'" '-1 1 -1 1' 20
f "grep 'Curve equation'" '0 -2 0 1' 20
f "grep 'Curve equation'" '1 0 0 0 0 -1' 20
f "grep 'Curve equation'" '-12 0 0 7' 20
f "grep 'Curve equation'" '1 -1' 20
f "grep 'Curve equation'" '0 1' 20

echo '---- errors ----'
e
e '1 2 3'
e '1 2 3' abc
e '1 2 3' 0
e '1 2 3' -5
e '1 2 3' 10 -w
e '1 2 3' 10 -n
e '1 2 3' 10 -n x
e '1 2 3' 10 -N
e '1 2 3' 10 -p
e '1 2 3' 10 -F
e '1 2 3' 10 -f
e '1 2 3' 10 -fs
e '1 2 3' 10 -fm
e '1 2 3' 10 -fe
e '1 2 3' 10 -fx
e '1 2 3' 10 -dl
e '1 2 3' 10 -du
e '1 2 3' 10 -dx 1
e '1 2 3' 10 -l
e '1 2 3' 10 -u
e '1 2 3' 10 -l 1 -l 2
e '1 2 3' 10 -u 1 -u 2
e '1 2 3' 10 -l 2 -u 1
e '1 2 3' 10 -l 1 -u 2 -l 1.5 -u 3
e '1 2 3' 10 -l x
e '1 2 3' 10 -u x
e '1 2 3' 10 -S x
e '1 2 3' 10 -p x
e '1 2 3' 10 -F x
e '1 2 3' 10 -N x
e '1 2 3' 10 -dl x
e '1 2 3' 10 -du x
t '1 2 3' 10 -q -S -k
e '1 2 3' 10 -q q
e 'x' 10
e '1 2 x' 10
e '' 10
e '   ' 10
e '5' 10
e '5 0' 10
e '1 2 0 0' 20 -q
e '1 2 0 0' 20 -q -k
e '1 2 0 0' 20 -q -s
# a hundred and one intervals
set --
i=0
while [ $i -le 100 ]; do set -- "$@" -l $((3*i)) -u $((3*i+1)); i=$((i+1)); done
e '1 2 3' 1000 "$@"
set --

echo '==== part 2: the report -v prints ===='
v '1 0 126 0 441' 100 -v -n 3 -N 5
v '1 0 126 0 441' 100 -v -s -k -j -x
v '1 0 126 0 441' 100 -v -dl 2 -du 50
v '1 0 126 0 441' 100 -v -l 0 -u 1
v '0 1 2' 40 -v
v '2 3 1' 50 -v
v '0 1 0 2' 50 -v
v '1 0 0 1' 50 -v
v '3 2 -1 12' 50 -v
v '5 0 0 1065023' 50 -v
v '1 0 0 3710369067405' 100000 -v -l 0 -u 0.00005
v '-1 0 -1' 50 -v
v '-1000000000000 0 0 0 1' 100 -v
v '-1000000 0 0 0 3' 100 -v
v '2 0 3' 200 -v
v '2 0 1 0 2' 100 -v
v '2 131 4 0 2' 100 -v -n 5 -N 10
v '2 63 1 0 2' 100 -v
v '2007238469666518094547220599513022568322942623866 2 1 2 2 0 1' 100 -v -n 3 -N 5
v '1 0 0 1' 100 -v -dl 5 -du 8
v '5 3 7' 100 -v -n 3 -N 5
v '5 3 7' 100 -v -n 3 -N 5 -j
v '1 1 9' 100 -v -n 3 -N 5
v '1 1 27' 100 -v -n 3 -N 5
v '1 3 27' 100 -v -n 3 -N 5
v '1 1 18' 100 -v -n 3 -N 5
v '1 2 3' 100 -v -n 3 -N 5
v '-36 -54 -5' 100 -v -n 3 -N 5
v '44 -2 -12' 100 -v -n 3 -N 5
v '32 48 55 -44' 100 -v -n 3 -N 5
v '1 2 1' 20 -v
v '0 0 1' 20 -v
v '10 10 5 -7 0 3 -2' 300 -v -n 4 -N 6
# the bits set per word, the one number of the report that depends on the
# curve alone
f "grep 'bits set per word' | sed 's/,.*//'" '1 0 126 0 441' 100 -v
f "grep 'bits set per word' | sed 's/,.*//'" '3 2 -1 12' 100 -v
f "grep 'bits set per word' | sed 's/,.*//'" '4 28 -53' 100 -v
f 'grep -c .' '1 0 126 0 441' 100 -v -q
f 'grep -c .' '1 0 126 0 441' 100 -q -v
