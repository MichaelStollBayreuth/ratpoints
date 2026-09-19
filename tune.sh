#!/bin/sh
# Find good values for the five machine-dependent constants that decide how
# many primes each sieving stage uses and which ones --
# RATPOINTS_SURVIVORS_PER_WORD, RATPOINTS_SP2_EXTRA, RATPOINTS_SP2_U0,
# RATPOINTS_COST_TABLE and RATPOINTS_SP3_PER_DENOM in ratpoints.h -- and write
# them to tuning.mk, which the Makefile includes.  Run it as "make tune"; it needs ./rptest and
# ./rptest-many, and it modifies no source file.
#
# Both tests are used throughout: they cover the two regimes that occur in
# practice (random curves, and curves with many rational points), and a value
# that suits one of them can be a poor choice for the other.
#
# Which pair of tests, and at what height bound, is set by TUNE_TESTS and
# TUNE_HEIGHT below.  The default is the pair "make test" uses, at their own
# height of 16383.  That is a short run in which two fifths of the time goes
# into work other than the two sieving phases (choosing the primes, the set-up
# per denominator, the exact checks), so a setting is judged partly on work it
# does not affect; "make tunehigh" measures the same thing on the
# large-height suites instead, where the sieve is nearly all of it.  Use that
# one if the runs that matter are long.  The two write the same tuning.mk and
# each starts from what the other left, so the cheap sweep can be run first and
# the expensive one asked only whether it wants to move.
#
# Measuring is the delicate part, and naive timing does not work.  The cost
# surface is flat -- anything within a factor of two of a good threshold costs
# under 3% -- while a laptop under sustained load drifts by 25% as it heats
# up, so a straight sweep measures the clock, not the settings.  Therefore:
#
#   * a warm-up brings the machine to a steady state before anything counts;
#   * every candidate is timed back to back with the current settings and
#     only the ratio of the two is kept, so drift cancels;
#   * the median over rounds is used, not a single reading;
#   * the current settings are themselves among the candidates, so their
#     measured ratio to themselves says how much noise is left; if that is
#     more than NOISE, the run is declared inconclusive;
#   * nothing is written unless the winner beats the current settings by more
#     than MARGIN, which keeps a noisy run from making the build slower.
#
# Takes a few minutes.  Do not run it under "make -j" or on a busy machine.
set -u

ROUNDS=${ROUNDS:-3}
MARGIN=${MARGIN:-0.97}     # accept only if the best ratio is below this
NOISE=${NOISE:-0.02}       # ... and the baseline's self-ratio is within this
WARMUP=${WARMUP:-20}       # seconds of load before measuring

# The candidates for each constant.  There are two ways to say what they are.
# R_VALUES, E_VALUES, U_VALUES, C_VALUES and Q_VALUES are an absolute ladder,
# bracketing the compiled-in values either way; a factor of two in the
# threshold is worth about one prime in the first phase.  R_FACTORS, E_DELTAS,
# U_FACTORS, C_FACTORS and Q_FACTORS instead describe a neighbourhood of the
# settings being measured against -- multiples of the threshold, of the run
# length, of the table cost and of the third stage's per-denominator cost, and
# offsets added to the number of primes -- and take precedence when they are
# set.
#
# Which to use depends on what a run costs.  "make tune" sweeps the ladder,
# since one timing there is three seconds.  "make tunehigh" costs two minutes a
# timing, so it starts from what "make tune" found and only asks whether a step
# either way is better: fifteen settings a round rather than twenty-four, which
# is more than half an hour off a run of two.  The assumption is that the two
# regimes do not want
# wildly different values; if a neighbourhood run
# moves a value, it has not finished looking, and should be run again from
# there.
#
# The offset is the third constant's business too.  Since version 3.0.0 it is
# the number of extra primes an arbitrarily long run wants, and a run of U
# numerator words gets sp2_extra/(1 + U0/U) of them, so the ladder for it is
# wider than it used to be and the value that suits a short run is not the
# value that suits a long one.  That is the point: U0 is what reconciles the
# two, and it is pinned by running "make tune" and "make tunehigh" one after
# the other, since they see very different U.
#
# The fourth constant is what one row of a sieve table costs, in units of one
# first-phase AND per word.  It decides how strongly the ranking of the primes
# prefers a small prime to a larger one that says a little more, which matters
# at a small height bound and hardly at all at a large one, where the tables
# are a negligible share of the run.  Its basin is flat -- on the machine the
# compiled-in value was measured on, anything within a factor of two of it was
# within half a per cent -- but it is the ratio of two different kinds of
# work, a scalar table build against a vector AND, and that ratio is the one
# most likely to differ between machines; it also depends on the register
# width, since the AND handles more words at once as the registers grow.
#
# The fifth constant is what carrying one third-stage prime costs per
# denominator, as a fraction of one exact check.  It decides whether the third
# stage runs at all on a curve (ratpoints.h, RATPOINTS_SP3_PER_DENOM); since
# 3.0.0 the stage's set-up is done on demand, for the denominators that bring a
# survivor that far, so the cost is smaller than it was when the constant was
# first fitted, and the ladder reaches down accordingly.
R_VALUES=${R_VALUES:-"0.0015 0.002 0.0045 0.0075"}
E_VALUES=${E_VALUES:-"4 6 9 13 18"}
U_VALUES=${U_VALUES:-"3e5 6e5 2.4e6 5e6"}
C_VALUES=${C_VALUES:-"10 20 70 140"}
Q_VALUES=${Q_VALUES:-"0.003 0.006 0.025 0.05"}
R_FACTORS=${R_FACTORS:-}
E_DELTAS=${E_DELTAS:-}
U_FACTORS=${U_FACTORS:-}
C_FACTORS=${C_FACTORS:-}
Q_FACTORS=${Q_FACTORS:-}

# The suites to tune on, as "program:reference" pairs, and the height bound to
# run them at (empty: each test's own default).  Set by "make tune" and
# "make tunehigh"; see the Makefile.
TUNE_TESTS=${TUNE_TESTS:-"./rptest:testbase ./rptest-many:testbase-many"}
TUNE_HEIGHT=${TUNE_HEIGHT:-}
[ -n "$TUNE_HEIGHT" ] && HFLAG="-h $TUNE_HEIGHT" || HFLAG=""

for t in $TUNE_TESTS; do
  for f in "${t%%:*}" "${t#*:}"; do
    [ -e "$f" ] || { echo "tune.sh: $f is missing; build it first" >&2; exit 1; }
  done
done
[ -e ratpoints.h ] || { echo "tune.sh: ratpoints.h is missing" >&2; exit 1; }

# The settings to measure against.  These are passed explicitly to every run,
# including the baseline, so that the comparison never depends on what happens
# to be compiled into the tests -- which matters when tuning.mk is already in
# effect from an earlier run, since then the compiled-in values are its
# values, not the ones in ratpoints.h.
DEF_R=`sed -n 's/^# *define  *RATPOINTS_SURVIVORS_PER_WORD  *\([0-9.eE+-]*\).*/\1/p' ratpoints.h`
DEF_E=`sed -n 's/^# *define  *RATPOINTS_SP2_EXTRA  *\([0-9]*\).*/\1/p' ratpoints.h`
DEF_U=`sed -n 's/^# *define  *RATPOINTS_SP2_U0  *\([0-9.eE+-]*\).*/\1/p' ratpoints.h`
DEF_C=`sed -n 's/^# *define  *RATPOINTS_COST_TABLE  *\([0-9.eE+-]*\).*/\1/p' ratpoints.h`
DEF_Q=`sed -n 's/^# *define  *RATPOINTS_SP3_PER_DENOM  *\([0-9.eE+-]*\).*/\1/p' ratpoints.h`
[ -n "$DEF_R" ] && [ -n "$DEF_E" ] && [ -n "$DEF_U" ] && [ -n "$DEF_C" ] && [ -n "$DEF_Q" ] \
  || { echo "tune.sh: cannot read the defaults from ratpoints.h" >&2; exit 1; }
if [ -f tuning.mk ] && [ "`sed -n 's/^TUNED_FOR *= *//p' tuning.mk`" = "${TUNE_CONFIG:-}" ]
then
  v=`sed -n 's/.*RATPOINTS_SURVIVORS_PER_WORD=\([^ 	]*\).*/\1/p' tuning.mk`
  [ -n "$v" ] && DEF_R=$v
  v=`sed -n 's/.*RATPOINTS_SP2_EXTRA=\([^ 	]*\).*/\1/p' tuning.mk`
  [ -n "$v" ] && DEF_E=$v
  v=`sed -n 's/.*RATPOINTS_SP2_U0=\([^ 	]*\).*/\1/p' tuning.mk`
  [ -n "$v" ] && DEF_U=$v
  v=`sed -n 's/.*RATPOINTS_COST_TABLE=\([^ 	]*\).*/\1/p' tuning.mk`
  [ -n "$v" ] && DEF_C=$v
  v=`sed -n 's/.*RATPOINTS_SP3_PER_DENOM=\([^ 	]*\).*/\1/p' tuning.mk`
  [ -n "$v" ] && DEF_Q=$v
  echo "tuning.mk is already in effect; measuring against its $DEF_R / $DEF_E / $DEF_U / $DEF_C / $DEF_Q"
fi
BASE="-r $DEF_R -R $DEF_E -U $DEF_U -C $DEF_C -Q $DEF_Q"

# a neighbourhood of those, if that is what was asked for
if [ -n "$R_FACTORS" ]; then
  R_VALUES=`awk -v r="$DEF_R" -v f="$R_FACTORS" \
    'BEGIN { n = split(f, a, " ")
             for (i = 1; i <= n; i++) printf "%.4g ", r*a[i] }'`
fi
if [ -n "$E_DELTAS" ]; then
  E_VALUES=`awk -v e="$DEF_E" -v d="$E_DELTAS" \
    'BEGIN { n = split(d, a, " ")
             for (i = 1; i <= n; i++) { v = e + a[i]
                                        if (v >= 0) printf "%d ", v } }'`
fi
if [ -n "$U_FACTORS" ]; then
  U_VALUES=`awk -v u="$DEF_U" -v f="$U_FACTORS" \
    'BEGIN { n = split(f, a, " ")
             for (i = 1; i <= n; i++) printf "%.4g ", u*a[i] }'`
fi
if [ -n "$C_FACTORS" ]; then
  C_VALUES=`awk -v c="$DEF_C" -v f="$C_FACTORS" \
    'BEGIN { n = split(f, a, " ")
             for (i = 1; i <= n; i++) printf "%.4g ", c*a[i] }'`
fi
if [ -n "$Q_FACTORS" ]; then
  Q_VALUES=`awk -v q="$DEF_Q" -v f="$Q_FACTORS" \
    'BEGIN { n = split(f, a, " ")
             for (i = 1; i <= n; i++) printf "%.4g ", q*a[i] }'`
fi

# Each stage carries the winners of the stages before it into every one of
# its candidates, so the constant it sweeps must have its current value among
# them as well: otherwise the combination "earlier winners, this constant
# unchanged" is never timed and can never win, and a threshold that won
# stage 1 by a clear margin could be lost again in stage 4 for want of a
# candidate.  (Stage 1 needs nothing: its current value is "current".)
case " $E_VALUES " in *" $DEF_E "*) ;; *) E_VALUES="$DEF_E $E_VALUES" ;; esac
case " $U_VALUES " in *" $DEF_U "*) ;; *) U_VALUES="$DEF_U $U_VALUES" ;; esac
case " $C_VALUES " in *" $DEF_C "*) ;; *) C_VALUES="$DEF_C $C_VALUES" ;; esac
case " $Q_VALUES " in *" $DEF_Q "*) ;; *) Q_VALUES="$DEF_Q $Q_VALUES" ;; esac
[ -n "$R_FACTORS$E_DELTAS$U_FACTORS$C_FACTORS$Q_FACTORS" ] \
  && echo "candidates: $R_VALUES/ $E_VALUES/ $U_VALUES/ $C_VALUES/ $Q_VALUES"

if command -v taskset >/dev/null 2>&1; then PIN="taskset -c 0"; else PIN=""; fi

TMP=`mktemp -d` || exit 1
trap 'rm -rf "$TMP"' EXIT INT TERM

echo "checking that the tests still pass ..."
for t in $TUNE_TESTS; do
  prog=${t%%:*}; base=${t#*:}
  $prog $HFLAG | cmp -s - "$base" \
    || { echo "tune.sh: $prog disagrees with $base" >&2; exit 1; }
done

# one timing = every test in TUNE_TESTS, so that a setting is judged on all of
# the regimes at once and not just on the one it happens to suit
time_tests() {
  tot=0
  for t in $TUNE_TESTS; do
    tt=`$PIN ${t%%:*} $HFLAG $1 -z -T` || exit 1
    tot=`awk -v a="$tot" -v b="$tt" 'BEGIN{printf "%.6f", a+b}'`
  done
  echo "$tot"
}

printf 'warming up (%ss) ' "$WARMUP"
end=`expr \`date +%s\` + $WARMUP`
while [ `date +%s` -lt $end ]; do time_tests "$BASE" > /dev/null; printf '.'; done
echo

# measure(): "label<TAB>args" lines from $1 -> "label median_ratio" in $2
measure() {
  : > "$TMP/raw"
  n=`wc -l < "$1"`
  r=1
  while [ $r -le $ROUNDS ]; do
    printf '  round %d/%d:' $r $ROUNDS
    off=`expr \( $r - 1 \) % $n`
    { tail -n +`expr $off + 1` "$1"
      [ $off -gt 0 ] && head -n $off "$1"
      true
    } > "$TMP/order"
    while IFS='	' read -r label args; do
      tc=`time_tests "$args"`
      tb=`time_tests "$BASE"`     # the current settings, right next to it
      echo "$label `awk -v c="$tc" -v b="$tb" 'BEGIN{printf "%.5f", c/b}'`" >> "$TMP/raw"
      printf ' .'
    done < "$TMP/order"
    echo
    r=`expr $r + 1`
  done
  sort "$TMP/raw" | awk '
    { v[$1] = v[$1] " " $2; k[$1] = 1 }
    END { for (a in k) { m = split(v[a], x, " ")
                         for (i = 1; i < m; i++) for (j = i+1; j <= m; j++)
                           if (x[j]+0 < x[i]+0) { t = x[i]; x[i] = x[j]; x[j] = t }
                         # true median: for an even number of rounds take the
                         # mean of the two central values, not the lower one --
                         # with ROUNDS=2 the lower one is the minimum, which
                         # picks up exactly the outliers this is meant to reject
                         med = (m % 2) ? x[int((m+1)/2)] \
                                       : (x[m/2] + x[m/2+1]) / 2
                         printf "%s %.5f\n", a, med } }' | sort > "$2"
}

report() { awk '{ printf "    %-14s %+7.1f%%\n", $1, 100*($2-1) }' "$1" | sort; }
best_of() { awk -v skip="${2:-}" '$1 != skip { if (m == "" || $2 < m) { m = $2; k = $1 } }
                                  END { print k }' "$1"; }

echo
echo "stage 1: the threshold, against the current $DEF_R (offset stays at $DEF_E)"
: > "$TMP/c1"
printf 'current\t%s\n' "$BASE" >> "$TMP/c1"
for v in $R_VALUES; do
  [ "$v" = "$DEF_R" ] && continue
  printf 'r=%s\t-r %s -R %s -U %s -C %s -Q %s\n' "$v" "$v" "$DEF_E" "$DEF_U" "$DEF_C" "$DEF_Q" >> "$TMP/c1"
done
measure "$TMP/c1" "$TMP/r1"
report "$TMP/r1"

BEST_R=`best_of "$TMP/r1"`
case $BEST_R in current) BEST_R=$DEF_R ;; r=*) BEST_R=`echo "$BEST_R" | sed 's/^r=//'` ;; esac

echo
echo "stage 2: the offset, with the threshold at $BEST_R"
: > "$TMP/c2"
printf 'current\t%s\n' "$BASE" >> "$TMP/c2"
for v in $E_VALUES; do
  printf 'e=%s\t-r %s -R %s -U %s -C %s -Q %s\n' "$v" "$BEST_R" "$v" "$DEF_U" "$DEF_C" "$DEF_Q" >> "$TMP/c2"
done
measure "$TMP/c2" "$TMP/r2"
report "$TMP/r2"

BEST_E=`best_of "$TMP/r2"`
case $BEST_E in current) BEST_E=$DEF_E ;; e=*) BEST_E=`echo "$BEST_E" | sed 's/^e=//'` ;; esac

echo
echo "stage 3: the run length at which a second-stage prime pays for itself,"
echo "         with the threshold at $BEST_R and the offset at $BEST_E"
: > "$TMP/c3"
printf 'current\t%s\n' "$BASE" >> "$TMP/c3"
for v in $U_VALUES; do
  printf 'u=%s\t-r %s -R %s -U %s -C %s -Q %s\n' "$v" "$BEST_R" "$BEST_E" "$v" "$DEF_C" "$DEF_Q" >> "$TMP/c3"
done
measure "$TMP/c3" "$TMP/r3"
report "$TMP/r3"

BEST_U=`best_of "$TMP/r3"`
case $BEST_U in current) BEST_U=$DEF_U ;; u=*) BEST_U=`echo "$BEST_U" | sed 's/^u=//'` ;; esac

echo
echo "stage 4: what a row of a sieve table costs, with the threshold at $BEST_R,"
echo "         the offset at $BEST_E and the run length at $BEST_U"
: > "$TMP/c4"
printf 'current\t%s\n' "$BASE" >> "$TMP/c4"
for v in $C_VALUES; do
  printf 'c=%s\t-r %s -R %s -U %s -C %s -Q %s\n' "$v" "$BEST_R" "$BEST_E" "$BEST_U" "$v" "$DEF_Q" >> "$TMP/c4"
done
measure "$TMP/c4" "$TMP/r4"
report "$TMP/r4"

BEST_C=`best_of "$TMP/r4"`
case $BEST_C in current) BEST_C=$DEF_C ;; c=*) BEST_C=`echo "$BEST_C" | sed 's/^c=//'` ;; esac

echo
echo "stage 5: what carrying a third-stage prime costs per denominator, with the"
echo "         threshold at $BEST_R, the offset at $BEST_E, the run length at $BEST_U"
echo "         and the table cost at $BEST_C"
: > "$TMP/c5"
printf 'current\t%s\n' "$BASE" >> "$TMP/c5"
for v in $Q_VALUES; do
  printf 'q=%s\t-r %s -R %s -U %s -C %s -Q %s\n' "$v" "$BEST_R" "$BEST_E" "$BEST_U" "$BEST_C" "$v" >> "$TMP/c5"
done
measure "$TMP/c5" "$TMP/r5"
report "$TMP/r5"

BEST_LBL=`best_of "$TMP/r5" current`
BEST_Q=`echo "$BEST_LBL" | sed 's/^q=//'`
BEST=`awk -v k="$BEST_LBL" '$1==k{print $2}' "$TMP/r5"`
# how far the current settings measured from themselves, in any stage: all
# are the same comparison of a binary with itself, so any one being far from
# 1 means the machine could not be measured on
SELF=`awk '$1=="current"{ d = $2 - 1; if (d < 0) d = -d; if (d > w) w = d }
           END { printf "%.5f", 1 + w }' "$TMP/r1" "$TMP/r2" "$TMP/r3" "$TMP/r4" "$TMP/r5"`

echo
verdict=`awk -v b="$BEST" -v s="$SELF" -v m="$MARGIN" -v z="$NOISE" \
  'BEGIN { d = s - 1; if (d < 0) d = -d
           if (d > z) print "noisy"; else if (b < m) print "accept"; else print "keep" }'`

# A candidate that carries every constant at its current value is the current
# settings measured a second time; if noise puts that row past the margin, it
# is not an improvement and must not be announced as one.
[ "$verdict" = accept ] && [ "$BEST_R" = "$DEF_R" ] && [ "$BEST_E" = "$DEF_E" ] \
  && [ "$BEST_U" = "$DEF_U" ] && [ "$BEST_C" = "$DEF_C" ] && [ "$BEST_Q" = "$DEF_Q" ] \
  && verdict=keep

case $verdict in
  noisy)
    echo "The current settings measured `awk -v s=$SELF 'BEGIN{printf "%.1f", 100*(s-1)}'`% away from themselves, so this"
    echo "machine is too noisy right now to tell the settings apart.  Nothing"
    echo "written.  Try again when it is idle, or with 'ROUNDS=6 make tune'."
    ;;
  keep)
    echo "Nothing beat the current settings ($DEF_R, $DEF_E, $DEF_U, $DEF_C, $DEF_Q) by the required"
    echo "`awk -v m=$MARGIN 'BEGIN{printf "%.0f", 100*(1-m)}'`%, so they are kept and nothing is written."
    ;;
  accept)
    cat > tuning.mk <<EOF
# Machine-dependent tuning, written by "make tune" on `date +%Y-%m-%d`.
# Delete this file to go back to the values compiled into ratpoints.h.
# TUNED_FOR records the configuration it was measured for; the Makefile
# ignores this file if the configuration has changed since.
# Measured on $TUNE_TESTS${TUNE_HEIGHT:+ at height $TUNE_HEIGHT}, starting
# from $DEF_R / $DEF_E / $DEF_U / $DEF_C / $DEF_Q.  That is a note to the reader,
# not something the Makefile looks at: "make tune" and "make tunehigh" write
# the same file and each takes the other's result as its starting point.
TUNED_FOR = ${TUNE_CONFIG:-unknown}
TUNEFLAGS = -DRATPOINTS_SURVIVORS_PER_WORD=$BEST_R -DRATPOINTS_SP2_EXTRA=$BEST_E -DRATPOINTS_SP2_U0=$BEST_U -DRATPOINTS_COST_TABLE=$BEST_C -DRATPOINTS_SP3_PER_DENOM=$BEST_Q
EOF
    echo "Wrote tuning.mk: threshold $BEST_R, offset $BEST_E, run length $BEST_U,"
    echo "table cost $BEST_C, third-stage cost $BEST_Q: `awk -v b=$BEST 'BEGIN{printf "%.1f", 100*(1-b)}'`% better than the current"
    echo "$DEF_R / $DEF_E / $DEF_U / $DEF_C / $DEF_Q."
    echo "Run 'make all' to rebuild with it; delete tuning.mk to discard it."
    ;;
esac
