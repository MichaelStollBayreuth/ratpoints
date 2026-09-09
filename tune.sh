#!/bin/sh
# Find good values for the two machine-dependent constants that decide how
# many primes each sieving stage uses -- RATPOINTS_SURVIVORS_PER_WORD and
# RATPOINTS_SP2_EXTRA in ratpoints.h -- and write them to tuning.mk, which
# the Makefile includes.  Run it as "make tune"; it needs ./rptest and
# ./rptest-many, and it modifies no source file.
#
# Both tests are used throughout: they cover the two regimes that occur in
# practice (random curves, and curves with many rational points), and a value
# that suits one of them can be a poor choice for the other.
#
# Which pair of tests, and at what height bound, is set by TUNE_TESTS and
# TUNE_HEIGHT below.  The default is the pair "make test" uses, at their own
# height of 16383.  That is a short run in which a fifth of the time goes into
# building sieve tables rather than into sieving, so a setting is judged partly
# on work it does not affect; "make tunehigh" measures the same thing on the
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
# R_VALUES and E_VALUES are an absolute ladder, bracketing the compiled-in
# 0.0075 and 5 either way; a factor of two in the threshold is worth about one
# prime in the first phase.  R_FACTORS and E_DELTAS instead describe a
# neighbourhood of the settings being measured against -- multiples of the
# threshold and offsets added to the other constant -- and take precedence when
# they are set.
#
# Which to use depends on what a run costs.  "make tune" sweeps the ladder,
# since one timing there is three seconds.  "make tunehigh" costs two minutes a
# timing, so it starts from what "make tune" found and only asks whether a step
# either way is better, which is most of an hour saved.  The assumption is that
# the two regimes do not want wildly different values; if a neighbourhood run
# moves a value, it has not finished looking, and should be run again from
# there.
R_VALUES=${R_VALUES:-"0.003 0.005 0.012 0.02"}
E_VALUES=${E_VALUES:-"3 5 7 10"}
R_FACTORS=${R_FACTORS:-}
E_DELTAS=${E_DELTAS:-}

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
[ -n "$DEF_R" ] && [ -n "$DEF_E" ] || { echo "tune.sh: cannot read the defaults from ratpoints.h" >&2; exit 1; }
if [ -f tuning.mk ] && [ "`sed -n 's/^TUNED_FOR *= *//p' tuning.mk`" = "${TUNE_CONFIG:-}" ]
then
  v=`sed -n 's/.*RATPOINTS_SURVIVORS_PER_WORD=\([^ 	]*\).*/\1/p' tuning.mk`
  [ -n "$v" ] && DEF_R=$v
  v=`sed -n 's/.*RATPOINTS_SP2_EXTRA=\([^ 	]*\).*/\1/p' tuning.mk`
  [ -n "$v" ] && DEF_E=$v
  echo "tuning.mk is already in effect; measuring against its $DEF_R / $DEF_E"
fi
BASE="-r $DEF_R -R $DEF_E"

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
[ -n "$R_FACTORS$E_DELTAS" ] && echo "candidates: $R_VALUES/ $E_VALUES"

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
  printf 'r=%s\t-r %s -R %s\n' "$v" "$v" "$DEF_E" >> "$TMP/c1"
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
  printf 'e=%s\t-r %s -R %s\n' "$v" "$BEST_R" "$v" >> "$TMP/c2"
done
measure "$TMP/c2" "$TMP/r2"
report "$TMP/r2"

BEST_LBL=`best_of "$TMP/r2" current`
BEST_E=`echo "$BEST_LBL" | sed 's/^e=//'`
BEST=`awk -v k="$BEST_LBL" '$1==k{print $2}' "$TMP/r2"`
# how far the current settings measured from themselves, in either stage: both
# are the same comparison of a binary with itself, so either one being far from
# 1 means the machine could not be measured on
SELF=`awk '$1=="current"{ d = $2 - 1; if (d < 0) d = -d; if (d > w) w = d }
           END { printf "%.5f", 1 + w }' "$TMP/r1" "$TMP/r2"`

echo
verdict=`awk -v b="$BEST" -v s="$SELF" -v m="$MARGIN" -v z="$NOISE" \
  'BEGIN { d = s - 1; if (d < 0) d = -d
           if (d > z) print "noisy"; else if (b < m) print "accept"; else print "keep" }'`

case $verdict in
  noisy)
    echo "The current settings measured `awk -v s=$SELF 'BEGIN{printf "%.1f", 100*(s-1)}'`% away from themselves, so this"
    echo "machine is too noisy right now to tell the settings apart.  Nothing"
    echo "written.  Try again when it is idle, or with 'ROUNDS=6 make tune'."
    ;;
  keep)
    echo "Nothing beat the current settings ($DEF_R, $DEF_E) by the required"
    echo "`awk -v m=$MARGIN 'BEGIN{printf "%.0f", 100*(1-m)}'`%, so they are kept and nothing is written."
    ;;
  accept)
    cat > tuning.mk <<EOF
# Machine-dependent tuning, written by "make tune" on `date +%Y-%m-%d`.
# Delete this file to go back to the values compiled into ratpoints.h.
# TUNED_FOR records the configuration it was measured for; the Makefile
# ignores this file if the configuration has changed since.
# Measured on $TUNE_TESTS${TUNE_HEIGHT:+ at height $TUNE_HEIGHT}, starting
# from $DEF_R / $DEF_E.  That is a note to the reader, not something the
# Makefile looks at: "make tune" and "make tunehigh" write the same file and
# each takes the other's result as its starting point.
TUNED_FOR = ${TUNE_CONFIG:-unknown}
TUNEFLAGS = -DRATPOINTS_SURVIVORS_PER_WORD=$BEST_R -DRATPOINTS_SP2_EXTRA=$BEST_E
EOF
    echo "Wrote tuning.mk: threshold $BEST_R, offset $BEST_E --"
    echo "`awk -v b=$BEST 'BEGIN{printf "%.1f", 100*(1-b)}'`% better than the current $DEF_R / $DEF_E."
    echo "Run 'make all' to rebuild with it; delete tuning.mk to discard it."
    ;;
esac
