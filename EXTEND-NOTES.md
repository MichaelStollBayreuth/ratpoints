# The prime-extension rule and the run-length estimate

Branch `extend-rule` off `v2.3` at 6503217, 2026-09-16 evening.  TODO item
27, which came out of the skeptic's review of item 26 (the mask modulo 64),
together with the run_shape weighting left over from item 18 (review F6).
Measured on the i7-1355U, pinned to P-core 0, in core cycles, paired and
interleaved with the base, medians of the per-round ratios, outputs
compared with the references before anything was timed (`pair.sh`).

## What was wrong

Two things in `sieving_info`'s model of the run, both older than item 26
and both exposed by it.

**The third stage could be starved.**  `sieving_info` looks at
`num_primes` (30) primes, drops the uninformative ones (f a square modulo
every residue: on a point-rich curve that is most of the small primes),
and looks at more only in two cases: in the main loop, while the primes
in hand are fewer than the first two phases want (`pnp < s1 + sp2_extra`);
and in the third stage's selection, when the pool is *exhausted*
(`sp3 >= pnp`).  The third stage takes primes from the pool in order of
their density r while a prime pays for itself, `S (1 - per_surv - r) >
per_denom`, and *breaks at the first that does not* -- without asking for
more.  So when the pool still holds a few poor primes (r near 1) the stage
stops, although further primes beyond `pn_lim` would have paid.  The mask
of item 26 lowered `bits_per_word`, hence `s1` by one on the curve
`456976 -448032 -255200 208380 61033 -12834 81`, hence `want` by one, so
the main loop stopped at 35 primes instead of 43; the pool left for the
third stage held one poor prime, the stage took nothing (29 primes instead
of 37 in all) and did 254116 exact checks instead of 5573 at height 200000,
1.85 per cent more instructions on that curve.  The stage's appetite was
there; the rule had no way to feed it.

**`U` over-counted the numerators of even denominators.**  `run_shape`
estimates the words of numerators the run sweeps as the sum over the
denominators kept of their numerator intervals, halved when only one
parity of numerator is in play (`which_bits != num_all`).  With `num_all`
the even denominators still sweep only the odd numerators (`sift()` forces
`num_odd` for them), which run_shape did not model: `U` was over-estimated
by up to a third on such curves -- and `U` is what the per-word costs of a
prime are spread over, and what `S`, the survivors per denominator, is
derived from.

## What was done

**The third stage's selection** (`sieving_info`, find_points.c).  The
loop that takes primes for the stage picks the best prime not yet spoken
for; when there is none, or the best does not pay for itself, it now looks
at a further prime -- as long as the stage still has appetite for one: a
prime of density `r_typ`, the mean density of the informative primes seen
so far on this curve, must pay by the stage's own rule, `S (1 - per_surv -
r_typ) > per_denom`.  Then it chooses again with the new prime in the
pool.  The two old cases -- pool exhausted, extend; best prime not worth
it, break -- are the two halves of this; what is new is that a poor prime
in hand no longer ends the search when better ones are to be had, and
that an empty pool is not refilled when even a typical prime would not
pay.  The appetite test is what keeps the cost in check: `S`, the
survivors per denominator when the stage begins, is small at small height
bounds, so nothing extra is looked at there; at height 200000 on a
point-rich curve it can look at every prime up to the bound (53 with
`PRIME_SIZE` 8; `examine_prime` costs some thousand instructions each,
nothing against such a run).  When the caller fixes `sp3` the stage takes
what it is told and looks further only when the pool is empty, as before.
On the curve that started this, `456976 -448032 -255200 208380 61033
-12834 81` at height 200000, the base looks at 35 primes and the third
stage takes none; now it looks at 41 and takes 6 (before the mask of
item 26 it was 43 and 7).

**run_shape's weighting.**  The three branches that count the classes of
the denominator kept (`b = k^2`, `b = d k^2`, all `b`) also count the even
ones, and the numerator total is scaled per class: an odd denominator
sweeps all numerators under `num_all` and half of them otherwise, an even
one always half.  The `[runshape]` instrumentation (`RP_PRIME_STATS` with
`RP_PHASE_TIMING`) prints the packing, so the prediction can be checked
per packing.

## What it is worth

Chain `chain-extend.sh` (18:03-18:20): the counter builds, `alt.sh`, then
`pair.sh`, d0681c4 against 6503217, 3 rounds of eight suites at the
default code placement -- no alignment variants, since the change is in
the set-up code and the sieve loops do not move.  Cycles new/base, medians
(raw reports in `review-2026-09-12/measurements-2026-09-16/extend/`):

| suite | cycles | [min, max] | instructions |
|---|---|---|---|
| test1 | 0.994 | [0.968, 1.029] | 0.998 |
| test1many | 0.987 | [0.977, 1.008] | 0.983 |
| testhigh | 0.996 | [0.993, 1.003] | 1.000 |
| testhighmany | 0.995 | [0.993, 1.016] | 0.996 |
| height 200 | 1.004 | [0.996, 1.004] | 1.000 |
| height 1000 | 1.002 | [0.999, 1.003] | 1.000 |
| height 4000 | 0.980 | [0.971, 1.005] | 0.999 |
| point-rich at 1000 | 1.005 | [0.982, 1.006] | 1.002 |

**Worth about 1% of `make test1many` and half a per cent of `make
testhighmany`; everything else within the noise of three rounds** (a point
and a half on the short suites).  Nothing gets slower: the small-height
suites move by less than half a per cent in either direction with
instruction counts within 0.2% of the base.  This is a correction of the
model, not an optimization, and its worth is in the counters.

The counters (`prelim-extend.txt`; `runshape.py` pairs the per-curve
lines of the statistics builds):

| suite | exact checks base -> new | primes looked at (mean) | third-stage primes (mean) |
|---|---|---|---|
| test1 | 48819 -> 51836 | 30.0 -> 30.0 | 16.0 -> 15.9 (74 curves lose one, 10 gain one) |
| test1many | 95953 -> 45029 | 33.0 -> 33.3 | 26.4 -> 26.8 (38 lose one or two, 34 gain one to six) |
| testdegrees | 2428 -> 2495 | 30.0 -> 30.0 | 16.6 -> 16.5 |
| height 1000 | 3576 -> 3765 | 30.0 -> 30.0 | 14.1 -> 14.0 |
| point-rich at 1000 | 7632 -> 8287 | 31.4 -> 30.7 | 23.6 -> 23.1 |
| testhigh | 299898 -> 296607 | 30.00 -> 30.01 | 20.6 -> 20.6 |
| testhighmany | 496629 -> 107529 | 44.3 -> 43.0 | 37.8 -> 37.5 (22 curves lose one, two gain six and seven) |

Two effects are mixed here.  The extension looks further where the stage
was starved: the two curves of testhighmany that gain six and seven
third-stage primes are the one that started this (35 -> 41 primes looked
at, 0 -> 6 taken; 254116 -> 5000-odd checks) and its neighbour, and they
are the whole of the suite's fall in checks -- 4.6 times fewer than
before, and well below the 250903 of before item 26.  The corrected `U`
is smaller on curves with even denominators, so `S`, the survivors per
denominator the third stage reckons with, is smaller too, and the stage
takes one prime less on a twelfth of the random curves: test1's checks
rise by 6% (3000 checks over 879 curves), which the cycles do not see.
Whether the third-stage constants, fitted with the inflated `U`, want
moving is for the tuning session after the sieve group.

**The prediction of `U`.**  `log(Uact/Upred)` over the curves, mean and
standard deviation, from the `[runshape]` lines (the base binary does not
print the packing; its per-packing means below are inferred from the new
binary's counts, since the change leaves `U` untouched for the packings
other than `num_all`):

| suite | base, all | new, all | new, num_all | new, others |
|---|---|---|---|---|
| test1 | +0.006 (sd 0.146) | +0.071 (sd 0.077) | +0.062 (n=454) | +0.075 to +0.095 (n=425) |
| testhigh | -0.016 (sd 0.141) | +0.049 (sd 0.070) | +0.048 | +0.045 to +0.058 |
| test1many | -0.214 (sd 0.171) | +0.009 (sd 0.160) | +0.029 (n=92) | -0.296 (n=6) |
| testhighmany | -0.181 (sd 0.092) | +0.081 (sd 0.070) | +0.081 (n=30) | -- |
| testdegrees | -0.154 (sd 0.240) | -0.044 (sd 0.155) | -0.022 (n=41) | -0.02 to -0.32 |

The spread halves on the random curves (0.146 -> 0.077 in the log, at
both heights), and the `num_all` curves, which the base over-predicted by
6% at 16383 while the other packings ran 8% above their prediction, now
sit where the others do: the base's mean near zero was two biases
cancelling.  On the point-rich suites, all `num_all`, the over-prediction
of 20% becomes an under-prediction of 1 to 8%.  What remains is common to
every packing: `U` is under-predicted by 5 to 8% on random curves at
both heights, which is the Jacobi factor -- run_shape says the symbol
lets half the denominators through, the exact count (item 24's review)
is 0.53 -- and at height 1000 by 60%, where each denominator's few
arrays are rounded up to whole ones and the estimate counts words.  Both
are older than this branch and are noted in TODO.md's left-overs.

`alt.sh` (`alt-extend.log`): all ten build configurations reproduce the
references (the statistics build's test3 differs by its report, as
always), test2, testhigh and testhighmany reproduce theirs, valgrind is
silent on the debug build and two optimised runs.
