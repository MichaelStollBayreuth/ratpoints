# Divisions and reductions modulo p

Notes on TODO item 11, taken broadly: every place the program divides, and
whether it can be done better.  Work done on branch `divisions-and-reductions`.

The short answer is that the divisions themselves cost almost nothing, and
that the helper written to avoid them cost 2 to 3 per cent.

## The divisions cost nothing

The hardware integer divider is busy for this fraction of all cycles
(`perf stat -e cpu_core/arith.idiv_active/u`, against `cpu_core/cycles/u`):

| suite | cycles | divider busy | share |
|---|---|---|---|
| `test1`        | 6.29e9  | 2.93e7 | 0.47% |
| `test1many`    | 6.28e9  | 7.1e6  | 0.11% |
| `testhigh`     | 3.28e11 | 1.95e9 | 0.60% |
| `testhighmany` | 1.89e11 | 6.3e8  | 0.33% |

Floating-point division is under 0.001% everywhere.  So removing *every*
division in the program could not have won more than 0.6%.

Two whole files turn out to contain no division at all in the generated code:
`init.c`, because the primes are compile-time constants in the generated
`sieve_init_<p>` functions, and `sturm.c`, whose only divide instruction was a
floating-point one (removed, see below).  `gen_find_points_h.c` and
`gen_init_sieve_h.c` run at build time and do not count.

## The helper written to avoid them cost 2 to 3 per cent

`mod(a, b)` -- there is a copy in `find_points.c` and another in `sift.c` --
reduces `a` modulo `b` by subtracting multiples of `b`, and divides only when
`a` is outside `[-16b, 16b)`.  Its comment says "try to avoid divisions", and
it does: measured over a whole run of the busiest caller,

| suite | calls | of which divide | mean \|a/b\| when it does |
|---|---|---|---|
| `test1`    | 47,570,362    | 72,132 (0.15%)      | 20.5 |
| `testhigh` | 1,339,607,011 | 174,303,916 (13.0%) | 39.5 |

The 174 million divisions on `testhigh` account for essentially the whole
divider occupancy in the table above, so the accounting closes.

The busiest caller is the one that sets the `start` field of each sieve before
the first phase (`sift.c`, in `_ratpoints_sift0`).  It runs `sp2` times per
call, and a call handles at most `RATPOINTS_ARRAY_SIZE` bit arrays, so at a
large height bound there are of the order of a billion of them.

Note what the second table says about the cutoff: at the default height the
helper practically never divides, and at the large one it divides on an eighth
of the calls, because the word number grows with the height bound while the
primes do not.  So the cutoff of 16 is not the problem either way.

**The problem is the other 87%.**  Replacing the helper at that one call site
with a multiply-high reduction is worth, measured in one binary with the two
forms selected at run time so that code alignment cannot confound the
comparison (see the note in the Makefile):

| suite | conditional subtraction | multiply-high | change |
|---|---|---|---|
| `test1`        | 6663923188   | 6478525477   | -2.78% |
| `test1many`    | 7529462076   | 7310254287   | -2.91% |
| `testhigh`     | 336456271864 | 329627920589 | -2.03% |
| `testhighmany` | 200398331603 | 195979944981 | -2.20% |

## Why, given that the divider was only 0.6% busy

Not the divisions, and not branch mispredictions either.  On `test1`, where
the helper divides on 0.15% of calls:

| | cycles | instructions | branch misses | divider | IPC |
|---|---|---|---|---|---|
| conditional subtraction | 6.52e9  | 1.272e10 | 6.21e7 | 2.90e7 | 1.95 |
| multiply-high           | 6.311e9 | 1.168e10 | 6.12e7 | 2.83e7 | 1.85 |

Branch misses are unchanged: the chain's branches are well predicted, because
consecutive calls see slowly varying word numbers and the same primes.  What
changes is the instruction count -- 1.04e9 fewer over 47.6 million calls, or
about 22 instructions per call.  The chain is five conditional subtractions
with a serial dependency between them; the replacement is two multiplies and a
shift.  The saved instructions were cheap and parallel, which is why IPC falls
while the cycle count falls too.

This is worth remembering as a general point: an "avoid the expensive
instruction" trick that costs twenty cheap instructions is a bad trade at any
frequency where the expensive instruction is rare.

## The replacement

With `m = 2^64/p` rounded up, `u mod p` is the top half of `(m*u mod 2^64)*p`,
exact for every `u` below `2^32`.  That reciprocal already existed: the third
stage (TODO item 10) put it in the sieve entry, one division per prime and
curve.  `sieve_spec` now carries it too, copied per denominator.

The sign is put back after reducing `|a|` rather than removed beforehand by
adding a multiple of `p`, because the caller cannot bound the word number
tightly enough for a shifted value to stay inside 32 bits.  Above that bound
the old helper is used, which needs the guard to leave room for the offset
added to the word number -- hence the `RATPOINTS_MAX_PRIME` slack in it.

Checked against `%` for all 53 compiled-in primes over the whole range of
arguments, including the largest, in both signs.

`-DRP_MULMOD_DIVIDE` puts the divisions back, and `-DRP_MOD_CHOICE` builds
both forms into one binary, chosen by the environment variable `RP_MOD_MUL`,
which is how the table above was measured.

## Everything else, and why it stays as it is

A fan-out over the rest of the source catalogued 38 more sites.  None is worth
changing.  The ones worth recording:

* **The forbidden-divisor test** (`find_points.c`, in the denominator loop)
  looks like a division per denominator -- 44 million of them on `testhigh` --
  and performs none, because the list it walks is almost always empty.  See
  the TODO for what that says about the code that fills it.
* **`jacobi1`** is called once per denominator and contains a `mod(f, b)`
  where `f` is the leading coefficient.  It divides only when the coefficient
  exceeds sixteen denominators, which is rare for the curves in the suites.
  The rest of it is already a binary algorithm with no division.
* **`relprime`** contains no division: it is a binary gcd.
* **The `bp_list` update** steps `b` past each prime with a `while` loop that
  runs 0.17 times per step on average.  Nothing beats that.
* **The Horner loop in `examine_prime`** reduces modulo `p` once per residue
  per prime per curve -- about 1700 per curve, so far too rare to matter.  Its
  other branch, taken for degree 8 and above, carries a test in every Horner
  step and costs 1.31 times the branch-free one at equal degree.  That belongs
  to the degree item, not this one.
* **`valuation1`** in the square-denominator path does about six divisions per
  denominator it is called for, but that path sees under 1% of denominators.
  The call could be removed outright, since the valuation of `d*b^2` is known
  from the valuation of `b`; not worth it at that frequency.
* **`pointer_align`** divides twice per curve, and the divisor is always a
  power of two, so a mask would do. Twice per curve.

## One latent bug, fixed on the way

`sturm.c` scaled an interval endpoint by dividing by `1 << del`, where `del`
is the bisection depth.  That shift is on an `int`, and the depth is bounded
by `args->sturm`, which the program caps at `LONG_LENGTH - 2`, so any Sturm
iteration count above 31 was undefined behaviour.  It is now `ldexp`, which is
exact, defined at every depth, and removes the division as well.  No test
output changes, at any setting of `-S`.
