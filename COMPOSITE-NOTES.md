# Composite sieving moduli

Branch `composite` off `v2.3` at cb0f187, 2026-09-17.  TODO item 21 (the
review's P9: verdict F3 on the products of primes, the reviewer's own
finding (b) on the prime powers).  Measured on the i7-1355U, pinned to
P-core 0, in core cycles, paired and interleaved with the base, medians of
the per-round ratios, outputs compared with the references before anything
was timed (`pair.sh`).

## What was wrong

A sieve table row only has to be periodic in the bit index with the
modulus as its period and be selected by the denominator's residue;
nothing in the first two phases uses primality.  The program sieved with
primes only, and small primes say little for what they cost: modulo 3 a
random curve admits two residues in three, so an AND with the row for 3
removes a third of the candidates, where an AND with a large prime's row
removes half.  Two things were on the table.

**Prime powers.**  Modulo 9 the admissible residues are those with f(x) a
square mod 9, and squares are 4 of 9 residues against 2 of 3: over the
thousand random test curves the mean density is 0.451 mod 9 against 0.666
mod 3, that is 1.19 large primes' worth of information for one AND and a
nine-row table, where 3 carries 0.61.  Mod 27: 1.37 (27 rows); 25: 1.41;
125: 1.51; 49: 1.28; 121 and 169: 1.2.  Nothing about the sieve changes
for them: for a unit denominator the row is f(a b^-1) square mod p^e, for
p | b it is frev(b a^-1) with frev(t) = t^D f(1/t) over the units a --
exactly what item 26 does for 2^6.

**Products.**  The row of a product m of coprime moduli is the AND of the
factors' rows (by the Chinese remainder theorem the residues are
independent, so its density is the product of the factors'): one AND per
word carries them all.  105 = 3*5*7 carries 2.26 primes' worth (the
verdict's wheel), 45 = 9*5 1.98 for 45 rows, 63 = 9*7 2.06, 225 = 9*25
2.6, and pairs of mid-size primes like 221 = 13*17 two primes' worth for
221 rows.  On the point-rich curves every small modulus is silent (the
mean density of 105 is 0.89, of 15 0.96), as the verdict said.

## What was done

**Candidates.**  After the primes have been looked at and the run's
length is known, `add_moduli()` (find_points.c) walks every odd m from 9
to `RATPOINTS_COMPOSITE_MAX`, factors it over the primes looked at, and
offers it as a candidate when it is not a prime, every factor is
informative, every prime involved is among the first 64, and its table
alone would not cost more than `RP_MODULUS_TABLE_MAX` (2) first-phase ANDs
per word -- which at height 100 rules out every modulus before any work is
done, so a short run pays nothing for this.  A factor p keeps the prime's
density; a factor p^e is examined the first time a modulus needs it
(`examine_power()`: f mod p^e at every residue for fsq, frev at every
multiple of p for gsq, the inverses of the units, and the density r as the
mean over the classes of b mod p^e, the classes with p | b counted only
when denominators divisible by p occur, as for a prime).  The candidate's
density is the product of its factors' -- exact, by the Chinese remainder
theorem -- and its key is what `prime_key` gives a prime of that size, so
the table cost m*min(D, m) is charged as for a prime.

**Selection.**  The candidates are ranked with the primes by key.
`take_entries()` walks the ranked pool for each stage: an entry sharing a
prime with one already taken is dropped from the pool (a mask over the
prime indices, `entry.mask`), the ones taken stay in front, and the stage's
own stopping rule applies -- the survivors-per-word target for the first
phase, the offset for the second, or the count the caller fixed.  What is
left after the second phase loses its composite entries, so the third
stage and `adapt_primes()` see primes only; `adapt_primes()` never demotes
a composite modulus out of the second phase either.  The `entry` type
carries the modulus, the mask and the factor codes; a composite's sieve
entry is made only when it is taken (`make_modulus()`), a prime power's
once (`make_power()`), so the ranking costs no entries.

**Rows.**  A prime power's row for the residue b is its pattern -- fsq at
x b^-1 for a unit b, gsq at b x^-1 over the units x for p | b -- laid out
over m bit arrays by `lay_out_pattern()` (init.c: the pattern repeated,
then 64 bits read at (x*RBA_LENGTH + 64 l) mod m for word l of bit array x,
plus the wrap-around copies).  A product's row is the AND of its factors'
rows, laid down in blocks of the factor's length (the factors divide m), a
factor's row being built first when it is not there yet
(`_ratpoints_sieve_init_product`); with 256-bit arrays that is one vector
AND per array and factor.  Neither touches sift.c: a composite entry has
p = m, a reciprocal, a bias, RBA_LENGTH^-1 mod m for the row shifts of the
numerator classes (`rbainv`, which replaces the look-up in the prime's
inverse table in `class_offsets`), and 2^-k mod m for the packing.

**Buffers.**  `se_buffer` holds an entry for every prime looked at, every
composite taken and every prime power that is a factor; `sieve_list` and
`magics` twice RATPOINTS_NUM_PRIMES; a new `pw_buffer` what the powers
record; the table buffer is enlarged after the selection for the tables of
the composites taken (`ensure_ba_buffer`, with `ba_buffer_arrays` recording
its size), before any table is built.  The verbose report says "moduli"
where composites may appear.

## Correctness

Every row a composite modulus builds can be checked against the
definition with `-DRP_VERIFY_MODULI` (init.c): bit N of the row for the
residue b is set exactly when F(N, b) is a square mod m and no prime of m
divides both N and b.  test1, test1many and testdegrees and the stride-8
curve at height 200000 ran under it with no mismatch.

## What it is worth

**Uncapped (171cb63, every odd composite below 256 offered).**  Cycles
new/base, 3 rounds: test1 0.955, test1many 1.019, testhigh 0.903,
testhighmany 1.068, test200 1.010, test1000 0.998, test4000 0.989,
point-rich at 1000 0.978; instructions test1 0.915, testhigh 0.777,
test1many 0.984, testhighmany 0.988.  The ranking takes composites on
every random curve: at 16383 the powers and small products (25 on 261 of
879 curves, 49 on 153, 9 on 134, 45 on 130, 27 on 125, 35, 77, 55, 63, 33,
39, 21 ...; mean sp1 11.7 -> 9.8, phase-1 rows 378 arrays in all), at
200000 the products of two mid-size primes (247 on 234 curves, 221 on 222,
253 on 215, 209, 203, 187, 225, 217 ...; sp1 7.6, phase-1 rows 895
arrays).  On the point-rich curves composites enter the second phase on 73
of 98 and 29 of 30 curves (171, 153, 187, 247, 125 ...), and there they
lose.

**Why the cycles do not follow the instructions.**  On the stride-8 curve
at 200000 the first phase's ANDs fall 45% (43.3M -> 23.6M) and its cycles
not at all (69.9M -> 69.2M): 1.6 cycles per AND with the base's small
primes (13, 31, 5, 17, 67, 79, 37, 89, 127, 19, 71; 555 arrays = 18 KB of
rows), 2.9 with 225, 221, 217, 209, 67, 79 (1018 arrays = 33 KB of rows,
and tables of 1.7 MB each against 180 KB for a prime near 67).  A small
prime's row lives in the first-level cache and is cycled through many
times per denominator; a row of 220 arrays does not, and a table of 1.7 MB
does not fit the second-level cache beside the others, so at a large
height bound, where most denominators sweep fewer arrays than the modulus
has, every line of the row comes from further away.  The cost model
charges every AND alike (the verdict's skeptic 1 warned that it has no
footprint term, F27), so it trades cheap ANDs for dear ones.  The
survivors of the first phase also rise on that curve (57837 -> 102649):
both choices meet the target of 0.0075 per word, the base's eleventh
prime overshot it, which is the rule's discrete step and not a fault.

**The cap.**  `RATPOINTS_COMPOSITE_MAX` bounds the moduli offered.
Cycles new/base, 3 rounds each, against the same base:

| suite | no cap (255) | cap 128 | cap 64 | cap 32 |
|---|---|---|---|---|
| test1 | 0.955 | 0.960 | 0.941 | 0.947 |
| test1many | 1.019 | 1.006 | 1.008 | 0.996 |
| testhigh | 0.903 | 0.895 | 0.890 | 0.938 |
| testhighmany | 1.068 | 1.005 | 0.973 | 0.997 |
| test1000 | 0.998 | 0.982 | 1.012 | 0.999 |

Instructions at cap 64: test1 0.927, test1many 0.994, testhigh 0.896,
testhighmany 0.998 -- and the cycles follow them, which they do not for
the larger moduli (cap 128: instructions 0.917 and 0.842 on test1 and
testhigh, cycles 0.960 and 0.895).  The moduli below 64 are the prime
powers 9, 25, 27, 49 and the products 15, 21, 33, 35, 39, 45, 51, 55, 57,
63: rows of at most 64 bit arrays, 2 KB at 256 bits, and tables of at
most 64 rows, which is where the model's assumption that an AND costs the
same whatever the modulus holds.

**Counters at cap 64** (new/base; prelim-c64.txt in the reports):

| suite | phase-1 ANDs | survivors of phase 1 | phase-2 ANDs | exact checks | table rows | mean sp1 base -> new |
|---|---|---|---|---|---|---|
| test1 | 0.864 | 0.968 | 0.906 | 1.100 | 0.854 | 11.68 -> 9.97 |
| test1many | 0.999 | 0.985 | 0.972 | 1.100 | 0.999 | 17.03 -> 16.97 |
| testdegrees | 0.837 | 1.006 | 0.964 | 1.127 | 0.848 | 13.12 -> 11.76 |
| test1 at 1000 | 0.967 | 0.972 | | 1.101 | 0.934 | 13.61 -> 13.34 |
| testhigh | 0.866 | 1.030 | 0.973 | 0.893 | 0.913 | 10.65 -> 9.01 |
| testhighmany | 1.000 | 1.000 | 0.992 | 1.226 | 1.003 | 19.13 -> 19.13 |

What is taken (cap 64): on test1 a composite in the first phase on 874 of
879 curves -- 25 on 286, 49 on 202, 45 on 137, 9 on 134, 27 on 128, 35 on
127, 55 on 121, 63 on 105, 33, 39, 21, 51, 57, 15 -- and in the second on
49; at 200000 55 on 294, 63 on 293, 49 on 202, 25 on 178, 45 on 171; at
1000 only 9 (on 286 curves), the table gate leaving nothing else; on the
point-rich suites the first phase takes a composite on 18 of 98 and 0 of
30 curves, the second phase 9, 25, 27, 49 on 42 of 98 and 26 of 30.  The
exact checks rise where the composites take primes out of the pool: the
third stage has fewer informative primes left and takes two fewer on test1
(sp3 16.2 -> 14.3); the cycles absorb it.  On the point-rich curves a
composite in the second phase lowers the model's S more than the actual
survivors fall (they sit on the floor of non-reduced points), so the stage
takes fewer primes there too (checks +23% on testhighmany, 2.7% faster
regardless).

## Not done, and why

* A footprint term in the cost model (a per-AND cost rising with the
  modulus, or a budget on the first phase's rows), which would let the
  large products in where they pay -- at 200000 they save 45% of the ANDs
  on some curves -- and keep them out where they do not.  It wants a
  cache-size constant and it concerns the primes too (the point-rich
  curves run the first phase on 20 primes' rows, 67 KB): TODO item 30, for
  the tuning session.
* Prime powers above 64 (81, 125, 121, 169, 243) and products above it
  (77, 91, 99, 105, 117 ...): excluded by the cap; 105 was the verdict's
  wheel.  Cap 128, which admits them, measured worse than 64 on every
  suite.
* Composite moduli in the third stage: it indexes is_f_square by
  (a b^-1) mod p, and a composite has no such table; the stage keeps to
  primes.
* Mod 3 information in the 2-adic style for every curve (the analogue of
  item 26 for 3-adic): that is what the prime power 9 or 27 in the
  ranking does, when it pays.
