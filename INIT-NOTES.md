# Where the time goes in `sieve_init`, and what item 4 was worth

Branch `sieve-init-byte-rotations`, 2026-09-09.  TODO item 4 proposed replacing
the bit-rotation recurrence in `CODE_INIT_SIEVE2` by byte rotations out of eight
precomputed copies, on the grounds that the recurrence is a serial chain and the
byte form would make the rows independent.  It was tried, and on its own it is
worth nothing; with the same trick applied to it that the rest of this file
turned out to need, it is worth 2.5% of the set-up, which is not enough to pay
for what it costs.  But looking for it found two other things in the same
function that are worth a great deal: the set-up is now **3.7 times faster** on
the primes that matter and `make test1` **13% faster**.

## How it was measured

`bench_init` builds every table for every prime and denominator class, checks
each one word for word against a reference computed straight from the
definition, and then times the same sweep.  It gained two things here:

* a `hot` argument.  Without it each table is written after the last, so a full
  sweep streams 5.3 MB and the largest primes measure how fast this machine
  writes to memory: 891 ns a table against 697 hot.  A real run rebuilds the
  tables for one curve, some tens of kilobytes, and then reuses them, so hot is
  much the closer of the two.  All the figures below are hot.
* `pmin`/`pmax` now restrict the correctness sweep as well as the timing, so
  that `perf` on a single prime is not swamped by the check.

**A warning about this benchmark.**  It builds the same tables from the same
`is_f_square` data thousands of times over, and the branch predictor learns
them.  That flatters any version with a data-dependent branch in it by an
amount that has no counterpart in a real run, where every curve brings new
data.  It is why the small primes below appear to *lose* from the change and do
not.  Where the two disagree, `make test1` is the arbiter.

## Where the time actually goes

`perf annotate` on `sieve_init_127`, and then two probes -- `RP_INIT_NOACC`
leaves out the pattern accumulation, `RP_INIT_NOREP` the copies that replicate
the table over a bit array -- give the split.  Before any change, three quarters
of the function is not the rotation at all but the loop that turns
`is_f_square` into the *p*-bit pattern:

| part | before | after both changes |
|---|---|---|
| accumulating the pattern from `is_f_square` | ~75% | 51% |
| the rotation recurrence | | 29% |
| replicating over `RBA_PACK` and wrapping the chunk | | 20% |

The 75% is the `perf annotate` share of the accumulation loop in the original
`sieve_init_127`; the right-hand column is the probes, on the build as it now
stands.  Either way the rotation the item aimed at is the smaller half of what
is left, and was never the larger part of anything.

## What was wrong with the accumulation

Two things, one behind the other.

**A branch on random data.**  `if(isfs[ab]) { work |= test; }` is a branch whose
outcome is a quadratic residue test: it has no pattern.  At p = 127 it is a
third of all branches executed and misses about half the time.  Shifting the
value into place instead -- `work |= (unsigned long)isfs[ab] << i` -- has
nothing to mispredict.  `is_f_square` is 0 or 1 by construction, so no
normalisation is needed.

**A three-cycle chain.**  With the branch gone the loop runs at 3.2 cycles a
row, which is exactly `ab += d; cmp; cmov` -- one row cannot start until the
one before has finished.  Since `a` runs over consecutive integers, four rows
can go at once: keep the four residues `a*d`, `(a+1)*d`, `(a+2)*d`, `(a+3)*d`
and step each by `4d mod p`.  The four chains are independent and overlap.
Groups of four never straddle a word, `LONG_LENGTH` being a multiple of four,
so only the tail of the last word needs a plain loop.

Both are in `init.c` and both have an escape hatch for comparison,
`RP_INIT_BRANCH` and `RP_INIT_ONEWAY`.

## What it is worth

`bench_init`, hot, ns per table:

| primes | original | branchless | + four-way | |
|---|---|---|---|---|
| 3..61 | 48.5 | 55.6 | 43.0 | 1.13x |
| 67..127 | 365.1 | 133.9 | 99.6 | 3.67x |
| 67..251 | 733.2 | 268.6 | 214.3 | 3.42x |

(The 3..61 column is the predictor artifact described above: there is nothing
real in the middle entry.)

End to end, medians of seven interleaved runs, seconds:

| | original | branchless | + four-way |
|---|---|---|---|
| `make test1` | 1.622 | 1.431 | 1.408 |
| `make test1many` | 1.611 | 1.572 | 1.574 |

`test1` is 0.868 at gcc's default code layout and 0.859 at
`-falign-loops=32 -falign-functions=32`, so the win is not a layout accident.
`test1many` is 0.977: the set-up is only 4.6% of it to begin with.  The share of
`make test1` spent in `sieve_init_*` falls from 22.3% to 10.1%.

The four-way step is only 1.5% of `make test1` on top of the branchless change,
against 27% of the set-up in isolation, because by then the set-up is under a
tenth of the run.  It earns its place on workloads of many curves at small
height bounds, which is exactly where the set-up matters at all.

## Item 4 itself: byte rotations

Implemented as `RP_INIT_BYTE` and correct on all 6026 tables.  The extended
pattern is held in eight copies, copy *r* shifted along by *r* bits, so that the
window at bit offset *s* is eight bytes of copy *s mod 8* at byte offset
*s/8* -- one unaligned load a row, no shifting, no rotation state.

**One row at a time it is worth nothing**, and the reason is the same thing that
was wrong with the accumulation.  The rotation loop it replaces runs at 2.86
cycles a row.  The byte loop has no rotation chain, but it acquires one of its
own in the offset: `s += LONG_LENGTH; if(s >= p) s -= p;` is again add, compare,
conditional move.  One three-cycle chain has been swapped for another.

**Four rows at a time it does win**, using four offsets a step of
`4*LONG_LENGTH` apart, exactly as in the accumulation.  The byte form is what
makes that convenient: an offset that happens to be a multiple of
`LONG_LENGTH` is an ordinary case for it, whereas a double-word shift would
have to branch round a shift by the word length.  Hot, ns a table:

| primes | rotate | byte, one row | byte, four rows |
|---|---|---|---|
| 67..127 | 99.9 | 107.6 | 97.4 |
| 67..251 | 212.1 | 205.3 | 187.3 |

**It is still not worth adopting.**  End to end it is 0.998 and 0.989 at two
code layouts -- nothing, or nearly.  Two and a half per cent of a set-up that
is now a tenth of the run cannot be anything else.  Against that it wants eight
scratch buffers, unaligned access through `memcpy`, and it is a rotation only
on a little-endian machine, which the rest of the file does not assume.  It
would be worth revisiting only if the set-up became a large share again, or at
a much larger `PRIME_SIZE`, where the 67..251 column says it is worth 12%.

A third form was tried on the way, `RP_INIT_WINDOW`: one copy of the extended
pattern, each row a double-word shift at bit offset
`s(a) = a*LONG_LENGTH mod p`, no eight copies and no byte addressing.  It is
20% worse hot and level end to end -- the same offset chain, plus two extra
shifts a row, and the multiple-of-`LONG_LENGTH` rows have to be peeled off
first.

## What is left

The set-up is now about a tenth of `make test1` and half a per cent of a run at
a large height bound, so there is not much more to be had here.  For the record,
the split at p in 67..127 after the two changes is 51% accumulation, 29%
rotation, 20% replicating over the bit array and wrapping the chunk.

The replication reloads every word it copies, though the rotation loop had it in
a register a moment earlier; folding the two together would save p loads.  That
is perhaps 2% of `make test1`, and it overlaps with TODO item 5, which would
change the layout being copied into, so it belongs there rather than here.

The earlier finding that interleaving eight independent `(p, b)` constructions
gave -24% (`PERFORMANCE-NOTES.md`) was measured against the branchy
accumulation.  It was overlapping *that* chain, and it has been overtaken: the
four-way step gets the same effect inside one table and needs no restructuring
of the caller.

The general lesson, since it caught two of the three loops in this file: on this
machine a loop that carries `x += d; if(x >= p) x -= p;` from one iteration to
the next costs three cycles a row whatever else it does, and the fix is not to
make the step cheaper but to run several of them at once.
