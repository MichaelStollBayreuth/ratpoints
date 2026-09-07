# Things to try

Ideas that have not been investigated yet, kept next to `PERFORMANCE-NOTES.md`,
which records what *has* been measured (including the dead ends).  Roughly in
order of expected payoff.

## 1. Multi-threading

By far the largest potential gain, and the one thing that would change the
scale of what the program can do.  Drew Sutherland's data point: the curve

    247747600 -985905640 567207969 2396040466 52485681 -470135160 82342800

to height 10^6 takes 92 s single-threaded on a Zen 5 (16 cores), and about
**5 s using all 16 cores** -- an 18-fold speedup, against the 13-14% that the
512-bit registers buy.  Everything else in this file is second-order compared
with this.

The structure of the program is already close to what is needed.  The natural
unit of work is one denominator `b`: the loops at `find_points.c:1889` and
`find_points.c:1969` (and the `use_squares` variants above them) run over `b`
independently, and `sift(b, ...)` for different `b` share nothing that is
written.  What has to be checked:

* All per-run scratch space already lives in `ratpoints_args` and is allocated
  in `ratpoints_init` (`work`, `se_buffer`, `ba_buffer_na`, `int_buffer`,
  `sieve_list`, `den_info`, `divisors`, `forb_ba`, `forbidden`).  Giving each
  thread its own `ratpoints_args` should therefore be enough; the generated
  tables (`sieves0`, `squares`, `inverses`) are read-only and can be shared.
  Memory cost: one `ba_buffer` per thread, a few MB each.
* `se->sieve[b]` is filled in lazily by `sieve_init_<p>` and cached in the
  `args`-owned buffer.  Per-thread `args` keeps that private, at the price of
  recomputing tables that a shared cache would build once.  Since a thread
  handles a contiguous range of `b`, and the table depends on `b mod p`,
  the caching is mostly wasted across threads anyway -- worth measuring
  before designing anything clever.
* The `process` callback is user code and is called from the sieve; it would
  become the user's synchronisation problem.  The library interface has to say
  so, or serialise the calls itself.  Note that exact checking is nearly free
  (see `PERFORMANCE-NOTES.md`), so serialising the callback costs little.
* `*quit` is a shared early-exit flag; it needs to become atomic (relaxed is
  fine -- it is only a hint to stop).
* Load balance: the work per `b` is not uniform (`which_bits`, the Sturm
  bounds, and the `use_squares` shortcut all vary with `b`), so a dynamic
  schedule over `b` is better than a static split.

Open design question: whether to thread inside the library (pthreads or
OpenMP, with a `num_threads` field in `ratpoints_args`) or to leave it to the
caller and only document that separate `args` are independent.  The second is
much less invasive and is probably the right first step -- the command-line
program can then fork threads over `b` ranges itself.

## 2. Phases 1 and 2 separately, for each register width

We have never separated the two phases when comparing 64/128/256/512-bit
builds; only the total.  The hunch is that **phase 2 may well be fastest at
64 bits**, so that the ideal build is not uniform in the register width.

The reasoning, which is what the measurement should test:

* Phase 1 is a stream of `sp1` loads and `AND`s per bit-array of survivors.
  Its cost per numerator is proportional to `sp1/W`, so it scales cleanly with
  the width `W`.
* Phase 2 costs (a) a scan over `N/W` bit-arrays, which also improves with `W`,
  plus (b) the `sp2-sp1` loop and the bit extraction for each *surviving*
  bit-array.  With a per-bit survival probability `q` of order `2^-sp1`
  (measured: 4.17% of 256-bit arrays survive at `sp1 = 11`, i.e.
  `q ~ 1.6e-4`), the number of surviving arrays per numerator is
  `(1-(1-q)^W)/W ~ q`, which is *independent of the width*.  So the count of
  phase-2 `AND`s does not fall as `W` grows, while each one is wider.
* Worse, the loop `for(n = sp2-sp1; n && TEST(nums); n--)` exits as soon as
  `nums` is all zero.  A wider `nums` holds more bits and therefore takes
  *more* primes to clear, so the early exit fires later at 512 bits than at 64.
  This is the effect that could make phase 2 genuinely slower on wide
  registers.

If that is confirmed, two consequences follow.  First, the optimal `sp1`
should *rise* with the register width, because the width only helps the phase
that scales.  Second, a mixed build -- phase 1 wide, phase 2 on 64-bit words
reinterpreting the same memory -- would be possible, since the bit-arrays are
just packed words (the `EXT0`/`EXT` extraction already reads them word by
word).  Measure first.

Instrumentation: the phase split was previously obtained by timing the two
loops in `_ratpoints_sift0` separately; do the same under `perf stat -e
cpu_core/cycles/` for each of `CCFLAGS64`/`128`/`256`/`512`, at several
heights.

## 3. Tie the grouping (and `sp1`) to the survivor rate, not to `sp1`

`PERFORMANCE-NOTES.md` tabulates the phase-2 group size against `sp1` and finds
a ridge crossing at about `sp1 = 13`.  But `sp1` is only a proxy: what the
group-skip actually responds to is **the fraction of bit-arrays surviving phase
1**, which depends on the curve (how many primes really cut down the residues)
as well as on `sp1`, and on the width `W`.

So:

* Instrument `_ratpoints_sift0` to report the phase-1 survival rate per curve,
  and re-plot the group-size table against the *measured* rate instead of
  against `sp1`.  If the curves collapse onto one, we have the right variable.
* Run the cycle-counting sweep over a *range of curves* (different genus,
  different numbers of points, different local behaviour), not just the single
  large-height curve used so far, and over `sp1` for each.
* If the dependence on the survivor rate is clean, then `sp1` itself should be
  chosen from the *expected* survivor rate rather than from the height bound
  alone.  That rate is cheap to predict: `sieving_info` already computes, for
  each prime, how many residues survive, so the product over the first `sp1`
  primes gives the expected rate before any sieving is done.  Choosing `sp1`
  (and the group size) from that product would adapt automatically to curves
  where the small primes happen to be unusually effective or unusually weak.

This is the concrete way to attack the three-dimensional `(sp1, sp2, group)`
problem: reduce it to a function of one measurable quantity.

## 4. Byte rotations instead of bit rotations in the set-up phase

In `CODE_INIT_SIEVE2` (`init.c`, primes above `LONG_LENGTH`) the table is built
by a serial recurrence that rotates the `p`-bit pattern by `LONG_LENGTH mod p`
bits for each successive word:

    help[a1] |= help[t] << diff_shift;
    si[a] = help[a1];
    a1 = t;
    help[a1] >>= diff;

`PERFORMANCE-NOTES.md` records that this loop is latency-bound on exactly this
recurrence -- every attempt to make it cheaper failed, and only *interleaving
independent `(p, b)` pairs* helped, because that breaks the dependence chain.
The idea here breaks the same chain a different way, without restructuring the
caller.

Any rotation by `r` bits is a rotation by `r mod 8` bits followed by a rotation
by whole bytes.  So precompute the pattern rotated by 0, 1, ..., 7 bits, each
stored *doubled* (the `p`-bit pattern concatenated with itself, `2p` bits).
Then every row of the table is an **unaligned copy out of one of those eight
buffers at a byte offset** -- no shifting, no recurrence, and all rows
independent of each other.  The eight buffers cost `8 * 2p` bits, i.e. 2 KB at
the largest prime (1021), so they stay in L1.

Why it should help: it converts a latency-bound serial chain into `p`
independent unaligned loads and stores, which is throughput-bound, and x86
unaligned loads are nearly free.  It should compose with the interleaving idea
rather than compete with it.

Things to watch: the byte offset advances by `(LONG_LENGTH mod p) / 8` per row
only when `LONG_LENGTH mod p` is a multiple of 8, which it is not in general --
so the *bit* offset mod 8 cycles through the eight buffers in a fixed pattern
of period 8 (or a divisor), and the byte offset advances irregularly.  That is
still just index arithmetic with no data dependence.  Also check the endianness
assumption: reading bytes out of a doubled `unsigned long` buffer at a byte
offset is only a rotation on a little-endian machine.

## 5. The memory layout of the sieves

Never looked at.  Some specific questions:

* Each table for `p > LONG_LENGTH` stores `RBA_PACK * p` words -- the pattern
  repeated `RBA_PACK` times, so that reading `RBA_PACK` consecutive words at
  any offset gives the right rows.  That means the footprint of the sieve
  buffer **grows linearly with the register width**: the same tables that take
  `x` MB at 64 bits take `8x` MB at 512.  This may well be part of why wide
  registers pay off less than the width suggests, and it is worth measuring
  the L2/L3 miss rate against the width.  Is there a layout that avoids the
  replication, e.g. by over-allocating one pattern by `RBA_PACK-1` words and
  reading unaligned?
* Phase 1 reads `sp1` tables per bit-array -- `sp1` concurrent streams with
  different strides, hitting `sp1` different pages.  Check the dTLB miss rate;
  transparent huge pages for `ba_buffer` might be free money.
* The tables are laid out in `ba_buffer` in the order they happen to be built
  (lazily, per `(p, b)`).  The phase-1 primes are the hot ones and are read
  together; allocating them contiguously, or separating the phase-1 and
  phase-2 tables into different regions, might improve locality at no cost.
* The small primes (`p < LONG_LENGTH`, in `sieves0`) are tiny and read on every
  pass; they should be permanently L1-resident.  Check that they are, and that
  the alignment padding does not push them out.
