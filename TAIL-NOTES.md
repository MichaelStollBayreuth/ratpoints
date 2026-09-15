# The tail leg: sieving a range's last bit arrays without padding them to a whole chunk

Branch `tail-leg` off `v2.3` at 9ca70c8, 2026-09-15.  The review's P2
(`review-2026-09-12/REVIEW-REPORT.md`; verdict F9 in `verdicts-final.md`,
both skeptics confirmed and corrected it), the first item of the sieve
group; TODO item 25.  Measured on the i7-1355U, pinned to P-core 0, in core
cycles, paired and interleaved with the base, medians of the per-round
ratios, outputs compared with the references before anything was timed
(`pair.sh` in `review-2026-09-12/`).

## What was wrong

`sift()` (find_points.c) cut the numerator interval of a denominator into
blocks of `array_size` bit arrays and rounded the last block up to a
multiple of `RATPOINTS_CHUNK` (16): `_ratpoints_sift0` then sieved the
padding with every prime of the first phase, the scan of the second phase
walked it, and a loop after the first phase zeroed it.  The verdicts
counted the padding at 14.2 per cent of all bit arrays swept in `make
test1` (height 16383; 11.2 in test1many, 13.3 in testdegrees), 80 per cent
at height 1000, 90 at height 200, 1.4 at 200000: a block at 16383 is 41
bit arrays long on average, so the last chunk is mostly padding, and at
height 1000 the whole interval of a denominator is 3 bit arrays in a chunk
of 16.

## What was done (655e10b)

The range is not padded any more.  In `_ratpoints_sift0` the chunk loop --
the hand-written ladder of sixteen registers, untouched -- runs over the
whole chunks, and the bit arrays left over, fewer than sixteen, are sieved
by `tail_leg()`, a helper inlined at the constant widths 8, 4, 2 and 1: one
leg for each set bit of the tail's length, each starting where the wider
ones stopped.  A leg is the first phase in miniature: the first prime
writes `w` registers from its table row ANDed with the 2-adic pattern, the
other `sp1-1` primes AND their rows in, the registers go to the survivors
array.  Two things it does not do.  It does not wrap the table pointers
around: the legs of a tail together advance by at most fifteen rows past
the pointers the chunk loop left, and every table carries
`RATPOINTS_CHUNK-1` wrap-around copies of its first rows at the end for
exactly this reason (init.c, gen_find_points_h.c), so a leg indexes the
pointers with its offset and reads within the table.  And it does not
store the pointers back: nothing reads `sieves[n].start` after the first
phase (the second phase has computed its own rows since the trio, 4bfc676),
and the next call recomputes them all on entry.  That makes a leg cheaper
than the ladder per prime by the pointer bookkeeping -- the load of `.end`,
the compare, the wrap-around loop and the store, about nine instructions
per prime in the compiled code -- which is what the verdict's fixed cost
`F` (105 instructions per leg at `sp1 = 11`) consisted of.

Why legs by the bits of the length and not the verdict's other variant, one
leg of the tail rounded up to a power of two: the verdict's arithmetic used
the ladder's `F` for a leg and concluded that a second leg only pays for
itself when it saves more than four bit arrays.  Without the bookkeeping a
leg's fixed cost is the loop over the primes, four instructions per prime,
roughly a third of that, and a leg pays when it saves more than one and a
half arrays; then the two variants come out within twenty instructions per
call of each other, with no padding at all and no zeroing loop to keep the
simpler one.  The companion change the verdict found, dropping the same
bookkeeping on the last *whole* chunk of a call, is now free whenever the
call has a tail (its last leg is a tail leg); for a call whose range is a
multiple of sixteen it would need a second copy of the ladder, and the
saving -- some nine instructions per prime and call, on the calls without
a tail -- is not worth one.

`n_pad` left `sift0`'s interface, the zeroing loop after the first phase
went with it, and `args->n_words` now counts the bit arrays that were
sieved: it counted the padded range, so the survival rate `adapt_primes`
works from was low by the padding share (14 per cent at 16383; the rule
fires from a million words on, so on `make test1` it did not fire at all,
and at 200000 the bias was 1.4 per cent).

What gcc made of it.  The helper's array of `w` registers must become `w`
vector registers, which is scalar replacement after the loops over it are
completely unrolled; at `-O2` gcc's early complete unrolling only unrolls
what does not grow the code, so the 4- and 8-iteration loops were still
loops when scalar replacement ran, the two wide legs kept their registers
on the stack (eight stores after the first prime, eight at the exit of the
prime loop, eight loads before the survivors are written; the loop body
itself was in registers), and `-fdump-tree-sra-details` said so
("Disqualifying reg - No scalar replacements to be created" for two of the
four instances).  `_Pragma("GCC unroll 8")` on the three loops over the
registers makes the early pass unroll them, and the compiled legs then
touch the stack only to reload the 2-adic pattern once per leg.
`_ratpoints_sift0` grows from 6097 to 8411 bytes, the loop over the primes
of each leg unrolled by four as `-funroll-loops` asks; with that loop kept
rolled (`RP_TAIL_NOUNROLL`, a switch for the measurement) it is 6852 bytes.
All the suites and `test3` byte-identical.

## Measurements

### What the sieve sweeps

`RP_PHASE_TIMING` builds of the base (9ca70c8) and of 655e10b, single runs,
the `arrays` counter (bit arrays swept by the first phase) and the number
of exact checks:

| suite | arrays, base | arrays, tail | ratio | exact checks, base -> tail |
|---|---|---|---|---|
| test1 (rptest) | 138092960 | 118456561 | 0.858 | 49643 -> 49659 |
| test1many (rptest-many) | 110958928 | 98509996 | 0.888 | 104969 -> 104901 |
| testdegrees | 5573216 | 4833977 | 0.867 | 2558 -> 2558 |
| rptest -h 1000 | 3281392 | 644657 | 0.197 | 3766 -> 3766 |
| rptest -h 200000 | 16083016336 | 15863072568 | 0.986 | 318996 -> 318914 |

The arrays fall by exactly the padding share the verdicts had counted (14.2,
11.2, 13.3, 80.4 and 1.37 per cent).  The number of exact checks moves by a
few dozen where it moves at all: the sieve's result on the real bit arrays
is the same to the bit (the legs read the very rows the padded chunk read),
so the only thing that can differ is what `adapt_primes` decides, and it
decides from `n_words`, which no longer counts the padding.  (Confirmed
below with the rule switched off.)

### Cycles, paired (655e10b against 9ca70c8)

Core cycles on P-core 0, interleaved, medians of the per-round ratios with
their range; instructions and branch misses as median ratios.  test200,
test1000 and test4000 are `rptest -h`; testmany1000 is `rptest-many -h
1000`.

| suite | default alignment (5 rounds) | -falign-loops=32 (3) | -falign-loops=64 (3) |
|---|---|---|---|
| test1 | **0.9104** [0.9068, 0.9206] (5); instr. 0.9307; misses 0.9553 | **0.9080** [0.9015, 0.9347] (3); instr. 0.9317; misses 0.9451 | **0.9057** [0.8826, 0.9231] (3); instr. 0.9242; misses 0.9488 |
| test1many | **0.9365** [0.8560, 0.9469] (5); instr. 0.9482; misses 1.0042 | **0.9340** [0.9125, 0.9394] (3); instr. 0.9485; misses 0.9961 | **0.9388** [0.9343, 0.9863] (3); instr. 0.9450; misses 0.9917 |
| testhigh | **0.9886** [0.9796, 1.0159] (5); instr. 0.9969; misses 0.9535 | **0.9805** [0.9767, 0.9840] (3); instr. 1.0012; misses 0.9540 | **0.9857** [0.9841, 0.9868] (3); instr. 0.9922; misses 0.9345 |
| testhighmany | **0.9871** [0.9828, 1.0059] (5); instr. 0.9974; misses 0.9788 | **1.0036** [1.0001, 1.0255] (3); instr. 0.9988; misses 0.9607 | **0.9731** [0.9582, 1.0193] (3); instr. 0.9952; misses 0.9797 |
| test200 | **0.9251** [0.9184, 0.9311] (5); instr. 0.9386; misses 0.9934 |  |  |
| test1000 | **0.8046** [0.7998, 0.8211] (5); instr. 0.8632; misses 0.9986 | **0.8118** [0.7960, 0.8477] (3); instr. 0.8632; misses 0.9974 | **0.8047** [0.7975, 0.8140] (3); instr. 0.8565; misses 0.9876 |
| test4000 | **0.8199** [0.8169, 0.8459] (5); instr. 0.8770; misses 0.9915 |  |  |
| testmany1000 | **0.7677** [0.7596, 0.7996] (5); instr. 0.8345; misses 1.0343 |  |  |

Nine per cent of `make test1`, six and a half of `make test1many`, one to
two of `make testhigh`, one or two of `make testhighmany` (the noisiest of
the four: 0.987, 1.004, 0.973), and a fifth of the run at the height bounds
1000 and 4000, a quarter of the point-rich set at 1000.  The verdicts had
said 3.5-4.5 per cent of test1 from instruction counts with a cycle factor
of two thirds; the instructions fell by 7 per cent, more than the 6.3 of
the prototype (no padding at all, no bookkeeping in the legs), and the
cycles fell by more than the instructions.  Both are consistent across the
three placements of the code.

### The loop over the primes of a leg: rolled or unrolled

`-funroll-loops` unrolls it by four (sift0 8411 bytes); with
`_Pragma("GCC unroll 1")` (`RP_TAIL_NOUNROLL`) it stays a loop (6852 bytes):

| suite | rolled against unrolled (3 rounds) |
|---|---|
| test1 | **0.9834** [0.9658, 1.0033] (3); instr. 1.0028; misses 1.0050 |
| test1many | **0.9998** [0.9995, 1.0055] (3); instr. 1.0055; misses 1.0022 |
| test1000 | **0.9995** [0.9977, 1.0079] (3); instr. 1.0041; misses 1.0075 |
| test4000 | **0.9987** [0.9981, 1.0168] (3); instr. 1.0046; misses 1.0176 |
| testmany1000 | **0.9996** [0.9940, 1.0062] (3); instr. 1.0131; misses 1.0109 |

Nothing in cycles, a quarter to one per cent more instructions; the verdict
had expected the rolled loop to give back a fifth of the gain.  The
decision is taken on the final form of the leg (next section).

### The form of the leg: indexing with an offset, or advancing the pointers

The leg as committed first (655e10b, "A") indexes each prime's pointer with
the offset at which the leg starts and leaves the pointers alone.  The
disassembly of the 8-wide leg showed gcc forming eight index registers for
`start + off + i` (so that the per-prime work is one load and eight ANDs
with base-plus-index addressing) and spilling five of them: about fourteen
instructions per prime instead of ten.  Two ways out were built and
measured, instruction counts exact (single pinned runs), cycles paired over
three rounds:

* "B", the prototype's form: each leg stores `start + w` back into the
  pointer (`sieves[n].start += w`; no wrap), so the next leg needs no
  offset and every access is base-plus-displacement -- one load, `w` ANDs,
  an add and a store per prime.  Result: 0.4 per cent *more* instructions
  than A on test1, 0.6 on test1many, 0.8 at height 1000, 1.4 on the
  point-rich set at 1000, and 1.00-1.01 in cycles.  The store costs the
  narrow legs, which are as frequent as the wide one, more than the spills
  cost the wide one.
* "A2", what was kept (f0e28a5): the 8-wide leg's offset is always zero and
  is written as the constant it is; the 4-, 2- and 1-wide legs get `t & 8`,
  `t & 12` and `t & 14`.  The spills go (scalar stack moves in `sift0`: 130
  in A, 99 here, 90 in the base), 0.9 per cent fewer instructions than A
  on test1 (7.310e9 -> 7.247e9), 0.7 on test1many, 1.0 at height 4000;
  three rounds of the one-second suites could not resolve that in cycles
  (test1 1.017 with a range of 0.96-1.03; test4000 0.991, testmany1000
  0.996), which is why the final chain below pairs A2 against A over seven
  rounds and against the base at the three placements.

`sift0`: 6097 bytes in the base, 8411 in A, 8246 in A2.

### The final form (f0e28a5) against A over seven rounds, and against the base

A2 against A, seven rounds (the suites at 16383 run in about a second, so
three rounds could not resolve a per cent): test1 **0.9885** [0.960,
1.031], test1many **1.0012** [0.978, 1.045], test4000 **0.9905** [0.990,
1.006]; instructions 0.9913, 0.9930, 0.9900.  A point on test1 and at
4000, nothing on the point-rich set, as the instruction counts said.

A2 (f0e28a5) against the base (9ca70c8), the same protocol as above:

| suite | default alignment (5 rounds) | -falign-loops=32 (3) | -falign-loops=64 (3) |
|---|---|---|---|
| test1 | **0.9047** [0.8844, 0.9099] (5); instr. 0.9226; misses 0.9558 | **0.8979** [0.8679, 0.9172] (3); instr. 0.9236; misses 0.9456 | **0.8862** [0.8819, 0.9300] (3); instr. 0.9163; misses 0.9517 |
| test1many | **0.9388** [0.8721, 0.9426] (5); instr. 0.9416; misses 1.0050 | **0.9335** [0.8850, 0.9379] (3); instr. 0.9419; misses 0.9960 | **0.9343** [0.9317, 0.9350] (3); instr. 0.9385; misses 0.9917 |
| testhigh | **0.9833** [0.9826, 0.9978] (5); instr. 0.9949; misses 0.9553 | **0.9842** [0.9686, 0.9867] (3); instr. 0.9992; misses 0.9556 | **0.9773** [0.9733, 0.9793] (3); instr. 0.9904; misses 0.9376 |
| testhighmany | **0.9981** [0.9830, 1.0062] (5); instr. 0.9957; misses 0.9807 | **0.9900** [0.9783, 1.0165] (3); instr. 0.9972; misses 0.9680 | **0.9749** [0.9541, 0.9771] (3); instr. 0.9936; misses 0.9751 |
| test200 | **0.9169** [0.9084, 0.9251] (5); instr. 0.9385; misses 0.9990 |  |  |
| test1000 | **0.8066** [0.7966, 0.8145] (5); instr. 0.8618; misses 0.9992 | **0.8020** [0.7992, 0.8231] (3); instr. 0.8618; misses 0.9909 | **0.8030** [0.7985, 0.8218] (3); instr. 0.8551; misses 0.9897 |
| test4000 | **0.8172** [0.7944, 0.8384] (5); instr. 0.8683; misses 0.9970 |  |  |
| testmany1000 | **0.7636** [0.7603, 0.7695] (5); instr. 0.8319; misses 1.0363 |  |  |

**Ten per cent of `make test1`** (9.5, 10.2 and 11.4 at the three
placements), **6.5 per cent of `make test1many`**, **2 per cent of `make
testhigh`** (1.6-2.3), between nothing and 2.5 per cent of `make
testhighmany` (0.998, 0.990, 0.975: its noise), **a fifth of the run at the
height bounds 1000 and 4000**, 8 per cent at 200, a quarter of the
point-rich set at 1000.  Instructions: 7.6-8.4 per cent of test1, 6 of
test1many, 0.1-1 of testhigh; the cycles fall by more than the
instructions on test1 and by less at 1000 (where the legs are narrow and
the loop over the primes is most of a leg's cost).  The verdicts had
predicted 3.5-4.5 per cent of test1 (from -6.3 per cent of instructions and
a cycle factor of two thirds), 3-5 of test1many, 0.2-0.7 of testhigh and
7-13 at 200-4000: right about the shape, low by a factor of two on the
random curves at 16383, because the prototype's legs still carried the
pointer bookkeeping and because the cycle factor turned out to be above
one rather than two thirds -- the padding's work was not the high-IPC kind
the verdict assumed.  These are the figures quoted in README.md, the manual
and ARCHITECTURE-TRIAGE.md.

### The checks

* Every suite and `test3` byte-identical after each commit; the build
  variants (plain 64-bit, AVX128, SSE, AVX-512 emulated, AVX without AVX2,
  `RATPOINTS_CHUNK=1` -- the generic arm, untouched --, `RATPOINTS_CHUNK`
  3, 8 and 12 -- other leg sets, 12 not a power of two --, the 64-bit
  second phase, `RP_MULMOD_DIVIDE`, phase timing, the debug build) all
  reproduce the references with no new warning; `make test1once`, test2,
  testhigh and testhighmany (the point-rich set at 200000, run separately
  because `alt.sh` had not built its binary) identical; valgrind clean on
  the debug build and on two optimised runs.
* The exact-check counts: with `RP_ADAPT_WORDS` raised beyond any run (the
  adaptation rule never fires) base and new tree report the same counts to
  the unit on test1 (49693 checks, 48430 `compute_bc` calls) and test1many
  (103433, 88172).  So the differences seen with the rule on (+16 and -68
  checks) are the rule seeing `n_words` without the padding, as it should.
* Branch misses in `sift0` (perf record on the base): the scan's exit, the
  second phase's data-dependent loop exits and the return account for
  nearly all of them; the per-prime wrap-around loop of the chunk loop is
  about four per cent.  The bookkeeping the legs drop was an instruction
  cost, not a misprediction cost.
