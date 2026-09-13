# Retiring the 32-bit architectures

Branch `retire-32bit` off `v2.3` at 2a305a2, 2026-09-13.  TODO item 19.
Michael, 2026-09-13: "I think we can retire the ability to work on 32-bit
architectures (whoever wants to use ratpoints on such a machine can use
v2.2.4)."  This file records what was removed, what was kept and why, and
how it was checked.  It stays on the branch; the merge strips it.

## What the code used to allow for

`rp-private.h` contemplated a word other than 64 bits in two places:

- `LONG_LENGTH` was `8*sizeof(long)`, `LONG_SHIFT` a conditional on it
  (16, 32 or 64 bits), `LONG_MASK` derived from the shift.  Being a
  `sizeof`, `LONG_LENGTH` could not appear in a preprocessor expression,
  which is why the two generated headers carry the comments they do.
- `#if __WORDSIZE != 64` switched the two 128-bit variants off, since a
  bit array there is two `unsigned long`s.  (`__WORDSIZE` is a glibc
  macro.  With a libc that does not define it, the test read `0 != 64` and
  switched the two variants off on a 64-bit machine as well, silently: a
  build asked for with `-DUSE_AVX128` came out as a 64-bit build.  So
  removing the block is not quite a no-op; on such a platform it gives the
  128-bit variants back.  The default build, `-DUSE_AVX`, was never
  affected.)

Everything else uses `LONG_LENGTH` as the name of a constant -- loop
bounds, the replication of 16-bit patterns to a word, the cap on the Sturm
depth, the mask macros -- and does not care what it is, as long as it is a
power of two.  The one other place that mentions a 32-bit target is the
cycle counter of `bench_check.c`, which accepted `__i386__`.

## What was done

`rp-private.h` now says

    #if (((ULONG_MAX >> 31) >> 31) >> 1) != 1
    # error "ratpoints needs a 64-bit long since version 2.3; ..."
    #endif
    typedef char rp_long_has_64_bits[(sizeof(unsigned long)*CHAR_BIT == 64) ? 1 : -1];
    #define LONG_LENGTH 64
    #define LONG_SHIFT 6
    #define LONG_MASK (LONG_LENGTH - 1)

and the `__WORDSIZE` block is gone.  The check is on `ULONG_MAX` from
`<limits.h>` (standard C) rather than on `__WORDSIZE` (glibc) or
`__SIZEOF_LONG__` (gcc/clang), and it is three shifts by at most 31 rather
than a comparison with `18446744073709551615UL`, so that it is right in
every preprocessor: a C89 preprocessor that does its arithmetic in a
32-bit `unsigned long` would truncate that literal to `ULONG_MAX` and
pass.  It tests for exactly 64 bits: a wider `long`, should one ever
appear, would get the wrong `LONG_LENGTH` and must not pass either.
Checked with `gcc -std=c89 -pedantic -E` on the four forms of `ULONG_MAX`
(a literal and glibc's `LONG_MAX * 2UL + 1UL`, for 32 and 64 bits): the
two 32-bit forms hit the `#error`, the two 64-bit forms pass.  (The
128-bit case cannot be tried this way: gcc's preprocessor computes in 64
bits and truncates such a constant, with a warning.  On a platform whose
`long` had 128 bits its preprocessor would compute in 128 bits, and the
value `ULONG_MAX >> 63` would be far above 1.)

The `typedef` is the same test in C proper, and it is there because the
preprocessor test rests on `<limits.h>` being conforming (C99 5.2.4.2.1:
the limits must be usable in `#if`).  The reviewer named two
non-conforming shapes: `((unsigned long)~0UL)` makes the `#if` itself a
syntax error (a confusing diagnostic, but still no build), and `(~0UL)`
passes on a 32-bit `long`, since every `#if` computes in its widest
unsigned type.  An array of negative size fails in every C compiler,
whatever the header says, and `sizeof` also settles the 128-bit case.
The preprocessor test stays for its message.
There is no 32-bit multilib on this machine, so a real `-m32` build was not
tried; the test above is the preprocessor's own arithmetic, which is what
the real build would run.

Where the check lives: `rp-private.h`, which every file of the library
includes except `sturm.c`, which includes only `ratpoints.h` and does
not reference `LONG_LENGTH` (it derives `LONG_MAX` as
`(unsigned long)(-1) >> 1`, which is right at any width; its shifts
`1L << (del-der)` depend on the word size only through the cap on
`args->sturm` in `find_points.c`, which is `LONG_LENGTH - 2`).  It is not in the public header on purpose: a
program on a 32-bit machine could not link the library anyway, and the
public header should not pull in `<limits.h>` for a test that cannot fail
where the library builds.

Note that "64-bit long" is the requirement, not "64-bit machine": under
the Windows ABI `long` has 32 bits on x86-64 as well.  The manual says so.

Seven casts that existed only because `LONG_LENGTH` was a `size_t` are
gone: the `(int)` in `WIDTH`, the `(long)` on the loop bounds and on the
cap of the Sturm depth in `find_points.c`, and `(long)(RBA_LENGTH-1)` in
the interval code (there it kept `high + RBA_LENGTH-1` signed for a
negative `high`; with an `int` constant it is signed by itself).

The rest is comments and documentation:

- `gen_init_sieve_h.c`, `gen_find_points_h.c`: the comments that explained
  why `LONG_LENGTH` / `RBA_LENGTH` cannot be used in a preprocessor
  expression.  Both generators stay: the generated files still depend on
  the list of primes (and the second on the register width, which
  `RBA_PACK` pins down as before).
- `bench_check.c`: `__i386__` dropped from the `rdtsc` condition.
- `Makefile`: the comment on `CCFLAGS64`, and the manual's copy of it in
  the section on building.
- Manual: the Requirements subsection states the 64-bit `long`; two
  descriptions of a bit array no longer list 32 bits; a paragraph in the
  2.3 change log.  The 2.1 and 2.1.1 entries, which mention 32-bit words
  and `__WORDSIZE`, are history and stay as they are.
- `README.md`: one paragraph.

## What was kept, and why

- `LONG_LENGTH`, `LONG_SHIFT` and `LONG_MASK` keep their names.  They are
  the vocabulary the code is written in ("a word"), and 64, 6 and 63
  scattered over the sources would say less.  What changes is that they are
  now literals, so they can be used in `#if` and in array bounds without a
  cast, and the review items that want a period-64 mask (P6) can rely on
  it.
- The loops that replicate a 16-, 8- or 4-bit pattern across a word
  (`for(i = 16; i < LONG_LENGTH; i <<= 1) { db |= db << i; }`) stay as
  loops: they run once per curve, and the loop says what it does.
- The `#ifndef __SSE2__` block stays: it is about the instruction set, not
  the word size, and it is what makes `-DUSE_AVX128` fall back to plain
  words on a machine without SSE2.
- The two comments that say a value has to "fit in 32 bits"
  (`RP_STAGE3_LIMIT`, `mod_mul`) are about the multiply-high reduction, not
  about the architecture.
- The `__SIZEOF_INT128__` fall-backs in `sift.c` and `find_points.c` stay:
  `__int128` is a compiler feature, not a word size, and the fall-back is
  also what `-DRP_MULMOD_DIVIDE` measures against.
- `(c & (0xff >> LONG_SHIFT)) == 0`, which breaks the lines of a DEBUG
  print in `find_points.c`, is the constant 3 now; its three siblings use
  `RBA_SHIFT`, which still varies with the register width, so it stays in
  the same form.

Type note: `LONG_LENGTH` used to have type `size_t` and `LONG_MASK` type
`unsigned long`; both are `int` now (`LONG_SHIFT` was an `int` already).
Every use compares, masks, multiplies or divides a quantity that is either
non-negative or, in the plain 64-bit build, a `long` that the old code
put through unsigned arithmetic and back (`mask_low = low - RBA_LENGTH *
w_low` with a negative `low`; `a0 = i * RBA_LENGTH` with a negative word
index): the values are the same on every two's-complement machine, and it
is the old code that needed the wrap-around.  One use was formally wrong
before: the DEBUG print `printf("... %ld ...", LONG_LENGTH*k + t, a)` in
`sift.c` passed an `unsigned long` to `%ld` (undefined by the letter of
the standard, harmless on every real machine, and diagnosed by gcc only
with `-Wformat-signedness`); it is a `long` now.  See the checks below for what
the type change did to the generated code.

## Checks

- Disassembly.  `objdump -d` of the four objects of the default (256-bit)
  build, before and after.  `init.o` and `sturm.o` are identical.  `sift.o`
  differs in six conditional jumps (`jbe` to `jle`: the comparisons of a
  non-negative `long` with `LONG_LENGTH` in the mask macros are signed now)
  and `find_points.o` in one shift (`shr` to `sar`, from
  `CEIL(height, LONG_LENGTH)`).  Same values on every input, since the
  quantities compared and shifted are non-negative.
- Tests.  `make test1 test1once test1many testdegrees test2 test3` at the
  default width; `test1 test1many testdegrees test3` at 64 bits
  (`CCFLAGS1=`), at 128 bits both ways (`-DUSE_AVX128`, `-DUSE_SSE`) and at
  512 bits emulated (`-DUSE_AVX512` without `-mavx512f`).  All pass.
- `-Wextra`, over the four library objects and the two generators: 42
  warnings before, 11 after.  The 31 that went were
  `-Wsign-compare` on comparisons with `LONG_LENGTH`, which was unsigned;
  none is new.
- Plain 64-bit build (`CCFLAGS1=`), the parent against the branch after
  the second commit: `sift.o`, `init.o` and `sturm.o` identical,
  `find_points.o` the same one `shr` to `sar`.  This is the build where
  the expressions that see a negative operand live (`RBA_LENGTH` is an
  `int` literal in the vector builds anyway).
- DEBUG builds at 64 bits, the parent against the branch, on
  `'1 0 126 0 441' 300` and `'-3 0 7 1 0 2 -5' 300`: output identical
  (36698 and 8579 lines).  The reviewer did six curves at two widths.
- After the second commit, the four objects of the default build are
  byte-identical to those of the first: the cast removals changed nothing.
- The `#error`: the preprocessor test described above.  The `typedef`:
  with `sizeof(unsigned int)` in place of `sizeof(unsigned long)`, a
  32-bit stand-in, gcc says "size of array is negative".
- Review: see the end of this file.

## Review

Two Opus agents, each in its own worktree: a sweep for anything still
contemplating another word size, and a skeptic set to refute the claim
that nothing behaves differently.

The sweep found the one real miss -- the manual's copy of the Makefile
sentence on `CCFLAGS64`, "word-size registers, which works on essentially
any machine" -- and the polish that went into the second commit: the seven
casts, the "32-bit machines" headlines in README and change log where the
requirement is a 64-bit `long`, a self-contradiction about `sturm.c` in
these notes, a dangling reference, and that the `#error` as first written
tested for at least 64 bits rather than exactly 64.

The skeptic upheld all four claims (no behavioural change; the
preprocessor test; nothing lost with the `__WORDSIZE` block; no other
dependence on the old types).  It enumerated every site where the
signedness changed by diffing `-Wsign-conversion -Wconversion` output
before and after, compared the old and new forms of the four expressions
that see negative operands in the plain build over 3.8 million inputs,
and diffed the DEBUG output of six curves at two widths (identical).  What
it added: the two non-conforming shapes of `ULONG_MAX` that the
preprocessor test cannot handle, answered by the `typedef`; the effect of
the old `__WORDSIZE` block on a libc without that macro, which these notes
had backwards; the DEBUG print that passed an `unsigned long` to `%ld`;
that `sturm.c` does depend on the word size through the cap on
`args->sturm`; and that the plain build had not been disassembled.  (Its
claim that the format mismatch is diagnosed by plain `-Wall` is not right:
gcc needs `-Wformat-signedness` for it.  Checked.)
