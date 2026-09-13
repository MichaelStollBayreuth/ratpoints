# Retiring the 32-bit architectures

Branch `retire-32bit` off `v2.3` at 2a305a2, 2026-09-13.  TODO item 19.
Michael, 2026-09-13: "I think we can retire the ability to work on 32-bit
architectures (whoever wants to use ratpoints on such a machine can use
v2.2.4)."  This file records what was removed, what was kept and why, and
how it was checked.  It stays on the branch; the merge strips it.

## What the code used to allow for

Only `rp-private.h` contemplated a word other than 64 bits, in two places:

- `LONG_LENGTH` was `8*sizeof(long)`, `LONG_SHIFT` a conditional on it
  (16, 32 or 64 bits), `LONG_MASK` derived from the shift.  Being a
  `sizeof`, `LONG_LENGTH` could not appear in a preprocessor expression,
  which is why the two generated headers carry the comments they do.
- `#if __WORDSIZE != 64` switched the two 128-bit variants off, since a
  bit array there is two `unsigned long`s.  (`__WORDSIZE` is a glibc
  macro; where it is not defined the test read `0 != 64` and switched them
  off as well, which was the right thing by accident.)

Everything else uses `LONG_LENGTH` as the name of a constant -- loop
bounds, the replication of 16-bit patterns to a word, the cap on the Sturm
depth, the mask macros -- and does not care what it is, as long as it is a
power of two.  The one other place that mentions a 32-bit target is the
cycle counter of `bench_check.c`, which accepted `__i386__`.

## What was done

`rp-private.h` now says

    #if ((ULONG_MAX >> 31) >> 1) == 0
    # error "ratpoints needs a 64-bit long since version 2.3; ..."
    #endif
    #define LONG_LENGTH 64
    #define LONG_SHIFT 6
    #define LONG_MASK (LONG_LENGTH - 1)

and the `__WORDSIZE` block is gone.  The check is on `ULONG_MAX` from
`<limits.h>` (standard C) rather than on `__WORDSIZE` (glibc) or
`__SIZEOF_LONG__` (gcc/clang), and it is two shifts by less than 32 rather
than a comparison with `18446744073709551615UL`, so that it is right in
every preprocessor: a C89 preprocessor that does its arithmetic in a
32-bit `unsigned long` would truncate that literal to `ULONG_MAX` and
pass.  Checked with `gcc -std=c89 -pedantic -E` on the four forms of
`ULONG_MAX` (a literal and glibc's `LONG_MAX * 2UL + 1UL`, for 32 and 64
bits): the two 32-bit forms hit the `#error`, the two 64-bit forms pass.
There is no 32-bit multilib on this machine, so a real `-m32` build was not
tried; the test above is the preprocessor's own arithmetic, which is what
the real build would run.

Where the check lives: `rp-private.h`, which every file of the library
includes (`sturm.c` includes only `ratpoints.h`, and has nothing that
depends on the word size).  It is not in the public header on purpose: a
program on a 32-bit machine could not link the library anyway, and the
public header should not pull in `<limits.h>` for a test that cannot fail
where the library builds.

Note that "64-bit long" is the requirement, not "64-bit machine": under
the Windows ABI `long` has 32 bits on x86-64 as well.  The manual says so.

The rest is comments and documentation:

- `gen_init_sieve_h.c`, `gen_find_points_h.c`: the comments that explained
  why `LONG_LENGTH` / `RBA_LENGTH` cannot be used in a preprocessor
  expression.  Both generators stay: the generated files still depend on
  the list of primes (and the second on the register width, which
  `RBA_PACK` pins down as before).
- `bench_check.c`: `__i386__` dropped from the `rdtsc` condition.
- `Makefile`: the comment on `CCFLAGS64`.
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

Type note: `LONG_LENGTH` used to have type `size_t` and `LONG_MASK` type
`unsigned long`; both are `int` now.  Every use compares or masks a
non-negative `long`, so the values are the same; only the signedness of a
few comparisons changes.  See the checks below for what that did to the
generated code.

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
- `-Wextra`: 42 warnings before, 11 after.  The 31 that went were
  `-Wsign-compare` on comparisons with `LONG_LENGTH`, which was unsigned;
  none is new.
- The `#error`: the preprocessor test described above.
- Review: see the end of this file.
