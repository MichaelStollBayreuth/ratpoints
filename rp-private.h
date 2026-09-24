/***********************************************************************
 * ratpoints-3.1.0                                                     *
 *  - A program to find rational points on hyperelliptic curves        *
 * Copyright (C) 2008, 2009, 2022, 2026  Michael Stoll                 *
 *                                                                     *
 * This program is free software: you can redistribute it and/or       *
 * modify it under the terms of the GNU General Public License         *
 * as published by the Free Software Foundation, either version 2 of   *
 * the License, or (at your option) any later version.                 *
 *                                                                     *
 * This program is distributed in the hope that it will be useful,     *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of      *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       *
 * GNU General Public License for more details.                        *
 *                                                                     *
 * You should have received a copy of version 2 of the GNU General     *
 * Public License along with this program.                             *
 * If not, see <http://www.gnu.org/licenses/>.                         *
 ***********************************************************************/

/***********************************************************************
 * rp-private.h                                                        *
 *                                                                     *
 * Header file with information local to the ratpoints code            *
 *                                                                     *
 * Michael Stoll, Apr 14, 2009; Jan 7-18, 2022; Sep 6, 2026            *
 * with changes by Bill Allombert, Dec 29, 2021                        *
 ***********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <string.h>

/* The code assumes that an unsigned long has 64 bits (version 2.2.4 is the
 * last one that accommodates a 32-bit long).  The test is for exactly 64
 * bits, and it is written so that it is right in every mode of preprocessor
 * arithmetic. */
#if (((ULONG_MAX >> 31) >> 31) >> 1) != 1
# error "ratpoints needs a 64-bit long since version 3.0.0; use version 2.2.4 on this machine"
#endif
/* The same test in C proper, for a <limits.h> that defines ULONG_MAX in a
 * form no #if can evaluate (a cast, or ~0UL, which every #if computes in its
 * widest type): an array of negative size does not compile. */
typedef char rp_long_has_64_bits[(sizeof(unsigned long)*CHAR_BIT == 64) ? 1 : -1];

#define LONG_LENGTH 64  /* number of bits in an unsigned long */
#define LONG_SHIFT 6    /* 2^LONG_SHIFT == LONG_LENGTH */
#define LONG_MASK (LONG_LENGTH - 1)

/* Check if SSE instructions can be used.  Both 128-bit variants need them:
 * USE_AVX128 uses only SSE2 intrinsics, in spite of its name. */
#ifndef __SSE2__
#undef USE_SSE
#undef USE_AVX128
#endif

#include "ratpoints.h"

/* Threads.  The library sieves on several threads when asked to (the
 * num_threads field): POSIX threads, and the gcc/clang atomic builtins for
 * the few words that one thread writes and another reads without a lock --
 * a pointer to a sieve table that a thread has just built, the count of the
 * correction's decisions, the flag that stops a search.  RP_PUBLISH writes
 * such a word after everything it points to is complete, RP_ACQUIRE reads
 * it so that what it points to is seen complete, and RP_RELAXED reads a
 * flag.  (Not a fence after a loop of relaxed loads: ThreadSanitizer does
 * not see through fences, and the acquire loads cost nothing measurable.)
 * With -DRATPOINTS_NO_THREADS they are the plain operations, and the
 * library sieves on the calling thread whatever num_threads says. */
#ifndef RATPOINTS_NO_THREADS
# define RP_ACQUIRE(p)    __atomic_load_n((p), __ATOMIC_ACQUIRE)
# define RP_RELAXED(p)    __atomic_load_n((p), __ATOMIC_RELAXED)
# define RP_PUBLISH(p, v) __atomic_store_n((p), (v), __ATOMIC_RELEASE)
#else
# define RP_ACQUIRE(p)    (*(p))
# define RP_RELAXED(p)    (*(p))
# define RP_PUBLISH(p, v) (*(p) = (v))
#endif

#define FLOOR(a,b) (((a) < 0) ? -(1 + (-(a)-1) / (b)) : (a) / (b))
#define CEIL(a,b) (((a) <= 0) ? -(-(a) / (b)) : 1 + ((a)-1) / (b))

/* Define interface for ratpoints_bit_array datatype:
 * RBA_LENGTH : number of bits
 * RBA_SHIFT  : 2^RBA_SHIFT == RBA_LENGTH
 * RBA_PACK   : number of words in a ratpoints_bit_array == RBA_LENGTH/LONG_LENGTH
 * RBA(a)     : fill a ratpoints_bit_array with copies of the word a
 * zero       : all bits zero == RBA(0UL)
 * AND(a,b)   : bit-wise and operation: a &= b
 *              (used in stage 2 to test several bit-arrays at once)
 * EXT0(a)    : extract first word (as unsigned long)
 * EXT(a,i)   : extract word with index i (as unsigned long)
 * TEST(a)    : tests if a is zero: TEST(a) == 0 <==> a == zero
 *              TEST should be fast if possible; it is used frequently
 *              in stage 2 of the sieve.
 * MASKL(a,s) : set lower s bits of a to zero
 * MASKU(a,s) : set upper s bits of a to zero
 *              MASKL and MASKU don't have to be terribly efficient;
 *              they are each executed once per denominator and interval.
 *              Both may assume that 0 <= s < RBA_LENGTH; the two call sites,
 *              in _ratpoints_sift0, pass mask_low and mask_high, which sift()
 *              computes as residues mod RBA_LENGTH.  (This matters: for
 *              s == RBA_LENGTH the versions below would either shift an
 *              unsigned long by LONG_LENGTH or address a word outside a.)
 */

#ifdef USE_AVX512
/* Use 512 bit AVX registers for the bit arrays */
/* The tests run this variant on a machine with AVX2 only, by compiling
 * with -DUSE_AVX512 but without -mavx512f: gcc then lowers the 64-byte
 * vectors to pairs of 256-bit operations, so that the whole code path
 * (RBA_PACK == 8, the mask macros, the second stage) is exercised.  What
 * such a run cannot check is the genuine 512-bit instructions, i.e., the
 * TEST macro below; with those, the test suite passes on AVX-512 hardware
 * (Zen 5). */

#include <immintrin.h>

#define RBA_LENGTH (512)
#define RBA_SHIFT (9)
#define RBA_PACK (8)
typedef unsigned long ratpoints_bit_array __attribute__ ((vector_size (64)));
#define RBA(a) ((ratpoints_bit_array){((unsigned long) a), ((unsigned long) a), \
                                      ((unsigned long) a), ((unsigned long) a), \
                                      ((unsigned long) a), ((unsigned long) a), \
                                      ((unsigned long) a), ((unsigned long) a)})
#define zero (RBA(0LL))
#define AND(a,b) ((a) = (a)&(b))
#define EXT0(a) ((unsigned long)(a)[0])
#define EXT(a,i) ((unsigned long)(a)[i])
#ifdef __AVX512F__
/* vptestmq sets bit i of the mask register if and only if word i is non-zero;
 * testing the resulting 8-bit mask against zero then compiles to a kortest.
 * This is the 512-bit analogue of the AVX2 version below.
 * The obvious fall-back (see the #else branch) is quite bad here: gcc spills
 * the whole 64-byte vector to the stack and reads it back in 8-byte pieces,
 * which defeats store-to-load forwarding in the innermost loop of stage 2. */
# define TEST(a) ( _mm512_test_epi64_mask((__m512i)(a), (__m512i)(a)) != 0 )
/* and the test the scan over the survivors uses; see TESTZ below */
# define TESTZ(a) ( _mm512_test_epi64_mask((__m512i)(a), (__m512i)(a)) == 0 )
#else
/* Fall-back version; also used when gcc lowers the 64-byte vectors itself
 * (i.e., when compiling with -DUSE_AVX512, but without -mavx512f). */
# define TEST(a) (EXT(a,0) || EXT(a,1) || EXT(a,2) || EXT(a,3) \
                   || EXT(a,4) || EXT(a,5) || EXT(a,6) || EXT(a,7))
#endif
#define MASKL(a,s) { unsigned long *survl = (unsigned long *)(a); long sh = (s); \
                     long l, qsh = sh>>LONG_SHIFT, rsh = sh & (LONG_LENGTH-1); \
                     for(l = 0; l < qsh; l++) { *survl++ = 0UL; }; *survl &= (~0UL)<<rsh; }
#define MASKU(a,s) { unsigned long *survl = (unsigned long *)(a); long sh = (s); \
                     long l, qsh = RBA_PACK-1 - (sh>>LONG_SHIFT), rsh = sh & (LONG_LENGTH-1); \
                     survl += qsh; *survl++ &= (~0UL)>>rsh; \
                     for(l = qsh+1; l < RBA_PACK; l++) { *survl++ = 0UL; } }
#ifndef RATPOINTS_CHUNK
/* Number of registers used in stage 1 of sieving.
 * One could use 32 here (there are as many ZMM registers),
 * but this would require extending the code in sift.c . */
# define RATPOINTS_CHUNK 16
#endif

#elif defined(USE_AVX)
/* Use 256 bit AVX registers for the bit arrays */

#include <immintrin.h>

#define RBA_LENGTH (256)
#define RBA_SHIFT (8)
#define RBA_PACK (4)
typedef unsigned long ratpoints_bit_array __attribute__ ((vector_size (32)));
#define AND(a,b) ((a) = (a)&(b))
#define EXT0(a) ((unsigned long)(a)[0])
#define EXT(a,i) ((unsigned long)(a)[i])
#ifdef __AVX2__
/* The following seems to be about the fastest way to test for zero,
 * see https://coderedirect.com/questions/445277/comparing-2-vectors-in-avx-avx2-c .
 * Note that this requires avx2, not just avx. */
# define TEST(a) ( _mm256_movemask_epi8(_mm256_cmpeq_epi8((__m256i)(a), (__m256i)zero)) != 0xffffffffU )
#elif defined(__AVX__)
/* This compiles to a vptest instruction */
# define TEST(a) ( !_mm256_testz_si256((__m256i)(a), (__m256i)(a)) )
#else
/* Fall-back version */
# define TEST(a) (EXT(a,0) || EXT(a,1) || EXT(a,2) || EXT(a,3))
#endif
#ifdef __AVX__
/* Whether a bit array is empty, for the loop in sift.c that scans the
 * survivors of the first stage.  That loop runs into a sentinel, so this
 * test is the whole of its body: a load, one vptest and the branch per bit
 * array.  (vptest against an all-ones register could take the memory
 * operand directly, but gcc 14 loads it anyway and then rebuilds the
 * constant inside the loop, which is one instruction more, not less.)
 * TEST above is left as it is: it tests a value that is already in a
 * register, and there the two forms cost the same.  This needs AVX, not
 * AVX2. */
# define TESTZ(a) ( _mm256_testz_si256((__m256i)(a), (__m256i)(a)) )
#endif
#define RBA(a) ((ratpoints_bit_array){((unsigned long) a), ((unsigned long) a), \
                                      ((unsigned long) a), ((unsigned long) a)})
#define zero (RBA(0LL))
#define MASKL(a,s) { unsigned long *survl = (unsigned long *)(a); long sh = (s); \
                     if(sh >= 2*LONG_LENGTH) \
                     { sh -= 2*LONG_LENGTH; survl[0] = 0UL; survl[1] = 0UL; \
                       if(sh >= LONG_LENGTH) \
                       { survl[2] = 0UL; survl[3] &= (~0UL)<<(sh - LONG_LENGTH); } \
                       else { survl[2] &= ~(0UL)<<sh; } } \
                     else if(sh >= LONG_LENGTH) { survl[0] = 0UL; survl[1] &= (~0UL)<<(sh - LONG_LENGTH); } \
                     else { survl[0] &= ~(0UL)<<sh; } }
#define MASKU(a,s) { unsigned long *survl = (unsigned long *)(a); long sh = (s); \
                     if(sh >= 2*LONG_LENGTH) \
                     { sh -= 2*LONG_LENGTH; survl[3] = 0UL; survl[2] = 0UL; \
                       if(sh >= LONG_LENGTH) \
                       { survl[0] &= ~(0UL)>>(sh - LONG_LENGTH); survl[1] = 0UL; } \
                       else { survl[1] &= ~(0UL)>>sh; } } \
                     else if(sh >= LONG_LENGTH) { survl[2] &= ~(0UL)>>(sh - LONG_LENGTH); survl[3] = 0UL; } \
                     else { survl[3] &= ~(0UL)>>sh; } }
#ifndef RATPOINTS_CHUNK
# define RATPOINTS_CHUNK 16  /* Number of registers used in stage 1 of sieving, max. 16. */
#endif

#elif defined(USE_AVX128)
/* Use 128 bit registers for the bit arrays */

#include <immintrin.h>

#define RBA_LENGTH (128)
#define RBA_SHIFT (7)
#define RBA_PACK (2)
typedef unsigned long ratpoints_bit_array __attribute__ ((vector_size (16)));
#define RBA(a) ((ratpoints_bit_array){((unsigned long) a), ((unsigned long) a)})
#define zero (RBA(0LL))
#define AND(a,b) ((a) = (a)&(b))
#define EXT0(a) ((unsigned long)(a)[0])
#define EXT(a,i) ((unsigned long)(a)[i])
/* See above for this definition of TEST(a) */
#define TEST(a) ( _mm_movemask_epi8(_mm_cmpeq_epi8((__m128i)(a), (__m128i)zero)) != 0xffffU )
#define MASKL(a,s) { unsigned long *survl = (unsigned long *)(a); long sh = (s); \
                     if(sh >= LONG_LENGTH) { survl[0] = 0UL; survl[1] &= (~0UL)<<(sh - LONG_LENGTH); } \
                     else { survl[0] &= ~(0UL)<<sh; } }
#define MASKU(a,s) { unsigned long *survl = (unsigned long *)(a); long sh = (s); \
                     if(sh >= LONG_LENGTH) { survl[0] &= ~(0UL)>>(sh - LONG_LENGTH); survl[1] = 0UL; } \
                     else { survl[1] &= ~(0UL)>>sh; } }
#ifndef RATPOINTS_CHUNK
# define RATPOINTS_CHUNK 16  /* Number of registers used in stage 1 of sieving, max. 16. */
#endif

#elif defined(USE_SSE)
/* Use SSE 128 bit SSE registers for the bit arrays */

#include <emmintrin.h>

#define RBA_LENGTH (128)
#define RBA_SHIFT (7)
#define RBA_PACK (2)
typedef __v2di ratpoints_bit_array;
#define RBA(a) ((__v2di){(a), (a)})
#define zero (RBA(0LL))
#define AND(a,b) ((a) = (ratpoints_bit_array)__builtin_ia32_andps((__v4sf)(a), (__v4sf)(b)))
#define EXT0(a) ((unsigned long)__builtin_ia32_vec_ext_v2di((__v2di)(a), 0))
/* It is a vector type, so subscript it like the other variants do. */
#define EXT(a,i) ((unsigned long)(a)[i])
#define TEST(a) (EXT0(a) || EXT(a,1))
#define MASKL(a,s) { unsigned long *survl = (unsigned long *)(a); long sh = (s); \
                     if(sh >= LONG_LENGTH) { survl[0] = 0UL; survl[1] &= (~0UL)<<(sh - LONG_LENGTH); } \
                     else { survl[0] &= ~(0UL)<<sh; } }
#define MASKU(a,s) { unsigned long *survl = (unsigned long *)(a); long sh = (s); \
                     if(sh >= LONG_LENGTH) { survl[0] &= ~(0UL)>>(sh - LONG_LENGTH); survl[1] = 0UL; } \
                     else { survl[1] &= ~(0UL)>>sh; } }
#ifndef RATPOINTS_CHUNK
# define RATPOINTS_CHUNK 16  /* Number of registers used in stage 1 of sieving, max. 16. */
#endif

#else
/* Use unsigned long for the bit arrays */

#define RBA_LENGTH LONG_LENGTH
#define RBA_SHIFT LONG_SHIFT
#define RBA_PACK (1)
typedef unsigned long ratpoints_bit_array;
#define RBA(a) ((ratpoints_bit_array)(a))
#define zero ((ratpoints_bit_array)0UL)
#define AND(a,b) ((a) &= (b))
#define EXT0(a) (a)
#define EXT(a,i) (a) /* just in case... */
#define TEST(a) (a)
#define MASKL(a,s) { *(a) &= ~(0UL)<<(s); }
#define MASKU(a,s) { *(a) &= ~(0UL)>>(s); }
#ifndef RATPOINTS_CHUNK
# define RATPOINTS_CHUNK 1  /* Leave optimization to the compiler... */
#endif
/* USE_LONG_IN_PHASE_2 makes the second stage sieve the survivors one 64-bit
 * word at a time instead of a whole bit array at a time; with RBA_PACK == 1
 * the two are the same code, so there is nothing to set here. */

#endif /* various register lengths */

/* Whether a bit array is all zero.  The loop that scans the survivors of
 * the first stage in sift.c uses this on the bit arrays in memory; the AVX
 * variants above have a form of their own, and everywhere else it is TEST
 * negated. */
#ifndef TESTZ
# define TESTZ(a) (!TEST(a))
#endif

/* The following is used for printing bit-arrays. */
#define WIDTH (LONG_LENGTH/4)

/* macro that prints a ratpoints_bit_array in hexadecimal */
#define PRINT_RBA(a) \
{ long i_; \
  for(i_ = RBA_PACK-1; i_; i_--) \
  { printf("%*.*lx", WIDTH, WIDTH, EXT((a), i_)); } \
  printf("%*.*lx ", WIDTH, WIDTH, EXT0(a)); }

/* set up data related to the set of primes considered */
#ifndef RATPOINTS_MAX_BITS_IN_PRIME
# define RATPOINTS_MAX_BITS_IN_PRIME 7
#endif

#if (RATPOINTS_MAX_BITS_IN_PRIME == 10)
# define RATPOINTS_NUM_PRIMES 171
# define RATPOINTS_MAX_PRIME 1021
# define RATPOINTS_MAX_PRIME_EVEN 1024

#elif (RATPOINTS_MAX_BITS_IN_PRIME == 9)
# define RATPOINTS_NUM_PRIMES 96
# define RATPOINTS_MAX_PRIME 509
# define RATPOINTS_MAX_PRIME_EVEN 512

#elif (RATPOINTS_MAX_BITS_IN_PRIME == 8)
# define RATPOINTS_NUM_PRIMES 53
# define RATPOINTS_MAX_PRIME 251
# define RATPOINTS_MAX_PRIME_EVEN 256

#elif (RATPOINTS_MAX_BITS_IN_PRIME == 7)
# define RATPOINTS_NUM_PRIMES 30
# define RATPOINTS_MAX_PRIME 127
# define RATPOINTS_MAX_PRIME_EVEN 128

#elif (RATPOINTS_MAX_BITS_IN_PRIME == 6)
# define RATPOINTS_NUM_PRIMES 17
# define RATPOINTS_MAX_PRIME 61
# define RATPOINTS_MAX_PRIME_EVEN 64

#elif (RATPOINTS_MAX_BITS_IN_PRIME == 5)
# define RATPOINTS_NUM_PRIMES 10
# define RATPOINTS_MAX_PRIME 31
# define RATPOINTS_MAX_PRIME_EVEN 32

#else
# define RATPOINTS_MAX_BITS_IN_PRIME 7
# define RATPOINTS_NUM_PRIMES 30
# define RATPOINTS_MAX_PRIME 127
# define RATPOINTS_MAX_PRIME_EVEN 128

#endif
/* so that RATPOINTS_MAX_PRIME < RATPOINTS_MAX_PRIME_EVEN
                               = 2^RATPOINTS_MAX_BITS_IN_PRIME */

/* define some datatypes */

/* This is used to hold the preliminary sieving information for one modulus
 * p, a prime or a composite modulus (see ratpoints_sieve_entry below).
 * The table at ptr has p bit arrays (and a few more repeating the first
 * ones, for the first stage to run past the end), and the one for word
 * number i is at index (i + offset) mod p.  Besides the shift that the
 * packing of the numerators needs (see rp_num_class), offset carries a
 * multiple of p of at least RP_ROW_BIAS, so that i + offset is never
 * negative for a word number the program handles and the reduction needs
 * no sign fix; the per-class table it is copied from in sift() (find_points.c)
 * has it built in.  start and end serve the first stage, which walks the
 * table with a pointer; the second stage computes the row from i directly. */
typedef struct { long p; long offset; ratpoints_bit_array *ptr;
                 ratpoints_bit_array *start; ratpoints_bit_array *end; } sieve_spec;

/* The multiple of p in the offset lies in [RP_ROW_BIAS, RP_ROW_BIAS + p).
 * The reductions in sift.c multiply by the reciprocal, which is exact below
 * 2^32, so with this value they are exact for every word number in
 * [-2^31, 2^31 - 2*RATPOINTS_MAX_PRIME_EVEN] (a composite modulus can
 * exceed the largest prime, by less than RATPOINTS_MAX_PRIME_EVEN);
 * _ratpoints_sift0 tests for that and divides otherwise. */
#define RP_ROW_BIAS 2147483648L

/* Reducing modulo a prime by multiplying rather than dividing.  With
 * m = 2^64/p rounded up (ULONG_MAX/p + 1; the sieve entry keeps it as
 * magic), the remainder of u modulo p is the top half of (m*u mod 2^64) * p,
 * and the quotient floor(u/p) is the top half of m*u; both are exact for
 * every u below 2^32, which every caller checks for in its own way.  The
 * callers are the third stage and its set-up, the start of the first stage
 * and the row look-up of the second (sift.c), and the reduction of the
 * denominator modulo each sieving prime and the Jacobi symbol test on the
 * denominators (find_points.c).
 *
 * Build with -DRP_MULMOD_DIVIDE to use the division everywhere instead,
 * which is what those callers cost without this. */
#if defined(__SIZEOF_INT128__) && !defined(RP_MULMOD_DIVIDE)
# define RP_MULMOD(u, p, m) \
    ((long)(unsigned long)(((__uint128_t)((m)*(unsigned long)(u)) \
                             * (unsigned long)(p)) >> 64))
# define RP_MULDIV(u, p, m) \
    ((long)(unsigned long)(((__uint128_t)(unsigned long)(u)*(m)) >> 64))
#else
/* the reciprocal is named so that the variables holding it stay used */
# define RP_MULMOD(u, p, m) \
    ((void)(m), (long)((unsigned long)(u) % (unsigned long)(p)))
# define RP_MULDIV(u, p, m) \
    ((void)(m), (long)((unsigned long)(u) / (unsigned long)(p)))
#endif

/* The largest value that is reduced that way; above it the callers fall
 * back on the division.  See mod_mul() and stage3() in sift.c. */
#define RP_MULMOD_LIMIT 4294967295L
/* The second stage keeps the value it reduces, a word number plus an offset
 * that carries RP_ROW_BIAS, below 2*RP_ROW_BIAS by the test at the head of
 * _ratpoints_sift0; that is only exact if the two limits agree.  An array
 * of negative size does not compile. */
typedef char rp_row_bias_within_mulmod_limit[
  (2*RP_ROW_BIAS - 1 <= RP_MULMOD_LIMIT) ? 1 : -1];

/* The position of the lowest set bit of a word: sift.c walks the set bits
 * of the survivors with it, find_points.c takes the odd part of a
 * denominator. */
#if defined(__GNUC__) || defined(__clang__)
# define RP_CTZL(w) ((long)__builtin_ctzl(w))
#else
/* Only reached on a compiler without the builtin; the rest of the program
 * needs gcc anyway once bit-arrays are used, but the plain unsigned long
 * build does not, so keep it buildable. */
static inline long RP_CTZL(unsigned long w)
{ long t = 0;

  while(!(w & 1UL)) { w >>= 1; t++; }
  return(t);
}
#endif

/* The inlining attributes sift.c relies on: accepted() and what it calls
 * must be inlined into the extraction sites, and the on-demand fill of the
 * third stage's data must not be, or gcc stops inlining accepted() (measured
 * at a per cent).  Empty on a compiler without them, for the reason above. */
#if defined(__GNUC__) || defined(__clang__)
# define RP_ALWAYS_INLINE __attribute__((always_inline))
# define RP_NOINLINE __attribute__((noinline))
#else
# define RP_ALWAYS_INLINE
# define RP_NOINLINE
#endif

/* What the third stage needs to test one numerator against one prime: the
 * prime, the inverse of the denominator modulo it, and the table saying
 * which residues admit points.  There is no sieve table, which is the point
 * of that stage: a prime costs it nothing per denominator, and the inverse
 * is looked up only for a denominator that brings a numerator this far.
 * binv == 0 means the prime divides the denominator, in which case it says
 * the same thing about every numerator and is skipped. */
typedef struct { long p; long binv; long bias; unsigned long magic;
                 const int *is_f_square; } check_spec;

/* The third stage reduces modulo p by multiplying instead of dividing, which
 * asks the value being reduced to fit in 32 bits; above this height bound it
 * would not, and the stage divides after all.  With the largest prime that
 * can be compiled in, that is a height of several million. */
#define RP_STAGE3_LIMIT 4294967296.0

/* The strides the numerators can be packed with are the powers of two up
 * to 64, the modulus of the 2-adic information: 2^k for 0 <= k < this. */
#define RP_NUM_STRIDES 7

/* What the sieve knows about the denominators of one residue class b mod
 * 64.  Their admissible numerators (those for which b^D f(a/b) is a square
 * mod 64, see get_2adic_info in find_points.c) lie in a single class a0 mod
 * 2^k, with k as large as that allows, so the bit arrays hold only those:
 * bit t of a bit array with word number 0 stands for the numerator
 * a0 + 2^k t, and a bit array sweeps 2^k times as many numerators as it has
 * bits.  bits is the 2-adic pattern in that packing -- the admissible ones
 * among a0 + 2^k t, of period 64/2^k in t and so one word repeated -- which
 * the first modulus ANDs into every bit array.  offset[n], for the n-th
 * modulus of sieve_list, is the shift of the table row that the packing needs, with
 * the multiple of p of sieve_spec built in: bit t wants the pattern for the
 * residue (a0 + 2^k t) b^-1 = (t + a0 2^-k) (b 2^-k)^-1 mod p, so the row for
 * the denominator is the one for b 2^-k mod p, read a0 2^-k bits further on,
 * that is a0 (2^k RBA_LENGTH)^-1 mod p bit arrays further on.  Classes with
 * the same k and a0 share one row of offsets.  A class without admissible
 * numerators has bits == 0 and is never sieved. */
typedef struct { ratpoints_bit_array bits; long k; long a0; const long *offset; }
        rp_num_class;

/* the type of the functions used for initializing the sieve */
typedef ratpoints_bit_array* (*ratpoints_init_fun)(void*, long, void*);

/* A sieving modulus need not be a prime.  A table row only has to be
 * periodic in the bit index with the modulus as its period, selected by the
 * denominator's residue, and the first two stages read it the same way
 * whatever the modulus is.  Three kinds of modulus are sieved with: a prime
 * p; a prime power p^e,
 * whose row for the residue b is the pattern of the a with F(a,b) a square
 * mod p^e -- f(a b^-1) for a unit b, frev(b a^-1) with frev(t) = t^D f(1/t)
 * for p | b and a unit a, exactly as get_2adic_info decides mod 2^6 --
 * which carries far more than the prime does (mod 9 about twice what mod 3
 * says, for nine rows); and a product of primes and prime powers, whose row
 * is the AND of its factors' rows and carries the information of all of
 * them for one AND per word.  Which moduli are used is decided by the same
 * ranking that chooses the primes, with the modulus's density the product
 * of its factors' and its own table cost; moduli sharing a prime exclude
 * one another.  The third stage, which indexes is_f_square by (a b^-1) mod
 * p, uses primes only.  The strides the numerators can be packed with
 * (rp_num_class) work with any odd modulus. */

/* the largest number of prime-power factors a modulus below 2^10 can have
 * (3*5*7*11 > 1024) */
#define RP_MAX_FACTORS 3
/* words holding a pattern of RATPOINTS_MAX_PRIME_EVEN bits (at least one:
 * the smallest prime size has 32) */
#define RP_MODWORDS ((RATPOINTS_MAX_PRIME_EVEN + LONG_LENGTH - 1)/LONG_LENGTH)

/* What a prime power p^e = m says about the curve: bit x of fsq is set when
 * f(x) is a square mod m, bit t of gsq (p | t) when frev(t) is one; inv[x]
 * is x^-1 mod m for a unit x; r the mean density of admissible numerators
 * over the classes of the denominator mod m, computed the same way as a
 * prime's (examine_prime), so that the two compete in one ranking.  Filled
 * by examine_power() in find_points.c. */
typedef struct { long p; long e; long m; double r; int np;
                 unsigned long fsq[RP_MODWORDS]; unsigned long gsq[RP_MODWORDS];
                 unsigned short inv[RATPOINTS_MAX_PRIME_EVEN]; }
        rp_power;

/* This is used to hold the sieving information for one modulus p: a prime
 * (nf == 0), a prime power (nf == 1 and pw set), or a product of nf >= 2
 * primes and prime powers, whose entries factor[] points to.  is_f_square
 * and inverses are the prime's tables and NULL for a composite modulus;
 * rbainv is RBA_LENGTH^-1 mod p, which the row shifts of the numerator
 * classes are computed from (class_offsets). */
typedef struct ratpoints_sieve_entry_s
        { ratpoints_init_fun init; long p; int *is_f_square;
          const long *inverses; unsigned long magic; double r;
          long bias;    /* the multiple of p in a row shift; see RP_ROW_BIAS */
          long rbainv;
          long dinv[RP_NUM_STRIDES]; /* 2^-k mod p: the denominator is
                                      * reduced to b 2^-k mod p for the
                                      * table row, see rp_num_class; filled
                                      * by class_offsets() for the strides
                                      * in use */
          long nf; struct ratpoints_sieve_entry_s *factor[RP_MAX_FACTORS];
          const rp_power *pw;
          ratpoints_bit_array* sieve[RATPOINTS_MAX_PRIME_EVEN]; }
        ratpoints_sieve_entry;

/* The state of one sieving thread.  Everything the path from a denominator
 * to its points writes is here, so that several of them can sieve the same
 * curve at once from one ratpoints_args that they only read: the survivors
 * array (array_size bit arrays, one spare for the alignment and one for the
 * sentinel the scan of the second stage runs into), the gmp temporaries
 * (work[0..2], and from work[3] on the coefficients times the powers of the
 * current denominator), the third stage's check_spec array, the numbers of
 * moduli in force, the counters of the block being sieved (which the run
 * folds into its own in block order, for the run-time correction; dec is
 * the correction's decision the block is sieved with), the two flags that
 * say what has been computed for the current denominator, and the flag
 * process() sets to stop the search.  A thread of a pool (find_points.c)
 * also knows its pool, the slot it puts its block's points into, and the
 * flag that tells it to abandon the block; on the calling thread the three
 * are NULL. */
typedef struct { ratpoints_args *args;
                 ratpoints_bit_array *survivors; void *survivors_na;
                 mpz_t *work;
                 check_spec *checks;
                 long sp1; long sp2; long sp3;
                 unsigned long n_words; unsigned long n_arrays;
                 unsigned long n_bits; unsigned long n_coprime;
                 unsigned long n_checks; unsigned long n_sifts;
                 long dec;
                 int compute_bc; int stage3_filled; int quit;
                 void *pool; void *slot; const int *stop; }
        rp_worker;

/* The following two functions are provided in init.c: the table row of a
 * prime power and of a product of moduli for the residue b, built into the
 * table buffer like the rows of the primes (the CODE_INIT_SIEVE functions
 * there), with -DRP_VERIFY_MODULI every row is checked against a direct
 * evaluation of F(a,b) mod m as it is built. */
ratpoints_bit_array *_ratpoints_sieve_init_power(void *se1, long b1, void *args1);
ratpoints_bit_array *_ratpoints_sieve_init_product(void *se1, long b1, void *args1);

/* The following function is provided in find_points.c : */
long _ratpoints_check_point(long a, long b, rp_worker *wk,
                 int process(long, long, const mpz_t, void*, int*), void *info);

/* The following function is provided in sift.c : */
/* cls is the numerator class of the denominator: its packing, which says
 * what numerator a bit stands for, and the 2-adic pattern every bit array
 * starts from, which the first stage ANDs in as it sieves rather than having
 * it written into the array beforehand.  mask_low and mask_high say how many
 * bits to clear at the two ends of the numerator interval (zero for an end
 * that is not a boundary); both are applied after the first stage, which
 * gives the same result because AND is commutative.  The range
 * w_high - w_low can be any length: the first stage takes it in chunks of
 * RATPOINTS_CHUNK bit arrays and sieves what is left over in narrower legs
 * (the arm for RATPOINTS_CHUNK 1 takes any length as well). */
long _ratpoints_sift0(long b, long w_low, long w_high,
           rp_worker *wk, const rp_num_class *cls,
           long mask_low, long mask_high, sieve_spec *sieves,
           check_spec *checks,
           int process(long, long, const mpz_t, void*, int*), void *info);

/* The following function is provided in sturm.c : */
long _ratpoints_compute_sturm(ratpoints_args*);
