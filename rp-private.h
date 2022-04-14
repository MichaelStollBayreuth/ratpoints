/***********************************************************************
 * ratpoints-2.2.1                                                     *
 *  - A program to find rational points on hyperelliptic curves        *
 * Copyright (C) 2008, 2009, 2022  Michael Stoll                       *
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
 * Michael Stoll, Apr 14, 2009; January 7-18, 2022                     *
 * with changes by Bill Allombert, Dec 29, 2021                        *
 ***********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define LONG_LENGTH (8*sizeof(long))
   /* number of bits in an unsigned long */
#define LONG_SHIFT ((LONG_LENGTH == 16) ? 4 : \
                    (LONG_LENGTH == 32) ? 5 : \
		    (LONG_LENGTH == 64) ? 6 : 0)
#define LONG_MASK (~(-(1UL<<LONG_SHIFT)))

/* Check if SSE instructions can be used.
 * We assume that one SSE word of 128 bit is two long's,
 * so check that a long is 64 bit = 8 byte. */
#ifndef __SSE2__
#undef USE_SSE
#endif
#if __WORDSIZE != 64
#undef USE_SSE
#endif

#include "ratpoints.h"

#define FLOOR(a,b) (((a) < 0) ? -(1 + (-(a)-1) / (b)) : (a) / (b))
#define CEIL(a,b) (((a) <= 0) ? -(-(a) / (b)) : 1 + ((a)-1) / (b))

/* Define interface for ratpoints_bit_array datatype:
 * RBA_LENGTH : number of bits
 * RBA_SHIFT  : 2^RBA_SHIFT == RBA_LENGTH
 * RBA_PACK   : number of words in a ratpoints_bit_array == RBA_LENGTH/LONG_LENGTH
 * RBA(a)     : fill a ratpoints_bit_array with copies of the word a
 * zero       : all bits zero == RBA(0UL)
 * AND(a,b)   : bit-wise and operation: a &= b
 * EXT0(a)    : extract first word (as unsigned long)
 * EXT(a,i)   : extract word with index i (as unsigned long)
 * TEST(a)    : tests if a is zero: TEST(a) == 0 <==> a == zero
 *              TEST should be fast if possible; it is used frequently
 *              in phase 2 of the sieve.
 * MASKL(a,s) : set lower s bits of a to zero
 * MASKU(a,s) : set upper s bits of a to zero
 *              MASKL and MASKU don't have to be terribly efficient;
 *              they are each executed once per denominator and interval.
 */

#ifdef USE_AVX512
/* Use 512 bit AVX registers for the bit arrays */
/* So far not used, since no suitable CPU available for testing... */

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
#define EXT0(a) ((unsigned long)a[0])
#define EXT(a,i) ((unsigned long)a[i])
/* there should be a faster way of doing the following; compare below */
#define TEST(a) (EXT(a,0) || EXT(a,1) || EXT(a,2) || EXT(a,3) \
                  || EXT(a,4) || EXT(a,5) || EXT(a,6) || EXT(a,7))
#define MASKL(a,s) { unsigned long *survl = (unsigned long *)(a); long sh = (s); \
                     long l, qsh = sh>>LONG_SHIFT, rsh = sh & (LONG_LENGTH-1); \
                     for(l = 0; l < qsh; l++) { *survl++ = 0UL; }; *survl &= (~0UL)<<rsh; }
#define MASKU(a,s) { unsigned long *survl = (unsigned long *)(a); long sh = (s); \
                     long l, qsh = RBA_PACK-1 - (sh>>LONG_SHIFT), rsh = sh & (LONG_LENGTH-1); \
                     survl += qsh; *survl++ &= (~0UL)>>rsh; \
                     for(l = qsh+1; l < RBA_PACK; l++) { *survl++ = 0UL; } }
#ifndef RATPOINTS_CHUNK
/* Number of registers used in phase 1 of sieving.
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
#define EXT0(a) ((unsigned long)a[0])
#define EXT(a,i) ((unsigned long)a[i])
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
# define RATPOINTS_CHUNK 16  /* Number of registers used in phase 1 of sieving, max. 16. */
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
#define EXT0(a) ((unsigned long)a[0])
#define EXT(a,i) ((unsigned long)a[i])
/* See above for this definition of TEST(a) */
#define TEST(a) ( _mm_movemask_epi8(_mm_cmpeq_epi8((__m128i)(a), (__m128i)zero)) != 0xffffU )
#define MASKL(a,s) { unsigned long *survl = (unsigned long *)(a); long sh = (s); \
                     if(sh >= LONG_LENGTH) { survl[0] = 0UL; survl[1] &= (~0UL)<<(sh - LONG_LENGTH); } \
                     else { survl[0] &= ~(0UL)<<sh; } }
#define MASKU(a,s) { unsigned long *survl = (unsigned long *)(a); long sh = (s); \
                     if(sh >= LONG_LENGTH) { survl[0] &= ~(0UL)>>(sh - LONG_LENGTH); survl[1] = 0UL; } \
                     else { survl[1] &= ~(0UL)>>sh; } }
#ifndef RATPOINTS_CHUNK
# define RATPOINTS_CHUNK 16  /* Number of registers used in phase 1 of sieving, max. 16. */
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
/* the following is a hack (here i is always 1) */
#define EXT(a,i) ((unsigned long)__builtin_ia32_vec_ext_v2di((__v2di)(a), 1))
#define TEST(a) (EXT0(a) || EXT(a,1))
#define MASKL(a,s) { unsigned long *survl = (unsigned long *)(a); long sh = (s); \
                     if(sh >= LONG_LENGTH) { survl[0] = 0UL; survl[1] &= (~0UL)<<(sh - LONG_LENGTH); } \
                     else { survl[0] &= ~(0UL)<<sh; } }
#define MASKU(a,s) { unsigned long *survl = (unsigned long *)(a); long sh = (s); \
                     if(sh >= LONG_LENGTH) { survl[0] &= ~(0UL)>>(sh - LONG_LENGTH); survl[1] = 0UL; } \
                     else { survl[1] &= ~(0UL)>>sh; } }
#ifndef RATPOINTS_CHUNK
# define RATPOINTS_CHUNK 16  /* Number of registers used in phase 1 of sieving, max. 16. */
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
#ifndef USE_LONG_IN_PHASE_2 /* Use simpler code in sift.c */
# define USE_LONG_IN_PHASE_2
#endif

#endif /* various register lengths */

/* The following is used for printing bit-arrays. */
#define WIDTH (int)(LONG_LENGTH/4)

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

/* this is used to hold the preliminary sieving information for one prime p */
typedef struct { long p; long offset; ratpoints_bit_array *ptr;
                 ratpoints_bit_array *start; ratpoints_bit_array *end; } sieve_spec;

/* this is used to record whether all / only even / only odd / no numerators
 * need to be considered */
typedef enum { num_all, num_even, num_odd, num_none } bit_selection;

/* the type of the functions used for initializing the sieve */
typedef ratpoints_bit_array* (*ratpoints_init_fun)(void*, long, void*);

/* this is used to hold the sieving information for one prime p */
typedef struct
        { ratpoints_init_fun init; long p; int *is_f_square;
          const long *inverses;
          long offset; ratpoints_bit_array* sieve[RATPOINTS_MAX_PRIME]; }
        ratpoints_sieve_entry;

/* The following function is provided in find_points.c : */
long _ratpoints_check_point(long a, long b, ratpoints_args *args, int *quit,
                 int process(long, long, const mpz_t, void*, int*), void *info);

/* The following function is provided in sift.c : */
long _ratpoints_sift0(long b, long w_low, long w_high,
           ratpoints_args *args, bit_selection which_bits,
           ratpoints_bit_array *survivors, sieve_spec *sieves, int *quit,
           int process(long, long, const mpz_t, void*, int*), void *info);

/* The following function is provided in sturm.c : */
long _ratpoints_compute_sturm(ratpoints_args*);
