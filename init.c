/***********************************************************************
 * ratpoints-2.2                                                       *
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
 * init.c                                                              *
 *                                                                     *
 * Macro definitions for the sieve_init functions                      *
 *                                                                     *
 * Michael Stoll, Apr 14, 2009, January 7, 2022                        *
 * with changes by Bill Allombert, Dec 29, 2021                        *
 ***********************************************************************/

#include "rp-private.h"

/* Define functions that initialize the sieve
 * for a given prime p and denominator b1 mod p. */

/* This is a bit different depending on whether p is smaller
 * or larger than the number LONG_LENGTH of bits in a word.
 * This is because we have to repeat a pattern of p bits. */

/* The following is for primes < LONG_LENGTH */
#define CODE_INIT_SIEVE1(prime) \
static ratpoints_bit_array *sieve_init_##prime(void *se1, long b1, void *args1) \
{ \
  ratpoints_sieve_entry *se = se1; \
  ratpoints_args *args = args1; \
  int *isfs = se->is_f_square; \
  long b = b1; \
  long lmp = LONG_LENGTH % (prime); \
  long ldp = LONG_LENGTH / (prime); \
  long p1 = (ldp + 1) * (prime); \
  long diff_shift = p1 & LONG_MASK; \
  long diff = LONG_LENGTH - diff_shift; \
  unsigned long help0;\
  { long a; \
    long d = se->inverses[b]; \
    long ab = 0; /* a/b mod p */ \
    unsigned long test = 1UL; \
    unsigned long he0 = 0UL; \
    for(a = 0; a < (prime); a++) \
    { if(isfs[ab]) { he0 |= test; } \
      ab += d; \
      if(ab >= (prime)) { ab -= (prime); } \
      test <<= 1; \
    } \
    help0 = he0; \
  } \
  \
  { unsigned long help1; \
    { /* repeat bit pattern floor(LONG_LENGTH/p) times */ \
      unsigned long pattern = help0; \
      long i; \
      /* the p * (floor(LONG_LENGTH/p) + 1) - LONG_LENGTH \
              = p - (LONG_LENGTH mod p) \
         upper bits into help[b][1] : \
         shift away the  LONG_LENGTH mod p  lower bits */ \
      help1 = pattern >> lmp; \
      for(i = (prime); i < LONG_LENGTH; i <<= 1) \
      { help0 |= help0 << i; } \
      /* \
      for(i = ldp; i; i--) \
      { pattern <<= (prime); help0 |= pattern; } \ */ \
    } \
    \
    { /* fill the bit pattern from help0/help1 into sieve[b][]. \
          sieve[b][a0] has the same semantics as help0/help1, \
          but here, a0 runs from 0 to p-1 and all bits are filled. */ \
      long a, k; \
      unsigned long *si = (unsigned long *)args->ba_next; \
      \
      args->ba_next += ((prime) + RATPOINTS_CHUNK-1)*sizeof(ratpoints_bit_array); \
      /* copy the first chunk into sieve[b][] */ \
      si[0] = help0; \
      /* now keep repeating the bit pattern, \
         rotating it in help0/help1 */ \
      for(a = 1 ; a < (prime); a++) \
      { unsigned long temp = help0 >> diff; \
        help0 = help1 | (help0 << diff_shift); \
        si[a] = help0; \
        help1 = temp; \
      } \
      /* copy into the next p*(RBA_PACK-1) long words \
       * (the compiler will eliminate this loop when RBA_PACK == 1) */ \
      for (a = 0; a < (prime); a++) \
      { for(k = 1; k < RBA_PACK; k++) \
        { si[a+k*(prime)] = si[a]; } \
      } \
      /* append a copy of the first (RATPOINTS_CHUNK-1)*RBA_PACK words at the end */ \
      for(k = 0; k < (RATPOINTS_CHUNK-1)*RBA_PACK; k++) \
      { si[(prime)*RBA_PACK + k] = si[k]; } \
      /* set sieve array and return the pointer */ \
      se->sieve[b] = (ratpoints_bit_array *)si; \
      return((ratpoints_bit_array *)si); \
  } } \
}

/* This is for p > LONG_LENGTH */
#define CODE_INIT_SIEVE2(prime) \
static ratpoints_bit_array *sieve_init_##prime(void *se1, long b1, void *args1) \
{ \
  ratpoints_sieve_entry *se = se1; \
  ratpoints_args *args = args1; \
  long p = (prime); \
  int *isfs = se->is_f_square; \
  long b = b1; \
  long wp = p >> LONG_SHIFT; \
  long diff_shift = p & LONG_MASK; \
  long diff = LONG_LENGTH - diff_shift; \
  unsigned long help[(p>>LONG_SHIFT) + 2]; \
  \
  /* initialize help */ \
  { unsigned long *he = &help[0]; \
    unsigned long *he1 = &he[(p>>LONG_SHIFT) + 2]; \
    while(he1 != he) { he1--; *he1 = 0UL; } \
  } \
  { unsigned long work = 0UL; \
    long a; \
    long ab = 0; /* a/b mod p */ \
    long d = se->inverses[b]; \
    long n = 0; \
    unsigned long test = 1UL;  \
    for(a = 0; a < p; ) \
    { if(isfs[ab]) { work |= test; } \
      ab += d; \
      if(ab >= p) { ab -= p; } \
      test <<= 1; \
      a++; \
      if((a & LONG_MASK) == 0) \
      { help[n] = work; n++; work = 0UL; test = 1UL; } \
    } \
    help[n] = work; \
  } \
  \
  { /* fill the bit pattern from help[] into sieve[b][]. \
       sieve[b][a0] has the same semantics as help[b][a0], \
       but here, a0 runs from 0 to p-1 and all bits are filled. */ \
    unsigned long *si = (unsigned long *)args->ba_next; \
    long a1; \
    long a, k; \
    \
    args->ba_next += (p + RATPOINTS_CHUNK-1)*sizeof(ratpoints_bit_array); \
    /* copy the first chunk from help[] into sieve[num][b][] */ \
    for(a = 0; a < wp; a++) { si[a] = help[a]; } \
    /* now keep repeating the bit pattern, rotating it in help */ \
    for(a1 = a ; a < p; a++) \
    { long t = (a1 == wp) ? 0 : a1+1; \
      help[a1] |= help[t]<<diff_shift; \
      si[a] = help[a1]; \
      a1 = t; \
      help[a1] >>= diff; \
    } \
    /* copy into the next p*(RBA_PACK-1) long words \
     * (the compiler will eliminate this loop when RBA_PACK == 1) */ \
    for (a = 0; a < (prime); a++) \
    { for(k = 1; k < RBA_PACK; k++) \
      { si[a+k*(prime)] = si[a]; } \
    } \
    /* append a copy of the first (RATPOINTS_CHUNK-1)*RBA_PACK words at the end */ \
    for(k = 0; k < (RATPOINTS_CHUNK-1)*RBA_PACK; k++) \
    { si[(prime)*RBA_PACK + k] = si[k]; } \
    /* set sieve array and return the pointer */ \
    se->sieve[b] = (ratpoints_bit_array *)si; \
    return((ratpoints_bit_array *)si); \
  } \
}

#include "init_sieve.h"
