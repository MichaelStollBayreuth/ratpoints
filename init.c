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
#include <string.h>

/* Define functions that initialize the sieve
 * for a given prime p and denominator b1 mod p. */

/* Setting bit i of w from is_f_square, which is 0 or 1.  The obvious
 *   if(isfs[ab]) { w |= 1UL << i; }
 * is a branch on data with no pattern in it: at p = 127 it accounts for a
 * third of all the branches executed and misses about half the time, and the
 * loop it sits in is three quarters of the whole set-up.  Shifting the value
 * into place instead has no branch to miss.  RP_INIT_BRANCH restores the old
 * form for comparison. */
#ifdef RP_INIT_BRANCH
# define RP_INIT_BIT(w, v, i) do { if(v) { (w) |= 1UL << (i); } } while(0)
#else
# define RP_INIT_BIT(w, v, i) ((w) |= (unsigned long)(v) << (i))
#endif
/* Four bits at a time.  What is left in the loop above once the branch has
 * gone is the step  ab += d; if(ab >= p) ab -= p;  and that is a chain: add,
 * compare, conditional move, three cycles, and the next row cannot start until
 * it finishes.  Measured at p = 127 the loop runs at 3.2 cycles a row, which
 * is exactly that.
 *
 * Since a runs over consecutive integers, four rows can be taken at once by
 * keeping the four residues a*d, (a+1)*d, (a+2)*d, (a+3)*d and stepping each
 * by 4d mod p.  The four chains are independent, so they overlap, and the loop
 * becomes bound by the work rather than by the latency.  Groups of four never
 * straddle a word, since LONG_LENGTH is a multiple of four; the caller runs a
 * plain loop over whatever is left of the last word.
 *
 * ab is the residue of the first of the four, so the caller can carry on from
 * it.  RP_INIT_ONEWAY does without all of this, for comparison. */
#ifdef RP_INIT_ONEWAY
# define RP_INIT_FOUR_DECL(prime)              long ab4_unused = 0; (void)ab4_unused;
# define RP_INIT_FOUR_LOOP(prime, w, i, lo, hi) (i) = (lo);
#else
# define RP_INIT_FOUR_DECL(prime) \
    long ab1 = d, ab2, ab3, d4; \
    ab2 = ab1 + d; if(ab2 >= (prime)) { ab2 -= (prime); } \
    ab3 = ab2 + d; if(ab3 >= (prime)) { ab3 -= (prime); } \
    d4  = ab3 + d; if(d4  >= (prime)) { d4  -= (prime); }

# define RP_INIT_FOUR_LOOP(prime, w, i, lo, hi) \
    for((i) = (lo); (i) + 4 <= (hi); (i) += 4) \
    { (w) |= ((unsigned long)isfs[ab]  << (i)) \
           | ((unsigned long)isfs[ab1] << ((i)+1)) \
           | ((unsigned long)isfs[ab2] << ((i)+2)) \
           | ((unsigned long)isfs[ab3] << ((i)+3)); \
      ab  += d4; if(ab  >= (prime)) { ab  -= (prime); } \
      ab1 += d4; if(ab1 >= (prime)) { ab1 -= (prime); } \
      ab2 += d4; if(ab2 >= (prime)) { ab2 -= (prime); } \
      ab3 += d4; if(ab3 >= (prime)) { ab3 -= (prime); } \
    }
#endif

/* The rows of the byte fill.  One at a time it is no better than the rotation
 * it replaces, and for the same reason the accumulation was slow: the offset
 * step  s += LONG_LENGTH; if(s >= p) s -= p;  is a three-cycle chain, so one
 * chain has merely been swapped for another.  Four at a time breaks it the way
 * RP_INIT_FOUR_LOOP does, with four offsets a step of 4*LONG_LENGTH apart.
 * The byte form is what makes this convenient: an offset that is a multiple of
 * LONG_LENGTH is not a special case for it, as it would be for a double-word
 * shift.  RP_INIT_BYTE1 keeps the one-at-a-time version for comparison. */
#define RP_INIT_BYTE_ROW(a, s) \
  memcpy(&si[a], (unsigned char *)&bufs[(s) & 7][0] + ((s) >> 3), \
         sizeof(unsigned long))
#ifdef RP_INIT_BYTE1
# define RP_INIT_BYTE_ROWS(prime) \
    { long s = 0; \
      for(a = 0; a < (prime); a++) \
      { RP_INIT_BYTE_ROW(a, s); \
        s += LONG_LENGTH; if(s >= (prime)) { s -= (prime); } } }
#else
# define RP_INIT_BYTE_ROWS(prime) \
    { long s0 = 0, s1 = LONG_LENGTH % (prime), s2 = (2*LONG_LENGTH) % (prime), \
           s3 = (3*LONG_LENGTH) % (prime), st = (4*LONG_LENGTH) % (prime); \
      for(a = 0; a + 4 <= (prime); a += 4) \
      { RP_INIT_BYTE_ROW(a,   s0); \
        RP_INIT_BYTE_ROW(a+1, s1); \
        RP_INIT_BYTE_ROW(a+2, s2); \
        RP_INIT_BYTE_ROW(a+3, s3); \
        s0 += st; if(s0 >= (prime)) { s0 -= (prime); } \
        s1 += st; if(s1 >= (prime)) { s1 -= (prime); } \
        s2 += st; if(s2 >= (prime)) { s2 -= (prime); } \
        s3 += st; if(s3 >= (prime)) { s3 -= (prime); } \
      } \
      for(; a < (prime); a++) \
      { RP_INIT_BYTE_ROW(a, s0); \
        s0 += LONG_LENGTH; if(s0 >= (prime)) { s0 -= (prime); } } }
#endif

/* A probe, not a variant: with RP_INIT_NOACC the pattern is never computed,
 * so what is timed is the filling in of the table and nothing else.  The
 * tables it produces are wrong; bench_init will say so. */
#ifdef RP_INIT_NOACC
# undef  RP_INIT_BIT
# undef  RP_INIT_FOUR_LOOP
# define RP_INIT_BIT(w, v, i) ((void)0)
# define RP_INIT_FOUR_LOOP(prime, w, i, lo, hi) (i) = (lo);
#endif
/* A second probe: RP_INIT_NOREP leaves out the copies that replicate the
 * table over the RBA_PACK words of a bit array and wrap the first chunk round
 * to the end.  Same caveat -- the tables are then wrong. */
#ifdef RP_INIT_NOREP
# define RP_INIT_REPLICATE(prime) /* nothing */
#else
# define RP_INIT_REPLICATE(prime) \
    /* copy into the next p*(RBA_PACK-1) long words \
     * (the compiler will eliminate this loop when RBA_PACK == 1) */ \
    for (a = 0; a < (prime); a++) \
    { for(k = 1; k < RBA_PACK; k++) \
      { si[a+k*(prime)] = si[a]; } \
    } \
    /* append a copy of the first (RATPOINTS_CHUNK-1)*RBA_PACK words at the end */ \
    for(k = 0; k < (RATPOINTS_CHUNK-1)*RBA_PACK; k++) \
    { si[(prime)*RBA_PACK + k] = si[k]; } \

#endif

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
    unsigned long he0 = 0UL; \
    RP_INIT_FOUR_DECL((prime)) \
    RP_INIT_FOUR_LOOP((prime), he0, a, 0, (prime)) \
    for(; a < (prime); a++) \
    { RP_INIT_BIT(he0, isfs[ab], a); \
      ab += d; \
      if(ab >= (prime)) { ab -= (prime); } \
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
      RP_INIT_REPLICATE(prime) \
      /* set sieve array and return the pointer */ \
      se->sieve[b] = (ratpoints_bit_array *)si; \
      return((ratpoints_bit_array *)si); \
  } } \
}

/* Filling in sieve[b][] for p > LONG_LENGTH, in two ways.
 *
 * What is wanted is  si[a] = the LONG_LENGTH bits of the p-periodic pattern
 * starting at bit  a*LONG_LENGTH,  which by periodicity is the window at bit
 * offset  s(a) = a*LONG_LENGTH mod p.
 *
 * The old way (RP_INIT_ROTATE) advances from one row to the next by rotating
 * the pattern in place.  That is two operations a row, but they form a chain:
 * row a cannot start until row a-1 has finished, for p-1 rows.
 *
 * The new way computes each row from the pattern directly.  Extend the pattern
 * past bit p far enough that every window wanted lies inside it -- one extra
 * word does it -- and then row a is a double-word shift at offset s(a), with
 * nothing carried from the row before but s itself.  The rows are independent,
 * so the processor can work on several at once.
 */
#define RP_INIT_FILL2_ROTATE(prime) \
    /* copy the first chunk from help[] into sieve[num][b][] */ \
    for(a = 0; a < wp; a++) { si[a] = help[a]; } \
    /* now keep repeating the bit pattern, rotating it in help */ \
    { long a1; \
    for(a1 = a ; a < p; a++) \
    { long t = (a1 == wp) ? 0 : a1+1; \
      help[a1] |= help[t]<<diff_shift; \
      si[a] = help[a1]; \
      a1 = t; \
      help[a1] >>= diff; \
    } } \

#define RP_INIT_FILL2_WINDOW(prime) \
    /* continue the pattern past its own end, so that the last window wanted, \
     * the one at bit p-1, still lies inside help[].  Bit p+k of the extended \
     * pattern is bit k of it, so this shifts help[] up by p bits and ors it \
     * in; only the top two words can be reached by a window, and help[wp+1] \
     * has to be written before help[wp] since the two coincide when wp == 1. */ \
    { help[wp+1] = (help[0] >> diff) | (help[1] << diff_shift); \
      help[wp]  |= help[0] << diff_shift; } \
    /* a window at a multiple of LONG_LENGTH is a word of help[] as it stands, \
     * and those are exactly the rows a = 0, ..., wp */ \
    for(a = 0; a <= wp; a++) { si[a] = help[a]; } \
    /* every other row straddles two words.  s is the bit offset of row a; \
     * it is the only thing that passes from one iteration to the next */ \
    { long s = LONG_LENGTH*(wp+1) - p; \
      for(a = wp+1; a < p; a++) \
      { long q = s >> LONG_SHIFT, r = s & LONG_MASK; \
        si[a] = (help[q] >> r) | (help[q+1] << (LONG_LENGTH - r)); \
        s += LONG_LENGTH; \
        if(s >= p) { s -= p; } \
      } } \

/* The third way, the one TODO item 4 proposes.  A rotation by r bits is a
 * rotation by r mod 8 bits followed by whole bytes, so hold the extended
 * pattern in eight copies, copy r shifted along by r bits, and every row is
 * then eight bytes read out of one of them at a byte offset: one unaligned
 * load, no shifting at all.  The eight copies cost eight passes over
 * wp+2 words to build, which for the primes that matter is a fifth of the
 * rows they save work on. */
#define RP_INIT_FILL2_BYTE(prime) \
    { help[wp+1] = (help[0] >> diff) | (help[1] << diff_shift); \
      help[wp]  |= help[0] << diff_shift; } \
    { long r, i; \
      unsigned long bufs[8][((prime) >> LONG_SHIFT) + 2]; \
      for(i = 0; i <= wp+1; i++) { bufs[0][i] = help[i]; } \
      for(r = 1; r < 8; r++) \
      { for(i = 0; i <= wp; i++) \
        { bufs[r][i] = (help[i] >> r) | (help[i+1] << (LONG_LENGTH - r)); } \
        bufs[r][wp+1] = help[wp+1] >> r; \
      } \
      RP_INIT_BYTE_ROWS(prime) \
    } \

#if defined(RP_INIT_BYTE)
# define RP_INIT_FILL2(prime) RP_INIT_FILL2_BYTE(prime)
#elif defined(RP_INIT_WINDOW)
# define RP_INIT_FILL2(prime) RP_INIT_FILL2_WINDOW(prime)
#else
# define RP_INIT_FILL2(prime) RP_INIT_FILL2_ROTATE(prime)
#endif

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
  { long n, i; \
    long ab = 0; /* a/b mod p */ \
    long d = se->inverses[b]; \
    RP_INIT_FOUR_DECL(p) \
    for(n = 0; n <= wp; n++) \
    { unsigned long work = 0UL; \
      long m = (n < wp) ? LONG_LENGTH : (p & LONG_MASK); \
      RP_INIT_FOUR_LOOP(p, work, i, 0, m) \
      for(; i < m; i++) \
      { RP_INIT_BIT(work, isfs[ab], i); \
        ab += d; \
        if(ab >= p) { ab -= p; } \
      } \
      help[n] = work; \
    } \
  } \
  \
  { /* fill the bit pattern from help[] into sieve[b][]. \
       sieve[b][a0] has the same semantics as help[b][a0], \
       but here, a0 runs from 0 to p-1 and all bits are filled. */ \
    unsigned long *si = (unsigned long *)args->ba_next; \
    long a, k; \
    \
    args->ba_next += (p + RATPOINTS_CHUNK-1)*sizeof(ratpoints_bit_array); \
    RP_INIT_FILL2(prime) \
    RP_INIT_REPLICATE(prime) \
    /* set sieve array and return the pointer */ \
    se->sieve[b] = (ratpoints_bit_array *)si; \
    return((ratpoints_bit_array *)si); \
  } \
}

#include "init_sieve.h"
