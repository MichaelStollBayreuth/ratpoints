/***********************************************************************
 * ratpoints-3.0.0                                                     *
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
 * init.c                                                              *
 *                                                                     *
 * Macro definitions for the sieve_init functions                      *
 *                                                                     *
 * Michael Stoll, Apr 14, 2009; Jan 7, 2022; Sep 21, 2026              *
 * with changes by Bill Allombert, Dec 29, 2021                        *
 ***********************************************************************/

#include "rp-private.h"

/* Define functions that initialize the sieve
 * for a given prime p and denominator b1 mod p. */

/* Setting bit i of w from is_f_square, which is 0 or 1.  The obvious
 *   if(isfs[ab]) { w |= 1UL << i; }
 * is a branch on data with no pattern in it: at p = 127 it accounts for a
 * third of all the branches executed and misses about half the time, and the
 * loop it sits in is three quarters of the whole set-up.  Shifting the value
 * into place instead has no branch to miss.  RP_INIT_BRANCH selects the
 * branching form for comparison. */
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
    RP_INIT_REPLICATE(prime) \
    /* set sieve array and return the pointer */ \
    se->sieve[b] = (ratpoints_bit_array *)si; \
    return((ratpoints_bit_array *)si); \
  } \
}

#include "init_sieve.h"

/************************************************************************
 * The rows of a composite modulus                                      *
 ************************************************************************/

/* Lay an m-periodic pattern of admissible bit indices out as a table row of
 * m bit arrays (plus the RATPOINTS_CHUNK-1 wrap-around copies), taken from
 * the table buffer like the rows of the primes.  Bit j of bit array x
 * stands for the bit index x*RBA_LENGTH + j, so word l of bit array x is
 * the 64 bits of the pattern from (x*RBA_LENGTH + 64 l) mod m on.  The
 * pattern is first repeated so that those 64 bits can be read from two
 * words wherever they start. */
static ratpoints_bit_array *lay_out_pattern(const unsigned long *pat, long m,
                                            ratpoints_args *args)
{ unsigned long rep[RP_MODWORDS + 2];
  unsigned long *si = (unsigned long *)args->ba_next;
  long words = (m + LONG_LENGTH + LONG_LENGTH - 1) >> LONG_SHIFT;
  long i, x, o, step;

  args->ba_next += (m + RATPOINTS_CHUNK-1)*sizeof(ratpoints_bit_array);

  /* the pattern, repeated over words*64 >= m + 64 bits */
  for(i = 0; i < words; i++) { rep[i] = 0UL; }
  for(o = 0; o < words*LONG_LENGTH; o += m)
  { /* OR the pattern in at bit offset o: word by word, shifted */
    long q = o >> LONG_SHIFT, s = o & LONG_MASK;
    long pw = (m + LONG_LENGTH - 1) >> LONG_SHIFT;

    for(i = 0; i < pw && q + i < words; i++)
    { rep[q + i] |= pat[i] << s;
      if(s && q + i + 1 < words) { rep[q + i + 1] |= pat[i] >> (LONG_LENGTH - s); }
    }
  }
  /* the pattern has period m: bits beyond the last full copy were ORed in
   * as far as they went; a copy that ran past words*64 was cut, which is
   * what is wanted */

  /* the rows: word l of bit array x starts at offset (x*RBA_LENGTH + 64 l)
   * mod m; the offset of the next word is 64 further on, mod m */
  step = LONG_LENGTH % m;
  o = 0;
  for(x = 0; x < m*RBA_PACK; x++)
  { long q = o >> LONG_SHIFT, s = o & LONG_MASK;

    si[x] = s ? (rep[q] >> s) | (rep[q + 1] << (LONG_LENGTH - s)) : rep[q];
    o += step; if(o >= m) { o -= m; }
  }
  /* the wrap-around copies */
  for(i = 0; i < (RATPOINTS_CHUNK-1)*RBA_PACK; i++) { si[m*RBA_PACK + i] = si[i]; }
  return((ratpoints_bit_array *)si);
}

#ifdef RP_VERIFY_MODULI
/* Check a row against what the program promises, factor by factor: for a
 * prime-power factor q = p^e of m, bit index N of the row for the residue
 * b is admissible when F(N, b) = sum c_j N^j b^(D-j) is a square mod q and
 * not both p | N and p | b; for a prime factor p with p | b the promise is
 * only that p does not divide N (the row is sieves0, whatever the curve:
 * when the leading coefficient is a non-square mod p those denominators
 * are excluded by the forbidden-divisor test, unless -F switched it off).
 * By the Chinese remainder theorem F is a square mod m exactly when it is
 * one mod every factor, so the row must be the AND of these.  The
 * wrap-around copies are checked too (the pattern is periodic in N).
 * Slow, for -DRP_VERIFY_MODULI builds only; aborts on the first mismatch. */
#include <stdlib.h>
static void verify_row(const ratpoints_bit_array *row, long m, long b,
                       long nf, const long *q, const long *p, const int *power,
                       ratpoints_args *args)
{ mpz_t *c = args->cof;
  long degree = args->degree, D = degree + (degree & 1);
  unsigned long cm[RP_MAX_FACTORS][D + 1], sq[RP_MAX_FACTORS][RP_MODWORDS];
  const unsigned long *w = (const unsigned long *)row;
  long N, k, i;

  for(i = 0; i < nf; i++)
  { for(k = 0; k <= degree; k++) { cm[i][k] = mpz_fdiv_ui(c[k], q[i]); }
    if(degree & 1) { cm[i][D] = 0UL; }
    for(k = 0; k < RP_MODWORDS; k++) { sq[i][k] = 0UL; }
    for(k = 0; k < q[i]; k++)
    { long s = (k*k) % q[i]; sq[i][s >> LONG_SHIFT] |= 1UL << (s & LONG_MASK); }
  }
  for(N = 0; N < (m + RATPOINTS_CHUNK-1)*RBA_LENGTH; N++)
  { long want = 1, bit;

    for(i = 0; i < nf && want; i++)
    { if(b % p[i] == 0 && !power[i]) { want = (N % p[i] != 0); }
      else
      { unsigned long F = 0UL, bpow = 1UL, Nq = (unsigned long)(N % q[i]);
        long j;

        for(j = D; j >= 0; j--)
        { F = (F*Nq + cm[i][j]*bpow) % (unsigned long)q[i];
          bpow = (bpow*(unsigned long)b) % (unsigned long)q[i];
        }
        want = (sq[i][F >> LONG_SHIFT] >> (F & LONG_MASK)) & 1UL;
        if(b % p[i] == 0 && N % p[i] == 0) { want = 0; }
      }
    }
    bit = (w[N >> LONG_SHIFT] >> (N & LONG_MASK)) & 1UL;
    if(bit != want)
    { fprintf(stderr, "RP_VERIFY_MODULI: modulus %ld, residue %ld, bit index %ld:"
              " row says %ld, F says %ld\n", m, b, N, bit, want);
      abort();
    }
  }
}
#endif

/* The row of a prime power m = p^e for the residue b of the denominator.
 * For a unit b the admissible bit indices x are those with f(x b^-1) a
 * square mod m, bit x b^-1 of fsq; for p | b they are the units x with
 * frev(b x^-1) a square, bit b x^-1 of gsq (see rp_power). */
ratpoints_bit_array *_ratpoints_sieve_init_power(void *se1, long b1, void *args1)
{ ratpoints_sieve_entry *se = se1;
  ratpoints_args *args = args1;
  const rp_power *pw = se->pw;
  long m = pw->m, p = pw->p, b = b1, x;
  unsigned long pat[RP_MODWORDS];
  ratpoints_bit_array *row;

  for(x = 0; x < RP_MODWORDS; x++) { pat[x] = 0UL; }
  if(b % p)
  { long binv = pw->inv[b], r = 0;

    for(x = 0; x < m; x++)
    { if((pw->fsq[r >> LONG_SHIFT] >> (r & LONG_MASK)) & 1UL)
      { pat[x >> LONG_SHIFT] |= 1UL << (x & LONG_MASK); }
      r += binv; if(r >= m) { r -= m; }
    }
  }
  else
  { for(x = 1; x < m; x++)
    { long t;

      if(x % p == 0) { continue; }
      t = (long)(((unsigned long)b*(unsigned long)pw->inv[x]) % (unsigned long)m);
      if((pw->gsq[t >> LONG_SHIFT] >> (t & LONG_MASK)) & 1UL)
      { pat[x >> LONG_SHIFT] |= 1UL << (x & LONG_MASK); }
    }
  }
  row = lay_out_pattern(pat, m, args);
#ifdef RP_VERIFY_MODULI
  { int power = 1; verify_row(row, m, b, 1, &m, &p, &power, args); }
#endif
  se->sieve[b] = row;
  return(row);
}

/* The row of a product of moduli for the residue b: the AND of the factors'
 * rows for b mod each factor, laid down in blocks of the factor's length
 * (the factors divide m).  A factor's row is built first when it is not
 * there yet. */
ratpoints_bit_array *_ratpoints_sieve_init_product(void *se1, long b1, void *args1)
{ ratpoints_sieve_entry *se = se1;
  ratpoints_args *args = args1;
  long m = se->p, b = b1, i, x, j;
  const ratpoints_bit_array *frow[RP_MAX_FACTORS] = {0}; /* nf >= 2 fills
                                                          * what is read */
  ratpoints_bit_array *row;

  for(i = 0; i < se->nf; i++)
  { ratpoints_sieve_entry *fe = se->factor[i];
    long bf = b % fe->p;

    frow[i] = fe->sieve[bf] ? fe->sieve[bf] : (*(fe->init))(fe, bf, args);
  }
  /* the factors' rows are read before the row is taken from the buffer, so
   * that a factor built just now does not land where this row goes */
  row = (ratpoints_bit_array *)args->ba_next;
  args->ba_next += (m + RATPOINTS_CHUNK-1)*sizeof(ratpoints_bit_array);
  { long q = se->factor[0]->p;

    for(x = 0; x < m; x += q)
    { for(j = 0; j < q; j++) { row[x + j] = frow[0][j]; } }
  }
  for(i = 1; i < se->nf; i++)
  { long q = se->factor[i]->p;

    for(x = 0; x < m; x += q)
    { for(j = 0; j < q; j++) { AND(row[x + j], frow[i][j]); } }
  }
  for(j = 0; j < RATPOINTS_CHUNK-1; j++) { row[m + j] = row[j]; }
#ifdef RP_VERIFY_MODULI
  { long q[RP_MAX_FACTORS], p[RP_MAX_FACTORS];
    int power[RP_MAX_FACTORS];

    for(i = 0; i < se->nf; i++)
    { q[i] = se->factor[i]->p;
      p[i] = se->factor[i]->pw ? se->factor[i]->pw->p : se->factor[i]->p;
      power[i] = (se->factor[i]->pw != NULL);
    }
    verify_row(row, m, b, se->nf, q, p, power, args);
  }
#endif
  se->sieve[b] = row;
  return(row);
}
