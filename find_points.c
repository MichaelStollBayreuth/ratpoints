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
 * find_points.c                                                       *
 *                                                                     *
 * Core program file for ratpoints                                     *
 *                                                                     *
 * Michael Stoll, September 21, 2009, January 7, 2022                  *
 * with changes by Bill Allombert, Dec 29, 2021                        *
 ***********************************************************************/

#include "rp-private.h"

#include "primes.h"
/* defines

   long prime[PRIMES1000]; */

#include "find_points.h"
/* defines

   static const int squares[RATPOINTS_NUM_PRIMES+1][RATPOINTS_MAX_PRIME];
     squares[n][x] = 1 if x is a square mod prime[n], 0 if not

   static const long offsets[RATPOINTS_NUM_PRIMES];
     offset[n] = (2*LONG_LENGTH)^(-1) mod prime[n]

   static const long inverses[RATPOINTS_NUM_PRIMES][RATPOINTS_MAX_PRIME];
     inverses[n][x] = x^(-1) mod prime[n] for x != 0 mod prime[n]

   ratpoints_bit_array sieves0[RATPOINTS_NUM_PRIMES][RATPOINTS_MAX_PRIME_EVEN]
     sieves0[n][x] has bit i set (0 <= x < prime[n])
       <==> x*LONG_LENGTH + i is not divisible by prime[n]
 */


#define MAX_DIVISORS 512
 /* Maximal length of array for squarefree divisors of leading coefficient */


extern ratpoints_init_fun sieve_init[RATPOINTS_NUM_PRIMES];

typedef struct { double r; ratpoints_sieve_entry *ssp; } entry;

typedef struct { int p; int val; int slope; } use_squares1_info;

typedef struct { long p;
                 unsigned long *start;
                 unsigned long *end;
                 unsigned long *curr; }
               forbidden_entry;

static const int squares16[16] = {1,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0};
 /* Says if a is a square mod 16, for a = 0..15 */

/**************************************************************************
 * Initialization and cleanup of ratpoints_args structure                 *
 **************************************************************************/

/* The following is needed to obtain the correct memory alignment
 * when using 256-bit or 512-bit "words".
 * (Added by Bill Allombert)
 */
void *pointer_align(void *xx, long m)
{
  /* return the smallest address >= xx that is divisible by m */
  unsigned long x = (unsigned long) xx;
  long r = x % m;
  if (r == 0) return xx;
  return (void *) (x + m - r);
}

/* NOTE: args->degree must be set */
void find_points_init(ratpoints_args *args)
{
  long work_len = 3 + (args->degree + 1);
  /* allocate space for work[] */
  mpz_t *work = malloc(work_len*sizeof(mpz_t));

#ifdef DEBUG
  printf("\nfind_points: initialize..."); fflush(NULL);
#endif

  /* and initialize the mpz_t's in it */
  { long i;
    for(i = 0; i < work_len; i++) mpz_init(work[i]);
  }

  /* insert in args */
  args->work = work;
  args->work_length = work_len;

  /* allocate space for se_buffer */
  args->se_buffer
    = (ratpoints_sieve_entry *) malloc(RATPOINTS_NUM_PRIMES
                                        * sizeof(ratpoints_sieve_entry));
  args->se_next = args->se_buffer;

  /* allocate space for ba_buffer */
  { long need = 0;
    long n;

    /* Figure out how much space is needed for the sieving information:
     * For each prime p, we need p arrays (one for each denominator mod p)
     * of length p + RATPOINTS_CHUNK-1.
     * (We add RATPOINTS_CHUNK-1 so that we can avoid a wrap-around in _ratpoints_sift0.)
     */
    for(n = 0; n < RATPOINTS_NUM_PRIMES; n++) { need += prime[n]*(prime[n] + RATPOINTS_CHUNK-1); }
    /* args->ba_buffer_na saves the start address of the block reserved by malloc,
     * so that it can be freed later.
     * Use  need+1  to have the necessary leeway for the alignment.
     */
    args->ba_buffer_na = malloc((need+1)*sizeof(ratpoints_bit_array));
    args->ba_buffer = pointer_align(args->ba_buffer_na, sizeof(ratpoints_bit_array));
    args->ba_next = args->ba_buffer;
  }

  /* allocate space for int_buffer */
  args->int_buffer
    = malloc(RATPOINTS_NUM_PRIMES*(RATPOINTS_MAX_PRIME+1)*sizeof(int));
  args->int_next = args->int_buffer;

  /* allocate sieve_list */
  args->sieve_list = malloc(RATPOINTS_NUM_PRIMES
                             * sizeof(ratpoints_sieve_entry*));

  /* allocate remaining data structures */
  args->den_info = malloc((PRIMES1000+2)*sizeof(use_squares1_info));
  args->divisors = malloc((MAX_DIVISORS+1)*sizeof(long));
  args->forb_ba = malloc((RATPOINTS_NUM_PRIMES + 1)*sizeof(forbidden_entry));
  args->forbidden = malloc((RATPOINTS_NUM_PRIMES + 1)*sizeof(long));

#ifdef DEBUG
  printf("done.\n"); fflush(NULL);
#endif
  return;
}

void find_points_clear(ratpoints_args *args)
{

#ifdef DEBUG
  printf("\nfind_points: clean up..."); fflush(NULL);
#endif

  /* clear mpz_t's in work[] */
  { long i;
    mpz_t *work = args->work;

    for(i = 0; i < args->work_length; i++) mpz_clear(work[i]);
  }

  /* free memory */
  free(args->work);
  free(args->se_buffer);
  free(args->ba_buffer_na);
  free(args->int_buffer);
  free(args->sieve_list);
  free(args->den_info);
  free(args->divisors);
  free(args->forb_ba);
  free(args->forbidden);

  /* clear pointer in args */
  args->work = NULL; args->work_length = 0;
  args->se_buffer = NULL; args->se_next = NULL;
  args->ba_buffer_na = NULL;
  args->ba_buffer = NULL; args->ba_next = NULL;
  args->int_buffer = NULL; args->int_next = NULL;
  args->sieve_list = NULL;
  args->den_info = NULL; args->divisors = NULL;
  args->forb_ba = NULL; args->forbidden = NULL;

#ifdef DEBUG
  printf("done.\n"); fflush(NULL);
#endif

  return;
}

/**************************************************************************
 * Helper function: valuation of gmp-integer at a prime                   *
 **************************************************************************/

#define VERY_BIG 1000

static long valuation(const mpz_t n, long p, long *r, mpz_t vvv)
{
  long v = 0;
  unsigned long rem;

  mpz_abs(vvv, n);
  if(mpz_cmp_ui(vvv, 0) == 0) { *r = 0; return(VERY_BIG); }
  rem = mpz_fdiv_q_ui(vvv, vvv, p);
  while(rem == 0)
  { v++;
    rem = mpz_fdiv_q_ui(vvv, vvv, p);
  }
  *r = rem;
  return(v);
}

/* Same for a long integer */
static long valuation1(long n, long p)
{
  long v = 0;
  unsigned long rem;
  unsigned long qn = abs(n);
  if(n == 0) { return(VERY_BIG); }
  rem = qn % p;
  while(rem == 0)
  { v++;
    qn = qn/p;
    rem = qn % p;
  }
  return(v);
}

/**************************************************************************
 * Try to avoid divisions                                                 *
 **************************************************************************/

static inline long mod(long a, long b)
{
  long b1 = b << 4; /* b1 = 16*b */

  if(a < -b1) { a %= b; if(a < 0) { a += b; } return(a); }
  if(a < 0) { a += b1; }
  else { if(a >= b1) { return(a % b); } }
  b1 >>= 1; /* b1 = 8*b */
  if(a >= b1) { a -= b1; }
  b1 >>= 1; /* b1 = 4*b */
  if(a >= b1) { a -= b1; }
  b1 >>= 1; /* b1 = 2*b */
  if(a >= b1) { a -= b1; }
  if(a >= b) { a -= b; }
  return(a);
}

/**************************************************************************
 * Helper function: Jacobi symbol                                         *
 **************************************************************************/

static inline int jacobi(long b, mpz_t tmp, const mpz_t lcf)
{ /* Jacobi symbol (leading coeff/b) */
  long f;

  /* avoid divisions as far as possible! */

  /* remove 2's from b */
  while((b & 1) == 0) b >>= 1;
  f = mpz_fdiv_r_ui(tmp, lcf, (unsigned long)b);
  if(f == 0) return(1);

  while(1)
  { long s = 1;
    long n = f;
    long m = b; /* m is odd, n is positive and < m */

    /* looking at (n/m) */
    while(!(n & 1))
    { if(m & 2) s = -s; /* change sign iff m = 3 or 5 mod 8 */
      if(m & 4) s = -s;
      n >>= 1;
    }
    while(1)
    { /* switch roles */
      if(n & m & 2) s = -s; /* change sign iff m, n = 3 mod 4 */
      /* now we are looking at (m/n) */
      while(m > n)
      { m -= n;
        do
        { if(n & 2) s = -s; /* change sign iff n = 3 or 5 mod 8 */
          if(n & 4) s = -s;
          m >>= 1;
        }
        while(!(m & 1));
      }
      if(m == n)
      { if(m == 1) return(s);
        /* otherwise, m is the gcd of f and b; remove it from b */
        b /= m; if(f >= b) f %= b;
        if(f == 0) return(1);
        break;
      }
      /* here m < n */
      /* switch roles */
      if(n & m & 2) s = -s; /* change sign iff m, n = 3 mod 4 */
      /* now we are looking at (n/m) */
      while(n > m)
      { n -= m;
        do
        { if(m & 2) s = -s; /* change sign iff m = 3 or 5 mod 8 */
          if(m & 4) s = -s;
          n >>= 1;
        }
        while(!(n & 1));
      }
      if(m == n)
      { if(m == 1) return(s);
        /* otherwise, m is the gcd of f and b; remove it from b */
        b /= m; if(f >= b) f %= b;
        if(f == 0) return(1);
        break;
      }
    }
  }
}

static inline int jacobi1(long b, const long lcf)
{ /* Jacobi symbol (leading coeff/b) */
  long f;
  int neg = 0;

  /* avoid divisions as far as possible! */

  /* remove 2's from b */
  while((b & 1) == 0) b >>= 1;
  f = lcf;
  if(f < 0) { f = -f; neg = 1; }
  if(b < 1UL<<(LONG_LENGTH - 4)) f = mod(f, b);
  if(f == 0) return(1);

  while(1)
  { long s = (neg && (b & 2)) ? -1 : 1;
    long n = f;
    long m = b; /* m is odd, n is positive */

    /* looking at (n/m) */
    while(!(n & 1))
    { if(m & 2) s = -s; /* change sign iff m = 3 or 5 mod 8 */
      if(m & 4) s = -s;
      n >>= 1;
    }
    while(1)
    { /* switch roles */
      if(n & m & 2) s = -s; /* change sign iff m, n = 3 mod 4 */
      /* now we are looking at (m/n) */
      while(m > n)
      { m -= n;
        do
        { if(n & 2) s = -s; /* change sign iff n = 3 or 5 mod 8 */
          if(n & 4) s = -s;
          m >>= 1;
        }
        while(!(m & 1));
      }
      if(m == n)
      { if(m == 1) return(s);
        /* otherwise, m is the gcd of f and b; remove it from b */
        b /= m; /* if(f >= b) f %= b; */
        if(f == 0) return(1);
        break;
      }
      /* here m < n */
      /* switch roles */
      if(n & m & 2) s = -s; /* change sign iff m, n = 3 mod 4 */
      /* now we are looking at (n/m) */
      while(n > m)
      { n -= m;
        do
        { if(m & 2) s = -s; /* change sign iff m = 3 or 5 mod 8 */
          if(m & 4) s = -s;
          n >>= 1;
        }
        while(!(n & 1));
      }
      if(m == n)
      { if(m == 1) return(s);
        /* otherwise, m is the gcd of f and b; remove it from b */
        b /= m; /* if(f >= b) f %= b; */
        if(f == 0) return(1);
        break;
      }
    }
  }
}

/************************************************************************
 * Set up information on possible denominators                          *
 * when polynomial is of odd degree with leading coefficient != +-1     *
 ************************************************************************/

static void setup_us1(ratpoints_args *args)
{
  mpz_t *work = args->work; /* abs. value of leading coeff. in work[0] */
  long count = 0;
  unsigned long i, v;
  unsigned long rem;

  /* typedef struct { int p; int val; int slope; } use_squares1_info; */
  use_squares1_info *den_info = (use_squares1_info *)args->den_info;
  long *divisors = (long *)args->divisors;

  /* find prime divisors of leading coefficient*/
  /* first p = 2 */

#ifdef DEBUG
  printf("\nsetup_us1: find v_2(lcf)..."); fflush(NULL);
#endif

  v = mpz_scan1(work[0], 0); /* find first 1-bit ==> 2-adic valuation */

#ifdef DEBUG
  printf(" = %ld\n", v); fflush(NULL);
#endif

  if(v > 0)
  { /* prime divisor found; divide it off */
    den_info[count].p = 2;
    mpz_fdiv_q_2exp(work[0], work[0], v); /* remove power of 2 */
    den_info[count].val = v;
    count++;
  }
  for(i = 0; i < PRIMES1000 && mpz_cmp_si(work[0], 1); i++)
  { int p = prime[i];

    if(mpz_cmp_si(work[0], p*p) < 0)
    { /* remaining part must be prime */

#ifdef DEBUG
      printf("\nsetup_us1: remaining factor");
      fflush(NULL);
#endif

      if(mpz_fits_slong_p(work[0]))
      {
        den_info[count].p = mpz_get_si(work[0]);
        den_info[count].val = 1;

#ifdef DEBUG
        printf(" = %d ==> fits into a long\n", den_info[count].p);
        fflush(NULL);
#endif

        count++;
        mpz_set_si(work[0], 1); /* divide it off */
      }

#ifdef DEBUG
      else
      { printf(" is too large\n"); fflush(NULL); }
#endif

      break;
    }
    else
    {

#ifdef DEBUG
      printf("\nsetup_us1: find v_%d(lcf)...", p); fflush(NULL);
#endif

      v = 0;
      rem = mpz_fdiv_q_ui(work[1], work[0], p);
      if(rem == 0)
      { /* prime divisor found; divide it off */
        den_info[count].p = p;
        while(rem == 0)
        { v++;
          mpz_set(work[0], work[1]);
          rem = mpz_fdiv_q_ui(work[1], work[0], p);
        }
        den_info[count].val = v;
        count++;
      }

#ifdef DEBUG
      printf(" = %ld\n", v); fflush(NULL);
#endif

  } }

#ifdef DEBUG
  printf("\nsetup_us1: %ld entries in den_info\n", count); fflush(NULL);
#endif

  den_info[count].p =  0; /* terminate  array */

  /* check if factorization is complete */
  if(mpz_cmp_si(work[0], 1) == 0)
  { /* set up array of squarefree divisors */
    long *div = &divisors[1];

    divisors[0] = 1;
    for(i = 0; i < count; i++)
    { /* multiply all divisors known so far by next prime */
      long *div0 = &divisors[0];
      long *div1 = div;

      for( ; div0 != div1; div0++)
      { long t = *div0 * (long)den_info[i].p;
        if(t <= args->b_high) { *div++ = t; }
        if(div >= &divisors[MAX_DIVISORS]) { break; }
      }
      if(div >= &divisors[MAX_DIVISORS]) { break; }
    }
    if(div < &divisors[MAX_DIVISORS])
    { *div = 0; /* terminate divisors array */

      /* note that we can use the information */
      args->flags |= RATPOINTS_USE_SQUARES1;

      /* set slopes in den_info */

#ifdef DEBUG
      printf("\nsetup_us1: compute slopes...\n"); fflush(NULL);
#endif

      for(i = 0; i < count; i++)
      { /* compute min{n : (d-k)*n > v_p(f_d) - v_p(f_k), k = 0,...,d-1} */
        int p = den_info[i].p;
        int v = den_info[i].val;
        int n = 1;
        int k;
        mpz_t *c = args->cof;
        long degree = args->degree;

        for(k = degree - 1; k >= 0; k--)
        { long dummy;
          int t = 1 + v - valuation(c[k], p, &dummy, work[0]);
          int m = CEIL(t, (degree - k));

          if(m > n) { n = m; }
        }

#ifdef DEBUG
        printf("  i = %ld (p = %d): slope = %d\n", i, p, n); fflush(NULL);
#endif

        den_info[i].slope = n;
      }
    }
    else
    {

#ifdef DEBUG
      printf("\nsetup_us1: too many divisors\n"); fflush(NULL);
#endif

    }
  }
  else
  {

#ifdef DEBUG
    printf("\nsetup_us1: no complete factorization\n"); fflush(NULL);
#endif

  }
  return;
}

/************************************************************************
 * Consider 2-adic information                                          *
 ************************************************************************/

static bit_selection get_2adic_info(ratpoints_args *args,
                                    unsigned long *den_bits,
                                    ratpoints_bit_array *num_bits)
{
  mpz_t *c = args->cof;
  long degree = args->degree;
  int is_f_square16[24];
  long cmp[degree+1]; /* The coefficients of f reduced modulo 16 */
  long npe = 0, npo = 0;
  bit_selection result;

#ifdef DEBUG
  printf("\nget_2adic_info: start...\n"); fflush(NULL);
#endif

  /* compute coefficients mod 16 */
  { long n;

    for(n = 0; n <= degree; n++) { cmp[n] = mpz_get_si(c[n]) & 0xf; }
  }

  /* determine if f(a) is a square mod 16, for a = 0..15 */
  { long a;

    for(a = 0 ; a < 16; a++)
    { unsigned long s = cmp[degree];
      long n;

      for(n = degree - 1 ; n >= 0 ; n--)
      { s *= a;
        s += cmp[n];
      }
      s &= 0xf;
      if((is_f_square16[a] = squares16[s]))
      { if(a & 1) { npo++; } else { npe++; } }
  } }

  /* even denominators:
     is_f_square16[16+k] says if f((2k+1)/2) is a square, k = 0..3
     is_f_square16[20+k] says if f((2k+1)/4) is a square, k = 0,1
     is_f_square16[22]   says if f(odd/8) is a square
     is_f_square16[23]   says if f(odd/2^n), n >= 4, can be a square */
  { long np1 = 0, np2 = 0, np3 = 0, np4 = 0;

    if(degree & 1)
    { long cf = 4*cmp[degree-1];
      long a;

      if(degree >= 2) { cf += 8*cmp[degree-2]; }
      for(a = 0; a < 4; a++)
      { /* Compute  2 c[d] k^d + 4 c[d-1] k^(d-1) + 8 c[d-2] k^(d-2), k = 2a+1.
           Note that k^d = k mod 8, k^(d-1) = 1 mod 8. */
        long k = 2*a+1;
        long s = (2*k*cmp[degree] + cf) & 0xf;

        if((is_f_square16[16+a] = squares16[s])) { np1++; }
      }
      if((is_f_square16[20] = squares16[(4*cmp[degree]) & 0xf])) { np2++; }
      if((is_f_square16[21] = squares16[(12*cmp[degree]) & 0xf])) { np2++; }
      if((is_f_square16[22] = squares16[(8*cmp[degree]) & 0xf])) { np3++; }
      is_f_square16[23] = 1; np4++;
    }
    else
    { long cf = (degree >= 2) ? 4*cmp[degree-2] : 0;
      long a;

      if(degree >= 3) { cf += 8*cmp[degree-3]; }
      for(a = 0; a < 4; a++)
      { /* compute c[d] k^d + 2 c[d-1] k^(d-1) + ... + 8 c[d-3] k^(d-3),
           k = 2a+1.
           Note that k^d = k^2 mod 16, k^(d-1) = k mod 8. */
        long k = 2*a+1;
        long s = ((cmp[degree]*k + 2*cmp[degree-1])*k + cf) & 0xf;

        if((is_f_square16[16+a] = squares16[s])) { np1++; }
      }
      if((is_f_square16[20] = squares16[(cmp[degree]+4*cmp[degree-1]) & 0xf]))
      { np2++; }
      if((is_f_square16[21] = squares16[(cmp[degree]+12*cmp[degree-1]) & 0xf]))
      {np2++; }
      if((is_f_square16[22] = squares16[(cmp[degree]+8*cmp[degree-1]) & 0xf]))
      {np3++; }
      if((is_f_square16[23] = squares16[cmp[degree]]))
      { np4++; }
    }

#ifdef DEBUG
    printf("\nis_f_square16 :\n[");
    { long a;

      for(a = 0; a < 23; a++) { printf("%d,", is_f_square16[a]); }
      printf("%d]\n", is_f_square16[23]);
    }
    fflush(NULL);
#endif

    /* set den_bits */
    { unsigned long db = 0;
      long i;

      if(npe + npo > 0) { db |= 0xaaaaUL; }
         /* odd denominators */
      if(np1 > 0)       { db |= 0x4444UL; }
         /* v_2(den) = 1 */
      if(np2 > 0)       { db |= 0x1010UL; }
         /* v_2(den) = 2 */
      if(np3 > 0)       { db |= 0x0100UL; }
         /* v_2(den) = 3 */
      if(np4 > 0)       { db |= 0x0001UL; }
         /* v_2(den) >= 4 */

      if(db == 0) { *den_bits = 0UL; return(num_none); }

      for(i = 16; i < LONG_LENGTH; i <<= 1) { db |= db << i; }

#ifdef DEBUG
      printf("\nden_bits: %*.*lx\n", WIDTH, WIDTH, db);
      fflush(NULL);
#endif

      *den_bits = db;
    }

    /* determine result */
    result = (npe == 0) ? ((npo == 0) ? num_none : num_odd)
                        : ((npo == 0) ? num_even : num_all);
  }

  { /* set up num_bits[16] */
    long b;

    /* odd denominators */
    switch(result)
    { case num_all:
        for(b = 1; b < 16; b += 2)
        { unsigned long work = 0;
          unsigned long bit = 1;
          long i;
          long invb = b; /* inverse of b mod 16 */

          if(b & 2) invb ^= 8;
          if(b & 4) invb ^= 8;
          for(i = 0; i < 16; i++)
          { if(is_f_square16[(invb*i) & 0xf]) { work |= bit; }
            bit <<= 1;
          }
          /* now repeat the 16 bits */
          for(i = 16; i < LONG_LENGTH; i <<= 1) { work |= work << i; }
          num_bits[b] = RBA(work);
        }
        break;

      case num_odd:
        for(b = 1; b < 16; b += 2)
        { unsigned long work = 0;
          unsigned long bit = 1;
          long i;
          long invb = b; /* inverse of b mod 16 */

          if(b & 2) invb ^= 8;
          if(b & 4) invb ^= 8;
          for(i = 1; i < 16; i += 2)
          { if(is_f_square16[(invb*i) & 0xf]) { work |= bit; }
            bit <<= 1;
          }
          /* now repeat the 8 bits */
          for(i = 8; i < LONG_LENGTH; i <<= 1) { work |= work << i; }
          num_bits[b] = RBA(work);
        }
        break;

      case num_even:
        for(b = 1; b < 16; b += 2)
        { unsigned long work = 0;
          unsigned long bit = 1;
          long i;
          long invb = b; /* inverse of b mod 16 */

          if(b & 2) invb ^= 8;
          if(b & 4) invb ^= 8;
          for(i = 0; i < 16; i += 2)
          { if(is_f_square16[(invb*i) & 0xf]) { work |= bit; }
            bit <<= 1;
          }
          /* now repeat the 8 bits */
          for(i = 8; i < LONG_LENGTH; i <<= 1) { work |= work << i; }
          num_bits[b] = RBA(work);
        }
        break;

      case num_none:
        for(b = 1; b < 16; b += 2) { num_bits[b] = zero; }
    }

    /* v_2(den) = 1 : only odd numerators */
    for(b = 1; b < 8; b += 2)
    { unsigned long work;
      unsigned long bit;
      long i;

      work = 0; bit = 1;
      for(i = 1; i < 16; i += 2)
      { if(is_f_square16[16 + (((b*i)>>1) & 0x3)]) { work |= bit; }
        bit <<= 1;
      }
      /* now repeat the 8 bits */
      for(i = 8; i < LONG_LENGTH; i <<= 1) { work |= work << i; }
      num_bits[2*b] = RBA(work);
    }

    /* v_2(den) = 2 : only odd numerators */
    for(b = 1; b < 4; b += 2)
    { unsigned long work = 0;
      unsigned long bit = 1;
      long i;

      work = 0; bit = 1;
      for(i = 1; i < 8; i += 2)
      { if(is_f_square16[20 + (((b*i)>>1) & 0x1)]) { work |= bit; }
        bit <<= 1;
      }
      /* now repeat the 4 bits */
      for(i = 4; i < LONG_LENGTH; i <<= 1) { work |= work << i; }
      num_bits[4*b] = RBA(work);
    }

    /* v_2(den) = 3, >= 4 : only odd numerators */
    num_bits[8] = (is_f_square16[22]) ? RBA(~(0UL)) : zero;
    num_bits[0] = (is_f_square16[23]) ? RBA(~(0UL)) : zero;
  }

#ifdef DEBUG
  printf("\nget_2adic_info: done.\n"); fflush(NULL);
#endif

  return(result);
}

/**************************************************************************
 * This is a comparison function needed for sorting in order to determine *
 * the `best' primes for sieving.                                         *
 **************************************************************************/

static int compare_entries(const void *a, const void *b)
{
  double diff = (((entry *)a)->r - ((entry *)b)->r);
  return (diff > 0) ? 1 : (diff < 0) ? -1 : 0;
}

/************************************************************************
 * Collect the sieving information                                      *
 ************************************************************************/

static long sieving_info(ratpoints_args *args,
                         int use_c_long, long *c_long,
                         ratpoints_sieve_entry **sieve_list)
/* This function either returns a prime p;
 * in this case, the curve has no points mod p, hence no rational points;
 * or else returns 0. */
{
  mpz_t *c = args->cof;
  long degree = args->degree;
  long fba = 0;
  long fdc = 0;
  long pn;
  long pnp = 0;
  entry prec[RATPOINTS_NUM_PRIMES];
    /* This array is used for sorting in order to
       determine the `best' sieving primes. */

  forbidden_entry *forb_ba = (forbidden_entry *)args->forb_ba;
  long *forbidden = (long *)args->forbidden;

  /* initialize sieve in se_buffer */
  for(pn = 0; pn < args->num_primes; pn++)
  { long coeffs_mod_p[degree+1];
           /* The coefficients of f reduced modulo p */
    long p = prime[pn];
    long n, a, np; /* np counts the x-coordinates that give points mod p */
    int *is_f_square = args->int_next;

    args->int_next += p + 1; /* need space for (p+1) int's */

#ifdef DEBUG
    printf("\nsieving_info: p = %ld\n", p);
    fflush(NULL);
#endif

    /* compute coefficients mod p */
    if(use_c_long)
    { for(n = 0; n <= degree; n++)
      { coeffs_mod_p[n] = mod(c_long[n], p); }
    }
    else
    { for(n = 0; n <= degree; n++)
      { coeffs_mod_p[n] = mpz_fdiv_r_ui(args->work[0], c[n], p); }
    }

    /* Determine the x-coords a mod p such that f(a) is a square mod p. */
    np = squares[pn][coeffs_mod_p[0]]; /* for a = 0, f(a) = constant term */
    is_f_square[0] = np;
    for(a = 1 ; a < p; a++)
    { unsigned long s = coeffs_mod_p[degree];
      /* try to avoid divisions (by p) */
      if((degree+1)*RATPOINTS_MAX_BITS_IN_PRIME <= LONG_LENGTH)
      { for(n = degree - 1 ; n >= 0 ; n--)
        { s *= a; s += coeffs_mod_p[n]; }
        /* here, s < p^(degree+1) <= max. long */
        s %= p;
      }
      else
      { for(n = degree - 1 ; n >= 0 ; n--)
        { s *= a; s += coeffs_mod_p[n];
          if(s+1 >= (1UL)<<(LONG_LENGTH - RATPOINTS_MAX_BITS_IN_PRIME))
          { s %= p; }
        }
        s %= p;
      }
      if((is_f_square[a] = squares[pn][s])) { np++; }
    }
    /* last entry says if there are points at infinity mod p */
    is_f_square[p] = (degree & 1) || squares[pn][coeffs_mod_p[degree]];

#ifdef DEBUG
    printf("\nis_f_square(p = %ld) : \n[", p);
    { long a;

      for(a = 0; a < p; a++) { printf("%d,", is_f_square[a]); }
      printf("%d]\n", is_f_square[p]);
    }
    fflush(NULL);
#endif

    /* check if there are no solutions mod p */
    if(np == 0 && !is_f_square[p])
    {
      return(p); /* if yes, return p --> no rational points */
    }

    /* Fill arrays with info for p */
    if(np < p)
    { /* only when there is some information */
      { double r = is_f_square[p] ? ((double)(np*(p-1) + p))/((double)(p*p))
                                  : (double)np/(double)p;

        prec[pnp].r = r;
      }

      /* set up sieve_entry :
         typedef struct
           { ratpoints_init_fun init; long p; int *is_f_square; int *inverses;
             long offset; (ratpoints_bit_array *)sieve[RATPOINTS_MAX_PRIME]; }
           ratpoints_sieve_entry;
       */
      { ratpoints_sieve_entry *se = (ratpoints_sieve_entry *)args->se_next;
        long i;

        args->se_next += sizeof(ratpoints_sieve_entry);
          /* one entry must be stored - note that se_next is of type void* */
        se->init = sieve_init[pn];
        se->p = p;
        se->is_f_square = is_f_square;
        se->inverses = &inverses[pn][0];
        se->offset = offsets[pn];
        se->sieve[0] = (ratpoints_bit_array *)&sieves0[pn][0];
        for(i = 1; i < p; i++) { se->sieve[i] = NULL; }

        prec[pnp].ssp = se;
      }
      pnp++;
    }

    if((args->flags & RATPOINTS_CHECK_DENOM)
         && fba + fdc < args->max_forbidden
         && !is_f_square[p])
    { /* record forbidden divisors of the denominator */
      if(coeffs_mod_p[degree] == 0)
      { /* leading coeff. divisible by p */
        long r;
        long v = valuation(c[degree], p, &r, args->work[0]);

        if((v & 1) || !squares[pn][r])
        { /* Can only get something when valuation is odd
             or when valuation is even and lcf is not a p-adic square.
             Compute smallest n such that if v(den) >= n, the leading
             term determines the valuation. Then we must have v(den) < n. */
          long n = 1;
          long k, pp;

          for(k = degree-1; k >= 0; k--)
          { if(coeffs_mod_p[k] == 0)
            { long dummy;
              long t = 1 + v - valuation(c[k], p, &dummy, args->work[0]);
              long m = CEIL(t, (degree-k));

              if(m > n) { n = m; }
          } }
          if(n == 1)
          { forb_ba[fba].p     = p;
            forb_ba[fba].start = &sieves0[pn][0];
            forb_ba[fba].end   = &sieves0[pn][p];
            forb_ba[fba].curr  = forb_ba[fba].start;
            fba++;
            pp = p;
          }
          else
          { for(pp = 1; n; n--) { pp *= p; } /* p^n */
            forbidden[fdc] = pp; fdc++;
          }

#ifdef DEBUG
          printf("\nexcluding denominators divisible by %ld\n", pp);
          fflush(NULL);
#endif

        }
      }
      else /* leading coefficient is a non-square mod p */
      { /* denominator divisible by p is excluded */
        forb_ba[fba].p     = p;
        forb_ba[fba].start = &sieves0[pn][0];
        forb_ba[fba].end   = &sieves0[pn][p];
        forb_ba[fba].curr  = forb_ba[fba].start;
        fba++;

#ifdef DEBUG
        printf("\nexcluding denominators divisible by %ld\n", p);
        fflush(NULL);
#endif

      }
    }

  } /* end for pn */

  /* update sp2 and sp1 if necessary */
  if(args->sp2 > pnp) { args->sp2 = pnp; }
  if(args->sp1 > args->sp2) { args->sp1 = args->sp2; }

  /* sort the array to get at the best primes */
  qsort(prec, pnp, sizeof(entry), compare_entries);

  /* put the sorted entries into sieve_list */
  { long n;

    for(n = 0; n < args->sp2; n++)
    { sieve_list[n] = prec[n].ssp; }
  }

  /* terminate array of forbidden divisors */
  if(args->flags & RATPOINTS_CHECK_DENOM)
  { long n;

    for(n = args->num_primes;
        fba + fdc < args->max_forbidden && n < RATPOINTS_NUM_PRIMES;
        n++)
    { long p = prime[n];

      if(p*p > args->b_high) break;
      if(mpz_kronecker_si(c[degree], p) == -1)
      { forb_ba[fba].p     = p;
        forb_ba[fba].start = &sieves0[n][0];
        forb_ba[fba].end   = &sieves0[n][p];
        forb_ba[fba].curr  = forb_ba[fba].start;
        fba++;

#ifdef DEBUG
        printf("\nexcluding denominators divisible by %ld\n", p);
        fflush(NULL);
#endif

      }
    }
    forb_ba[fba].p = 0;        /* terminating zero */
    forbidden[fdc] = 0;        /* terminating zero */
    args->max_forbidden = fba + fdc; /* note actual number */
  }

  if(fba + fdc == 0)
  { args->flags &= ~RATPOINTS_CHECK_DENOM; }

#ifdef DEBUG
  printf("\nsieving_info: done.\n"); fflush(NULL);
#endif

  return(0);
}

/**************************************************************************
 * The sieving procedure itself                                           *
 **************************************************************************/

static
long sift(long b, ratpoints_bit_array *survivors, ratpoints_args *args,
          bit_selection which_bits, ratpoints_bit_array bits16,
          ratpoints_sieve_entry **sieve_list, long *bp_list, int *quit,
          int process(long, long, const mpz_t, void*, int*), void *info)
{
  long total = 0;
  /* typedef struct { long p; long offset; ratpoints_bit_array *ptr; }
             sieve_spec; */
  sieve_spec ssp[args->sp2];
  int do_setup = 1;

#ifdef DEBUG
  printf("\nsift(b = %ld): start...\n", b); fflush(NULL);
#endif

  if((b & 1) == 0) { which_bits = num_odd; } /* even denominator */

  /* Note that b is new */
  args->flags |= RATPOINTS_COMPUTE_BC;

  { long k;
    long height = args->height;

    for(k = 0; k < args->num_inter; k++)
    { long low, high;
      /* For each of the positivity intervals,
       * determine relevant interval [low, high] of numerators. */
      { ratpoints_interval inter = args->domain[k];

        if(b*inter.low <= -height)
        { low = -height; }
        else
        { if(b*inter.low > height) { return(total); } /* remaining numerator intervals are empty */
          low = ceil(b*inter.low);
        }
        if(b*inter.up >= height)
        { high = height; }
        else
        { if(b*inter.up < -height) { continue; } /* this numerator interval is empty */
          high = floor(b*inter.up);
        }
      }

#ifdef DEBUG
      printf("\nsift: numerator interval [%ld, %ld]\n", low, high);
      fflush(NULL);
#endif

      if(do_setup)
      { /* set up the sieve information */
        long n;

        do_setup = 0; /* only do it once for every b */

#ifdef DEBUG
        printf("\nsift: set up sieve...\n");
        fflush(NULL);
#endif

        for(n = 0; n < args->sp2; n++)
        { ratpoints_sieve_entry *se = sieve_list[n];
          long p = se->p;
          long bp = bp_list[n];
          ratpoints_bit_array *sptr;

          if(which_bits != num_all) /* divide by 2 mod p */
          { bp = (bp & 1) ? (bp+p) >> 1 : bp >> 1; }
          sptr = se->sieve[bp];

          ssp[n].p = p;
          ssp[n].offset = (which_bits == num_odd) ? se->offset : 0;

#ifdef DEBUG
          printf("\np = %ld, bp = %ld, offset = %ld\n", p, bp, ssp[n].offset);
          fflush(NULL);
#endif
          /* copy if already initialized, else initialize */
          ssp[n].ptr = sptr ? sptr : (*(se->init))(se, bp, args);
          /* put a meaningful value in the start field */
          ssp[n].start = ssp[n].ptr;
          /* set the end field */
          ssp[n].end = ssp[n].ptr + p;

#ifdef DEBUG
          if(!sptr)
          { long a, c = 0;

            printf("\nsieve(%ld, %ld) [high numerators to the left]:", p, bp);
            for(a = p-1; a >= 0; a--, c++)
            { if((c & (0xff >> RBA_SHIFT)) == 0) { printf("\n"); }
              PRINT_RBA(ssp[n].ptr[a]);
            }
            printf("\n");
            fflush(NULL);
          }
#endif

        }
      }

      switch(which_bits)
      { case num_all: break;
        case num_none: break;
        case num_odd: low >>= 1; high--; high >>= 1; break;
        case num_even: low++; low >>= 1; high >>= 1; break;
      }

      /* now turn the bit interval into [low, high[ */
      high++;

      if(low < high)
      { long w_low, w_high;
        long w_low0, w_high0;
        long range = args->array_size;

        /* Now the range of longwords (= bit_arrays) */
        w_low = low >> RBA_SHIFT; /* FLOOR(low, RBA_LENGTH); */
        w_high = (high + (long)(RBA_LENGTH-1)) >> RBA_SHIFT;
                                 /* CEIL(high, RBA_LENGTH); */
        w_low0 = w_low;
        w_high0 = w_low0 + range;
        for( ; w_low0 < w_high; w_low0 = w_high0, w_high0 += range)
        { if(w_high0 > w_high)
          { w_high0 = w_high; range = w_high0 - w_low0; }
          /* initialize the bits */
          { long i;

            for(i = range; i; i--) { survivors[i-1] = bits16; }
          }
          /* boundary words */
          if(w_low0 == w_low)
          /* set lower bits of the first bit array to zero */
          { MASKL(survivors, low - RBA_LENGTH * w_low); }

#ifdef DEBUG
          printf("\nsurvivors[0] = ");
          PRINT_RBA(survivors[0]);
          printf("\n");
          fflush(NULL);
#endif

          if(w_high0 == w_high)
          /* set upper bits of the last bit array to zero */
          { MASKU(&survivors[range-1], RBA_LENGTH * w_high - high); }

#ifdef DEBUG
          printf("survivors[%ld] = ", range-1);
          PRINT_RBA(survivors[range-1]);
          printf("\n");
          fflush(NULL);
#endif

#if (RATPOINTS_CHUNK > 1)
          /* if necessary, increase range to be a multiple of RATPOINTS_CHUNK
           * and fill extra words with zeros. */
          while(range%RATPOINTS_CHUNK != 0)
          { survivors[range] = zero; range++, w_high0++; }
#endif

          total += _ratpoints_sift0(b, w_low0, w_high0, args, which_bits,
                                    survivors, &ssp[0], quit, process, info);
          if(*quit) return(total);
      } }
  } }

  return(total);
}

/**************************************************************************
 * Find points by looping over the denominators and sieving numerators    *
 **************************************************************************/

/*
typedef struct {mpz_t *cof; long degree; long height;
                ratpoints_interval *domain; long num_inter;
                long b_low; long b_high; long sp1; long sp2;
                long array_size;
                long sturm; long num_primes; long max_forbidden;
                unsigned int flags;
        ** from here: private data **
                mpz_t *work; long work_length;
                void *se_buffer; void *se_next;
                void *ba_buffer; void *ba_next;
                int *int_buffer; int *int_next;
                void *den_info; void *divisors;
                void *forb_ba; void *forbidden;
               }
        ratpoints_args;
*/

/* The first three entries of work[] are temporary mpz_t storage,
   the remaining ones constitue an array bc[] that
   will hold the coefficents of the polynomial,
   multiplied by powers of the denominator b */

long find_points_work(ratpoints_args *args,
                 int process(long, long, const mpz_t, void*, int*), void *info)
{
  long total = 0;       /* total counts the points */
  int quit = 0;
  mpz_t *c = args->cof;
  long degree = args->degree;
  long height = args->height;
  mpz_t *work = args->work;

  int point_at_infty = 0; /* indicates if there are points at infinity */
  int lcfsq = mpz_perfect_square_p(c[degree]);

  forbidden_entry *forb_ba = (forbidden_entry *)args->forb_ba;
  long *forbidden = (long *)args->forbidden;
    /* The forbidden divisors, a zero-terminated array.
       Used when degree is even and leading coefficient is not a square */

  use_squares1_info *den_info = (use_squares1_info *)args->den_info;
  long *divisors = (long *)args->divisors;
    /* These are used when degree is odd and leading coeff. is not +-1 */

  long c_long[degree+1]; /* Stores the coefficients as longs if possible */
  int use_c_long = 0;    /* Flag that says if c_long[] is set */

  ratpoints_sieve_entry **sieve_list = (ratpoints_sieve_entry **)args->sieve_list;
  bit_selection which_bits = num_all;
  unsigned long den_bits;
  ratpoints_bit_array num_bits[16];

  args->flags &= RATPOINTS_FLAGS_INPUT_MASK;
  args->flags |= RATPOINTS_CHECK_DENOM;

  /* initialize memory management */
  args->se_next = args->se_buffer;
  args->ba_next = args->ba_buffer;
  args->int_next = args->int_buffer;

#ifdef DEBUG
  printf("\nfind_points_work: start...\n"); fflush(NULL);
#endif

  if(c == NULL) return(RATPOINTS_BAD_ARGS);
  if(args->work_length < 3 + degree+1) return(RATPOINTS_WORK_LENGTH_TOO_SMALL);
  /* Eliminate leading zero coefficients */
  { long old_degree = degree;

    while(degree > 0 && mpz_cmp_si(c[degree], 0) == 0) { degree--; }
    args->degree = degree;
    if((degree+1)>>1 < (old_degree+1)>>1)
    { /* Polynomial not squarefree as a binary form of even degree */
      return(RATPOINTS_NON_SQUAREFREE);
  } }
  if(degree <= 0) return(RATPOINTS_BAD_ARGS);

#ifdef DEBUG
  printf("\nfind_points_work: sanity checks...\n"); fflush(NULL);
#endif

  /* Some sanity checks */
  if(args->num_inter < 0) { args->num_inter = 0; }

  if(args->num_primes < 0)
  { args->num_primes = RATPOINTS_DEFAULT_NUM_PRIMES; }
  if(args->sp1 < 0) { args->sp1 = RATPOINTS_DEFAULT_SP1; }
  if(args->sp2 < 0) { args->sp2 = RATPOINTS_DEFAULT_SP2; }

  if(args->num_primes > RATPOINTS_NUM_PRIMES)
  { args->num_primes = RATPOINTS_NUM_PRIMES; }
  if(args->sp2 > args->num_primes) { args->sp2 = args->num_primes; }
  if(args->sp1 > args->sp2) { args->sp1 = args->sp2; }

  if(height < 1) { return(RATPOINTS_BAD_ARGS); }
  if(args->b_low < 1) { args->b_low = 1; }
  if(args->b_high < 1) { args->b_high = height; }
  if(args->b_high > height) { args->b_high = height; }
  if(args->max_forbidden < 0)
  { args->max_forbidden = RATPOINTS_DEFAULT_MAX_FORBIDDEN; }
  if(args->max_forbidden > RATPOINTS_NUM_PRIMES)
  { args->max_forbidden = RATPOINTS_NUM_PRIMES; }
  if(args->array_size <= 0) { args->array_size = RATPOINTS_ARRAY_SIZE; }
  { long s = 2*CEIL(height, LONG_LENGTH);
    if(args->array_size > s) { args->array_size = s; }
  }
  /* make sure that array size is a multiple of RATPOINTS_CHUNK */
  args->array_size = CEIL(args->array_size, RATPOINTS_CHUNK)*RATPOINTS_CHUNK;
  if(args->sturm > (long)(LONG_LENGTH - 2))
  { args->sturm = (long)(LONG_LENGTH - 2); }

  /* Don't reverse if intervals are specified or limits for the denominator
     are given */
  if(args->num_inter > 0 || args->b_low > 1 || args->b_high < height)
  { args->flags |= RATPOINTS_NO_REVERSE; }

  if(args->flags & RATPOINTS_VERBOSE)
  { printf("\nfind_points:\n");
    printf("  degree:       %ld\n", args->degree);
    printf("  coefficients:");
    { long n;

      for(n = 0; n <= degree; n++)
      { printf(" "); mpz_out_str(NULL, 10, args->cof[n]); }
    }
    printf("\n");
    printf("  height bound: %ld\n", args->height);
    printf("  denominators from %ld to %ld\n", args->b_low, args->b_high);
    printf("  number of primes to consider:     %3ld\n", args->num_primes);
    printf("  number of primes for sieving:     %3ld\n", args->sp2);
    printf("  number of primes for first stage: %3ld\n", args->sp1);
    printf("  maximal number of `forbidden divisors': %ld\n",
           args->max_forbidden);
    if(args->sturm >= 0)
    { printf("  iterations for isolations of connected components: %ld\n",
             args->sturm);
    }
    else
    { printf("  no isolation of connected components to be done\n"); }
    if(args->flags & RATPOINTS_NO_CHECK)
    { printf("  do not verify the points\n"); }
    if(args->flags & RATPOINTS_NO_REVERSE)
    { printf("  do not reverse the polynomial\n"); }
    if(args->flags & RATPOINTS_NO_JACOBI)
    { printf("  do not perform Jacobi symbol test\n"); }
    printf("\n");
  }

#ifdef DEBUG
  printf("\nfind_points_work: check whether to reverse polynomial\n");
  fflush(NULL);
#endif

  /* Check if reversal of polynomial might be better:
    * case 1: degree is even, but trailing coefficient is zero
    * case 2: degree is even, leading coefficient is a square, but
              trailing coefficient is not
    * case 3: degree is odd, leading coefficient is not +-1,
              trailing coefficient is zero, coeff. of x is +-1
  */
  if(!((args->flags) & RATPOINTS_NO_REVERSE))
  { if(args->flags & RATPOINTS_VERBOSE)
    { printf("Check if polynomial should be reversed "
             "for better performance:\n");
    }
    if((degree & 1) == 0)
    { if(mpz_cmp_si(c[0], 0) == 0) /* case 1 */
      { long n;

        if(mpz_cmp_si(c[1], 0) == 0)
        { return(RATPOINTS_NON_SQUAREFREE); /* divisible by x^2 */ }
        args->flags |= RATPOINTS_REVERSED;
        for(n = 0; n < degree>>1; n++)
        { mpz_set(work[0], c[n]);
          mpz_set(c[n], c[degree-n]);
          mpz_set(c[degree-n], work[0]);
        }
        degree--; args->degree = degree;
        if(args->flags & RATPOINTS_VERBOSE)
        { printf("  even degree, zero constant term ==> reverse\n\n"); }
      }
      else
      { if(lcfsq && !mpz_perfect_square_p(c[0])) /* case 2 */
        { long n;

          args->flags |= RATPOINTS_REVERSED;
          for(n = 0; n < degree>>1; n++)
          { mpz_set(work[0], c[n]);
            mpz_set(c[n], c[degree-n]);
            mpz_set(c[degree-n], work[0]);
          }
          lcfsq = 0;
          if(args->flags & RATPOINTS_VERBOSE)
          { printf("  even degree, leading coefficient is a square, "
                   "constant term is not a square ==> reverse\n\n");
          }
      } }
    }
    else /* now degree is odd */
    { mpz_abs(work[0], c[degree]);
      mpz_abs(work[1], c[1]);
      if(mpz_cmp_si(work[0], 1) != 0
          && mpz_cmp_si(c[0], 0) == 0
          && mpz_cmp_si(work[1], 1) == 0) /* case 3*/
      { long n;

        args->flags |= RATPOINTS_REVERSED;
        for(n = 1; n < degree>>1; n++)
        { mpz_set(work[0], c[n]);
          mpz_set(c[n], c[degree+1-n]);
          mpz_set(c[degree+1-n], work[0]);
        }
        if(args->flags & RATPOINTS_VERBOSE)
        { printf("  odd degree, leading coefficient not +/-1, zero "
                 "constant term, coefficient of x is +/-1 ==> reverse\n\n");
        }
      }
  } }
  if(args->flags & RATPOINTS_VERBOSE)
  { if(!(args->flags & RATPOINTS_REVERSED))
    { printf("  criteria are not met ==> don't reverse\n\n"); }
  }

#ifdef DEBUG
  if(args->flags & RATPOINTS_REVERSED)
  { printf("\nfind_points_work: polynomial reversed.\n"); fflush(NULL); }
#endif

  /* Check is coefficients are small (i.e., fit into a long) */
  { long i;
    int flag = 1;

    for(i = 0; i <= degree; i++)
    { if(mpz_fits_slong_p(c[i])) { c_long[i] = mpz_get_si(c[i]); }
      else { flag = 0; break; }
    }
    use_c_long = flag;
  }

#ifdef DEBUG
  printf("\nfind_points_work: compute connected components\n"); fflush(NULL);
#endif

  /* Deal with the intervals */
  if(args->num_inter == 0)
  /* default interval (effectively ]-infty,infty[) if none is given */
  { if(args->domain == NULL)  return(RATPOINTS_BAD_ARGS);
    args->domain[0].low = -height; args->domain[0].up = height;
    args->num_inter = 1;
  }

  if(args->sturm >= 0)
  { long ret;

    if(args->flags & RATPOINTS_VERBOSE)
    { printf("Isolate the connected components:\n"); }
    ret = _ratpoints_compute_sturm(args);
    if(args->flags & RATPOINTS_VERBOSE)
    { if(ret < 0) { printf("  polynomial is not squarefree ==> stop\n\n"); }
      else
      if(ret == 0)
      { printf("  polynomial is always negative ==> no points\n\n"); }
      else
      { long n;

        printf("  can restrict to the following intervals:\n  ");
        for(n = 0; n < args->num_inter; n++)
        { printf("[%lf, %lf] ", args->domain[n].low, args->domain[n].up); }
        printf("\n\n");
      }
    }
    if(ret <= 0) /* not squarefree or no real points */
    { if(ret == 0) { return(0); }
      return(RATPOINTS_NON_SQUAREFREE);
    }
  }

  /* Point(s) at infinity? */
  if((degree & 1) || lcfsq)
  { args->flags &= ~RATPOINTS_CHECK_DENOM;
    point_at_infty = 1;
    if(args->flags & RATPOINTS_VERBOSE)
    { printf("There are points at infinity\n\n"); }
  }

  /* Can use only squares as denoms if degree is odd and poly is +-monic */
  if(degree & 1)
  { mpz_set(work[1], c[degree]);
    mpz_abs(work[0], work[1]);
    if(mpz_cmp_si(work[0], 1) == 0)
    { args->flags |= RATPOINTS_USE_SQUARES;
      if(args->flags & RATPOINTS_VERBOSE)
      { printf("Degree is odd, leading coefficient is +/-1\n");
        printf("  ==> can restrict to squares as denominators\n\n");
      }
    }
    else /* set up information on divisors of leading coefficient */
    { if(args->flags & RATPOINTS_VERBOSE)
      { printf("Degree is odd, leading coefficient is not +/-1\n");
        printf("  ==> can restrict denominators\n"
               "      to squares times certain "
               "divisors of the leading coefficient:\n");
      }
      setup_us1(args);
      if(args->flags & RATPOINTS_VERBOSE)
      { if(args->flags & RATPOINTS_USE_SQUARES1)
        { long n;

          printf("    divisors:");
          for(n = 0; divisors[n]; n++)
          { printf(" %ld", divisors[n]); }
          printf("\n\n");
        }
        else
        { printf("  no complete factorization obtained, or too many divisors\n"
                 "  ==> cannot use this feature\n\n");
        }
      }
    }
  }

  /* deal with f mod powers of 2 */
  if(args->flags & RATPOINTS_VERBOSE)
  { printf("Obtain information from the polynomial mod 16:\n"); }
  which_bits = get_2adic_info(args, &den_bits, &num_bits[0]);
  /* which_bits says whether to consider even and/or odd numerators
     when the denominator is odd.

     Bit k in den_bits is 0 if b congruent to k mod LONG_LENGTH need
     not be considered as a denominator.

     Bit k in num_bits[b] is 0 is numerators congruent to
     k (which_bits = den_all) / 2k (which_bits = den_even) /
     2k+1 (which_bits = den_odd)
     need not be considered for denominators congruent to b mod 16.
   */

#ifdef DEBUG
  { long i, c = 0;

    printf("\nusing %s numerators for odd denominators\n",
           (which_bits == num_none) ? "no"
            : (which_bits == num_even) ? "even"
            : (which_bits == num_odd) ? "odd"
            : "all");
    printf("\nden_bits: %*.*lx\n", WIDTH, WIDTH, den_bits);
    printf("\nnum_bits for b = 15, 14, ..., 0 mod 16 "
           "[high numerators to the left]:");
    for(i = 15; i >= 0; i--, c++)
    { if((c & (0xff >> LONG_SHIFT)) == 0) { printf("\n"); }
      printf(" %*.*lx", WIDTH, WIDTH, EXT0(num_bits[i]));
    }
    printf("\n\n");
    fflush(NULL);
  }
#else
  if(args->flags & RATPOINTS_VERBOSE)
  { printf("  use %s numerators for odd denominators\n\n",
           (which_bits == num_none) ? "no"
            : (which_bits == num_even) ? "even"
            : (which_bits == num_odd) ? "odd"
            : "all");
  }
#endif

  /* set up the sieve data structure */
  if(args->flags & RATPOINTS_VERBOSE)
  { printf("Find the points mod p for the first %ld odd primes p:\n",
           args->num_primes);
  }
  { long ret = sieving_info(args, use_c_long, &c_long[0], sieve_list);

    if(ret)
    {

#ifdef DEBUG
      printf("\nno points mod p = %ld ==> return(0)\n", ret);
#else
      if(args->flags & RATPOINTS_VERBOSE)
      { printf("  no points mod p = %ld ==> no rational points\n\n", ret); }
#endif

      return(0);
  } }

#ifdef DEBUG
  { long n;

    printf("\n%ld primes for first stage:\n", args->sp1);
    for(n = 0; n < args->sp1; n++)
    { printf(" %ld", sieve_list[n]->p); }
    printf("\n\n%ld primes for second stage:\n", args->sp2 - args->sp1);
    for( ; n < args->sp2; n++)
    { printf(" %ld", sieve_list[n]->p); }
    printf("\n");
    fflush(NULL);
  }
#else
  if(args->flags & RATPOINTS_VERBOSE)
  { long n;

    printf("  use %ld primes for first stage:\n   ", args->sp1);
    for(n = 0; n < args->sp1; n++)
    { printf(" %ld", sieve_list[n]->p); }
    printf("\n  use %ld primes for second stage:\n   ", args->sp2 - args->sp1);
    for( ; n < args->sp2; n++)
    { printf(" %ld", sieve_list[n]->p); }
    printf("\n\n");
  }
#endif

  /* deal with point(s) at infinity */
  if(point_at_infty)
  { long a = 1, b = 0;

#ifdef DEBUG
    printf("\nfind_points_work: points at infinity...\n"); fflush(NULL);
#else
    if(args->flags & RATPOINTS_VERBOSE)
    { printf("Points at infinity:\n"); }
#endif

    if(args->flags & RATPOINTS_REVERSED) { a = 0; b = 1; }

    if(args->flags & RATPOINTS_NO_CHECK)
    { mpz_set_si(work[0], 0);
      total += process(a, b, work[0], info, &quit);
    }
    else
    { if(degree & 1)
      { mpz_set_si(work[0], 0);
        total += process(a, b, work[0], info, &quit);
      }
      else
      { mpz_sqrt(work[0], c[degree]);
        total += process(a, b, work[0], info, &quit);
        if(!quit && !((args->flags) & RATPOINTS_NO_Y))
        { mpz_neg(work[0], work[0]);
          total += process(a, b, work[0], info, &quit);
        }
      }
    }

    if(quit)
    {
      return(total);
    }
    if(args->flags & RATPOINTS_VERBOSE) { printf("\n"); }
  }

#ifdef DEBUG
  printf("\nfind_points_work: start sieving...\n"); fflush(NULL);
#else
  if(args->flags & RATPOINTS_VERBOSE)
  { printf("Now start the sieving procedure...\n\n"); }
#endif

  /* now do the sieving */
  { ratpoints_bit_array *survivors;
    void *survivors_na;

#ifdef DEBUG
    printf("\nfind_points_work: allocating space for survivors...");
    fflush(NULL);
#endif

    /* allocate space for survivors array; make sure of correct alignment */
    survivors_na = malloc((args->array_size+1)*sizeof(ratpoints_bit_array));
    survivors = (ratpoints_bit_array *)
                pointer_align(survivors_na, sizeof(ratpoints_bit_array));
#ifdef DEBUG
    printf(" done\n");
    fflush(NULL);
#endif

    if(args->flags & (RATPOINTS_USE_SQUARES | RATPOINTS_USE_SQUARES1))
    { if(args->flags & RATPOINTS_USE_SQUARES)
      /* need only take squares as denoms */
      { long b, bb;
        long bp_list[args->sp2];
        long last_b = args->b_low;

#ifdef DEBUG
        printf("\n  using squares\n");
        fflush(NULL);
#endif

        { long n;

          for(n = 0; n < args->sp2; n++)
          { bp_list[n] = mod(args->b_low, sieve_list[n]->p); }
        }

        for(b = 1; bb = b*b, bb <= args->b_high; b++)
        { if(bb >= args->b_low)
          { ratpoints_bit_array bits = num_bits[bb & 0xf];

            if(TEST(bits))
            { long n;
              long d = bb - last_b;

              /* fill bp_list */
              for(n = 0; n < args->sp2; n++)
              { bp_list[n] = mod(bp_list[n] + d, sieve_list[n]->p); }
              last_b = bb;

              total += sift(bb, survivors, args, which_bits, bits,
                            sieve_list, &bp_list[0],
                            &quit, process, info);
              if(quit) { break; }
            }

#ifdef DEBUG
            else
            { printf("\nb = %ld: excluded mod 16\n", b);
              fflush(NULL);
            }
#endif

        } }
      }
      else /* args->flags & RATPOINTS_USE_SQUARES1 */
      { long *div = &divisors[0];
        long b, bb;
        long bp_list[args->sp2];

#ifdef DEBUG
        printf("\n  using squares times divisors of leading coefficient\n");
        fflush(NULL);
#endif

        for( ; *div; div++)
        { long last_b = *div;

#ifdef DEBUG
          printf("\n  divisor = %ld\n", *div);
          fflush(NULL);
#endif

          { long n;

            for(n = 0; n < args->sp2; n++)
            { bp_list[n] = mod(*div, sieve_list[n]->p); }
          }

          for(b = 1; bb = (*div)*b*b, bb <= args->b_high; b++)
          { if(bb >= args->b_low)
            { int flag = 1;
              ratpoints_bit_array bits = num_bits[bb & 0xf];

              if(EXT0(bits))
              { long i;
                long n;
                long d = bb - last_b;

                /* fill bp_list */
                for(n = 0; n < args->sp2; n++)
                { bp_list[n] = mod(bp_list[n] + d, sieve_list[n]->p); }
                last_b = bb;

                for(i = 0; den_info[i].p; i++)
                { int v = valuation1(bb, den_info[i].p);
                  if((v >= den_info[i].slope)
                       && ((v + (den_info[i].val)) & 1))
                  { flag = 0; break; }
                }
                if(flag)
                {
                  total += sift(bb, survivors, args, which_bits, bits,
                                sieve_list, &bp_list[0],
                                &quit, process, info);
                  if(quit) { break; }
                }
              }

#ifdef DEBUG
              else
              { printf("\nb = %ld: excluded mod 16\n", b);
                fflush(NULL);
              }
#endif

          } }
        if(quit) { break; }
        }
    } }
    else
    { if(args->flags & RATPOINTS_CHECK_DENOM)
      { long *forb;
        long b;
        long bp_list[args->sp2];
        long last_b = args->b_low;
        unsigned long b_bits;

#ifdef DEBUG
        printf("\n  taking account of forbidden divisors of the denominator\n");
        fflush(NULL);
#endif

        { long n;

          for(n = 0; n < args->sp2; n++)
          { bp_list[n] = mod(args->b_low, sieve_list[n]->p); }
        }

#ifdef DEBUG
        printf("\n  bp_list initialized\n");
        fflush(NULL);
#endif

        { forbidden_entry *fba = &forb_ba[0];
          long b_low = args->b_low;
          long w_low = (b_low-1) >> LONG_SHIFT;

          b_bits = den_bits;
          while(fba->p)
          { fba->curr = fba->start + mod(w_low, fba->p);
            b_bits &= *(fba->curr);
            fba++;
          }
          b_bits >>= (b_low-1) & LONG_MASK;
        }

#ifdef DEBUG
          printf("\n  initial b_bits = %*.*lx\n", WIDTH, WIDTH, b_bits);
          fflush(NULL);
#endif

        for(b = args->b_low; b <= args->b_high; b++)
        { ratpoints_bit_array bits = num_bits[b & 0xf];

          if((b & LONG_MASK) == 0)
          { /* next b_bits */
            forbidden_entry *fba = &forb_ba[0];

            b_bits = den_bits;
            while(fba->p)
            { fba->curr++;
              if(fba->curr == fba->end) { fba->curr = fba->start; }
              b_bits &= *(fba->curr);
              fba++;
            }
          }
          else
          { b_bits >>= 1; }

#ifdef DEBUG
          printf("\n  b_bits = %*.*lx\n", WIDTH, WIDTH, b_bits);
          fflush(NULL);
#endif

          if((b_bits & 1) && EXT0(bits))
          { /* check if denominator is excluded */
            for(forb = &forbidden[0] ; *forb && (b % (*forb)); forb++) {};

#ifdef DEBUG
            if(*forb)
            { printf("\nb = %ld: excluded mod %ld\n", b, *forb);
              fflush(NULL);
            }
#endif

            if(*forb == 0
                && ((args->flags & RATPOINTS_NO_JACOBI)
                      || (use_c_long
                           ? jacobi1(b, c_long[degree])
                           : jacobi(b, work[0], c[degree])) == 1))
            { long n;
              long d = b - last_b;

              /* fill bp_list */
              for(n = 0; n < args->sp2; n++)
              { long bp = bp_list[n] + d;
                long p = sieve_list[n]->p;

                while(bp >= p) { bp -= p; }
                bp_list[n] = bp;
              }
              last_b = b;

              total += sift(b, survivors, args, which_bits, bits,
                            sieve_list, &bp_list[0],
                            &quit, process, info);
              if(quit) { break; }
            }

#ifdef DEBUG
            else
            { if(*forb == 0)
              { printf("\nb = %ld: excluded by Jacobi symbol\n", b);
                fflush(NULL);
            } }
#endif

          }
        }
      } /* if(args->flags & RATPOINTS_CHECK_DENOM) */
      else
      { long b;
        long bp_list[args->sp2];
        long last_b = args->b_low;

        { long n;

          for(n = 0; n < args->sp2; n++)
          { bp_list[n] = mod(args->b_low, sieve_list[n]->p); }
        }

        for(b = args->b_low; b <= args->b_high; b++)
        { ratpoints_bit_array bits = num_bits[b & 0xf];

          if(EXT0(bits))
          { long n;
            long d = b - last_b;

            /* fill bp_list */
            for(n = 0; n < args->sp2; n++)
            { long bp = bp_list[n] + d;
              long p = sieve_list[n]->p;

              while(bp >= p) { bp -= p; }
              bp_list[n] = bp;
            }
            last_b = b;

            total += sift(b, survivors, args, which_bits, bits,
                          sieve_list, &bp_list[0],
                          &quit, process, info);
            if(quit) { break; }
          }

#ifdef DEBUG
          else
          { printf("\nb = %ld: excluded mod 16\n", b);
            fflush(NULL);
          }
#endif

      } }
    }
    /* de-allocate memory */
    free(survivors_na);
  }

#ifdef DEBUG
  printf("\nfind_points_work: done. total = %ld.\n", total); fflush(NULL);
#endif

  return(total);
}

/**************************************************************************
 * The wrapper function, doing init, work, and clear                      *
 **************************************************************************/

long find_points(ratpoints_args *args,
                 int process(long, long, const mpz_t, void*, int*), void *info)
{
  long result;

  /* first initialize */
  find_points_init(args);

  /* then do the work */
  result = find_points_work(args, process, info);

  /* now clean up */
  find_points_clear(args);

  /* and return the result */
  return(result);
}


/**************************************************************************
 * Check a `survivor' of the sieve if it really gives a point.            *
 * This function is called by _ratpoints_sift0(), see sift.c .            *
 **************************************************************************/

long _ratpoints_check_point(long a, long b, ratpoints_args *args, int *quit,
                 int process(long, long, const mpz_t, void*, int*), void *info)
{
  mpz_t *c = args->cof;
  long degree = args->degree;
  int reverse = args->flags & RATPOINTS_REVERSED;
  long total = 0;
  mpz_t *work = args->work;
  mpz_t *bc = &work[3];

  if(!(args->flags & RATPOINTS_NO_CHECK))
  { long k;

    /* Compute F(a, b), where F is the homogenized version of f
       of smallest possible even degree  */
    if(args->flags & RATPOINTS_COMPUTE_BC)
    { /* compute entries bc[k] = c[k] * b^(degree-k), k < degree */

#ifdef DEBUG
      printf("\ncheck_point: compute bc[] (b = %ld)\n", b);
      fflush(NULL);
#endif

      mpz_set_si(work[0], 1); /* work[0] contains the successive powers of b */
      for(k = degree-1; k >= 0; k--)
      { mpz_mul_ui(work[0], work[0], b);
        mpz_mul(bc[k], c[k], work[0]);
      }
      /* note that bc[] has been computed for the current b */
      args->flags &= ~RATPOINTS_COMPUTE_BC;
    }

#ifdef DEBUG
    printf("check_point: computing f(a = %ld, b = %ld)\n", a, b);
    fflush(NULL);
#endif

    mpz_set(work[2], c[degree]); /* use work[2] to accumulate the result */
    for(k = degree-1; k >= 0; k--)
    { mpz_mul_si(work[2], work[2], a);
      mpz_add(work[2], work[2], bc[k]);
    }
    /* if degree is odd, need to multiply again by b
     * to get value of binary form of even degree */
    if(degree & 1) mpz_mul_ui(work[2], work[2], b);
    /* check if f(x,z) is a square; if so, process the point(s) */
    if(mpz_cmp_si(work[2], 0) >= 0)
    { mpz_sqrtrem(work[0], work[1], work[2]);
      /* work[0] = isqrt(work[2]), work[1] = remainder,
       * so the y-coordinate is in work[0] */
      if(mpz_cmp_si(work[1], 0) == 0)
      {

#ifdef DEBUG
        printf("check_point: found point (a = %ld, b = %ld)\n", a, b);
        fflush(NULL);
#endif

        if(reverse)
        { if(a >= 0) { total += process(b, a, work[0], info, quit); }
          else { total += process(-b, -a, work[0], info, quit); }
        }
        else total += process(a, b, work[0], info, quit);
        /* process opposite point if necessary */
        if(!*quit && mpz_cmp_si(work[0], 0) != 0
                  && !((args->flags) & RATPOINTS_NO_Y))
        { mpz_neg(work[0], work[0]);
          if(reverse)
          { if(a >= 0) { total += process(b, a, work[0], info, quit); }
            else { total += process(-b, -a, work[0], info, quit); }
          }
          else { total += process(a, b, work[0], info, quit); }
        }
    } }
  } /* if(!no_check) */
  else /* arg->flags & RATPOINTS_NO_CHECK : no computation */
  { mpz_set_si(work[0], 0);
    if(reverse)
    { if(a >= 0) { total += process(b, a, work[0], info, quit); }
      else { total += process(-b, -a, work[0], info, quit); }
    }
    else { total += process(a, b, work[0], info, quit); }
  }
  return(total);
}
