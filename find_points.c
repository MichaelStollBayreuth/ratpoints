/***********************************************************************
 * ratpoints-2.2.3                                                     *
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
 * find_points.c                                                       *
 *                                                                     *
 * Core program file for ratpoints                                     *
 *                                                                     *
 * Michael Stoll, Sep 21, 2009; Jan 7, 2022; Sep 6, 2026               *
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
     offset[n] = (2*RBA_LENGTH)^(-1) mod prime[n]

   static const long inverses[RATPOINTS_NUM_PRIMES][RATPOINTS_MAX_PRIME];
     inverses[n][x] = x^(-1) mod prime[n] for x != 0 mod prime[n]

   unsigned long sieves0[RATPOINTS_NUM_PRIMES]
                        [RBA_PACK*(RATPOINTS_MAX_PRIME_EVEN + RATPOINTS_CHUNK-1)];
     the sieving information for denominator b == 0 mod prime[n]:
     sieves0[n][x] has bit i set (0 <= x < prime[n])
       <==> x*LONG_LENGTH + i is not divisible by prime[n],
     and this pattern of prime[n] words is then repeated, so that the array
     can be read as prime[n] + RATPOINTS_CHUNK-1 bit-arrays
     (the last RATPOINTS_CHUNK-1 of them being the wrap-around copies that
      _ratpoints_sift0 relies on in phase 1; compare init.c).
     The array carries an alignment attribute, since it is accessed through
     pointers of type  ratpoints_bit_array * ; see gen_find_points_h.c.
 */

/* Development instrumentation, see the head of sift.c .  The exact check
 * is timed there as a whole; this brackets the part of it that is done once
 * per denominator, so that the two can be told apart. */
#ifdef RP_PHASE_TIMING
#include <x86intrin.h>
extern unsigned long long _rp_bc_cycles, _rp_bc_calls;
# define RP_BC_TIC(t) unsigned long long t = __rdtsc()
# define RP_BC_TOC(t) do { _rp_bc_cycles += __rdtsc() - (t); _rp_bc_calls++; } \
                      while(0)
/* and the loop that steps b modulo each sieving prime, which is per
 * denominator and per prime, so that a third stage using further primes
 * would pay it whether or not a survivor turns up */
extern unsigned long long _rp_bp_cycles, _rp_bp_dens, _rp_bp_steps;
extern unsigned long long _rp_arrays_swept;
/* and building one sieve table, which is the fixed cost a prime has to earn
 * back over the run; see run_shape */
extern unsigned long long _rp_init_cycles, _rp_init_calls, _rp_init_rows;
extern unsigned long long _rp_setup_cycles, _rp_setup_dens;
# define RP_SETUP_TIC(t) unsigned long long t = __rdtsc()
# define RP_SETUP_TOC(t) do { _rp_setup_cycles += __rdtsc() - (t); \
                              _rp_setup_dens++; } while(0)
# define RP_INIT_TIC(t) unsigned long long t = __rdtsc()
# define RP_INIT_TOC(t, n) do { _rp_init_cycles += __rdtsc() - (t); \
                                _rp_init_calls++; _rp_init_rows += (n); } \
                           while(0)
# define RP_BP_TIC(t) unsigned long long t = __rdtsc()
# define RP_BP_TOC(t, n) do { _rp_bp_cycles += __rdtsc() - (t); _rp_bp_dens++; \
                              _rp_bp_steps += (n); } while(0)
/* and the whole of sift(), so that what is in none of the regions above can
 * be seen */
extern unsigned long long _rp_sift_cycles, _rp_sift_calls;
# define RP_SIFT_TIC(t) unsigned long long t = __rdtsc()
# define RP_SIFT_TOC(t) do { _rp_sift_cycles += __rdtsc() - (t); \
                             _rp_sift_calls++; } while(0)
#else
# define RP_BC_TIC(t)
# define RP_BC_TOC(t)
# define RP_BP_TIC(t)
# define RP_SIFT_TIC(t)
# define RP_SIFT_TOC(t)
# define RP_BP_TOC(t, n)
# define RP_INIT_TIC(t)
# define RP_INIT_TOC(t, n)
# define RP_SETUP_TIC(t)
# define RP_SETUP_TOC(t)
#endif


#define MAX_DIVISORS 512
 /* Maximal length of array for squarefree divisors of leading coefficient */


extern ratpoints_init_fun sieve_init[RATPOINTS_NUM_PRIMES];

typedef struct { double r; double key; ratpoints_sieve_entry *ssp; } entry;

typedef struct { int p; int val; int slope; } use_squares1_info;

typedef struct { long p;
                 unsigned long *start;
                 unsigned long *end;
                 unsigned long *curr; }
               forbidden_entry;
  /* a prime p no denominator may be divisible by; tested with a bit array */

typedef struct { long p; unsigned long mask; } forbidden_val;
  /* a prime p and the set of valuations v_p(b) no denominator b may have:
     bit m of mask is set <==> v_p(b) = m is excluded.  Tested by division.
     See forbidden_valuations() for where these come from. */

static const int squares16[16] = {1,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0};
 /* Says if a is a square mod 16, for a = 0..15 */

/**************************************************************************
 * Initialization and cleanup of ratpoints_args structure                 *
 **************************************************************************/

/* The following is needed to obtain the correct memory alignment
 * when using 256-bit or 512-bit "words": malloc() only guarantees 16 bytes.
 * The callers allocate one bit-array more than they need, so that there is
 * always enough room to move the start address up.
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
/* Reserve the space for the sieving information for the first n primes.
 * For each prime p we may need p arrays (one for each denominator mod p) of
 * length p + RATPOINTS_CHUNK-1, the CHUNK-1 so that _ratpoints_sift0 can
 * avoid a wrap-around.  The sum grows with the square of the largest prime,
 * which is why n matters: taking it to be RATPOINTS_NUM_PRIMES asks for 5 MB
 * when primes go up to 127, but 1.7 GB when they go up to 1021, whether or
 * not the large primes are ever looked at.
 * args->ba_buffer_na keeps the address malloc returned, so that it can be
 * freed later; the +1 leaves the leeway needed for the alignment.
 * args->ba_buffer_primes records what the block is good for.
 */
static void alloc_ba_buffer(ratpoints_args *args, long n)
{ long need = 0;
  long i;

  for(i = 0; i < n; i++) { need += prime[i]*(prime[i] + RATPOINTS_CHUNK-1); }
  args->ba_buffer_na = malloc((need+1)*sizeof(ratpoints_bit_array));
  args->ba_buffer = pointer_align(args->ba_buffer_na, sizeof(ratpoints_bit_array));
  args->ba_next = args->ba_buffer;
  args->ba_buffer_primes = n;
}

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

  /* allocate space for ba_buffer, for the first RATPOINTS_DEFAULT_NUM_PRIMES
   * primes; find_points_work enlarges it should a caller ask for more.  It
   * cannot be sized from args->num_primes here, because the documented way of
   * using the library sets that field between find_points_init and
   * find_points_work, not before.
   */
  alloc_ba_buffer(args, (RATPOINTS_DEFAULT_NUM_PRIMES < RATPOINTS_NUM_PRIMES)
                          ? RATPOINTS_DEFAULT_NUM_PRIMES
                          : RATPOINTS_NUM_PRIMES);

  /* allocate space for int_buffer */
  args->int_buffer
    = malloc(RATPOINTS_NUM_PRIMES*(RATPOINTS_MAX_PRIME+1)*sizeof(int));
  args->int_next = args->int_buffer;

  /* allocate sieve_list */
  args->sieve_list = malloc(RATPOINTS_NUM_PRIMES
                             * sizeof(ratpoints_sieve_entry*));

  /* and the third stage's working copy of what it needs per denominator.
   * It lives here rather than on sift()'s stack because that function is
   * entered once per denominator, and enlarging its frame by this much was
   * measured to cost several per cent all by itself. */
  args->stage3_list = malloc(RATPOINTS_NUM_PRIMES * sizeof(check_spec));

  /* the reciprocals _ratpoints_sift0 reduces word numbers with.  They belong
   * to the primes, not to the denominators, so they are filled in once per
   * curve; and they are kept out of sieve_spec because that structure is read
   * in the innermost loop of the first phase, where its size tells. */
  args->magics = malloc(RATPOINTS_NUM_PRIMES * sizeof(unsigned long));

  /* allocate remaining data structures */
  args->den_info = malloc((PRIMES1000+2)*sizeof(use_squares1_info));
  args->divisors = malloc((MAX_DIVISORS+1)*sizeof(long));
  args->forb_ba = malloc((RATPOINTS_NUM_PRIMES + 1)*sizeof(forbidden_entry));
  args->forbidden = malloc((RATPOINTS_NUM_PRIMES + 1)*sizeof(forbidden_val));

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
  free(args->stage3_list);
  free(args->magics);
  free(args->den_info);
  free(args->divisors);
  free(args->forb_ba);
  free(args->forbidden);

  /* clear pointer in args */
  args->work = NULL; args->work_length = 0;
  args->se_buffer = NULL; args->se_next = NULL;
  args->ba_buffer_na = NULL; args->ba_buffer_primes = 0;
  args->ba_buffer = NULL; args->ba_next = NULL;
  args->int_buffer = NULL; args->int_next = NULL;
  args->sieve_list = NULL; args->stage3_list = NULL;
  args->magics = NULL;
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
  unsigned long qn = labs(n);
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
  if(b < 1UL<<(LONG_LENGTH - 5)) f = mod(f, b); /* mod() forms 16*b */
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
      { long p = (long)den_info[i].p;
        /* the product may not fit when b_high is large, so compare first */
        if(*div0 <= args->b_high / p) { *div++ = *div0 * p; }
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

      if(db == 0)
      { /* No residue class of the denominator admits any numerator.  Return
         * early -- but fill num_bits[] first: the caller reads all sixteen
         * entries (bits_per_word, run_shape, the test on every denominator),
         * and the denominator loop that runs when RATPOINTS_CHECK_DENOM is
         * off looks at nothing else. */
        long i;

        for(i = 0; i < 16; i++) { num_bits[i] = zero; }
        *den_bits = 0UL;
        return(num_none);
      }

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

/* Primes are ranked by what they say per unit of what they cost, not by
 * what they say alone; key is set by prime_key() below.  With the cost of a
 * table switched off the key is monotone in r and this is the old order. */
static int compare_entries(const void *a, const void *b)
{
  double diff = (((entry *)a)->key - ((entry *)b)->key);
  return (diff > 0) ? 1 : (diff < 0) ? -1 : 0;
}

/* Beyond the second phase a prime builds no table, so all that separates
 * two of them is what they say: there the order is by density alone. */
static int compare_by_r(const void *a, const void *b)
{
  double diff = (((entry *)a)->r - ((entry *)b)->r);
  return (diff > 0) ? 1 : (diff < 0) ? -1 : 0;
}

/* What one more prime costs the sieve, per numerator word, in units of what
 * a first-phase prime costs there.
 *
 * per_word is the part that is paid for every word (or for every surviving
 * bit array, which comes to the same thing once multiplied by the survival
 * rate): 1 in the first phase, COST_PHASE2*rate in the second, and the
 * third stage's own cost per survivor in the third.  The other two terms are
 * paid once and spread over the run: the sieve table, which the third stage
 * does not build, and the step of bp_list, which every stage pays.
 */
static double prime_cost(long p, double per_word, int tabled,
                         double cost_table, double u_words, double n_denoms)
{ double cost = per_word + RATPOINTS_COST_BP*n_denoms/u_words;

  if(tabled)
  { /* a prime of the first two phases has a sieve_spec filled in for it as
     * well, once per denominator */
    cost += RATPOINTS_COST_SETUP*n_denoms/u_words;
    if(cost_table > 0.0)
    { double builds = (n_denoms < (double)p) ? n_denoms : (double)p;

      cost += cost_table*(double)p*builds/u_words;
    }
  }
  return(cost);
}

/* The rank of a prime: what it costs divided by what it says.  A prime
 * multiplies the survival rate by r, so what it says is -log(r), and the
 * best set of primes for a given total cost is found by taking them in
 * increasing order of this ratio. */
static double prime_key(double r, long p, double per_word, int tabled,
                        double cost_table, double u_words, double n_denoms)
{ double info = -log(r);

  if(info <= 0.0) { return(1.0e300); }
  return(prime_cost(p, per_word, tabled, cost_table, u_words, n_denoms)/info);
}

/* What one exact check costs for this curve, in the units of
 * RATPOINTS_CHECK_REFERENCE.  The two third-stage constants are fractions of
 * one check, so they have to be divided by this; see the comment on
 * RATPOINTS_CHECK_STEP in ratpoints.h for what the formula is counting.
 *
 * Everything it needs is known before the first prime is looked at: the
 * degree, the largest coefficient and the height bound fix the size of
 * F(a,b) = c[degree]*a^degree + ... + c[0]*b^degree, and with it the size of
 * every number the check touches.  The height bound is used for both a and
 * b, which is what they are bounded by; a denominator range narrower than
 * that makes the estimate a little high, and by less than the rounding to
 * whole limbs does.
 */
static double check_cost(const ratpoints_args *args)
{ mpz_t *c = args->cof;
  long degree = args->degree;
  long k;
  double hbits = log((double)args->height + 1.0)/log(2.0);
  double cbits = 1.0;
  double fbits, limbs, mid, root, cost;

  if(args->check_cost > 0.0) { return(args->check_cost); }

  for(k = 0; k <= degree; k++)
  { if(mpz_sgn(c[k]) != 0)
    { double b = (double)mpz_sizeinbase(c[k], 2);

      if(b > cbits) { cbits = b; }
  } }

  fbits = cbits + (double)degree*hbits;
  limbs = fbits/(double)LONG_LENGTH;          /* the size of F, in limbs */
  mid = 0.5*(cbits + fbits)/(double)LONG_LENGTH; /* the mean over the loop */
  root = ceil(0.5*limbs);                     /* limbs of the square root */
  if(root < 1.0) { root = 1.0; }

  cost = (double)degree*(RATPOINTS_CHECK_STEP + RATPOINTS_CHECK_LIMB*mid)
          + RATPOINTS_CHECK_CALL + RATPOINTS_CHECK_ROOT*(root - 1.0);
  /* An odd degree needs one more multiplication, by b, to make the form of
   * even degree.  The term is larger than that multiplication alone, because
   * it was fitted to what the sieve shows and an odd degree also restricts
   * the denominators to squares, which leaves fewer survivors per
   * denominator and so a colder check. */
  if(degree & 1)
  { cost += RATPOINTS_CHECK_STEP + RATPOINTS_CHECK_LIMB*limbs; }
  return(cost);
}

/* ----------------------------------------------------------------------
 * Correcting the number of primes from what the sieve is actually doing
 *
 * sieving_info picks sp1, sp2 and sp3 before anything has been sieved, from
 * R(n), the product of the densities.  That prediction is wrong in a way no
 * amount of care will fix.  The non-reduced representations (k*a, k*b) of a
 * rational point are the same rational number, so they give the same value
 * of f and pass every prime test there is: a floor of survivors outlives any
 * amount of sieving, and only the test for common factors removes it.  R(n)
 * cannot see that floor, so it always overshoots -- which is why sp2 could
 * never be predicted and ended up as a tuned offset.
 *
 * The sieve itself knows better.  The scan visits every bit array anyway, so
 * counting the non-empty ones costs an increment on a path taken half a per
 * cent of the time; the survivors of the second phase and of the test for
 * common factors are counted as cheaply.  Two such counts pin both terms of
 *
 *          S(n) = floor + chance * R(n)
 *
 * and the marginal rule can then be applied to the curve in hand rather than
 * to the predicted one.
 *
 * Changing the number of primes part-way through a run is safe by
 * construction: sieving only ever removes numerators that cannot be points,
 * so using more or fewer of them for later denominators changes the running
 * time and nothing else.  What it must not do is get ahead of bp_list, which
 * is why sp3_valid says how many of its entries are up to date.
 * ---------------------------------------------------------------------- */

/* how much data is wanted before the first correction, in numerator words;
 * after that the next one waits until twice as much has been seen */
#define RP_ADAPT_WORDS 1000000UL
#define RP_ADAPT_ARRAYS 1000UL   /* ...and this many non-empty bit arrays */
#define RP_ADAPT_BITS 200UL      /* ...and this many survivors of phase 2 */

static void adapt_primes(ratpoints_args *args)
{ ratpoints_sieve_entry **sieve_list
    = (ratpoints_sieve_entry **)args->sieve_list;
  double u = args->run_words, d = args->run_denoms;
  double cost_table = (args->cost_table >= 0.0) ? args->cost_table
                                                : RATPOINTS_COST_TABLE;
  double words = (double)args->n_words;
  double s1, s2, r1, r2, chance, level, s, rate;
  long n, sp1 = args->sp1, sp2 = args->sp2, max = args->sp3_max;
  /* 1 (the default) corrects the third stage only; 2 also corrects sp2 */
  long mode = (args->adapt < 0) ? 1 : args->adapt;

  /* next time, when twice as much has been seen */
  args->adapt_at = args->n_words + args->n_words;

  if(words <= 0.0 || sp2 <= sp1 || sp1 <= 0) { return; }
  if(args->n_arrays < RP_ADAPT_ARRAYS || args->n_bits < RP_ADAPT_BITS)
  { return; }

  /* The two rates the run has shown, per numerator word.  The first is
   * cumulative -- sp1 never moves, so every word swept measures the same
   * thing -- but the second is not: everything downstream of the first phase
   * was counted under whatever sp2 was in force, so those counters are reset
   * whenever sp2 changes and only the words since then divide into them. */
  { double words_2 = (double)(args->n_words - args->n_words_2);

    if(words_2 <= 0.0) { return; }
    s1 = (double)args->n_arrays/words;
    s2 = (double)args->n_bits/words_2;
  }

  r1 = 1.0;
  for(n = 0; n < sp1; n++) { r1 *= sieve_list[n]->r; }
  r2 = r1;
  for(n = sp1; n < sp2; n++) { r2 *= sieve_list[n]->r; }
  if(r1 - r2 <= 0.0 || s1 <= s2) { return; }

  chance = (s1 - s2)/(r1 - r2);
  level = s1 - chance*r1;          /* the floor */
  if(level < 0.0) { level = 0.0; }

  /* How many primes the second phase should use.  Adding the next one costs
   * what it does on every bit array still in play, plus the fixed costs it
   * has to earn back over the run, and saves the survivors it removes -- all
   * of which would otherwise be extracted, tested for common factors, run
   * through the third stage and sometimes checked exactly.
   *
   * This is the more adventurous half of the correction, and it is asked for
   * separately (adapt >= 2), because it is a second answer to a question the
   * scaled offset already answers: both decide sp2, one from a measurement
   * and one from a fit, and whichever is applied later wins.  Correcting the
   * third stage (below) is not like that -- there the measurement replaces an
   * estimate that nothing else supplies. */
  if(mode >= 2)
  { long want = sp2;
    /* what a survivor costs downstream.  The exact check is one term of it,
     * and the only one the degree moves; how many survivors reach the check
     * is measured rather than assumed, since the counters are here anyway. */
    double cost_surv = RATPOINTS_COST_SURVIVOR
                        + ((double)args->n_checks/(double)args->n_bits)
                           *RATPOINTS_COST_CHECK*(args->check_rel - 1.0);

    rate = r2;
    s = level + chance*rate;
    for(n = sp2; n < max; n++)
    { double r = sieve_list[n]->r;
      double next = level + chance*rate*r;
      double cost = prime_cost(sieve_list[n]->p, RATPOINTS_COST_PHASE2*s, 1,
                               cost_table, u, d);

      if((s - next)*cost_surv <= cost) { break; }
      rate *= r; s = next; want++;
    }
    if(want == sp2)
    { /* nothing to add: see whether the last one is still worth having */
      for(n = sp2 - 1; n > sp1; n--)
      { double r = sieve_list[n]->r;
        double prev = level + chance*rate/r;
        double cost = prime_cost(sieve_list[n]->p,
                                 RATPOINTS_COST_PHASE2*prev, 1,
                                 cost_table, u, d);

        if((prev - s)*cost_surv > cost) { break; }
        rate /= r; s = prev; want--;
      }
    }
    if(want != sp2)
    { /* what was counted downstream belongs to the old sp2 */
      args->n_bits = 0; args->n_coprime = 0; args->n_checks = 0;
      args->n_sifts = 0; args->n_words_2 = args->n_words;
    }
    args->sp2 = want;
  }
  else { s = level + chance*r2; }  /* sp2 stands; the rate is what it was */

  /* And how many the third stage should use.  Here the measurement is the
   * number of survivors a denominator brings to the stage -- after the test
   * for common factors, which is the one thing no prime can help with -- so
   * it replaces both the predicted rate and the fitted fraction that stood
   * for the coprimality test. */
  if(args->sp3_extra < 0 && args->n_sifts > 0)
  { /* both fractions are of one exact check, which is dearer at a high
     * degree or with large coefficients: see check_cost() */
    double per_denom = ((args->sp3_per_denom >= 0.0) ? args->sp3_per_denom
                                                     : RATPOINTS_SP3_PER_DENOM)
                         /args->check_rel;
    double per_surv = RATPOINTS_SP3_PER_SURVIVOR/args->check_rel;
    double S = (double)args->n_coprime/(double)args->n_sifts;
    double sp2_old = level + chance*r2;
    long sp3;

    /* the second phase may just have moved; carry the measurement across */
    if(sp2_old > 0.0) { S *= s/sp2_old; }

    for(sp3 = args->sp2; sp3 < max; sp3++)
    { double r = sieve_list[sp3]->r;

      if(S*(1.0 - per_surv - r) <= per_denom) { break; }
      S *= r;
    }
    args->sp3 = sp3;
  }
  else if(args->sp3 < args->sp2) { args->sp3 = args->sp2; }

#ifdef RP_PRIME_STATS
  fprintf(stderr, "[adapt] words=%lu s1=%.3g s2=%.3g floor=%.3g"
          " sp1=%ld sp2=%ld sp3=%ld\n", args->n_words, s1, s2, level,
          args->sp1, args->sp2, args->sp3);
#endif
}

/* How many primes the first phase needs: enough of them that the expected
 * number of surviving numerators per 64-bit word falls to the target.
 * prec[] must be sorted by increasing r.  If the target cannot be reached
 * with the primes available, all of them are used. */
static long primes_for_phase_1(entry *prec, long pnp,
                               double bits_per_word, double target)
{ double rate = 1.0;
  long n;

  for(n = 0; n < pnp; n++)
  { rate *= prec[n].r;
    if(bits_per_word*rate <= target) { return(n + 1); }
  }
  return(pnp > 0 ? pnp : 1);
}

/* How many primes the second phase adds to the first.  A phase-2 prime is
 * paid for once -- its sieve table, and its step in bp_list -- and then used
 * for the whole run, so how many are worth having depends on how long the
 * run is; see RATPOINTS_SP2_U0 in ratpoints.h .  With u0 = 0 this is a flat
 * offset, which is what every version before 2.3 used. */
static long phase_2_offset(long extra, double u0, double u_words)
{ double e;

  if(u0 <= 0.0 || u_words <= 0.0) { return(extra); }
  e = (double)extra/(1.0 + u0/u_words);
  return((long)(e + 0.5));
}

/* u modulo p for a full-width u, by Barrett reduction.  With
 * m = floor(2^64/p), the quotient floor(u*m/2^64) is floor(u/p) or one less,
 * so one conditional subtraction finishes the job: two multiplications in
 * place of a division.  The Horner loop below runs this once for every
 * residue and every prime, which is O(degree*p) per prime per curve.
 *
 * This is not the reduction sift.c uses.  That one is cheaper still, but it
 * is exact only below 2^32, and the accumulator here runs up to p^(degree+1).
 */
#ifdef __SIZEOF_INT128__
static inline unsigned long barrett(unsigned long u, unsigned long p,
                                    unsigned long m)
{ unsigned long r = u - (unsigned long)(((__uint128_t)u*m) >> 64)*p;

  return((r >= p) ? r - p : r);
}
#else
static inline unsigned long barrett(unsigned long u, unsigned long p,
                                    unsigned long m)
{ (void)m; return(u % p); }
#endif

/* How many Horner steps the accumulator survives without a reduction: after
 * k of them it is below p^(k+1), so k+1 must not exceed the number of
 * primes' worth of bits in a long. */
#define RP_HORNER_STEPS ((long)(LONG_LENGTH/RATPOINTS_MAX_BITS_IN_PRIME) - 1)

/* Look at one prime and record what it says about the curve.
 *
 * Fills in the table is_f_square[0..p], where entry a says whether f(a) is a
 * square modulo p and the last entry whether there are points at infinity,
 * and counts the residues that admit points.  When the prime carries any
 * information -- that is, when some residue does not -- a sieve entry is
 * built for it and *prec_entry is filled in with that entry and the density
 * r of the admissible residues.
 *
 * Returns 1 in that case, 0 when the prime says nothing, and -1 when the
 * curve has no points modulo p at all, so that it has no rational points.
 * coeffs_mod_p and is_f_square_p hand back what the caller needs for the
 * test on forbidden divisors of the denominator.
 */
static int examine_prime(ratpoints_args *args, long pn,
                         int use_c_long, long *c_long,
                         long *coeffs_mod_p, int **is_f_square_p,
                         entry *prec_entry)
{
  mpz_t *c = args->cof;
  long degree = args->degree;
  long p = prime[pn];
  long n, a, np; /* np counts the x-coordinates that give points mod p */
  unsigned long recip = ULONG_MAX/(unsigned long)p; /* = floor(2^64/p) */
  int *is_f_square = args->int_next;

  args->int_next += p + 1; /* need space for (p+1) int's */
  *is_f_square_p = is_f_square;

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
  /* Evaluate f at every residue by Horner, reducing on a fixed schedule.
   * Up to RP_HORNER_STEPS steps fit in a long without one, which at the
   * default PRIME_SIZE is every degree up to 7; beyond that the accumulator
   * is reduced every RP_HORNER_STEPS steps.  The schedule is fixed rather
   * than decided by testing the accumulator, which is what this used to do:
   * that test is a data-dependent branch in the innermost loop, and with the
   * reduction now two multiplications instead of a division it is cheaper to
   * reduce on a schedule than to work out whether to.
   * (It also puts right what raising PRIME_SIZE from 7 to 8 did to degree 8:
   * it moved the boundary of the division-free path from degree 8 to 7, so
   * genus 3 with an even model took a conditional division in every step.) */
  if(degree <= RP_HORNER_STEPS)
  { for(a = 1 ; a < p; a++)
    { unsigned long s = coeffs_mod_p[degree];

      for(n = degree - 1 ; n >= 0 ; n--)
      { s *= a; s += coeffs_mod_p[n]; }
      /* here, s < p^(degree+1) <= max. long */
      s = barrett(s, p, recip);
      if((is_f_square[a] = squares[pn][s])) { np++; }
    }
  }
  else
  { for(a = 1 ; a < p; a++)
    { unsigned long s = coeffs_mod_p[degree];
      long k = 0;

      for(n = degree - 1 ; n >= 0 ; n--)
      { s *= a; s += coeffs_mod_p[n];
        if(++k == RP_HORNER_STEPS) { s = barrett(s, p, recip); k = 0; }
      }
      s = barrett(s, p, recip);
      if((is_f_square[a] = squares[pn][s])) { np++; }
    }
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
  if(np == 0 && !is_f_square[p]) { return(-1); }

  if(np >= p) { return(0); } /* the prime carries no information */

  { double r = is_f_square[p] ? ((double)(np*(p-1) + p))/((double)(p*p))
                              : (double)np/(double)p;

    prec_entry->r = r;
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
    /* the reciprocal the third stage reduces with; see stage3() in sift.c .
     * One division per prime and curve, against one per survivor saved. */
    se->magic = ULONG_MAX/(unsigned long)p + 1;
    /* the entry keeps the density too, so that the choice of primes can be
     * revisited during the run, when prec[] is long gone */
    se->r = prec_entry->r;
    se->offset = offsets[pn];
    /* sieves0 is 64-bit words, but is read as bit-arrays; it is given
     * the alignment of ratpoints_bit_array in gen_find_points_h.c . */
    se->sieve[0] = (ratpoints_bit_array *)&sieves0[pn][0];
    for(i = 1; i < p; i++) { se->sieve[i] = NULL; }

    prec_entry->ssp = se;
  }
  return(1);
}

/************************************************************************
 * Collect the sieving information                                      *
 ************************************************************************/

/* The number of numerators denominator b has to consider: the part of
 * b*domain that lies within the height bound.  Used both to size the run
 * (run_shape below) and to estimate how many survivors a denominator brings
 * to the third stage. */
static double numerators_for(const ratpoints_args *args, double b, double H)
{ double sum = 0.0;
  long k;

  for(k = 0; k < args->num_inter; k++)
  { double lo = b*args->domain[k].low, up = b*args->domain[k].up;

    if(lo < -H) { lo = -H; }
    if(up > H) { up = H; }
    if(up > lo) { sum += up - lo; }
  }
  return(sum);
}

/* How big the run is: the number of denominators that will actually be
 * sifted, and the number of 64-bit words of numerators they sweep between
 * them.  Both are wanted by the rule that picks the sieving primes, because
 * two of the costs of a prime are paid once and then spread over the whole
 * run -- its sieve table, built for at most p denominator classes, and its
 * entry in bp_list, stepped once per denominator.  Per word of numerators
 * those come to k*p*min(D,p)/U and l*D/U, and they are the reason the best
 * number of primes at a height bound of 200000 is not the best number at
 * 16383.
 *
 * Nothing here needs any sieving.  The denominators that get sifted are
 * those that pass the 2-adic mask on b, have an admissible numerator at all,
 * are not divisible by a forbidden divisor and pass the Jacobi symbol test
 * where it applies; the first three are periodic and are counted exactly,
 * and the fourth lets through half of what is left.  The numerators of one
 * denominator are piecewise linear in b with a handful of breakpoints, so a
 * midpoint sample over the range of b is accurate to a fraction of a per
 * cent.
 *
 * The result is an estimate, and a biased one -- the Jacobi factor is an
 * average, and the valuation test of the use_squares1 path is not modelled
 * at all.  That is by design: it is used only to compare a fixed cost with a
 * per-word one, where being right to within a factor of about 1.5 moves the
 * chosen number of primes by less than one.
 */
#define RUN_SHAPE_SAMPLES 64

/* The fraction of all integers b with v_p(b) in the set the mask describes
 * (bit m set <==> v_p(b) = m); the density of v_p(b) = m is (p-1)/p^(m+1). */
static double forbidden_fraction(long p, unsigned long mask)
{
  double f = 0.0;
  double q = 1.0/(double)p; /* 1/p^m */
  long m;

  for(m = 1; m < (long)LONG_LENGTH && (mask >> m) != 0; m++)
  { if((mask >> m) & 1) { f += q*(1.0 - 1.0/(double)p); }
    q /= (double)p;
  }
  return(f);
}

static void run_shape(ratpoints_args *args, bit_selection which_bits,
                      unsigned long den_bits,
                      const ratpoints_bit_array *num_bits,
                      long fba, long fdc,
                      double *n_denom, double *u_words)
{ double H = (double)args->height;
  double keep = 1.0;    /* fraction of the candidates that reach sift() */
  double count = 0.0;   /* candidate denominators */
  double nums = 0.0;    /* numerators they sweep, before that fraction */
  long i, j;

  if(args->flags & RATPOINTS_USE_SQUARES)
  { /* the denominators are the squares in [b_low, b_high] */
    double klo = ceil(sqrt((double)args->b_low));
    double khi = floor(sqrt((double)args->b_high));
    long good = 0;

    if(khi >= klo)
    { count = khi - klo + 1.0;
      for(i = 0; i < RUN_SHAPE_SAMPLES; i++)
      { double k = klo + (khi - klo)*((double)i + 0.5)/RUN_SHAPE_SAMPLES;
        nums += numerators_for(args, k*k, H);
      }
      nums *= count/RUN_SHAPE_SAMPLES;
    }
    /* only the mask on b mod 16 applies, and b = k^2 mod 16 has period 8 */
    for(j = 0; j < 8; j++)
    { if(EXT0(num_bits[(j*j) & 0xf])) { good++; } }
    keep = (double)good/8.0;
  }
  else if(args->flags & RATPOINTS_USE_SQUARES1)
  { /* squares times the divisors of the leading coefficient */
    long *divisors = (long *)args->divisors;
    long n;
    long good = 0, tried = 0;

    for(n = 0; divisors[n]; n++)
    { double d = (double)divisors[n];
      double klo = ceil(sqrt((double)args->b_low/d));
      double khi = floor(sqrt((double)args->b_high/d));

      if(klo < 1.0) { klo = 1.0; }
      if(khi >= klo)
      { double c = khi - klo + 1.0;

        for(i = 0; i < RUN_SHAPE_SAMPLES; i++)
        { double k = klo + (khi - klo)*((double)i + 0.5)/RUN_SHAPE_SAMPLES;
          nums += c*numerators_for(args, d*k*k, H)/RUN_SHAPE_SAMPLES;
        }
        count += c;
      }
      for(j = 0; j < 8; j++, tried++)
      { if(EXT0(num_bits[(divisors[n]*j*j) & 0xf])) { good++; } }
    }
    if(tried) { keep = (double)good/(double)tried; }
  }
  else
  { /* every denominator in the range is a candidate */
    double blo = (double)args->b_low, bhi = (double)args->b_high;
    long good = 0;

    if(bhi >= blo)
    { count = bhi - blo + 1.0;
      for(i = 0; i < RUN_SHAPE_SAMPLES; i++)
      { double b = blo + (bhi - blo)*((double)i + 0.5)/RUN_SHAPE_SAMPLES;
        nums += numerators_for(args, b, H);
      }
      nums *= count/RUN_SHAPE_SAMPLES;
    }
    /* b congruent to j modulo 64 is tested against bit j of den_bits (the
     * loop shifts before it tests, which is what puts b and the bit index
     * in step) and against num_bits[b mod 16] */
    for(j = 0; j < 64; j++)
    { if(((den_bits >> j) & 1UL) && EXT0(num_bits[j & 0xf])) { good++; } }
    keep = (double)good/64.0;

    if(args->flags & RATPOINTS_CHECK_DENOM)
    { forbidden_entry *fb = (forbidden_entry *)args->forb_ba;
      forbidden_val *fd = (forbidden_val *)args->forbidden;

      for(i = 0; i < fba; i++) { keep *= 1.0 - 1.0/(double)fb[i].p; }
      for(i = 0; i < fdc; i++)
      { keep *= 1.0 - forbidden_fraction(fd[i].p, fd[i].mask); }
      /* the Jacobi symbol lets through half of the rest */
      if(args->flags & RATPOINTS_USE_JACOBI) { keep *= 0.5; }
    }
  }

  /* Only every other numerator is looked at unless both parities are in
   * play.  That includes num_none, which says only that no odd denominator
   * has an admissible numerator: the even denominators are still sieved,
   * with num_odd forced in sift(), so exactly half the numerators are looked
   * at for them.  (Setting nums to zero here put U on its floor of 1, which
   * switched the second and third stages off for those curves.)  Not
   * modelled: with num_all the even denominators also sweep only half, so U
   * is over-estimated by up to a third there; see TODO item 18. */
  if(which_bits != num_all) { nums *= 0.5; }

  *n_denom = keep*count;
  *u_words = keep*nums/(double)LONG_LENGTH;
  if(*n_denom < 1.0) { *n_denom = 1.0; }
  if(*u_words < 1.0) { *u_words = 1.0; }
}

/**************************************************************************
 * The p-adic test on the denominator, for a prime p that divides the
 * leading coefficient.
 *
 * If p does not divide the leading coefficient c[d] (d even), a denominator
 * divisible by p is excluded exactly when c[d] is a non-square mod p, since
 * then F(a,b) = c[d] a^d mod p is one too.  If p | c[d], then F(a,b) is 0
 * mod p for every such denominator, which is a square, and the mod-p test
 * cannot say anything; the question is p-adic.
 *
 * Write v_p(b) = m >= 1, so p does not divide a, and w_j = v_p(c[d-j]).
 * The term c[d-j] a^(d-j) b^j of F(a,b) has valuation w_j + j*m.  If one
 * term has strictly smaller valuation than all the others, it sets v_p(F)
 * and F/p^v mod p, and F is not a square if that valuation is odd, nor if j
 * is even and the unit part of c[d-j] is a non-square mod p (with d even,
 * a^(d-j) and (b/p^m)^j are both squares then).  With j odd, (b/p^m)^j runs
 * through non-squares as well as squares as b varies, so nothing follows;
 * and when two terms tie for the minimum they may cancel, which is exactly
 * what happens at a point of such a curve, so nothing follows either.
 *
 * For m large the leading term j = 0 wins on its own, so the excluded set
 * is a tail of ones from some n on when v_p(c[d]) is odd or the unit part of
 * c[d] is a non-square -- "no denominator is divisible by p^n" -- possibly
 * with a few isolated valuations below it.  The commonest cases: v_p(c[d])
 * = 1 and p | c[d-1] exclude p itself, and v_p(c[d]) = 1 with p not dividing
 * c[d-1] exclude p^2 (at m = 1 the two top terms tie).
 *
 * Returns the excluded valuations as a bit mask, bit m set <==> v_p(b) = m
 * is excluded, for the m with p^m <= b_high; 0 if there are none.  The
 * caller decides how to test for them.
 */
static unsigned long forbidden_valuations(ratpoints_args *args, long pn,
                                          const long *coeffs_mod_p)
{
  mpz_t *c = args->cof;
  long degree = args->degree;
  long p = prime[pn];
  long w[degree+1]; /* w[j] = v_p(c[degree-j]); a zero coefficient is absent */
  long r[degree+1]; /* the unit part of c[degree-j] mod p, in [1, p-1] */
  int present[degree+1];
  unsigned long mask = 0;
  long j, m, pm;

  for(j = 0; j <= degree; j++)
  { long k = degree - j;

    present[j] = (mpz_sgn(c[k]) != 0);
    if(!present[j]) { w[j] = 0; r[j] = 0; }
    else if(coeffs_mod_p[k] == 0)
    { w[j] = valuation(c[k], p, &r[j], args->work[0]);
      /* valuation() works on |c[k]|, and the sign matters here */
      if(mpz_sgn(c[k]) < 0) { r[j] = p - r[j]; }
    }
    else { w[j] = 0; r[j] = coeffs_mod_p[k]; }
  }

  for(m = 1, pm = p; pm <= args->b_high && m < (long)LONG_LENGTH - 1; m++)
  { long jmin = 0, vmin = w[0];
    int tie = 0;

    for(j = 1; j <= degree; j++)
    { if(present[j])
      { long v = w[j] + j*m;

        if(v < vmin) { vmin = v; jmin = j; tie = 0; }
        else if(v == vmin) { tie = 1; }
    } }
    if(!tie && ((vmin & 1) || (((jmin & 1) == 0) && !squares[pn][r[jmin]])))
    { mask |= 1UL << m; }
    if(pm > args->b_high/p) { break; } /* p^(m+1) > b_high */
    pm *= p;
  }
  return(mask);
}

static long sieving_info(ratpoints_args *args,
                         int use_c_long, long *c_long,
                         ratpoints_sieve_entry **sieve_list,
                         double bits_per_word, int may_extend,
                         bit_selection which_bits, unsigned long den_bits,
                         const ratpoints_bit_array *num_bits)
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
  forbidden_val *forbidden = (forbidden_val *)args->forbidden;

  /* How many primes to look at.  The loop below may raise this: see the
   * comment at its end. */
  long pn_lim = args->num_primes;
  double target = (args->survivors_per_word > 0.0) ? args->survivors_per_word
                                                   : RATPOINTS_SURVIVORS_PER_WORD;
  long sp2_extra = (args->sp2_extra >= 0) ? args->sp2_extra
                                          : RATPOINTS_SP2_EXTRA;
  double sp2_u0 = (args->sp2_u0 >= 0.0) ? args->sp2_u0 : RATPOINTS_SP2_U0;
  double cost_table = (args->cost_table >= 0.0) ? args->cost_table
                                                : RATPOINTS_COST_TABLE;

  /* Whether the number of primes is ours to choose, and so ours to correct
   * as the run goes on; see adapt_primes.  If the caller fixed sp1 or sp2,
   * they stay fixed. */
  args->adapt_at = (args->adapt != 0 && args->sp1 < 0 && args->sp2 < 0)
                     ? RP_ADAPT_WORDS : ULONG_MAX;
  args->sp3_valid = 0;

  /* What one exact check costs on this curve, relative to the curves the
   * third-stage constants were tuned on.  It is what those constants are
   * fractions of, so a curve whose check is dearer -- a high degree, or
   * large coefficients, or both -- is worth more third-stage primes. */
  args->check_rel = check_cost(args)/RATPOINTS_CHECK_REFERENCE;
  if(args->check_rel <= 0.0) { args->check_rel = 1.0; }

  /* How big the run is.  This is wanted before the first prime is looked at,
   * because the rule that decides whether to look past the primes we were
   * given uses the same offset as the final choice does; it is computed
   * again below, once the forbidden divisors are known and the estimate can
   * take them into account. */
  run_shape(args, which_bits, den_bits, num_bits, 0, 0,
            &args->run_denoms, &args->run_words);
  sp2_extra = phase_2_offset(sp2_extra, sp2_u0, args->run_words);

  /* initialize sieve in se_buffer */
  for(pn = 0; pn < pn_lim; pn++)
  { long coeffs_mod_p[degree+1];
           /* The coefficients of f reduced modulo p */
    long p = prime[pn];
    int *is_f_square;
    int info = examine_prime(args, pn, use_c_long, c_long,
                             &coeffs_mod_p[0], &is_f_square, &prec[pnp]);

    if(info < 0)
    { return(p); /* no points mod p, hence no rational points */ }
    if(info > 0)
    { prec[pnp].key = prime_key(prec[pnp].r, p, 1.0, 1, cost_table,
                                args->run_words, args->run_denoms);
      pnp++;
    }

    if((args->flags & RATPOINTS_CHECK_DENOM)
         && fba + fdc < args->max_forbidden)
    { /* record forbidden divisors of the denominator */
      if(coeffs_mod_p[degree] == 0)
      { /* p divides the leading coefficient: see forbidden_valuations */
        unsigned long mask = forbidden_valuations(args, pn, &coeffs_mod_p[0]);

        if(mask)
        { /* the valuations a denominator can have at all */
          unsigned long all = 0;
          long m, pm;

          for(m = 1, pm = p;
              pm <= args->b_high && m < (long)LONG_LENGTH - 1; m++)
          { all |= 1UL << m;
            if(pm > args->b_high/p) { break; }
            pm *= p;
          }
          if((mask & all) == all)
          { /* every one of them: no denominator is divisible by p, which
             * the bit arrays test for */
            forb_ba[fba].p     = p;
            forb_ba[fba].start = &sieves0[pn][0];
            forb_ba[fba].end   = &sieves0[pn][p];
            forb_ba[fba].curr  = forb_ba[fba].start;
            fba++;
          }
          else
          { forbidden[fdc].p = p; forbidden[fdc].mask = mask; fdc++; }

#ifdef DEBUG
          printf("\nexcluding denominators b with v_%ld(b) in %#lx"
                 " (bit m <==> v = m)\n", p, mask);
          fflush(NULL);
#endif

        }
      }
      else if(!is_f_square[p])
      { /* leading coefficient is a non-square mod p:
         * a denominator divisible by p is excluded */
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

    /* Once the primes we were told to look at are used up, look at more if
     * the choice below would otherwise be cramped.  The second phase sieves
     * the survivors of the first with sp2 - sp1 further primes, and there
     * have to be that many left over; a curve with very many rational points
     * makes f a square modulo every residue for the smallest primes, so those
     * carry no information and are dropped above, and without this the second
     * phase can end up with nothing to sieve with at all.
     * Adding a prime can only lower sp1 -- the n smallest of a larger set have
     * a smaller product -- and can only raise pnp, so the shortfall shrinks
     * with every prime added and this stops at the first one that is enough.
     */
    if(may_extend && pn + 1 == pn_lim && pn_lim < RATPOINTS_NUM_PRIMES)
    { long s1, want;

      qsort(prec, pnp, sizeof(entry), compare_entries);
      s1 = (args->sp1 >= 0) ? args->sp1
                            : primes_for_phase_1(prec, pnp, bits_per_word, target);
      want = (args->sp2 >= 0) ? args->sp2 : s1 + sp2_extra;
      if(pnp < want) { pn_lim++; }
    }

  } /* end for pn */

  /* Terminate the array of forbidden divisors, having first looked for
   * more of them among the primes the loop above did not reach.  This is
   * done here, before the primes are chosen, because the choice needs to
   * know how many denominators will survive these tests: see run_shape. */
  if((args->flags & RATPOINTS_CHECK_DENOM)
       && !mpz_perfect_square_p(c[degree]))
       /* the test below asks for a non-square residue, which a square
        * leading coefficient never has; such a curve gets here since the
        * valuation test above applies to it */
  { long n;

    for(n = pn_lim;
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
  }
  if(args->flags & RATPOINTS_CHECK_DENOM)
  { forb_ba[fba].p = 0;        /* terminating zero */
    forbidden[fdc].p = 0;      /* terminating zero */
    /* args->max_forbidden is an input and stays one: writing the number
     * found into it latched the test off for a caller that reuses args. */
  }

  /* nothing for the checked denominator loop to test? then use the plain one */
  if(fba + fdc == 0 && !(args->flags & RATPOINTS_USE_JACOBI))
  { args->flags &= ~RATPOINTS_CHECK_DENOM; }

  /* the sieve tables live in a block that was reserved for args->num_primes
   * primes; if the loop went further, that block has to grow.  Nothing has
   * been taken from it yet -- the tables are built lazily during the sieving
   * itself -- so it can simply be replaced. */
  if(pn_lim > args->ba_buffer_primes)
  { free(args->ba_buffer_na);
    alloc_ba_buffer(args, pn_lim);
  }

  /* The run shape again, now that the forbidden divisors are known and can
   * be taken off the denominator count; see run_shape.  The keys the primes
   * were given inside the loop used the first estimate, which does not know
   * about those divisors and so overstates the run, so they are computed
   * again here before anything is sorted for good. */
  { long e = (args->sp2_extra >= 0) ? args->sp2_extra : RATPOINTS_SP2_EXTRA;
    long n;

    run_shape(args, which_bits, den_bits, num_bits, fba, fdc,
              &args->run_denoms, &args->run_words);
    sp2_extra = phase_2_offset(e, sp2_u0, args->run_words);
    for(n = 0; n < pnp; n++)
    { prec[n].key = prime_key(prec[n].r, prec[n].ssp->p, 1.0, 1, cost_table,
                              args->run_words, args->run_denoms);
    }
  }

  /* sort the array to get at the best primes */
  qsort(prec, pnp, sizeof(entry), compare_entries);

  /* Choose sp1 and sp2 unless they were given.
   * prec[] is now sorted by increasing r, where r is the density of the
   * numerators that are admissible modulo the corresponding prime, so the
   * expected fraction of numerators surviving the first n primes is the
   * product of the first n values of r.  Multiplied by the number of bits
   * actually set in a bit-array to begin with, that is the expected number
   * of survivors per bit-array; see the comment on
   * RATPOINTS_SURVIVORS_PER_WORD in ratpoints.h . */
  if(args->sp1 < 0)
  { args->sp1 = primes_for_phase_1(prec, pnp, bits_per_word, target); }

  /* Rank what is left again, for the second phase.  There a prime is applied
   * only to the bit arrays that survived the first phase, so its per-word
   * cost is smaller by the survival rate -- which makes the fixed cost of
   * its table weigh far more heavily, and the size of the prime matter far
   * more than it does in the first phase. */
  if(args->sp1 >= 0 && args->sp1 < pnp)
  { long n;
    double rate = bits_per_word;

    for(n = 0; n < args->sp1; n++) { rate *= prec[n].r; }
    for(n = args->sp1; n < pnp; n++)
    { prec[n].key = prime_key(prec[n].r, prec[n].ssp->p,
                              RATPOINTS_COST_PHASE2*rate, 1, cost_table,
                              args->run_words, args->run_denoms);
    }
    qsort(&prec[args->sp1], pnp - args->sp1, sizeof(entry), compare_entries);
  }

  if(args->sp2 < 0) { args->sp2 = args->sp1 + sp2_extra; }

  /* update sp2 and sp1 if necessary */
  if(args->sp2 > pnp) { args->sp2 = pnp; }
  if(args->sp1 > args->sp2) { args->sp1 = args->sp2; }


  /* put the sorted entries into sieve_list */
  { long n;

    for(n = 0; n < args->sp2; n++)
    { sieve_list[n] = prec[n].ssp; }
  }

  /* Choose sp3, the number of primes the third stage adds to those two.
   * That stage tests one surviving numerator at a time and needs no sieve
   * table, so a prime costs it one test per survivor and one subtraction per
   * denominator, and nothing per curve beyond what has been done here.  A
   * prime is therefore worth adding as long as the survivors it removes are
   * worth more than the denominators it is carried through, which is the
   * rule below; see RATPOINTS_SP3_PER_SURVIVOR in ratpoints.h .  Both costs
   * are fractions of one exact check, and what one check costs depends on
   * the curve, so both are divided by check_rel: at a high degree or with
   * large coefficients the check is dearer and more primes are worth having,
   * and at degree 3 or 4 it is cheaper and fewer are.
   *
   * S is the expected number of survivors a denominator still has when the
   * stage begins.  It is the number of numerators the denominator considers,
   * thinned by the sixteen-fold pre-sieve and by the primes of the first two
   * phases, and thinned again by the test for common factors, which runs
   * before this stage and which no prime can help with: a numerator sharing
   * a factor with the denominator stands for a fraction that has already
   * been looked at with a smaller denominator, so it passes every prime.
   */
  { long sp3 = args->sp2;
    long sp3_want = (args->sp3_extra >= 0) ? args->sp2 + args->sp3_extra
                                           : RATPOINTS_NUM_PRIMES;
    double per_denom = ((args->sp3_per_denom >= 0.0) ? args->sp3_per_denom
                                                     : RATPOINTS_SP3_PER_DENOM)
                        /args->check_rel;
    double per_surv = RATPOINTS_SP3_PER_SURVIVOR/args->check_rel;
    double S;

    if(sp3_want > RATPOINTS_NUM_PRIMES) { sp3_want = RATPOINTS_NUM_PRIMES; }

    /* The average number of numerators a denominator has to consider.  The
     * run shape already has the total in words, so this is just the mean per
     * denominator; run_shape counts only the numerators that are actually
     * looked at, so no further halving is wanted here. */
    S = args->run_words*(double)LONG_LENGTH/args->run_denoms;
    { long n;

      S *= bits_per_word/(double)LONG_LENGTH;
      for(n = 0; n < args->sp2; n++) { S *= prec[n].r; }
      S *= RATPOINTS_SP3_COPRIME;
    }

    while(sp3 < sp3_want)
    { double r;

      if(sp3 >= pnp)
      { /* the primes looked at so far are used up: look at one more */
        long coeffs_mod_p[degree+1];
        int *is_f_square;
        int info;

        if(!may_extend || pn_lim >= RATPOINTS_NUM_PRIMES) { break; }
        info = examine_prime(args, pn_lim, use_c_long, c_long,
                             &coeffs_mod_p[0], &is_f_square, &prec[pnp]);
        pn_lim++;
        if(info < 0)
        { return(prime[pn_lim-1]); /* no points mod p */ }
        if(info == 0) { continue; } /* it says nothing; try the next one */
        /* the third stage builds no table, so its primes are ranked by what
         * they say alone, which is what the selection below does */
        prec[pnp].key = prec[pnp].r;
        pnp++;
      }

      /* the best of the primes not yet spoken for; only the ones this stage
       * takes need to be in order, so this is a selection sort that stops
       * as soon as the rule below does */
      { long m, best = sp3;

        for(m = sp3 + 1; m < pnp; m++)
        { if(prec[m].r < prec[best].r) { best = m; } }
        if(best != sp3)
        { entry t = prec[sp3]; prec[sp3] = prec[best]; prec[best] = t; }
      }

      r = prec[sp3].r;
      if(args->sp3_extra < 0 && S*(1.0 - per_surv - r) <= per_denom)
      { break; }
      S *= r;
      sieve_list[sp3] = prec[sp3].ssp;
      sp3++;
    }
    args->sp3 = sp3;

    /* Put the rest of the primes in sieve_list too, in the order the third
     * stage would take them.  They cost nothing to keep -- no table is built
     * and no bp_list entry stepped until a prime is actually used -- and
     * having them there is what lets adapt_primes() reach for one more
     * during the run. */
    { long n;

      if(pnp > args->sp2)
      { qsort(&prec[args->sp2], pnp - args->sp2, sizeof(entry),
              compare_by_r);
      }
      for(n = args->sp2; n < pnp; n++) { sieve_list[n] = prec[n].ssp; }
      args->sp3_max = pnp;
      if(args->sp3_max < args->sp3) { args->sp3_max = args->sp3; }
    }
  }

  /* The third stage may have looked at further primes, and those are now in
   * sieve_list, where adapt_primes() can promote one into the second phase
   * during the run -- at which point it does build a sieve table.  So the
   * buffer the tables come out of has to cover every prime looked at, not
   * just the ones the first two phases started with.  Nothing has been taken
   * from it yet: the tables are built lazily while sieving. */
  if(pn_lim > args->ba_buffer_primes)
  { free(args->ba_buffer_na);
    alloc_ba_buffer(args, pn_lim);
  }

  /* the reciprocals the first phase reduces word numbers with, in the order
   * the primes are used; see the note in find_points_init */
  { long n;
    unsigned long *magics = (unsigned long *)args->magics;

    for(n = 0; n < args->sp3_max; n++) { magics[n] = sieve_list[n]->magic; }
  }


#ifdef RP_PRIME_STATS
  /* Development instrumentation: one line per curve saying how many primes
   * carried information, how the three stages divide them up, and the
   * density r of each. */
  { long n;

    fprintf(stderr, "[primestats] pn_lim=%ld pnp=%ld sp1=%ld sp2=%ld sp3=%ld"
            " bpw=%.2f U=%.6g D=%.6g", pn_lim, pnp, args->sp1, args->sp2,
            args->sp3, bits_per_word, args->run_words, args->run_denoms);
    for(n = 0; n < pnp; n++)
    { fprintf(stderr, " %ld:%.4f", prec[n].ssp->p, prec[n].r); }
    fprintf(stderr, "\n");
  }
#endif

  if(args->flags & RATPOINTS_VERBOSE)
  { printf("  %.1f bits set per word, %ld primes looked at"
           " ==> use %ld primes in the first phase, %ld altogether,\n"
           "  and %ld more in the third stage\n",
           bits_per_word, pn_lim, args->sp1, args->sp2,
           args->sp3 - args->sp2);
  }


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
  sieve_spec ssp[args->sp2 > 0 ? args->sp2 : 1]; /* length 0 is undefined */
  /* what the third stage needs per denominator; see find_points_init on why
   * it is not an array here */
  check_spec *csp = (check_spec *)args->stage3_list;
  int do_setup = 1;
  RP_SIFT_TIC(t_sift);

  args->n_sifts++;

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
        { if(b*inter.low > height)
          { RP_SIFT_TOC(t_sift);
            return(total); /* remaining numerator intervals are empty */
          }
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
        RP_SETUP_TIC(t_setup);

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
          if(sptr) { ssp[n].ptr = sptr; }
          else
          { RP_INIT_TIC(t_init);
            ssp[n].ptr = (*(se->init))(se, bp, args);
            RP_INIT_TOC(t_init, p);
          }
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

        /* and the primes of the third stage, which need no table: only the
         * inverse of b modulo each of them.  It is a table lookup, not a
         * division, because the inverses modulo every prime that can be used
         * are compiled in (see gen_find_points_h.c).  Note that bp is the
         * denominator itself here, not halved as it is above: the third
         * stage tests the numerator, not the bit that stands for it. */
        for(n = args->sp2; n < args->sp3; n++)
        { ratpoints_sieve_entry *se = sieve_list[n];
          long bp = bp_list[n];
          long m = n - args->sp2;

          csp[m].p = se->p;
          csp[m].is_f_square = se->is_f_square;
          csp[m].binv = bp ? se->inverses[bp] : 0;
          csp[m].magic = se->magic;
          /* the numerator is shifted by this multiple of p to make it
           * non-negative, so that the reduction can be the cheap one; a zero
           * says the shifted value would not fit and the slow path is to be
           * taken (see stage3() in sift.c) */
          csp[m].bias = ((double)se->p*(double)(2*args->height)
                           < RP_STAGE3_LIMIT)
                          ? se->p*args->height : 0;
        }
        RP_SETUP_TOC(t_setup);
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
          /* The bit arrays are not written here.  The first phase's
           * first prime ANDs the 2-adic pattern in as it sieves, which
           * saves a store and a load on every one of them; what is left
           * for sift0 to do afterwards is the two boundary words and the
           * padding, and it is told about them like this. */
          { long mask_low = 0, mask_high = 0, n_pad = 0;

            if(w_low0 == w_low)
            /* lower bits of the first bit array are to be set to zero */
            { mask_low = low - RBA_LENGTH * w_low; }
            if(w_high0 == w_high)
            /* upper bits of the last bit array are to be set to zero */
            { mask_high = RBA_LENGTH * w_high - high; }

#if (RATPOINTS_CHUNK > 1)
            /* if necessary, increase the range to a multiple of
             * RATPOINTS_CHUNK; the extra bit arrays are zeroed there too */
            while((range + n_pad)%RATPOINTS_CHUNK != 0) { n_pad++; }
            range += n_pad; w_high0 += n_pad;
#endif

            total += _ratpoints_sift0(b, w_low0, w_high0, args, which_bits,
                                      survivors, bits16, mask_low, mask_high,
                                      n_pad, &ssp[0], &csp[0],
                                      quit, process, info);
            if(*quit) { RP_SIFT_TOC(t_sift); return(total); }
      } } }
  } }

  RP_SIFT_TOC(t_sift);
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

static long find_points_work_1(ratpoints_args *args,
                 int process(long, long, const mpz_t, void*, int*), void *info);

long find_points_work(ratpoints_args *args,
                 int process(long, long, const mpz_t, void*, int*), void *info)
{
  /* The input fields of args stay what the caller set them to.  The search
   * normalises them and, where they say "choose", used to write its choice
   * into them, which made a caller that fills args once and then loops over
   * curves run every curve after the first with the first one's choices.  So
   * they are saved here and put back on the way out; what was chosen is
   * reported in sp1_used, sp2_used and sp3_used.  Deliberately left as the
   * search made them: cof and degree, which describe the polynomial the
   * search worked with -- reversed (RATPOINTS_REVERSED says so) or with
   * leading zero coefficients dropped (no flag reports that; the caller
   * supplied them) --, num_inter and domain[], which on return describe the
   * region that was actually searched (the intervals given, intersected
   * with the positivity region of f; the program prints them), and the flag
   * bits that report on the run. */
  long b_low = args->b_low, b_high = args->b_high;
  long sp1 = args->sp1, sp2 = args->sp2;
  long array_size = args->array_size, sturm = args->sturm;
  long num_primes = args->num_primes, max_forbidden = args->max_forbidden;
  long result;

  /* sp3 is a working field the caller never sets, and the three _used
   * fields are outputs: they stay 0 when the search ends before it has
   * chosen its primes (no real points, nothing admissible mod 16, ...). */
  args->sp3 = 0;
  args->sp1_used = 0; args->sp2_used = 0; args->sp3_used = 0;
  result = find_points_work_1(args, process, info);
  args->b_low = b_low; args->b_high = b_high;
  args->sp1 = sp1; args->sp2 = sp2;
  args->array_size = array_size; args->sturm = sturm;
  args->num_primes = num_primes; args->max_forbidden = max_forbidden;
  return(result);
}

static long find_points_work_1(ratpoints_args *args,
                 int process(long, long, const mpz_t, void*, int*), void *info)
{
  long total = 0;       /* total counts the points */
  int quit = 0;
  /* Whether the caller left the number of primes to us.  If it did,
   * sieving_info may look past RATPOINTS_DEFAULT_NUM_PRIMES for the curves
   * that need it; an explicit num_primes is a hard limit. */
  int np_is_default = (args->num_primes < 0);
  mpz_t *c = args->cof;
  long degree = args->degree;
  long height = args->height;
  mpz_t *work = args->work;

  int point_at_infty = 0; /* indicates if there are points at infinity */
  int sturm_empty = 0;    /* the positivity region misses the search domain */
  int lcfsq = mpz_perfect_square_p(c[degree]);

  forbidden_entry *forb_ba = (forbidden_entry *)args->forb_ba;
  forbidden_val *forbidden = (forbidden_val *)args->forbidden;
    /* The forbidden divisors, two zero-terminated arrays: primes that are
       tested for with bit arrays, and primes with a set of valuations that
       are tested for by division.  Used when the degree is even. */

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

  /* the counts that say what the sieve actually did, which the choice of
   * primes is corrected from as the run goes on */
  args->n_words = 0; args->n_arrays = 0; args->n_bits = 0;
  args->n_coprime = 0; args->n_checks = 0; args->n_sifts = 0;
  args->n_words_2 = 0;

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
  /* a negative sp1 or sp2 means "choose it from the curve", which is done
   * in sieving_info() once the densities of the primes are known */

  if(args->num_primes > RATPOINTS_NUM_PRIMES)
  { args->num_primes = RATPOINTS_NUM_PRIMES; }
  /* find_points_init sized ba_buffer for the default number of primes; if this
   * call wants more, enlarge it.  Nothing points into the buffer at this
   * moment -- it is a bump allocator and ba_next was reset above -- so it can
   * simply be replaced. */
  if(args->num_primes > args->ba_buffer_primes)
  { free(args->ba_buffer_na);
    alloc_ba_buffer(args, args->num_primes);
  }
  if(args->sp2 > args->num_primes) { args->sp2 = args->num_primes; }
  if(args->sp2 >= 0 && args->sp1 > args->sp2) { args->sp1 = args->sp2; }

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
  { args->flags |= RATPOINTS_NO_REVERSE_AUTO; }
    /* a private bit, cleared on entry: setting RATPOINTS_NO_REVERSE itself
     * would latch the caller's input for the next call */

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
    if(args->sp2 < 0)
    { printf("  number of primes for sieving:     (chosen from the curve)\n"); }
    else { printf("  number of primes for sieving:     %3ld\n", args->sp2); }
    if(args->sp1 < 0)
    { printf("  number of primes for first stage: (chosen from the curve)\n"); }
    else { printf("  number of primes for first stage: %3ld\n", args->sp1); }
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
    if(args->flags & (RATPOINTS_NO_REVERSE | RATPOINTS_NO_REVERSE_AUTO))
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
  if(!((args->flags) & (RATPOINTS_NO_REVERSE | RATPOINTS_NO_REVERSE_AUTO)))
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
        for(n = 1; n <= degree>>1; n++)
        { mpz_set(work[0], c[n]);
          mpz_set(c[n], c[degree+1-n]);
          mpz_set(c[degree+1-n], work[0]);
        }
        if(args->flags & RATPOINTS_VERBOSE)
        { printf("  odd degree, leading coefficient not +/-1, zero "
                 "constant term, coefficient of x is +/-1 ==> reverse\n");
          printf("  coefficients now:");
          { long n;

            for(n = 0; n <= degree; n++)
            { printf(" "); mpz_out_str(NULL, 10, args->cof[n]); }
            printf("\n\n");
          }
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
  if(args->domain == NULL) { return(RATPOINTS_BAD_ARGS); }
  if(args->num_inter == 0)
  /* default interval (effectively ]-infty,infty[) if none is given */
  { args->domain[0].low = -height; args->domain[0].up = height;
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
      { printf("  polynomial is negative on the whole search region"
               " ==> no affine points\n\n"); }
      else
      { long n;

        printf("  can restrict to the following intervals:\n  ");
        for(n = 0; n < args->num_inter; n++)
        { printf("[%lf, %lf] ", args->domain[n].low, args->domain[n].up); }
        printf("\n\n");
      }
    }
    if(ret < 0) { return(RATPOINTS_NON_SQUAREFREE); }
    if(ret == 0)
    { /* No real point in the search region: either f is negative everywhere
       * (then the degree is even and the leading coefficient negative, so
       * there is no point at infinity either) or the positivity region does
       * not meet [-H, H].  In the second case the points at infinity are
       * still there -- an odd degree always has one, an even degree with a
       * square leading coefficient has two, and after a reversal that point
       * is an affine point of the caller's curve -- so only the sieve is
       * skipped, below, after those points have been dealt with. */
      if(!((degree & 1) || lcfsq)) { return(0); }
      sturm_empty = 1;
    }
  }

  /* Point(s) at infinity? */
  if((degree & 1) || lcfsq)
  { point_at_infty = 1;
    if(args->flags & RATPOINTS_VERBOSE)
    { printf("There are points at infinity\n\n"); }
  }
  /* The tests on the denominator.  Odd degree has its own (use_squares);
   * for even degree the Jacobi symbol test wants the leading coefficient
   * not to be a square, and the p-adic test of forbidden_valuations does
   * not mind, so a square leading coefficient keeps the checked loop as a
   * possibility.  sieving_info switches it off if nothing comes of it. */
  if(degree & 1) { args->flags &= ~RATPOINTS_CHECK_DENOM; }
  else if(!lcfsq && !(args->flags & RATPOINTS_NO_JACOBI))
  { args->flags |= RATPOINTS_USE_JACOBI; }

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

  if(den_bits == 0 && !point_at_infty)
  { /* No residue class of the denominator mod 16 admits any numerator, so
     * there is no affine point.  There is none at infinity either: an odd
     * degree always leaves the class v_2(b) >= 4 admissible, so the degree
     * is even, and the leading coefficient then is not a square mod 16,
     * hence not a square (the test on point_at_infty only says so). */
    if(args->flags & RATPOINTS_VERBOSE)
    { printf("  no denominator admits a numerator mod 16 ==> no points\n\n"); }
    return(total);
  }

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
  { /* The mean number of bits set in one word of a bit-array on entry to
     * the sieve.  num_bits[b] holds the admissible numerators for
     * denominators b mod 16, as one word repeated through the bit-array, so
     * the population count of that word is what is wanted; the denominators
     * with no admissible numerator at all are skipped, so the mean is taken
     * over the non-zero entries only.
     * Per word rather than per bit-array on purpose: measurements across
     * register widths show that the survivor rate at the best sp1 is
     * constant per word, not per bit-array (see RATPOINTS_SURVIVORS_PER_WORD
     * in ratpoints.h). */
    double bits_per_word = 0.0;
    { long i, nz = 0, tot = 0;

      for(i = 0; i < 16; i++)
      { long c = __builtin_popcountl(EXT0(num_bits[i]));

        if(c) { tot += c; nz++; }
      }
      if(nz) { bits_per_word = (double)tot/(double)nz; }
    }
    { long ret = sieving_info(args, use_c_long, &c_long[0], sieve_list,
                              bits_per_word, np_is_default,
                              which_bits, den_bits, &num_bits[0]);

    if(ret)
    {

#ifdef DEBUG
      printf("\nno points mod p = %ld ==> return(0)\n", ret);
#else
      if(args->flags & RATPOINTS_VERBOSE)
      { printf("  no points mod p = %ld ==> no rational points\n\n", ret); }
#endif

      return(0);
  } } }

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
    printf("\n  use %ld primes for third stage:\n   ", args->sp3 - args->sp2);
    for( ; n < args->sp3; n++)
    { printf(" %ld", sieve_list[n]->p); }
    printf("\n  one exact check is put at %.0f cycles, %.2f times what it"
           " costs\n    on the curves the third stage was tuned on\n",
           args->check_rel*RATPOINTS_CHECK_REFERENCE, args->check_rel);
    if(args->flags & RATPOINTS_CHECK_DENOM)
    { forbidden_entry *fb = forb_ba;
      forbidden_val *fd = forbidden;

      printf("  denominators excluded:\n   ");
      for( ; fb->p; fb++) { printf(" %ld|b", fb->p); }
      for( ; fd->p; fd++)
      { long p = fd->p, m, mmax, first, pm;

        /* mmax: the largest valuation a denominator can have at p;
         * first: where the tail of excluded valuations up to mmax begins */
        for(mmax = 1, pm = p; pm <= args->b_high/p; mmax++) { pm *= p; }
        for(first = mmax + 1; first > 1 && ((fd->mask >> (first-1)) & 1);
            first--) {}
        for(m = 1; m < first; m++)
        { if((fd->mask >> m) & 1) { printf(" v_%ld(b)=%ld", p, m); } }
        if(first <= mmax) { printf(" %ld^%ld|b", p, first); }
      }
      if(args->flags & RATPOINTS_USE_JACOBI)
      { printf(" (lcf/b) = -1"); }
      printf("\n");
    }
    printf("\n");
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
  if(sturm_empty) { return(total); } /* nothing left to sieve */

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
        long bp_list[args->sp3_max > 0 ? args->sp3_max : 1];
          /* sp3_max, not sp3: adapt_primes may reach for a
           * further prime as the run goes on */
        long last_b = args->b_low;

#ifdef DEBUG
        printf("\n  using squares\n");
        fflush(NULL);
#endif

        { long n;

          for(n = 0; n < args->sp3; n++)
          { bp_list[n] = mod(args->b_low, sieve_list[n]->p); }
          args->sp3_valid = args->sp3;
        }

        for(b = 1; bb = b*b, bb <= args->b_high; b++)
        { if(bb >= args->b_low)
          { ratpoints_bit_array bits = num_bits[bb & 0xf];

            if(TEST(bits))
            { long n;
              long d = bb - last_b;

              /* fill bp_list, after any correction to how many primes
               * the sieve is using (see adapt_primes): one just brought into
               * play has no entry yet and is set from the denominator. */
              if(args->n_words >= args->adapt_at) { adapt_primes(args); }
              RP_BP_TIC(t_bp);
              { long nv = (args->sp3_valid < args->sp3) ? args->sp3_valid
                                                        : args->sp3;

                for(n = 0; n < nv; n++)
                { bp_list[n] = mod(bp_list[n] + d, sieve_list[n]->p); }
                for(n = nv; n < args->sp3; n++)
                { bp_list[n] = mod(bb, sieve_list[n]->p); }
                args->sp3_valid = args->sp3;
              }
              RP_BP_TOC(t_bp, args->sp3);
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
        long bp_list[args->sp3_max > 0 ? args->sp3_max : 1];
          /* sp3_max, not sp3: adapt_primes may reach for a
           * further prime as the run goes on */

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

            for(n = 0; n < args->sp3; n++)
            { bp_list[n] = mod(*div, sieve_list[n]->p); }
            args->sp3_valid = args->sp3;
          }

          for(b = 1; bb = (*div)*b*b, bb <= args->b_high; b++)
          { if(bb >= args->b_low)
            { int flag = 1;
              ratpoints_bit_array bits = num_bits[bb & 0xf];

              if(EXT0(bits))
              { long i;
                long n;
                long d = bb - last_b;

                /* fill bp_list; see the note at the same place above */
                if(args->n_words >= args->adapt_at) { adapt_primes(args); }
                RP_BP_TIC(t_bp);
                { long nv = (args->sp3_valid < args->sp3) ? args->sp3_valid
                                                          : args->sp3;

                  for(n = 0; n < nv; n++)
                  { bp_list[n] = mod(bp_list[n] + d, sieve_list[n]->p); }
                  for(n = nv; n < args->sp3; n++)
                  { bp_list[n] = mod(bb, sieve_list[n]->p); }
                  args->sp3_valid = args->sp3;
                }
                RP_BP_TOC(t_bp, args->sp3);
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
      { forbidden_val *forb;
        long b;
        long bp_list[args->sp3_max > 0 ? args->sp3_max : 1];
          /* sp3_max, not sp3: adapt_primes may reach for a
           * further prime as the run goes on */
        long last_b = args->b_low;
        unsigned long b_bits;

#ifdef DEBUG
        printf("\n  taking account of forbidden divisors of the denominator\n");
        fflush(NULL);
#endif

        { long n;

          for(n = 0; n < args->sp3; n++)
          { bp_list[n] = mod(args->b_low, sieve_list[n]->p); }
          args->sp3_valid = args->sp3;
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
          { /* check if denominator is excluded: is v_p(b) one of the
             * valuations the entry for p forbids? */
            for(forb = &forbidden[0];
                forb->p && !((forb->mask >> valuation1(b, forb->p)) & 1);
                forb++) {};

#ifdef DEBUG
            if(forb->p)
            { printf("\nb = %ld: excluded, v_%ld(b) = %ld\n",
                     b, forb->p, valuation1(b, forb->p));
              fflush(NULL);
            }
#endif

            if(forb->p == 0
                && (!(args->flags & RATPOINTS_USE_JACOBI)
                      || (use_c_long
                           ? jacobi1(b, c_long[degree])
                           : jacobi(b, work[0], c[degree])) == 1))
            { long n;
              long d = b - last_b;

              /* fill bp_list; see the note at the same place above */
              if(args->n_words >= args->adapt_at) { adapt_primes(args); }
              RP_BP_TIC(t_bp);
              { long nv = (args->sp3_valid < args->sp3) ? args->sp3_valid
                                                        : args->sp3;

                for(n = 0; n < nv; n++)
                { long bp = bp_list[n] + d;
                  long p = sieve_list[n]->p;

                  while(bp >= p) { bp -= p; }
                  bp_list[n] = bp;
                }
                for(n = nv; n < args->sp3; n++)
                { bp_list[n] = mod(b, sieve_list[n]->p); }
                args->sp3_valid = args->sp3;
              }
              RP_BP_TOC(t_bp, args->sp3);
              last_b = b;

              total += sift(b, survivors, args, which_bits, bits,
                            sieve_list, &bp_list[0],
                            &quit, process, info);
              if(quit) { break; }
            }

#ifdef DEBUG
            else
            { if(forb->p == 0)
              { printf("\nb = %ld: excluded by Jacobi symbol\n", b);
                fflush(NULL);
            } }
#endif

          }
        }
      } /* if(args->flags & RATPOINTS_CHECK_DENOM) */
      else
      { long b;
        long bp_list[args->sp3_max > 0 ? args->sp3_max : 1];
          /* sp3_max, not sp3: adapt_primes may reach for a
           * further prime as the run goes on */
        long last_b = args->b_low;

        { long n;

          for(n = 0; n < args->sp3; n++)
          { bp_list[n] = mod(args->b_low, sieve_list[n]->p); }
          args->sp3_valid = args->sp3;
        }

        for(b = args->b_low; b <= args->b_high; b++)
        { ratpoints_bit_array bits = num_bits[b & 0xf];

          if(EXT0(bits))
          { long n;
            long d = b - last_b;

            /* fill bp_list; see the note at the same place above */
            if(args->n_words >= args->adapt_at) { adapt_primes(args); }
            RP_BP_TIC(t_bp);
            { long nv = (args->sp3_valid < args->sp3) ? args->sp3_valid
                                                      : args->sp3;

              for(n = 0; n < nv; n++)
              { long bp = bp_list[n] + d;
                long p = sieve_list[n]->p;

                while(bp >= p) { bp -= p; }
                bp_list[n] = bp;
              }
              for(n = nv; n < args->sp3; n++)
              { bp_list[n] = mod(b, sieve_list[n]->p); }
              args->sp3_valid = args->sp3;
            }
            RP_BP_TOC(t_bp, args->sp3);
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

#if defined(RP_PRIME_STATS) && defined(RP_PHASE_TIMING)
  /* Development instrumentation: what run_shape predicted for this curve
   * against what the run actually did.  Needs both switches, since the
   * counters it reads belong to the phase timing. */
  { static unsigned long long last_arrays = 0, last_dens = 0;

    fprintf(stderr, "[runshape] Upred=%.6g Uact=%.6g Dpred=%.6g Dact=%.6g"
            " words=%lu arrays=%lu bits=%lu coprime=%lu checks=%lu\n",
            args->run_words,
            (double)(_rp_arrays_swept - last_arrays)*(double)RBA_PACK,
            args->run_denoms, (double)(_rp_bp_dens - last_dens),
            args->n_words, args->n_arrays, args->n_bits,
            args->n_coprime, args->n_checks);
    last_arrays = _rp_arrays_swept; last_dens = _rp_bp_dens;
  }
#endif

  /* report the primes the sieve used (adapt_primes may have moved sp2 and
   * sp3 during the run; find_points_work zeroed these on entry) */
  args->sp1_used = args->sp1; args->sp2_used = args->sp2;
  args->sp3_used = args->sp3;

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
      RP_BC_TIC(t_bc);

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
      RP_BC_TOC(t_bc);
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
