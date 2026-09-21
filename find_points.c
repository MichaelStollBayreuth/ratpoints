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
/* and the loop that computes b modulo each modulus of the first two phases,
 * once per denominator and modulus; the primes of the third stage have no
 * such cost, since their set-up is done on demand (fill_checks in sift.c) */
extern unsigned long long _rp_bp_cycles, _rp_bp_dens, _rp_bp_steps;
extern unsigned long long _rp_arrays_swept;
/* and building one sieve table, which is the fixed cost a modulus has to
 * earn back over the run; see run_shape */
extern unsigned long long _rp_init_cycles, _rp_init_calls, _rp_init_rows;
extern unsigned long long _rp_setup_cycles, _rp_setup_dens;
/* and the phase counters themselves, for the per-curve line at the end of
 * find_points_work */
extern unsigned long long _rp_phase1_cycles, _rp_phase2_cycles;
extern unsigned long long _rp_check_cycles, _rp_sift0_calls;
#ifdef RP_PHASE_COUNTS
extern unsigned long long _rp_and2;
#endif
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

/* A candidate modulus in the ranking: its density r, its key and the cost
 * the key was made from (what the modulus costs per word swept in the
 * first phase, one AND plus its fixed costs spread over the run; the rule
 * that ends the first phase weighs it, see take_entries), the modulus p,
 * the primes it involves as a mask over their indices (a modulus involving
 * one of the primes beyond the first 64 is not offered), its entry -- a
 * prime's from examine_prime(), a composite's made by make_modulus() when
 * the ranking takes it and NULL until then -- and its prime-power factors
 * as codes: pn for the prime prime[pn], RATPOINTS_NUM_PRIMES + i for the
 * i-th prime power examined (rp_power).  nf is the number of factors of a
 * composite modulus, 0 for a prime, and -1 for a prime that says only that
 * numerator and denominator are not both divisible by it (examine_prime). */
typedef struct { double r; double key; double cost; long p; unsigned long mask;
                 ratpoints_sieve_entry *ssp; short nf; short fac[RP_MAX_FACTORS]; }
        entry;

/* a bound on the number of composite moduli offered: at most the odd
 * numbers below the limit */
#define RP_MAX_MODULI (RATPOINTS_MAX_PRIME_EVEN/2)

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
/* Reserve the space for the sieving information for the first n primes,
 * and for extra further bit arrays (the tables of the composite moduli the
 * ranking has taken, see ensure_ba_buffer).  For each prime p we may need p
 * arrays (one for each denominator mod p) of length p + RATPOINTS_CHUNK-1,
 * the CHUNK-1 so that _ratpoints_sift0 can avoid a wrap-around.  The sum
 * grows with the square of the largest prime, which is why n matters:
 * taking it to be RATPOINTS_NUM_PRIMES asks for 5 MB when primes go up to
 * 127, but 1.7 GB when they go up to 1021, whether or not the large primes
 * are ever looked at.
 * args->ba_buffer_na keeps the address malloc returned, so that it can be
 * freed later; the +1 leaves the leeway needed for the alignment.
 * args->ba_buffer_primes and args->ba_buffer_arrays record what the block
 * is good for.
 */
static long ba_buffer_need(long n, long extra)
{ long need = extra;
  long i;

  for(i = 0; i < n; i++) { need += prime[i]*(prime[i] + RATPOINTS_CHUNK-1); }
  return(need);
}

static void alloc_ba_buffer(ratpoints_args *args, long n, long extra)
{ long need = ba_buffer_need(n, extra);

  args->ba_buffer_na = malloc((need+1)*sizeof(ratpoints_bit_array));
  args->ba_buffer = pointer_align(args->ba_buffer_na, sizeof(ratpoints_bit_array));
  args->ba_next = args->ba_buffer;
  args->ba_buffer_primes = n;
  args->ba_buffer_arrays = need;
}

/* Make sure the block holds the tables of the first n primes and extra
 * further bit arrays (the tables of composite moduli, TODO item 21).  Only
 * to be called while nothing has been taken from the block: it is a bump
 * allocator, reset at the start of every curve, and the tables are built
 * lazily while sieving, so before the first denominator it can simply be
 * replaced. */
static void ensure_ba_buffer(ratpoints_args *args, long n, long extra)
{ if(ba_buffer_need(n, extra) > args->ba_buffer_arrays
     || n > args->ba_buffer_primes)
  { free(args->ba_buffer_na);
    alloc_ba_buffer(args, n, extra);
  }
}

/* How many prime powers p^e, e >= 2, lie below RATPOINTS_MAX_PRIME_EVEN: the
 * size of the buffer for what examine_power() records about them. */
static long num_powers(void)
{ long n = 0, pn;

  for(pn = 0; pn < RATPOINTS_NUM_PRIMES; pn++)
  { long p = prime[pn], m = p*p;

    if(m >= RATPOINTS_MAX_PRIME_EVEN) { break; }
    for( ; m < RATPOINTS_MAX_PRIME_EVEN; m *= p) { n++; }
  }
  return(n);
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

  /* allocate space for se_buffer: an entry for every prime looked at, for
   * every composite modulus the ranking takes (at most one per prime it
   * involves, so at most as many again) and for every prime power that is a
   * factor of one (TODO item 21) */
  args->se_buffer
    = (ratpoints_sieve_entry *) malloc((2*RATPOINTS_NUM_PRIMES + num_powers())
                                        * sizeof(ratpoints_sieve_entry));
  args->se_next = args->se_buffer;
  /* and for what the prime powers say about the curve */
  args->pw_buffer = malloc((num_powers() > 0 ? num_powers() : 1)
                             * sizeof(rp_power));

  /* allocate space for ba_buffer, for the first RATPOINTS_DEFAULT_NUM_PRIMES
   * primes; find_points_work enlarges it should a caller ask for more.  It
   * cannot be sized from args->num_primes here, because the documented way of
   * using the library sets that field between find_points_init and
   * find_points_work, not before.
   */
  alloc_ba_buffer(args, (RATPOINTS_DEFAULT_NUM_PRIMES < RATPOINTS_NUM_PRIMES)
                          ? RATPOINTS_DEFAULT_NUM_PRIMES
                          : RATPOINTS_NUM_PRIMES, 0);

  /* allocate space for int_buffer */
  args->int_buffer
    = malloc(RATPOINTS_NUM_PRIMES*(RATPOINTS_MAX_PRIME+1)*sizeof(int));
  args->int_next = args->int_buffer;

  /* allocate sieve_list: the moduli taken and the primes left over (at most
   * RATPOINTS_NUM_PRIMES entries in all, since each owns a prime no other
   * has; the doubling is headroom) */
  args->sieve_list = malloc(2*RATPOINTS_NUM_PRIMES
                             * sizeof(ratpoints_sieve_entry*));

  /* and the third stage's working copy of what it needs per denominator.
   * It lives here rather than on sift()'s stack because that function is
   * entered once per denominator, and enlarging its frame by this much was
   * measured to cost several per cent all by itself. */
  args->stage3_list = malloc(RATPOINTS_NUM_PRIMES * sizeof(check_spec));

  /* the reciprocals _ratpoints_sift0 reduces word numbers with.  They belong
   * to the moduli, not to the denominators, so they are filled in once per
   * curve; and they are kept out of sieve_spec because that structure is read
   * in the innermost loop of the first phase, where its size tells.  As many
   * as sieve_list has (RATPOINTS_NUM_PRIMES would do, since every entry in
   * the list owns a prime no other has; the doubling is headroom). */
  args->magics = malloc(2*RATPOINTS_NUM_PRIMES * sizeof(unsigned long));

  /* allocate remaining data structures */
  args->den_info = malloc((PRIMES1000+2)*sizeof(use_squares1_info));
  args->divisors = malloc((MAX_DIVISORS+1)*sizeof(long));
  args->forb_ba = malloc((PRIMES1000 + 1)*sizeof(forbidden_entry));
  args->forbidden = malloc((PRIMES1000 + 1)*sizeof(forbidden_val));
  /* the bit patterns for forbidden divisors beyond the compiled table of
   * primes are built per curve, in a buffer that grows as needed; see
   * sieving_info */
  args->forb_words = NULL; args->forb_words_len = 0;

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
  free(args->pw_buffer);
  free(args->den_info);
  free(args->divisors);
  free(args->forb_ba);
  free(args->forbidden);
  free(args->forb_words);

  /* clear pointer in args */
  args->work = NULL; args->work_length = 0;
  args->se_buffer = NULL; args->se_next = NULL;
  args->ba_buffer_na = NULL; args->ba_buffer_primes = 0;
  args->ba_buffer_arrays = 0; args->pw_buffer = NULL;
  args->ba_buffer = NULL; args->ba_next = NULL;
  args->int_buffer = NULL; args->int_next = NULL;
  args->sieve_list = NULL; args->stage3_list = NULL;
  args->magics = NULL;
  args->den_info = NULL; args->divisors = NULL;
  args->forb_ba = NULL; args->forbidden = NULL;
  args->forb_words = NULL; args->forb_words_len = 0;

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
 * Helper function: the least k >= 1 with k^2 >= n, for n >= 1 -- where   *
 * the loops over the square denominators start.  sqrt in double is exact *
 * only up to 2^53 and n may exceed that, so the result is corrected by a *
 * step either way; the comparisons avoid forming k^2, which could        *
 * overflow near LONG_MAX:  k^2 < n  <==>  k <= (n-1)/k.                  *
 *************************************************************************/

static long ceil_sqrt(long n)
{
  long k = (long)sqrt((double)n);

  while(k > 1 && k - 1 > (n - 1)/(k - 1)) { k--; }  /* (k-1)^2 >= n */
  while(k <= (n - 1)/k) { k++; }                    /* k^2 < n */
  return(k);
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

/* The Jacobi symbol test on the denominators without the Jacobi symbol.
 *
 * What the test asks of a denominator b is that (lcf/b*) = 1, where b* is b
 * with the prime factors of 2*lcf taken out; jacobi1 above computes exactly
 * that, by a binary gcd-like loop of some two hundred instructions with half
 * a dozen mispredicted branches, for about a quarter of all denominators.
 * Written out with lcf = +-2^v * prod q_i^(e_i), the symbol is
 *
 *   (-1/b*)^neg * (2/b*)^v * prod_i (q_i/b*)^(e_i) ,
 *
 * where the first two factors depend on b* mod 8 alone, and quadratic
 * reciprocity turns (q_i/b*) into the Legendre symbol (b* mod q_i / q_i)
 * times a sign that depends on b* mod 4.  So the test needs a table of the
 * non-squares modulo each odd prime factor of lcf with an odd exponent, and
 * a table of eight signs; the primes with an even exponent contribute
 * nothing, but still have to be taken out of b.  That makes the test one
 * multiplication, one table look-up and one exclusive or per prime.
 *
 * jacobi_setup prepares this for one curve.  It applies when every odd prime
 * factor of lcf is in prime[] (trial division finds them) and the
 * denominators stay below 2^32 (the reductions multiply by a reciprocal, see
 * RP_MULDIV); otherwise it says so, and the denominator loop calls jacobi1
 * or jacobi as before.  A leading coefficient that fits a long always fits
 * the tables: distinct odd primes below 1024 with a product below 2^63 sum
 * to at most 6057 (the six largest and a 7), and only those with an odd
 * exponent need a table. */
#define RP_JACOBI_PRIMES 16   /* odd prime factors of lcf, at most */
#define RP_JACOBI_TABLE 8192  /* bytes of non-square tables, at most */

typedef struct { long nq;                   /* the odd primes of lcf */
                 long q[RP_JACOBI_PRIMES];
                 unsigned long magic[RP_JACOBI_PRIMES];
                 const unsigned char *nonsq[RP_JACOBI_PRIMES];
                   /* nonsq[i][r] = 1 iff r is a non-square modulo q[i];
                      NULL when the exponent of q[i] is even */
                 unsigned char sign[8];
                   /* sign[b mod 8] = 1 iff the factors (-1/b), (2/b) and
                      the reciprocity signs multiply to -1 */
               } jacobi_info;

static int jacobi_setup(jacobi_info *ji, unsigned char *tab, long tab_len,
                        const mpz_t lcf, mpz_t tmp, long b_high)
{ long i, v, n3 = 0, used = 0;
  int neg = (mpz_sgn(lcf) < 0);

  if(b_high > RP_MULMOD_LIMIT) { return(0); }
  ji->nq = 0;
  mpz_abs(tmp, lcf);
  v = mpz_scan1(tmp, 0);
  mpz_fdiv_q_2exp(tmp, tmp, v);
  for(i = 0; i < PRIMES1000 && mpz_cmp_ui(tmp, 1) != 0; i++)
  { long q = prime[i], e = 0; /* prime[] holds the odd primes */

    if(!mpz_divisible_ui_p(tmp, q)) { continue; }
    do { mpz_divexact_ui(tmp, tmp, q); e++; } while(mpz_divisible_ui_p(tmp, q));
    if(ji->nq == RP_JACOBI_PRIMES) { return(0); }
    ji->q[ji->nq] = q;
    ji->magic[ji->nq] = ULONG_MAX/(unsigned long)q + 1;
    if(e & 1)
    { unsigned char *ns = tab + used;
      long x, s;

      if(used + q > tab_len) { return(0); }
      used += q;
      /* mark the squares x^2 mod q for 0 < x < q/2, which are all of them,
       * stepping from x^2 to (x+1)^2 by 2x+1; entry 0 is never looked at */
      for(x = 0; x < q; x++) { ns[x] = 1; }
      for(x = 1, s = 1; x <= q/2; x++)
      { ns[s] = 0;
        s += 2*x + 1; if(s >= q) { s -= q; }
      }
      ji->nonsq[ji->nq] = ns;
      if(q & 2) { n3++; } /* q = 3 mod 4: reciprocity brings a sign */
    }
    else { ji->nonsq[ji->nq] = NULL; }
    ji->nq++;
  }
  if(mpz_cmp_ui(tmp, 1) != 0) { return(0); } /* a prime beyond the table */
  for(i = 0; i < 8; i++)
  { int s = 0;

    /* (-1/b) = -1 and (q/b) = -(b/q) for q = 3 mod 4, both iff b = 3 mod 4 */
    if((i & 3) == 3) { s ^= (n3 + neg) & 1; }
    /* (2/b) = -1 iff b = 3 or 5 mod 8 */
    if(i == 3 || i == 5) { s ^= v & 1; }
    ji->sign[i] = s;
  }
  return(1);
}

/* Is (lcf/b*) = 1 ?  b > 0 is a denominator below 2^32. */
static inline int jacobi_test(long b, const jacobi_info *ji)
{ long i;
  unsigned char s;

  b >>= RP_CTZL((unsigned long)b); /* the odd part */
  for(i = 0; i < ji->nq; i++)      /* take the primes of lcf out */
  { long q = ji->q[i];
    unsigned long m = ji->magic[i];

    for(;;)
    { long k = RP_MULDIV(b, q, m); /* b/q rounded down */

      if(b - k*q) { break; }
      b = k;
    }
  }
  s = ji->sign[b & 7];
  for(i = 0; i < ji->nq; i++)
  { if(ji->nonsq[i])
    { long q = ji->q[i];

      s ^= ji->nonsq[i][b - RP_MULDIV(b, q, ji->magic[i])*q];
    }
  }
  return(s == 0);
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

/* Which numerators a denominator can have is decided modulo 64.  With D
 * the degree of f rounded up to an even number, F(a,b) = b^D f(a/b) is a
 * form of even degree D with integer coefficients, so F(a,b) mod 64 depends
 * on a and b mod 64 only, and a point needs F(a,b) to be a square: one of
 * the twelve residues
 *   0, 1, 4, 9, 16, 17, 25, 33, 36, 41, 49, 57
 * mod 64 -- bit s of RP_SQUARES64 is set for these s.  Modulo 16 the
 * squares are four residues of sixteen, modulo 64 twelve of sixty-four; the
 * finer modulus halves what is accepted when F = 0 mod 4 (an odd F is a
 * square mod 64 exactly when it is one mod 8, so nothing changes there).
 * The information costs nothing in the sieve -- a pattern of period 64 is
 * one repeated word in every packing, just as one of period 16 was -- and
 * 64 is where it stops paying: modulo 256 the squares are 44 of 256, which
 * would remove one admissible class in twelve more -- a quarter of what
 * going from 16 to 64 removes -- for sixteen times the set-up.
 *
 * For an odd denominator, b^D is the square of a unit, so F(a,b) is a
 * square mod 64 exactly when f(a b^-1 mod 64) is.  For an even one the
 * numerator is odd, a^D is the square of a unit, and F(a,b) = a^D frev(b/a)
 * with frev(t) = t^D f(1/t) = c_0 t^D + c_1 t^(D-1) + ... + c_D (c_D = 0
 * when the degree is odd), so F(a,b) is a square exactly when
 * frev(b a^-1 mod 64) is.  Two words thus decide everything: bit k of fsq
 * says that f(k) is a square mod 64, bit t of gsq that frev(t) is.  The
 * rule is exact for every class of the denominator; the hand-derived
 * congruences for f(odd/2), f(odd/4), ... of the mod-16 code are not
 * needed.
 *
 * What the two words say is then put in the form the sieve uses, per class
 * of the denominator mod 64 (rp_num_class in rp-private.h).  The admissible
 * numerators of a class are a union of residue classes mod 64, and where
 * they lie in a single class a0 mod 2^k the bit arrays hold only the
 * numerators a0 + 2^k t: the odd denominators of a random curve take their
 * numerators from one class mod 4 on a fifth of the curves and from one
 * class mod 8 on a sixth, and the two-fold packing that "only odd
 * numerators" used to be is the case k = 1.  k and a0 are read off the
 * pattern.  The theory says that k is the same for every odd class, and for
 * every class of the same 2-adic valuation of b -- the pattern of b = u b'
 * with u a unit mod 64 is u times that of b', and multiplication by a unit
 * keeps a set of residues within one class mod 2^k -- which the verbose
 * report relies on; the sieve does not. */
#define RP_SQUARES64 0x0202021202030213UL

static void get_2adic_info(ratpoints_args *args, unsigned long *den_bits,
                           rp_num_class *cls)
{
  mpz_t *c = args->cof;
  long degree = args->degree;
  long D = degree + (degree & 1); /* the degree of the form F(a,b) */
  unsigned long cm[D + 1];        /* the coefficients of f mod 64, cm[D] = 0
                                   * for odd degree */
  unsigned long fsq = 0UL;        /* bit k set: f(k) is a square mod 64 */
  unsigned long gsq = 0UL;        /* bit t set, t even: frev(t) is one */
  unsigned long nb[64];           /* bit a of nb[b]: the numerator a is
                                   * admissible for the denominators b mod 64 */
  unsigned long db = 0UL;         /* the denominator classes with a pattern */
  long b, k;

#ifdef DEBUG
  printf("\nget_2adic_info: start...\n"); fflush(NULL);
#endif

  /* the coefficients mod 64, as non-negative residues */
  for(k = 0; k <= degree; k++) { cm[k] = mpz_fdiv_ui(c[k], 64); }
  if(degree & 1) { cm[D] = 0UL; }

  /* the two tables by Horner's rule; the arithmetic of unsigned long is
   * exact modulo 2^64, so the residue mod 64 is taken at the end */
  for(k = 0; k < 64; k++)
  { unsigned long s = cm[D];
    long n;

    for(n = D - 1; n >= 0; n--) { s = s*(unsigned long)k + cm[n]; }
    if((RP_SQUARES64 >> (s & 0x3f)) & 1UL) { fsq |= 1UL << k; }
  }
  for(k = 0; k < 64; k += 2)
  { unsigned long s = cm[0];
    long n;

    for(n = 1; n <= D; n++) { s = s*(unsigned long)k + cm[n]; }
    if((RP_SQUARES64 >> (s & 0x3f)) & 1UL) { gsq |= 1UL << k; }
  }

  for(b = 0; b < 64; b++) { nb[b] = 0UL; }

  /* Odd denominators: the numerator a is admissible for b when f(a b^-1)
   * is a square, that is, when a = k b for a k with bit k of fsq set -- so
   * every such k is scattered to a = k b in the pattern of every odd b. */
  { unsigned long w = fsq;

    while(w)
    { k = RP_CTZL(w); w &= w - 1UL;
      for(b = 1; b < 64; b += 2) { nb[b] |= 1UL << ((k*b) & 0x3f); }
    }
  }

  /* Even denominators: the numerator is odd, and a is admissible for b when
   * frev(b a^-1) is a square, that is, when b = t a for a t with bit t of
   * gsq set -- so for every such t and every odd a, bit a of the pattern of
   * b = t a is set.  For t = 0 this is the class b = 0 mod 64, which admits
   * every odd numerator when frev(0) = c_D is a square (always when the
   * degree is odd). */
  { unsigned long w = gsq;

    while(w)
    { long a, t = RP_CTZL(w);

      w &= w - 1UL;
      for(a = 1; a < 64; a += 2) { nb[(t*a) & 0x3f] |= 1UL << a; }
    }
  }

  /* The packing of each class: the largest k with every admissible
   * numerator in one class a0 mod 2^k, and the pattern in that packing --
   * bit t for the numerator a0 + 2^k t, of period 64/2^k in t and so one
   * word repeated.  A class without a pattern gets k = 0, a0 = 0 and an
   * empty pattern; it is never sieved, but its entry is read (bits_per_word,
   * run_shape, the test on every denominator). */
  { /* bit a set for a = 0 mod 2^k, a = 0..63 */
    static const unsigned long class_mask[RP_NUM_STRIDES]
      = {~0UL, 0x5555555555555555UL, 0x1111111111111111UL,
         0x0101010101010101UL, 0x0001000100010001UL,
         0x0000000100000001UL, 0x0000000000000001UL};

    for(b = 0; b < 64; b++)
    { unsigned long w = nb[b];
      unsigned long packed = 0UL;
      long a0 = 0;

      k = 0;
      if(w)
      { long a = RP_CTZL(w); /* the least admissible numerator, which fixes
                              * the class mod 2^k for every k */

        /* the largest k for which every set bit is at a = a0 mod 2^k;
         * upwards, since most classes stop at the first test */
        for(k = 1; k < RP_NUM_STRIDES; k++)
        { if(w & ~(class_mask[k] << (a & ((1L << k) - 1)))) { break; } }
        k--;
        a0 = a & ((1L << k) - 1);
        if(k == 0) { packed = w; }
        else
        { /* the set bits -- at most 64/2^k of them -- each at its place
           * t = (a - a0)/2^k in the packing; a runs over 0..63 and the
           * pattern has period 64/2^k in t, so every set bit lands below
           * that and the word is then replicated */
          while(w)
          { long t;

            a = RP_CTZL(w); w &= w - 1UL;
            t = ((a - a0) & 0x3f) >> k;
            packed |= 1UL << t;
          }
          for(a = LONG_LENGTH >> k; a < LONG_LENGTH; a <<= 1)
          { packed |= packed << a; }
        }
        db |= 1UL << b;
      }
      cls[b].bits = RBA(packed);
      cls[b].k = k;
      cls[b].a0 = a0;
      cls[b].offset = NULL; /* set once the primes are known, class_offsets */
    }
  }
  *den_bits = db;

#ifdef DEBUG
  printf("\nfsq = %016lx, gsq = %016lx, den_bits = %016lx\n", fsq, gsq, db);
  for(b = 0; b < 64; b++)
  { if(nb[b])
    { printf("  b = %2ld mod 64: numerators %ld mod %ld, packed pattern %016lx\n",
             b, cls[b].a0, 1L << cls[b].k, EXT0(cls[b].bits));
  } }
  printf("\nget_2adic_info: done.\n"); fflush(NULL);
#endif
}

/* How many distinct packings (k, a0) the classes with a pattern have: the
 * number of rows class_offsets() builds. */
static long num_packings(const rp_num_class *cls)
{ unsigned long seen[2] = {0UL, 0UL}; /* the keys 2^k + a0, below 128 */
  long b, n = 0;

  for(b = 0; b < 64; b++)
  { long key = (1L << cls[b].k) + cls[b].a0;

    if(EXT0(cls[b].bits) && !((seen[key >> 6] >> (key & 63)) & 1UL))
    { seen[key >> 6] |= 1UL << (key & 63); n++; }
  }
  return(n > 0 ? n : 1);
}

/* The row shifts of the numerator classes (rp_num_class), once the moduli
 * are known: for a class with stride 2^k and offset a0 and the n-th modulus
 * p of sieve_list, a0 (2^k RBA_LENGTH)^-1 mod p, plus the multiple of p that
 * keeps the word number plus the shift non-negative (RP_ROW_BIAS).  Classes
 * with the same packing share a row.  offsets has room for num_packings()
 * rows of np entries, one per modulus that may come to be sieved with
 * (adapt_primes can promote a prime into the second phase during the run).
 * On the way the sieve entries get 2^-k mod p for the strides in use (dinv,
 * read by fill_bp_list), by halving from 1: 2^-1 mod p is (p+1)/2, which
 * works for any odd modulus. */
static void class_offsets(rp_num_class *cls, ratpoints_sieve_entry **sieve_list,
                          long np, long *offsets)
{ /* the row of the packing (k, a0), if built: a0 < 2^k, so 2^k + a0 is a
   * key below 2^RP_NUM_STRIDES */
  long *row_of[1L << RP_NUM_STRIDES];
  long rinv[np > 0 ? np : 1][RP_NUM_STRIDES]; /* (2^k RBA_LENGTH)^-1 mod p */
  long b, n, kmax = 0, rows = 0;

  for(b = 0; b < 64; b++)
  { if(EXT0(cls[b].bits) && cls[b].k > kmax) { kmax = cls[b].k; } }
  for(n = 0; n < np; n++)
  { ratpoints_sieve_entry *se = sieve_list[n];
    long p = se->p, k, d = 1;
    long e = se->rbainv;

    for(k = 0; k <= kmax; k++)
    { se->dinv[k] = d; rinv[n][k] = e;
      d = (d & 1) ? (d + p) >> 1 : d >> 1;
      e = (e & 1) ? (e + p) >> 1 : e >> 1;
    }
  }

  for(b = 0; b < (1L << RP_NUM_STRIDES); b++) { row_of[b] = NULL; }
  for(b = 0; b < 64; b++)
  { long key = (1L << cls[b].k) + cls[b].a0;

    if(!EXT0(cls[b].bits)) { continue; }
    if(row_of[key]) { cls[b].offset = row_of[key]; }
    else
    { long *row = offsets + rows*np;

      for(n = 0; n < np; n++)
      { ratpoints_sieve_entry *se = sieve_list[n];

        row[n] = RP_MULMOD(cls[b].a0*rinv[n][cls[b].k], se->p, se->magic)
                   + se->bias;
      }
      cls[b].offset = row_of[key] = row;
      rows++;
    }
  }
}

/**************************************************************************
 * This is a comparison function needed for sorting in order to determine *
 * the `best' primes for sieving.                                         *
 **************************************************************************/

/* Moduli are ranked by what they say per unit of what they cost, not by
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

/* What one more modulus -- a prime, or a composite modulus -- costs the
 * sieve, per word swept, in units of what a first-phase AND costs there.
 *
 * per_word is the part that is paid for every word (or for every surviving
 * bit array, which comes to the same thing once multiplied by the survival
 * rate): 1 in the first phase, COST_PHASE2*rate in the second, and the
 * third stage's own cost per survivor in the third.  The other two terms are
 * paid once and spread over the run: the sieve table, which the third stage
 * does not build, and the entry of bp_list, which the first two phases pay
 * for every denominator.  (The third stage's primes have neither: their
 * set-up is done on demand, and they are ranked by another rule, below.)
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

/* The rank of a modulus: what it costs divided by what it says.  A modulus
 * multiplies the survival rate by r, so what it says is -log(r), and the
 * best set of moduli for a given total cost is found by taking them in
 * increasing order of this ratio. */
static double prime_key(double r, long p, double per_word, int tabled,
                        double cost_table, double u_words, double n_denoms)
{ double info = -log(r);

  if(info <= 0.0) { return(1.0e300); }
  return(prime_cost(p, per_word, tabled, cost_table, u_words, n_denoms)/info);
}

/* Key a candidate for the first phase, and keep the cost the key was made
 * from: the rule that ends the phase compares it with what the candidate
 * would save (take_entries).  call_cost is what a first-phase modulus pays
 * per word for the calls of the sieve, RATPOINTS_COST_CALL per call spread
 * over the run (the row pointer's reduction at the head of every call and
 * the narrower legs past the last whole chunk; the second phase has
 * neither, it finds its rows directly). */
static void phase_1_key(entry *e, double cost_table, double u_words,
                        double n_denoms, double call_cost)
{ double info = -log(e->r);
  /* the bit arrays a denominator sweeps, and the part of the row it walks
   * and has to fetch from beyond the first-level cache (RATPOINTS_COST_LINE
   * per line of 64 bytes, once per denominator) */
  double arrays = u_words/(n_denoms*(double)RBA_PACK);
  double walk = ((double)e->p < arrays) ? (double)e->p : arrays;

  e->cost = prime_cost(e->p, 1.0, 1, cost_table, u_words, n_denoms)
            + call_cost
            + RATPOINTS_COST_LINE*walk
                *((double)sizeof(ratpoints_bit_array)/64.0)*n_denoms/u_words;
  e->key = (info <= 0.0) ? 1.0e300 : e->cost/info;
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
 * time and nothing else.  Nor can it get ahead of bp_list: fill_bp_list()
 * makes this correction first and then computes every entry the current
 * number of primes asks for.
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
    { /* nothing to add: see whether the last one is still worth having --
       * a composite modulus stays, since the third stage cannot take it,
       * and as this walks down from the top it stops there, so the primes
       * below a composite are not demoted either */
      for(n = sp2 - 1; n > sp1 && sieve_list[n]->nf == 0; n--)
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

/* What a word that survives the first phase costs from there on, relative
 * to what it costs in a long run.  It meets the extra moduli of the second
 * phase, one AND on a surviving word each (RATPOINTS_COST_PHASE2) for as
 * long as it survives them, and what survives those is extracted and
 * tested, RATPOINTS_COST_SURVIVOR per survivor (about one per surviving
 * word at these rates).  In a long run extra is eleven and the second phase
 * kills nearly every survivor cheaply, and the factor is one.  In a run of
 * a few thousand words extra is 0 or 1 and every surviving word reaches the
 * extraction, which costs several times more, so the first phase should
 * sieve harder than the target says.  RATPOINTS_SURVIVORS_PER_WORD is
 * fitted at 16383, where extra is about two: there this factor and the cost
 * factor of phase_1_wants are both about three at the modulus where the
 * phase stops and nearly cancel.  The moduli the second phase would take
 * are approximated by the next entries of the pool, in key order (the phase
 * ranks them again by a key of its own, which favours smaller moduli of
 * higher density: the survivors modelled here die a little faster than the
 * real ones will, so the factor errs low).  The long-run value is an
 * endless second phase of the mean density of the next RATPOINTS_SP2_EXTRA
 * entries -- the second phase of a long run -- whether or not this run has
 * them; when the pool has nothing beyond the modulus in hand, its own
 * density stands in.  The result is floored at one: a surviving word cannot
 * cost less than the long run's, whatever the densities say. */
static double downstream_factor(const entry *prec, long n, long pnp,
                                long extra)
{ double rate = 1.0, cost = 0.0, rho, k;
  long j, end, end_ref;

  /* the next extra entries, and the next RATPOINTS_SP2_EXTRA, as far as the
   * pool goes; written so that a huge extra (-R) cannot overflow */
  end = (extra >= pnp - n - 1) ? pnp : n + 1 + extra;
  end_ref = (RATPOINTS_SP2_EXTRA >= pnp - n - 1) ? pnp
                                                 : n + 1 + RATPOINTS_SP2_EXTRA;
  for(j = n + 1; j < end; j++)
  { cost += RATPOINTS_COST_PHASE2*rate; rate *= prec[j].r; }
  cost += rate*RATPOINTS_COST_SURVIVOR;
  if(end_ref > n + 1)
  { double ref = 1.0;

    for(j = n + 1; j < end_ref; j++) { ref *= prec[j].r; }
    rho = pow(ref, 1.0/(double)(end_ref - n - 1));
  }
  else { rho = prec[n].r; }
  k = cost*(1.0 - rho)/RATPOINTS_COST_PHASE2;
  return((k < 1.0) ? 1.0 : k);
}

/* The rule that ends the first phase.  Modulus n of the pool is worth
 * adding while the expected survivors per 64-bit word it removes --
 * bits_per_word times the product rate of the densities taken so far,
 * times 1 - r_n --, weighted by what they cost downstream (above), exceed
 * the target times what the modulus costs per word.  (Until the tuning
 * session of 3.0.0 the rule compared the survivors the modulus meets, not
 * the ones it removes, and the target had the typical 1 - r of a half
 * folded in.)  The target, RATPOINTS_SURVIVORS_PER_WORD, is
 * fitted for a long run and a modulus that costs one AND per word; the
 * fixed costs of a modulus -- its tables, its sieve_spec and bp_list
 * entries per denominator, its row pointer per call and the fetch of its
 * row per denominator -- are spread over the words of the run in
 * entry.cost and raise the bar for it.  Over a long run that is by a few
 * per cent; at a height bound of a few hundred, where the run is a few
 * dozen words and a table has more rows than that, it is by a factor of a
 * hundred for the small primes and several hundred at the modulus where the
 * phase stops, and the phase stops after some seven primes where it used to
 * take a dozen whose tables were a fifth of the run (TODO item 29).  The
 * first modulus is always taken: the chunked sieve writes the 2-adic
 * pattern on its pass, so its ANDs cost nothing beyond that (its tables do,
 * but a sieve without a first phase is no sieve). */
static int phase_1_wants(const entry *prec, long n, long pnp, long taken,
                         double bits_per_word, double rate, double target,
                         long extra)
{ return(taken == 0
         || bits_per_word*rate*(1.0 - prec[n].r)
              *downstream_factor(prec, n, pnp, extra)
              > target*prec[n].cost); }

/* How many primes the first phase would need under that rule.  prec[] must
 * be sorted by increasing key.  If the target cannot be reached with the
 * primes available, all of them are used.  Since 3.0.0 the choice itself is
 * made by take_entries(), with the composite moduli in the pool; this
 * estimate, over the primes alone, serves the rule that decides whether to
 * look at more primes.  (It used to err on the high side, a composite
 * modulus being able only to lower the count; with the cost in the rule a
 * composite's table can also raise it, so it is an estimate and no bound.) */
static long primes_for_phase_1(entry *prec, long pnp,
                               double bits_per_word, double target, long extra)
{ double rate = 1.0;
  long n;

  for(n = 0; n < pnp; n++)
  { if(!phase_1_wants(prec, n, pnp, n, bits_per_word, rate, target, extra))
    { return(n); }
    rate *= prec[n].r;
  }
  return(pnp > 0 ? pnp : 1);
}

/* How many primes the second phase adds to the first.  A phase-2 prime is
 * paid for once -- its sieve table, and its entry in bp_list -- and then used
 * for the whole run, so how many are worth having depends on how long the
 * run is; see RATPOINTS_SP2_U0 in ratpoints.h .  With u0 = 0 this is a flat
 * offset, which is what every version before 3.0.0 used. */
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
#define RP_HORNER_STEPS (LONG_LENGTH/RATPOINTS_MAX_BITS_IN_PRIME - 1)

/* Look at one prime and record what it says about the curve.
 *
 * Fills in the table is_f_square[0..p], where entry a says whether f(a) is a
 * square modulo p and the last entry whether there are points at infinity,
 * and counts the residues that admit points.  When the prime carries any
 * information a sieve entry is built for it and *prec_entry is filled in
 * with that entry and the density r of the admissible residues.
 *
 * Returns 1 when some residue admits no point.  Returns 2 when every residue
 * does but denominators divisible by p occur: such a prime still says one
 * thing, namely that numerator and denominator are not both divisible by p
 * -- the row of the denominators divisible by p (sieves0) admits only the
 * numerators that are not -- which is a density of 1 - 1/p^2.  That is
 * little for a modulus of its own, but it comes free as a factor of a
 * composite modulus (33 in place of 11, 15 when both 3 and 5 are like
 * this), and curves with very many rational points have several such
 * primes; it is worth 1.5% on them at a height bound of 16383 (TODO item
 * 31).  The entry is marked nf = -1: it is a candidate for the first two
 * stages and a factor for add_moduli, but it does not count as one of the
 * primes the look-further rule wants, and the third stage, which runs after
 * the test for common factors, has no use for it.  Returns 0 when the prime
 * says nothing at all, and -1 when the curve has no points modulo p, so
 * that it has no rational points.
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

  if(np >= p && !is_f_square[p])
  { return(0); } /* the prime carries no information */

  { /* The mean density of admissible numerators over the classes of the
     * denominator mod p: np/p for the p-1 unit classes, and (p-1)/p for the
     * class divisible by p, whose row (sieves0) admits the numerators not
     * divisible by p; that class counts only when such denominators occur.
     * A prime power's density (examine_power) is computed the same way, and
     * the two compete in one ranking.  Until the tuning session of 3.0.0 the
     * last class counted as 1, which is 1/p^2 too much -- the numerators
     * that row removes are the ones the test for common factors removes
     * anyway, and the constants had been fitted to the mixture (item 21
     * measured the exact convention at 5% of testhighmany under those
     * constants); it is one of the estimate corrections refitted as a
     * group. */
    double r = is_f_square[p] ? ((double)((np + 1)*(p-1)))/((double)(p*p))
                              : (double)np/(double)p;

    prec_entry->r = r;
    prec_entry->p = p;
    prec_entry->mask = (pn < LONG_LENGTH) ? 1UL << pn : 0UL;
    prec_entry->nf = (np >= p) ? -1 : 0; /* -1: coprimality alone */
  }

  /* set up sieve_entry :
     typedef struct
       { ratpoints_init_fun init; long p; int *is_f_square;
         const long *inverses; unsigned long magic; double r;
         long bias; long dinv[RP_NUM_STRIDES];
         ratpoints_bit_array* sieve[RATPOINTS_MAX_PRIME]; }
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
    /* the reciprocal the sieve reduces with: the third stage (see stage3()
     * in sift.c), and the first and second phases through the magics array
     * of sieving_info.  One division per prime and curve, against one per
     * survivor and one per call saved. */
    se->magic = ULONG_MAX/(unsigned long)p + 1;
    /* the entry keeps the density too, so that the choice of primes can be
     * revisited during the run, when prec[] is long gone */
    se->r = prec_entry->r;
    /* the multiple of p in the row shifts, so that the word number plus the
     * shift is never negative; see RP_ROW_BIAS and class_offsets(), which
     * also fills dinv[] once the strides in use are known */
    se->bias = p*((RP_ROW_BIAS + p - 1)/p);
    se->rbainv = inverses[pn][RP_MULMOD(RBA_LENGTH, p, se->magic)];
    se->nf = 0; se->pw = NULL; /* a prime */
    /* sieves0 is 64-bit words, but is read as bit-arrays; it is given
     * the alignment of ratpoints_bit_array in gen_find_points_h.c . */
    se->sieve[0] = (ratpoints_bit_array *)&sieves0[pn][0];
    for(i = 1; i < p; i++) { se->sieve[i] = NULL; }

    prec_entry->ssp = se;
  }
  return((np >= p) ? 2 : 1);
}

/************************************************************************
 * Collect the sieving information                                      *
 ************************************************************************/

/* The number of numerators denominator b has to consider: the part of
 * b*domain that lies within the height bound.  Used both to size the run
 * (run_shape below) and to estimate how many survivors a denominator brings
 * to the third stage. */
static double numerators_for(const ratpoints_args *args, double b, double H,
                             double *n_inter, double *n_clip_lo,
                             double *n_clip_hi)
{ double sum = 0.0;
  long k, n = 0, clo = 0, chi = 0;

  for(k = 0; k < args->num_inter; k++)
  { double lo = b*args->domain[k].low, up = b*args->domain[k].up;
    int cl = 0, ch = 0;

    if(lo <= -H) { lo = -H; cl = 1; }
    if(up >= H) { up = H; ch = 1; }
    if(up > lo) { sum += up - lo; n++; clo += cl; chi += ch; }
  }
  /* the intervals that were not empty, and how many of their ends the
   * height bound cut (find_points_work clips them the same way) */
  *n_inter = (double)n; *n_clip_lo = (double)clo; *n_clip_hi = (double)chi;
  return(sum);
}

/* What the rounding of an interval to whole bit arrays adds at an end the
 * height bound cut, for the numerator class cl: the bits of the class run
 * from ceil((-H - a0)/2^k) to floor((H - a0)/2^k) (find_points_work), and
 * the padding is what floor and ceil to a multiple of RBA_LENGTH add below
 * the first and above the last.  Fixed for the class, since every such end
 * sits at the same bit; an end the bound did not cut sits anywhere in its
 * bit array and pads half of one on average. */
static void clipped_padding(const rp_num_class *cl, long H,
                            double *pad_lo, double *pad_hi)
{ long k = cl->k, a0 = cl->a0;
  long lo_bit = (-H - a0 + (1L << k) - 1) >> k;
  long hi_bit = ((H - a0) >> k) + 1;

  /* x & (RBA_LENGTH - 1) is x mod RBA_LENGTH, non-negative, for negative
   * x as well (two's complement) */
  *pad_lo += (double)(lo_bit & (RBA_LENGTH - 1));
  *pad_hi += (double)((-hi_bit) & (RBA_LENGTH - 1));
}

/* How big the run is: the number of denominators that will actually be
 * sifted, and the number of 64-bit words they sweep between them, padding
 * included.  Both are wanted by the rule that picks the sieving primes, because
 * two of the costs of a prime are paid once and then spread over the whole
 * run -- its sieve table, built for at most p denominator classes, and its
 * entry in bp_list, computed once per denominator.  Per word of numerators
 * those come to k*p*min(D,p)/U and l*D/U, and they are the reason the best
 * number of primes at a height bound of 200000 is not the best number at
 * 16383.
 *
 * Nothing here needs any sieving.  The denominators that get sifted are
 * those whose class mod 64 admits a numerator, that are not divisible by a
 * forbidden divisor and that pass the Jacobi symbol test where it applies;
 * the first two are periodic and are counted exactly, and the third lets
 * through 0.53 of what is left.  The numerators of one
 * denominator are piecewise linear in b with a handful of breakpoints, so a
 * midpoint sample over the range of b is accurate to a fraction of a per
 * cent.  What the sieve sweeps is more than the numerators: each interval
 * of each denominator is rounded outwards to whole bit arrays
 * (find_points_work), whatever the packing -- a third of the words at a
 * height bound of 1000, where a denominator has a few bit arrays, a tenth
 * at 4000, two per cent at 16383.  An end the height bound cut sits at the
 * same bit for every denominator of a class, so its padding is exact
 * (clipped_padding); an end inside the bound pads half a bit array on
 * average.  (One whole bit array per interval, the first version of this,
 * was a fifth too much at 1000: most ends are cut, and at a height bound of
 * 2^n - 1 the cut ends pad nothing at all.)
 * u_words counts the words swept, padding included, since that is what the
 * per-word costs are spread over; u_pad says how many of them are padding,
 * for the estimates that want the numerators themselves; n_calls is the
 * number of calls of the sieve, which a first-phase modulus pays a fixed
 * cost for (RATPOINTS_COST_CALL).
 *
 * The result is an estimate, and a biased one -- the Jacobi factor is an
 * average, and the valuation test of the use_squares1 path is not modelled
 * at all.  That is by design: it is used only to compare a fixed cost with a
 * per-word one, where being right to within a factor of about 1.5 moves the
 * chosen number of primes by less than one.
 */
#define RUN_SHAPE_SAMPLES 64

/* the share of a class's numerators its bit arrays hold: one in 2^k, k the
 * stride of the class (rp_num_class) */
static const double rp_inv_stride[RP_NUM_STRIDES]
  = {1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625};

/* The fraction of all integers b with v_p(b) in the set the mask describes
 * (bit m set <==> v_p(b) = m); the density of v_p(b) = m is (p-1)/p^(m+1). */
static double forbidden_fraction(long p, unsigned long mask)
{
  double f = 0.0;
  double q = 1.0/(double)p; /* 1/p^m */
  long m;

  for(m = 1; m < LONG_LENGTH && (mask >> m) != 0; m++)
  { if((mask >> m) & 1) { f += q*(1.0 - 1.0/(double)p); }
    q /= (double)p;
  }
  return(f);
}

static void run_shape(ratpoints_args *args, unsigned long den_bits,
                      const rp_num_class *cls,
                      long fba, long fdc,
                      double *n_denom, double *u_words, double *u_pad,
                      double *n_calls)
{ double H = (double)args->height;
  double keep = 1.0;    /* fraction of the candidates that reach sift() */
  double count = 0.0;   /* candidate denominators */
  double nums = 0.0;    /* numerators they sweep, before that fraction */
  double inters = 0.0;  /* the non-empty intervals they sweep them in */
  double clo = 0.0, chi = 0.0; /* how many of their ends the bound cut */
  double ni, nlo, nhi;
  long good = 0;        /* classes of b (mod 64, or of k for b = k^2) kept */
  double packed = 0.0;  /* the sum over them of 1/stride: the share of their
                         * numerators the bit arrays hold */
  double pad_lo = 0.0, pad_hi = 0.0; /* the sum over them of what a cut end
                                      * pads (clipped_padding) */
  long Hl = args->height;
  long i, j;

  if(args->flags & RATPOINTS_USE_SQUARES)
  { /* the denominators are the squares in [b_low, b_high] */
    double klo = ceil(sqrt((double)args->b_low));
    double khi = floor(sqrt((double)args->b_high));

    if(khi >= klo)
    { count = khi - klo + 1.0;
      for(i = 0; i < RUN_SHAPE_SAMPLES; i++)
      { double k = klo + (khi - klo)*((double)i + 0.5)/RUN_SHAPE_SAMPLES;
        nums += numerators_for(args, k*k, H, &ni, &nlo, &nhi);
        inters += ni; clo += nlo; chi += nhi;
      }
      nums *= count/RUN_SHAPE_SAMPLES; inters *= count/RUN_SHAPE_SAMPLES;
      clo *= count/RUN_SHAPE_SAMPLES; chi *= count/RUN_SHAPE_SAMPLES;
    }
    /* only the pattern for b mod 64 applies, and k^2 mod 64 has period 32
     * in k */
    for(j = 0; j < 32; j++)
    { const rp_num_class *cl = &cls[(j*j) & 0x3f];

      if(EXT0(cl->bits))
      { good++; packed += rp_inv_stride[cl->k];
        clipped_padding(cl, Hl, &pad_lo, &pad_hi);
      }
    }
    keep = (double)good/32.0;
  }
  else if(args->flags & RATPOINTS_USE_SQUARES1)
  { /* squares times the divisors of the leading coefficient */
    long *divisors = (long *)args->divisors;
    long n;
    long tried = 0;

    for(n = 0; divisors[n]; n++)
    { double d = (double)divisors[n];
      double klo = ceil(sqrt((double)args->b_low/d));
      double khi = floor(sqrt((double)args->b_high/d));

      if(klo < 1.0) { klo = 1.0; }
      if(khi >= klo)
      { double c = khi - klo + 1.0;

        for(i = 0; i < RUN_SHAPE_SAMPLES; i++)
        { double k = klo + (khi - klo)*((double)i + 0.5)/RUN_SHAPE_SAMPLES;
          nums += c*numerators_for(args, d*k*k, H, &ni, &nlo, &nhi)
                   /RUN_SHAPE_SAMPLES;
          inters += c*ni/RUN_SHAPE_SAMPLES;
          clo += c*nlo/RUN_SHAPE_SAMPLES; chi += c*nhi/RUN_SHAPE_SAMPLES;
        }
        count += c;
      }
      /* every divisor's 32 classes weigh the same here, whatever its share
       * of the denominators -- a bias keep has always had */
      for(j = 0; j < 32; j++, tried++)
      { const rp_num_class *cl
          = &cls[((unsigned long)divisors[n]*(unsigned long)(j*j)) & 0x3f];

        if(EXT0(cl->bits))
        { good++; packed += rp_inv_stride[cl->k];
          clipped_padding(cl, Hl, &pad_lo, &pad_hi);
        }
      }
    }
    if(tried) { keep = (double)good/(double)tried; }
  }
  else
  { /* every denominator in the range is a candidate */
    double blo = (double)args->b_low, bhi = (double)args->b_high;

    if(bhi >= blo)
    { count = bhi - blo + 1.0;
      for(i = 0; i < RUN_SHAPE_SAMPLES; i++)
      { double b = blo + (bhi - blo)*((double)i + 0.5)/RUN_SHAPE_SAMPLES;
        nums += numerators_for(args, b, H, &ni, &nlo, &nhi);
        inters += ni; clo += nlo; chi += nhi;
      }
      nums *= count/RUN_SHAPE_SAMPLES; inters *= count/RUN_SHAPE_SAMPLES;
      clo *= count/RUN_SHAPE_SAMPLES; chi *= count/RUN_SHAPE_SAMPLES;
    }
    /* bit j of den_bits is set exactly when the denominators congruent to j
     * modulo 64 have a numerator pattern (the word for the denominators
     * 64w..64w+63 has b at bit b mod 64) */
    { unsigned long w = den_bits;

      while(w)
      { j = RP_CTZL(w); w &= w - 1UL;
        good++; packed += rp_inv_stride[cls[j].k];
        clipped_padding(&cls[j], Hl, &pad_lo, &pad_hi);
      }
    }
    keep = (double)good/64.0;

    if(args->flags & RATPOINTS_CHECK_DENOM)
    { forbidden_entry *fb = (forbidden_entry *)args->forb_ba;
      forbidden_val *fd = (forbidden_val *)args->forbidden;

      for(i = 0; i < fba; i++) { keep *= 1.0 - 1.0/(double)fb[i].p; }
      for(i = 0; i < fdc; i++)
      { keep *= 1.0 - forbidden_fraction(fd[i].p, fd[i].mask); }
      /* the Jacobi symbol lets through half of the rest -- a little more
       * than half, since the denominators whose odd part divides into the
       * leading coefficient pass unconditionally: the count is 0.53 (item
       * 24's review; the runs of 2026-09-18 put Dact/Dpred at 1.07 on nine
       * tenths of the random curves at 16383 and 200000, which is 0.535) */
      if(args->flags & RATPOINTS_USE_JACOBI) { keep *= 0.53; }
    }
  }

  /* The bit arrays of a denominator hold one numerator in 2^k, k the stride
   * of its class (rp_num_class), so a class sweeps that share of its
   * numerators; the classes counted above are equally frequent among the
   * denominators, so the mean over the classes kept is the factor.  (Until
   * 3.0.0 only the packing by parity existed, and counting the even
   * denominators of a curve with numerators of both parities at full width
   * over-estimated U by up to a third.) */
  if(good > 0) { nums *= packed/(double)good; }

  /* The rounding to whole bit arrays: the bits of an interval run from
   * floor(low/RBA_LENGTH) to ceil(high/RBA_LENGTH) bit arrays.  An end the
   * height bound cut pads what clipped_padding says for its class, the
   * mean over the classes kept; an end inside the bound pads half a bit
   * array on average.  (Until 3.0.0 the ranges were padded to whole chunks of
   * RATPOINTS_CHUNK bit arrays on top of that; the tail legs of item 25
   * took that away, this is what is left.) */
  { double pad = (good > 0)
                   ? clo*pad_lo/(double)good + chi*pad_hi/(double)good
                       + (2.0*inters - clo - chi)*0.5*(double)RBA_LENGTH
                   : 0.0;
    /* and the calls of the sieve: one per interval and denominator while
     * an interval fits array_size bit arrays, which it does below a height
     * bound of some 30000; beyond that as many as it takes (the mean
     * interval stands for all of them) */
    double asz = (args->array_size > 0) ? (double)args->array_size
                                        : (double)RATPOINTS_ARRAY_SIZE;
    double per = (inters > 0.0) ? (nums + pad)/inters/(double)RBA_LENGTH
                                : 0.0;  /* bit arrays per interval */

    *n_denom = keep*count;
    *u_words = keep*(nums + pad)/(double)LONG_LENGTH;
    *u_pad = keep*pad/(double)LONG_LENGTH;
    *n_calls = keep*inters*ceil(per/asz);
  }
  /* the floors, for a run too small to estimate: one denominator, one word
   * of numerators (no padding), one call */
  if(*n_denom < 1.0) { *n_denom = 1.0; }
  if(*u_words < 1.0) { *u_words = 1.0; *u_pad = 0.0; }
  if(*n_calls < 1.0) { *n_calls = 1.0; }
}

/* The mean number of bits set per 64-bit word of a bit array on entry to
 * the sieve, over the words the run sweeps.  cls[b].bits holds the
 * admissible numerators of the denominators b mod 64 in their packing, as
 * one word repeated through the bit array, so the population count of that
 * word is the bits per word of the class; the classes are those the run
 * visits (run_shape: every class with a pattern on the plain path, k^2 mod
 * 64 for k = 0..31 with squares as denominators, d k^2 for each divisor d
 * of the leading coefficient with squares times divisors), each as often
 * as it comes up, and each weighted by the words it sweeps -- one in 2^k of
 * its numerators.  Until 3.0.0 the mean was over all 64 classes unweighted,
 * which on the square paths counted the 52 classes never visited; a monic
 * curve of odd degree sieves the twelve square classes only, where the
 * unweighted mean is higher (item 26's review; by 16 to 61% on the curves
 * probed in the tuning session).  The weighting by words takes more away
 * than that, on the square paths and on the plain one alike: the classes
 * with the larger strides pack their admissible numerators densely and
 * sweep few words, so per word swept the mean comes out below the old
 * one, by 3 to 11% on the square paths.  What the value is for is the
 * survivors per word swept, and that is what it now is.
 * Per word rather than per bit-array on purpose: measurements across
 * register widths show that the survivor rate at the best sp1 is constant
 * per word, not per bit-array (see RATPOINTS_SURVIVORS_PER_WORD in
 * ratpoints.h). */
static void bpw_add(const rp_num_class *cl, double *tot, double *wsum)
{ long c = __builtin_popcountl(EXT0(cl->bits));

  if(c) { *tot += rp_inv_stride[cl->k]*(double)c; *wsum += rp_inv_stride[cl->k]; }
}

static double mean_bits_per_word(const ratpoints_args *args,
                                 const rp_num_class *cls,
                                 unsigned long den_bits)
{ double tot = 0.0, wsum = 0.0;
  long j;

  if(args->flags & RATPOINTS_USE_SQUARES)
  { for(j = 0; j < 32; j++) { bpw_add(&cls[(j*j) & 0x3f], &tot, &wsum); } }
  else if(args->flags & RATPOINTS_USE_SQUARES1)
  { long *divisors = (long *)args->divisors;
    long n;

    for(n = 0; divisors[n]; n++)
    { for(j = 0; j < 32; j++)
      { bpw_add(&cls[((unsigned long)divisors[n]*(unsigned long)(j*j)) & 0x3f],
                &tot, &wsum);
      }
    }
  }
  else
  { unsigned long w = den_bits;

    while(w) { j = RP_CTZL(w); w &= w - 1UL; bpw_add(&cls[j], &tot, &wsum); }
  }
  return((wsum > 0.0) ? tot/wsum : 0.0);
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

  for(m = 1, pm = p; pm <= args->b_high && m < LONG_LENGTH - 1; m++)
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

/************************************************************************
 * Composite sieving moduli (TODO item 21)                              *
 ************************************************************************/

/* x^-1 mod m for x coprime to m, by the extended Euclidean algorithm */
static long modinv(long x, long m)
{ long r0 = m, r1 = x % m, s0 = 0, s1 = 1; /* r_i = s_i x mod m */

  if(r1 < 0) { r1 += m; }
  while(r1)
  { long q = r0/r1, t;

    t = r0 - q*r1; r0 = r1; r1 = t;
    t = s0 - q*s1; s0 = s1; s1 = t;
  }
  return((s0 < 0) ? s0 + m : s0);
}

/* Look at the prime power m = p^e and record what it says (rp_power in
 * rp-private.h): f mod m at every residue, frev mod m at every multiple of
 * p, the inverses of the units, and the density r -- the mean over the
 * classes of the denominator mod m of the share of numerators its row
 * admits, exactly: np/m for a unit b; for p | b the map x -> b x^-1 takes
 * the units x onto the residues of the same valuation as b, each equally
 * often, so the row admits the share of those residues with frev a square,
 * times the share of units among the numerators.  (A prime's density,
 * examine_prime, is the case e = 1 of this: (np+1)(p-1)/p^2 with the class
 * it divides.)  pinf says whether denominators divisible
 * by p occur at all (is_f_square[p] of the prime): when they do not, r is
 * the mean over the unit classes alone, as it is for a prime.  Returns 1
 * when the power says more than nothing. */
static int examine_power(ratpoints_args *args, rp_power *pw, long p, long e,
                         int use_c_long, long *c_long, int pinf)
{ mpz_t *c = args->cof;
  long degree = args->degree, D = degree + (degree & 1);
  unsigned long cm[D + 1], sq[RP_MODWORDS];
  long m = 1, k, x, np = 0;
  double r;

  for(k = 0; k < e; k++) { m *= p; }
  pw->p = p; pw->e = e; pw->m = m;
  for(k = 0; k <= degree; k++)
  { cm[k] = use_c_long ? (unsigned long)mod(c_long[k], m)
                       : mpz_fdiv_ui(c[k], (unsigned long)m); }
  if(degree & 1) { cm[D] = 0UL; }

  for(k = 0; k < RP_MODWORDS; k++) { sq[k] = 0UL; pw->fsq[k] = 0UL; pw->gsq[k] = 0UL; }
  for(x = 0; x < m; x++)
  { long t = (x*x) % m; sq[t >> LONG_SHIFT] |= 1UL << (t & LONG_MASK); }
#define RP_SQ_M(v) ((sq[(v) >> LONG_SHIFT] >> ((v) & LONG_MASK)) & 1UL)
  /* f(x) mod m by Horner, reduced every step (m < 2^10, so the accumulator
   * stays below 2^21) */
  for(x = 0; x < m; x++)
  { unsigned long v = cm[degree];

    for(k = degree - 1; k >= 0; k--) { v = (v*(unsigned long)x + cm[k]) % (unsigned long)m; }
    if(RP_SQ_M(v)) { pw->fsq[x >> LONG_SHIFT] |= 1UL << (x & LONG_MASK); np++; }
  }
  /* frev(t) = c_0 t^D + ... + c_D at the multiples of p */
  for(x = 0; x < m; x += p)
  { unsigned long v = cm[0];

    for(k = 1; k <= D; k++) { v = (v*(unsigned long)x + cm[k]) % (unsigned long)m; }
    if(RP_SQ_M(v)) { pw->gsq[x >> LONG_SHIFT] |= 1UL << (x & LONG_MASK); }
  }
#undef RP_SQ_M
  pw->inv[0] = 0;
  for(x = 1; x < m; x++) { pw->inv[x] = (x % p) ? (unsigned short)modinv(x, m) : 0; }
  pw->np = (int)np;

  { double units = (double)(m - m/p);

    r = units*(double)np/(double)m; /* the unit classes, density np/m each */
    if(pinf)
    { long v, pv = p;

      for(v = 1; v <= e; v++, pv *= p)
      { /* the classes with v_p(b) = v: phi(m/p^v) of them, one for v = e;
         * as many residues t with v_p(t) = v, of which good have frev(t)
         * square; a row's density is that share times the share of units,
         * and the count of classes cancels against the count of t */
        long good = 0, t;

        for(t = 0; t < m; t += pv)
        { if(v < e && ((t/pv) % p) == 0) { continue; }
          if((pw->gsq[t >> LONG_SHIFT] >> (t & LONG_MASK)) & 1UL) { good++; }
        }
        r += (double)good*units/(double)m;
      }
      r /= (double)m;
    }
    else { r /= units; }
    pw->r = r;
  }
  return((r < 1.0 - 1.0e-9) ? 1 : 0);
}

/* the largest table cost per numerator word, in units of a first-phase
 * AND, at which a composite modulus is offered at all: above it the modulus
 * cannot pay whatever it says, and at a small height bound this keeps the
 * set-up from looking at the prime powers */
#define RP_MODULUS_TABLE_MAX 2.0

/* Offer the composite moduli below RATPOINTS_MAX_PRIME_EVEN as candidates:
 * every odd m that is not a prime, factored into prime powers, with r the
 * product of its factors' densities and the key the ranking uses for a
 * prime of that size, provided every factor says something and every prime
 * involved is among those looked at.  A prime power is examined the first
 * time a modulus needs it; pinf[pn] says whether denominators divisible by
 * prime[pn] occur.  Candidates are appended to prec[] from *pnp_p on. */
static void add_moduli(ratpoints_args *args, entry *prec, long *pnp_p,
                       ratpoints_sieve_entry **prime_se, const int *pinf,
                       long pn_lim, long *npw_p, int use_c_long, long *c_long,
                       double cost_table, double call_cost)
{ rp_power *pws = (rp_power *)args->pw_buffer;
  long npw_max = num_powers();
  /* the power p^e of prime[pn] as an index into pws: -1 not examined, -2
   * examined and useless; e < 10 since 3^10 > 1024 */
  long pw_idx[RATPOINTS_NUM_PRIMES][10];
  double u = args->run_words, d = args->run_denoms;
  long pnp = *pnp_p, m, pn, e;

  for(pn = 0; pn < pn_lim; pn++)
  { for(e = 0; e < 10; e++) { pw_idx[pn][e] = -1; } }

  for(m = 9; m < RATPOINTS_MAX_PRIME_EVEN && m <= RATPOINTS_COMPOSITE_MAX; m += 2)
  { long rest = m, nf = 0, fac[RP_MAX_FACTORS];
    double r = 1.0, builds = (d < (double)m) ? d : (double)m;
    unsigned long mask = 0UL;
    int ok = 1;

    /* the table alone too dear for the run: not a candidate */
    if(cost_table*(double)m*builds/u > RP_MODULUS_TABLE_MAX) { continue; }
    /* factor m over the primes looked at */
    for(pn = 0; ok && rest > 1 && pn < pn_lim; pn++)
    { long p = prime[pn];

      if(p*p > rest)
      { /* what is left is a prime: find it, or give up */
        for( ; pn < pn_lim && prime[pn] < rest; pn++) {}
        if(pn >= pn_lim || prime[pn] != rest) { ok = 0; break; }
        p = rest;
      }
      if(rest % p) { continue; }
      for(e = 0; rest % p == 0; e++) { rest /= p; }
      if(nf >= RP_MAX_FACTORS || pn >= LONG_LENGTH) { ok = 0; break; }
      if(e == 1)
      { if(prime_se[pn] == NULL) { ok = 0; break; } /* says nothing */
        r *= prime_se[pn]->r; fac[nf] = pn;
      }
      else
      { long idx = pw_idx[pn][e];

        if(idx == -1)
        { if(*npw_p >= npw_max) { ok = 0; break; }
          idx = examine_power(args, &pws[*npw_p], p, e, use_c_long, c_long,
                              pinf[pn]) ? (*npw_p)++ : -2;
          pw_idx[pn][e] = idx;
        }
        if(idx < 0) { ok = 0; break; }
        r *= pws[idx].r; fac[nf] = RATPOINTS_NUM_PRIMES + idx;
      }
      mask |= 1UL << pn; nf++;
    }
    if(!ok || rest > 1) { continue; }
    if(nf == 1 && fac[0] < RATPOINTS_NUM_PRIMES) { continue; } /* m prime */
    if(r >= 1.0) { continue; }
    if(pnp >= RATPOINTS_NUM_PRIMES + RP_MAX_MODULI) { break; }
    prec[pnp].r = r; prec[pnp].p = m; prec[pnp].mask = mask;
    prec[pnp].ssp = NULL; prec[pnp].nf = (short)nf;
    for(e = 0; e < nf; e++) { prec[pnp].fac[e] = (short)fac[e]; }
    phase_1_key(&prec[pnp], cost_table, u, d, call_cost);
    pnp++;
  }
  *pnp_p = pnp;
}

/* The sieve entry of a prime power, made once, when a modulus needs it */
static ratpoints_sieve_entry *make_power(ratpoints_args *args, const rp_power *pw,
                                         double r)
{ ratpoints_sieve_entry *se = (ratpoints_sieve_entry *)args->se_next;
  long m = pw->m, i;

  args->se_next += sizeof(ratpoints_sieve_entry);
  se->init = _ratpoints_sieve_init_power;
  se->p = m; se->is_f_square = NULL; se->inverses = NULL;
  se->magic = ULONG_MAX/(unsigned long)m + 1;
  se->r = r;
  se->bias = m*((RP_ROW_BIAS + m - 1)/m);
  se->rbainv = modinv(RBA_LENGTH % m, m);
  se->nf = 1; se->pw = pw;
  for(i = 0; i < m; i++) { se->sieve[i] = NULL; }
  return(se);
}

/* The sieve entry of a composite modulus the ranking has taken: a prime
 * power's, or a product's with its factors' entries -- the primes' from
 * examine_prime, the powers' made here as needed (power_se, by index). */
static ratpoints_sieve_entry *make_modulus(ratpoints_args *args, entry *en,
                                           ratpoints_sieve_entry **prime_se,
                                           ratpoints_sieve_entry **power_se)
{ rp_power *pws = (rp_power *)args->pw_buffer;
  ratpoints_sieve_entry *se;
  long i;

  if(en->nf == 1)
  { long idx = en->fac[0] - RATPOINTS_NUM_PRIMES;

    if(power_se[idx] == NULL) { power_se[idx] = make_power(args, &pws[idx], en->r); }
    se = power_se[idx];
  }
  else
  { long m = en->p;

    se = (ratpoints_sieve_entry *)args->se_next;
    args->se_next += sizeof(ratpoints_sieve_entry);
    se->init = _ratpoints_sieve_init_product;
    se->p = m; se->is_f_square = NULL; se->inverses = NULL;
    se->magic = ULONG_MAX/(unsigned long)m + 1;
    se->r = en->r;
    se->bias = m*((RP_ROW_BIAS + m - 1)/m);
    se->rbainv = modinv(RBA_LENGTH % m, m);
    se->nf = en->nf; se->pw = NULL;
    for(i = 0; i < en->nf; i++)
    { long code = en->fac[i];

      if(code < RATPOINTS_NUM_PRIMES) { se->factor[i] = prime_se[code]; }
      else
      { long idx = code - RATPOINTS_NUM_PRIMES;

        if(power_se[idx] == NULL)
        { power_se[idx] = make_power(args, &pws[idx], pws[idx].r); }
        se->factor[i] = power_se[idx];
      }
    }
    for(i = 0; i < m; i++) { se->sieve[i] = NULL; }
  }
  en->ssp = se;
  return(se);
}

/* Take candidates from prec[from..*pnp_p) in the order they are in, for one
 * stage of the sieve: an entry sharing a prime with one taken before (used,
 * a mask of primes) is dropped from the pool, the ones taken stay where
 * they are, from prec[from] on.  With want >= 0 that many are taken; with
 * want < 0 the first-phase rule decides (phase_1_wants, with extra the
 * number of moduli the second phase will add): entries are taken while the
 * expected survivors per word swept that the entry removes -- bits_per_word
 * (per word swept, see bpw_swept in sieving_info) times the product rate of
 * the densities so far, times 1 - r --, weighted by what they cost
 * downstream, exceed target times the entry's cost per word.  When the pool
 * runs out first, all of
 * it is taken.  What is left of the pool has nothing in common with what
 * was taken.  Returns the number taken. */
static long take_entries(entry *prec, long from, long *pnp_p,
                         unsigned long *used, long want, double *rate,
                         double bits_per_word, double target, long extra)
{ long n = from, taken = 0, k;

  while(n < *pnp_p)
  { /* an entry that is no candidate leaves before the rule looks at it:
     * the rule reads the entry's own cost and density */
    if(prec[n].mask & *used)
    { for(k = n; k + 1 < *pnp_p; k++) { prec[k] = prec[k+1]; }
      (*pnp_p)--;
      continue;
    }
    if(want >= 0 ? taken >= want
                 : !phase_1_wants(prec, n, *pnp_p, taken, bits_per_word, *rate,
                                  target, extra))
    { break; }
    *used |= prec[n].mask; *rate *= prec[n].r; taken++; n++;
  }
  k = n;
  while(k < *pnp_p)
  { if(prec[k].mask & *used)
    { long j;

      for(j = k; j + 1 < *pnp_p; j++) { prec[j] = prec[j+1]; }
      (*pnp_p)--;
    }
    else { k++; }
  }
  return(taken);
}

static long sieving_info(ratpoints_args *args,
                         int use_c_long, long *c_long,
                         ratpoints_sieve_entry **sieve_list,
                         double bits_per_word, int may_extend,
                         unsigned long den_bits, const rp_num_class *cls)
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
  entry prec[RATPOINTS_NUM_PRIMES + RP_MAX_MODULI];
    /* This array is used for sorting in order to
       determine the `best' sieving moduli: the primes, then the composite
       candidates (TODO item 21). */
  /* the entries of the informative primes, by index; whether denominators
   * divisible by each prime occur; the prime powers' entries, by index */
  ratpoints_sieve_entry *prime_se[RATPOINTS_NUM_PRIMES];
  ratpoints_sieve_entry *power_se[RATPOINTS_NUM_PRIMES];
  int pinf[RATPOINTS_NUM_PRIMES];
  long npw = 0; /* prime powers examined so far */
  unsigned long used = 0UL; /* the primes the moduli taken involve */
  double u_pad = 0.0; /* of run_words, the padding to whole bit arrays */
  double n_calls = 0.0; /* calls of the sieve the run will make */
  double call_cost = 0.0; /* what they cost a first-phase modulus, per word */
  /* The survivors per word the two rules of the first phase and the key of
   * the second compare with the costs per word: bits_per_word is per word
   * of numerators, the costs are spread over every word swept, and the
   * padding to whole bit arrays (u_pad) is swept and masked, so it
   * carries none.  (The third stage's S counts per denominator and takes
   * the padding off itself.) */
  double bpw_swept = bits_per_word;

  forbidden_entry *forb_ba = (forbidden_entry *)args->forb_ba;
  forbidden_val *forbidden = (forbidden_val *)args->forbidden;

  /* How many primes to look at.  The loop below may raise this: see the
   * comment at its end. */
  long pn_lim = args->num_primes;
  long n_weak = 0; /* the primes among the candidates that say only that
                    * numerator and denominator are coprime (nf = -1) */
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
  run_shape(args, den_bits, cls, 0, 0,
            &args->run_denoms, &args->run_words, &u_pad, &n_calls);
  call_cost = RATPOINTS_COST_CALL*n_calls/args->run_words;
  bpw_swept = bits_per_word*(args->run_words - u_pad)/args->run_words;
  sp2_extra = phase_2_offset(sp2_extra, sp2_u0, args->run_words);

  for(pn = 0; pn < RATPOINTS_NUM_PRIMES; pn++)
  { prime_se[pn] = NULL; power_se[pn] = NULL; pinf[pn] = 0; }

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
    pinf[pn] = is_f_square[p];
    if(info > 0)
    { phase_1_key(&prec[pnp], cost_table, args->run_words, args->run_denoms,
                  call_cost);
      prime_se[pn] = prec[pnp].ssp;
      if(info == 2) { n_weak++; }
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
              pm <= args->b_high && m < LONG_LENGTH - 1; m++)
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
     * say next to nothing (they are not counted here, see n_weak) or nothing
     * at all (they are dropped above), and without this the second phase can
     * end up with nothing to sieve with at all.
     * Adding a prime can only raise pnp, and with the old stop rule could
     * only lower sp1 (the n smallest of a larger set have a smaller
     * product); with the cost in the rule (phase_1_wants) a new entry can
     * also move the stop by one the other way, so the shortfall need not
     * shrink at every step.  It shrinks on the whole, and the loop is bounded
     * by RATPOINTS_NUM_PRIMES in any case.
     */
    if(may_extend && pn + 1 == pn_lim && pn_lim < RATPOINTS_NUM_PRIMES)
    { long s1, want;

      qsort(prec, pnp, sizeof(entry), compare_entries);
      s1 = (args->sp1 >= 0) ? args->sp1
                            : primes_for_phase_1(prec, pnp, bpw_swept,
                                                 target, sp2_extra);
      want = (args->sp2 >= 0) ? args->sp2 : s1 + sp2_extra;
      if(pnp - n_weak < want) { pn_lim++; }
    }

  } /* end for pn */

  /* Terminate the array of forbidden divisors, having first looked for
   * more of them among the primes the loop above did not reach.  This is
   * done here, before the primes are chosen, because the choice needs to
   * know how many denominators will survive these tests: see run_shape.
   *
   * The search goes up to the square root of the height bound, or to the
   * end of prime[] if that comes first, and so beyond the compiled table of
   * sieving primes, whose word patterns in sieves0 the arrays used to be
   * limited to; for a prime beyond it the patterns are built here, in a
   * buffer that stays with args.  Why the square root: the Jacobi symbol
   * test lets a denominator through when it has an even number of bad
   * primes -- those with (lcf/p) = -1 -- and with every bad prime up to
   * sqrt(b_high) in the arrays, what gets through both tests is the product
   * of two bad primes beyond the table (b = q1*q2 or 2*q1*q2), which at a
   * height bound of 200000 with the table ending at 251 was 1.2 per cent of
   * the denominators sifted (the review's count, item P13). */
  if((args->flags & RATPOINTS_CHECK_DENOM)
       && !mpz_perfect_square_p(c[degree]))
       /* the test below asks for a non-square residue, which a square
        * leading coefficient never has; such a curve gets here since the
        * valuation test above applies to it */
  { long n, first = fba, words = 0;

    for(n = pn_lim; fba + fdc < args->max_forbidden && n < PRIMES1000; n++)
    { long p = prime[n];

      if(p*p > args->b_high) break;
      if(mpz_kronecker_si(c[degree], p) == -1)
      { forb_ba[fba].p = p;
        if(n < RATPOINTS_NUM_PRIMES)
        { forb_ba[fba].start = &sieves0[n][0];
          forb_ba[fba].end   = &sieves0[n][p];
        }
        else
        { forb_ba[fba].start = NULL; words += p; } /* built below */
        fba++;

#ifdef DEBUG
        printf("\nexcluding denominators divisible by %ld\n", p);
        fflush(NULL);
#endif

      }
    }

    /* the patterns for the primes beyond the table: p words for the prime
     * p, word r of them for the word numbers congruent to r mod p, with bit
     * j clear iff p divides 64*r + j -- what sieves0 holds for the compiled
     * primes (gen_find_points_h.c), and what the denominator loop expects */
    if(words > args->forb_words_len)
    { free(args->forb_words);
      args->forb_words = malloc(words*sizeof(unsigned long));
      args->forb_words_len = (args->forb_words == NULL) ? 0 : words;
    }
    if(words > args->forb_words_len)
    { /* no memory for the patterns: do without the primes beyond the table,
       * which are the last entries added */
      while(fba > first && forb_ba[fba-1].start == NULL) { fba--; }
    }
    { unsigned long *row = (unsigned long *)args->forb_words;

      for(n = first; n < fba; n++)
      { if(forb_ba[n].start == NULL)
        { long p = forb_ba[n].p, r, m;

          for(r = 0; r < p; r++) { row[r] = ~0UL; }
          for(m = 0; m < LONG_LENGTH*p; m += p)
          { row[m >> LONG_SHIFT] &= ~(1UL << (m & LONG_MASK)); }
          forb_ba[n].start = row;
          forb_ba[n].end   = row + p;
          row += p;
        }
        forb_ba[n].curr = forb_ba[n].start;
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
  ensure_ba_buffer(args, pn_lim, 0);

  /* The run shape again, now that the forbidden divisors are known and can
   * be taken off the denominator count; see run_shape.  The keys the primes
   * were given inside the loop used the first estimate, which does not know
   * about those divisors and so overstates the run, so they are computed
   * again here before anything is sorted for good. */
  { long e = (args->sp2_extra >= 0) ? args->sp2_extra : RATPOINTS_SP2_EXTRA;
    long n;

    run_shape(args, den_bits, cls, fba, fdc,
              &args->run_denoms, &args->run_words, &u_pad, &n_calls);
    call_cost = RATPOINTS_COST_CALL*n_calls/args->run_words;
    bpw_swept = bits_per_word*(args->run_words - u_pad)/args->run_words;
    sp2_extra = phase_2_offset(e, sp2_u0, args->run_words);
    for(n = 0; n < pnp; n++)
    { phase_1_key(&prec[n], cost_table, args->run_words, args->run_denoms,
                  call_cost);
    }
    /* and the composite moduli join the candidates, keyed the same way */
    add_moduli(args, prec, &pnp, prime_se, pinf, pn_lim, &npw,
               use_c_long, c_long, cost_table, call_cost);
  }

  /* sort the array to get at the best moduli */
  qsort(prec, pnp, sizeof(entry), compare_entries);

  /* Choose sp1 and sp2 unless they were given.
   * prec[] is now sorted by increasing key, what a modulus costs per unit
   * of what it says; r is the density of the numerators that are admissible
   * modulo the modulus, so the expected fraction of numerators surviving
   * the first n moduli is the product of their r.  Multiplied by the number
   * of bits actually set in a bit-array to begin with, that is the expected
   * number of survivors per bit-array; see the comment on
   * RATPOINTS_SURVIVORS_PER_WORD in ratpoints.h .  A modulus is taken while
   * that, weighted by what a survivor costs downstream, exceeds the target
   * times what the modulus costs per word, its fixed costs included
   * (phase_1_wants). */
  /* The candidates are taken in the order of their keys, but a modulus
   * sharing a prime with one already taken is passed over: the composite
   * moduli carry their primes' information, so the prime (or another
   * product involving it) would add nothing (take_entries). */
  { double rate = 1.0;

    args->sp1 = take_entries(prec, 0, &pnp, &used, args->sp1, &rate,
                             bpw_swept, target, sp2_extra);

    /* Rank what is left again, for the second phase.  There a modulus is
     * applied only to the bit arrays that survived the first phase, so its
     * per-word cost is smaller by the survival rate -- which makes the fixed
     * cost of its table weigh far more heavily, and the size of the modulus
     * matter far more than it does in the first phase. */
    if(args->sp1 < pnp)
    { long n;

      rate *= bpw_swept;
      for(n = args->sp1; n < pnp; n++)
      { prec[n].key = prime_key(prec[n].r, prec[n].p,
                                RATPOINTS_COST_PHASE2*rate, 1, cost_table,
                                args->run_words, args->run_denoms);
      }
      qsort(&prec[args->sp1], pnp - args->sp1, sizeof(entry), compare_entries);
    }
    { long want = (args->sp2 >= 0) ? args->sp2 - args->sp1 : sp2_extra;

      if(want < 0) { want = 0; }
      args->sp2 = args->sp1
                  + take_entries(prec, args->sp1, &pnp, &used, want, &rate,
                                 0.0, 0.0, 0);
    }
  }

  /* The third stage and the correction during the run see primes only, and
   * only primes that exclude a residue: the composite moduli not taken leave
   * the pool (their primes stay, unless a modulus taken involves them, in
   * which case they left already), and so do the primes that say no more
   * than that numerator and denominator are coprime (nf = -1), which the
   * test for common factors has seen to by then. */
  { long n = args->sp2;

    while(n < pnp)
    { if(prec[n].nf != 0)
      { long k;

        for(k = n; k + 1 < pnp; k++) { prec[k] = prec[k+1]; }
        pnp--;
      }
      else { n++; }
    }
  }

  /* put the moduli taken into sieve_list, making the entries of the
   * composite ones */
  { long n;

    for(n = 0; n < args->sp2; n++)
    { if(prec[n].ssp == NULL) { make_modulus(args, &prec[n], prime_se, power_se); }
      sieve_list[n] = prec[n].ssp;
    }
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
   * thinned by the 2-adic pre-sieve and by the primes of the first two
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
     * denominator, less the padding to whole bit arrays, which holds no
     * numerators (the boundary masks clear it); run_shape counts only the
     * numerators that are actually looked at, so no further halving is
     * wanted here. */
    S = (args->run_words - u_pad)*(double)LONG_LENGTH/args->run_denoms;
    { long n;

      S *= bits_per_word/(double)LONG_LENGTH;
      for(n = 0; n < args->sp2; n++) { S *= prec[n].r; }
      S *= RATPOINTS_SP3_COPRIME;
    }

    while(sp3 < sp3_want)
    { long best = -1;

      /* the best of the primes not yet spoken for; only the ones this stage
       * takes need to be in order, so each round picks one, and the scan is
       * repeated when a prime has been added to the pool */
      if(sp3 < pnp)
      { long m;

        best = sp3;
        for(m = sp3 + 1; m < pnp; m++)
        { if(prec[m].r < prec[best].r) { best = m; } }
      }

      /* When the pool holds no prime that pays for itself -- none at all, or
       * none good enough -- look at a further prime, as long as the stage
       * still has appetite for a prime of the quality this curve has been
       * offering: r_typ, the mean density of the informative primes seen so
       * far, stands for the one about to be looked at.  On a curve with
       * very many rational points the small primes are useless and the
       * informative ones come late, and the stage used to stop at the first
       * poor prime in hand although the next ones would have paid (item 27).
       * When the caller fixed sp3 the stage takes what it is told to, and
       * looks further only when the pool is empty. */
      if(best < 0
         || (args->sp3_extra < 0
             && S*(1.0 - per_surv - prec[best].r) <= per_denom))
      { long coeffs_mod_p[degree+1];
        int *is_f_square;
        int info;

        if(!may_extend || pn_lim >= RATPOINTS_NUM_PRIMES) { break; }
        if(args->sp3_extra < 0)
        { double r_typ = 1.0; /* no informative prime seen: assume none */
          long n, cnt = 0;

          if(pnp > 0)
          { r_typ = 0.0;
            /* a modulus that says only "coprime" is not what is looked for */
            for(n = 0; n < pnp; n++)
            { if(prec[n].nf >= 0) { r_typ += prec[n].r; cnt++; } }
            r_typ = (cnt > 0) ? r_typ/(double)cnt : 1.0;
          }
          if(S*(1.0 - per_surv - r_typ) <= per_denom) { break; }
        }
        info = examine_prime(args, pn_lim, use_c_long, c_long,
                             &coeffs_mod_p[0], &is_f_square, &prec[pnp]);
        pn_lim++;
        if(info < 0)
        { return(prime[pn_lim-1]); /* no points mod p */ }
        if(info == 0 || info == 2)
        { continue; } /* it says nothing, or nothing this stage can use,
                       * which runs after the test for common factors; try
                       * the next one */
        /* the third stage builds no table, so its primes are ranked by what
         * they say alone, which is what the selection above does */
        prec[pnp].key = prec[pnp].r;
        prec[pnp].cost = 0.0; /* never a first-phase candidate */
        pnp++;
        continue; /* choose again with the new prime in the pool */
      }

      if(best != sp3)
      { entry t = prec[sp3]; prec[sp3] = prec[best]; prec[best] = t; }
      S *= prec[sp3].r;
      sieve_list[sp3] = prec[sp3].ssp;
      sp3++;
    }
    args->sp3 = sp3;

    /* Put the rest of the primes in sieve_list too, in the order the third
     * stage would take them.  They cost nothing to keep -- no table is built
     * and nothing computed for it until a prime is actually used -- and
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
  { long extra = 0, n;

    /* the tables of the composite moduli taken, and of the prime powers
     * that are factors of them (a prime factor's table is in the primes'
     * share) */
    for(n = 0; n < args->sp2; n++)
    { ratpoints_sieve_entry *se = sieve_list[n];

      if(se->nf > 0)
      { long i;

        extra += se->p*(se->p + RATPOINTS_CHUNK-1);
        for(i = 0; i < se->nf && se->pw == NULL; i++)
        { ratpoints_sieve_entry *fe = se->factor[i];

          if(fe->pw) { extra += fe->p*(fe->p + RATPOINTS_CHUNK-1); }
        }
      }
    }
    ensure_ba_buffer(args, pn_lim, extra);
  }

  /* the reciprocals the first two phases reduce word numbers with, in the order
   * the moduli are used; see the note in find_points_init */
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
            " bpw=%.2f bpws=%.2f U=%.6g D=%.6g pad=%.6g calls=%.6g", pn_lim,
            pnp, args->sp1, args->sp2, args->sp3, bits_per_word, bpw_swept,
            args->run_words, args->run_denoms, u_pad, n_calls);
    for(n = 0; n < pnp; n++)
    { fprintf(stderr, " %ld:%.4f", prec[n].p, prec[n].r); }
    fprintf(stderr, "\n");
  }
#endif

  if(args->flags & RATPOINTS_VERBOSE)
  { printf("  %.1f bits set per word, %ld primes looked at"
           " ==> use %ld moduli in the first phase, %ld altogether,\n"
           "  and %ld primes more in the third stage\n",
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
          const rp_num_class *cls,
          ratpoints_sieve_entry **sieve_list, long *bp_list, int *quit,
          int process(long, long, const mpz_t, void*, int*), void *info)
{
  long total = 0;
  /* typedef struct { long p; long offset; ratpoints_bit_array *ptr;
                     ratpoints_bit_array *start; ratpoints_bit_array *end; }
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

  /* Note that b is new */
  args->flags |= RATPOINTS_COMPUTE_BC;
  args->stage3_filled = 0; /* see fill_checks() in sift.c */

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
          long bp = bp_list[n]; /* b 2^-k mod p, see fill_bp_list */
          ratpoints_bit_array *sptr = se->sieve[bp];

          ssp[n].p = p;
          /* the shift of the row for the packing of the class, with the
           * multiple of p that keeps the row index non-negative built in
           * (see rp_num_class and sieve_spec in rp-private.h) */
          ssp[n].offset = cls->offset[n];

#ifdef DEBUG
          printf("\np = %ld, bp = %ld, offset = %ld (+ bias %ld)\n",
                 p, bp, ssp[n].offset - se->bias, se->bias);
          fflush(NULL);
#endif
          /* copy if already initialized, else initialize */
          if(sptr) { ssp[n].ptr = sptr; }
          else
          { RP_INIT_TIC(t_init);
            ssp[n].ptr = (*(se->init))(se, bp, args);
            RP_INIT_TOC(t_init, p);
          }
          /* the end of the table, which the first phase's wrap-around
           * compares against; the start field is set by sift0 at the head
           * of every call, for the first-phase moduli, and nothing reads
           * it before that */
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

        /* the primes of the third stage need no table, only the inverse of
         * b modulo each of them, and fill_checks() in sift.c looks that up
         * on the first numerator that reaches the stage: most denominators
         * bring none that far */
        RP_SETUP_TOC(t_setup);
      }

      /* From numerators to bits: bit t stands for a0 + 2^k t, so the bits
       * are ceil((low - a0)/2^k) .. floor((high - a0)/2^k); the shifts of
       * signed values round down, as they do throughout the program.  (Only
       * the two-fold packings existed until 3.0.0: low >>= 1 and so on.) */
      { long k = cls->k, a0 = cls->a0;

        low = (low - a0 + (1L << k) - 1) >> k;
        high = (high - a0) >> k;
      }

      /* The bit interval is [low, high], both ends included.  (It used to
       * be made half-open by high++, and the bit arrays counted as
       * CEIL(high, RBA_LENGTH) = (high + RBA_LENGTH - 1) >> RBA_SHIFT; both
       * overflow a long when the last bit is within RBA_LENGTH of LONG_MAX,
       * which a height bound that close to it reaches, and the interval
       * was then dropped without a word.) */
      if(low <= high)
      { long w_low, w_high;
        long w_low0, w_high0;
        long range = args->array_size;

        /* Now the range of longwords (= bit_arrays): the one holding the
         * first bit to the one past the last; the shifts round down for
         * negative values too */
        w_low = low >> RBA_SHIFT; /* FLOOR(low, RBA_LENGTH); */
        w_high = (high >> RBA_SHIFT) + 1;
        w_low0 = w_low;
        w_high0 = w_low0 + range;
        for( ; w_low0 < w_high; w_low0 = w_high0, w_high0 += range)
        { if(w_high0 > w_high)
          { w_high0 = w_high; range = w_high0 - w_low0; }
          /* The bit arrays are not written here.  The first phase's
           * first prime ANDs the 2-adic pattern in as it sieves, which
           * saves a store and a load on every one of them; what is left
           * for sift0 to do afterwards is the two boundary words, and it
           * is told about them like this.  The range is not padded to a
           * multiple of RATPOINTS_CHUNK either, since 3.0.0: sift0 sieves
           * the bit arrays past the last whole chunk in narrower legs.
           * (Until then the padding was sieved with every prime of the
           * first phase, walked by the scan of the second and zeroed --
           * a seventh of all bit arrays swept at height 16383, and four
           * fifths of them at height 1000.) */
          { long mask_low = 0, mask_high = 0;

            if(w_low0 == w_low)
            /* lower bits of the first bit array are to be set to zero */
            { mask_low = low - RBA_LENGTH * w_low; }
            if(w_high0 == w_high)
            /* upper bits of the last bit array are to be set to zero:
             * those above the last bit, high mod RBA_LENGTH */
            { mask_high = RBA_LENGTH - 1 - (high & (RBA_LENGTH - 1)); }

            total += _ratpoints_sift0(b, w_low0, w_high0, args, cls,
                                      survivors, mask_low, mask_high,
                                      &ssp[0], &csp[0], quit, process, info);
            if(*quit) { RP_SIFT_TOC(t_sift); return(total); }
      } } }
  } }

  RP_SIFT_TOC(t_sift);
  return(total);
}

/**************************************************************************
 * Find points by looping over the denominators and sieving numerators    *
 **************************************************************************/

/* The denominator, divided by 2^k modulo each prime of the first two phases
 * -- k the stride its numerators are packed with, so that this is the
 * residue whose sieve table the denominator reads (see rp_num_class) -- in
 * bp_list, computed afresh for every denominator by the multiply-high
 * reduction (RP_MULMOD; two of them for a denominator beyond 2^32 divided by
 * the largest prime that can be compiled in, a division beyond 2^32).  It
 * used to be
 * stepped from the previous denominator, bp += d followed by while(bp >= p)
 * bp -= p: one to three data-dependent branches per prime and denominator,
 * mispredicted a quarter of the time, an eighth of all the branch misses of
 * make test1.  Computing it afresh also does away with the bookkeeping of
 * which entries were up to date when adapt_primes had just brought another
 * prime into play, which is why that correction is made here first: it is
 * due once the sieve has swept as many words as adapt_at says.  The primes
 * of the third stage are not in the list any more: what that stage needs is
 * looked up when a numerator reaches it, see fill_checks() in sift.c. */
static inline void fill_bp_list(long b, long k, long *bp_list,
                                ratpoints_args *args,
                                ratpoints_sieve_entry **sieve_list)
{ long n, sp2;
  const unsigned long *magics = (const unsigned long *)args->magics;

  if(args->n_words >= args->adapt_at) { adapt_primes(args); }
  sp2 = args->sp2;
  RP_BP_TIC(t_bp);
  /* The reduction is exact below 2^32.  b times 2^-k mod p stays below that
   * for a denominator below 2^32 divided by the largest prime that can be
   * compiled in (some ten million at the default prime size), and one
   * reduction does; up to 2^32 itself b is reduced first and the product
   * then, still without a division; beyond, mod() divides, as it always
   * did. */
  if(b <= RP_MULMOD_LIMIT/RATPOINTS_MAX_PRIME_EVEN)
  { for(n = 0; n < sp2; n++)
    { ratpoints_sieve_entry *se = sieve_list[n];

      bp_list[n] = RP_MULMOD(b*se->dinv[k], se->p, magics[n]);
    }
  }
  else if(b <= RP_MULMOD_LIMIT)
  { for(n = 0; n < sp2; n++)
    { ratpoints_sieve_entry *se = sieve_list[n];
      long p = se->p;

      bp_list[n] = RP_MULMOD(RP_MULMOD(b, p, magics[n])*se->dinv[k], p,
                             magics[n]);
    }
  }
  else
  { for(n = 0; n < sp2; n++)
    { ratpoints_sieve_entry *se = sieve_list[n];
      long p = se->p;

      bp_list[n] = RP_MULMOD(mod(b, p)*se->dinv[k], p, magics[n]);
    }
  }
  RP_BP_TOC(t_bp, sp2);
}

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

  /* the degree sizes an array of find_points_work_1 before that function
   * can check anything, so a negative one is refused here */
  if(args->degree < 0) { return(RATPOINTS_BAD_ARGS); }
  /* sp3 is a working field the caller never sets, and the three _used
   * fields are outputs: they stay 0 when the search ends before it has
   * chosen its primes (no real points, nothing admissible mod 64, ...). */
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
  int lcfsq;             /* whether the leading coefficient is a square;
                           set once the degree is known, below */

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
  unsigned long den_bits;
  rp_num_class cls[64]; /* the numerator classes, by b mod 64 */

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
  lcfsq = mpz_perfect_square_p(c[degree]);

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
  ensure_ba_buffer(args, args->num_primes, 0);
  if(args->sp2 > args->num_primes) { args->sp2 = args->num_primes; }
  if(args->sp2 >= 0 && args->sp1 > args->sp2) { args->sp1 = args->sp2; }

  if(height < 1) { return(RATPOINTS_BAD_ARGS); }
  if(args->b_low < 1) { args->b_low = 1; }
  if(args->b_high < 1) { args->b_high = height; }
  if(args->b_high > height) { args->b_high = height; }
  if(args->max_forbidden < 0)
  { args->max_forbidden = RATPOINTS_DEFAULT_MAX_FORBIDDEN; }
  if(args->max_forbidden > PRIMES1000)
  { args->max_forbidden = PRIMES1000; }
  if(args->array_size <= 0) { args->array_size = RATPOINTS_ARRAY_SIZE; }
  { long s = 2*CEIL(height, LONG_LENGTH);
    if(args->array_size > s) { args->array_size = s; }
  }
  /* make sure that array size is a multiple of RATPOINTS_CHUNK, so that of
   * the blocks a numerator interval is cut into only the last one can end
   * in a partial chunk (see _ratpoints_sift0) */
  args->array_size = CEIL(args->array_size, RATPOINTS_CHUNK)*RATPOINTS_CHUNK;
  if(args->sturm > LONG_LENGTH - 2) { args->sturm = LONG_LENGTH - 2; }

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
    { printf("  number of moduli for sieving:     (chosen from the curve)\n"); }
    else { printf("  number of moduli for sieving:     %3ld\n", args->sp2); }
    if(args->sp1 < 0)
    { printf("  number of moduli for first stage: (chosen from the curve)\n"); }
    else { printf("  number of moduli for first stage: %3ld\n", args->sp1); }
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
  { printf("Obtain information from the polynomial mod 64:\n"); }
  get_2adic_info(args, &den_bits, &cls[0]);
  /* Bit k in den_bits is 0 if b congruent to k mod 64 need not be
     considered as a denominator; cls[b] says how the numerators of the
     denominators b mod 64 are packed into the bit arrays, and which of
     them are admissible (rp_num_class in rp-private.h). */

  if(den_bits == 0 && !point_at_infty)
  { /* No residue class of the denominator mod 64 admits any numerator, so
     * there is no affine point.  There is none at infinity either: an odd
     * degree always leaves the class b = 0 mod 64 admissible, so the degree
     * is even, and the leading coefficient then is not a square mod 64,
     * hence not a square (the test on point_at_infty only says so). */
    if(args->flags & RATPOINTS_VERBOSE)
    { printf("  no denominator admits a numerator mod 64 ==> no points\n\n"); }
    return(total);
  }

#ifdef DEBUG
  printf("\nden_bits: %*.*lx\n\n", WIDTH, WIDTH, den_bits);
  fflush(NULL);
#else
  if(args->flags & RATPOINTS_VERBOSE)
  { /* the stride is the same for every class of one 2-adic valuation of the
     * denominator (see get_2adic_info), so one class of each kind says it */
    /* one class of each kind of denominator: b odd, 2 mod 4, 4 mod 8, ...,
     * 32 mod 64, 0 mod 64 */
    static const long rep[] = {1, 2, 4, 8, 16, 32, 0};
    long v;

    printf("  the bit arrays hold one numerator in (for b odd, 2 mod 4,"
           " 4 mod 8, ..., 32 mod 64, 0 mod 64; - = none admissible):");
    for(v = 0; v < (long)(sizeof(rep)/sizeof(rep[0])); v++)
    { if(EXT0(cls[rep[v]].bits)) { printf(" %ld", 1L << cls[rep[v]].k); }
      else { printf(" -"); }
    }
    printf("\n\n");
  }
#endif

  /* set up the sieve data structure */
  if(args->flags & RATPOINTS_VERBOSE)
  { printf("Find the points mod p for the first %ld odd primes p:\n",
           args->num_primes);
  }
  { /* The mean number of bits set in one word of a bit-array on entry to
     * the sieve, over the classes of denominators the run visits and the
     * words they sweep; see mean_bits_per_word. */
    double bits_per_word = mean_bits_per_word(args, cls, den_bits);
    { long ret = sieving_info(args, use_c_long, &c_long[0], sieve_list,
                              bits_per_word, np_is_default,
                              den_bits, &cls[0]);

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

    printf("  use %ld moduli for first stage:\n   ", args->sp1);
    for(n = 0; n < args->sp1; n++)
    { printf(" %ld", sieve_list[n]->p); }
    printf("\n  use %ld moduli for second stage:\n   ", args->sp2 - args->sp1);
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
    /* the row shifts of the numerator classes: one row per distinct
     * packing (a dozen or so; at most 64), over every prime that may come to
     * be sieved with */
    long offsets[num_packings(&cls[0])*(args->sp3_max > 0 ? args->sp3_max : 1)];

#ifdef DEBUG
    printf("\nfind_points_work: allocating space for survivors...");
    fflush(NULL);
#endif

    /* allocate space for survivors array; make sure of correct alignment.
     * One spare bit array pays for the alignment, and one more for the
     * sentinel that the scan in _ratpoints_sift0 runs into: it sits just
     * past the range, and the range can be all of array_size. */
    survivors_na = malloc((args->array_size+2)*sizeof(ratpoints_bit_array));
    survivors = (ratpoints_bit_array *)
                pointer_align(survivors_na, sizeof(ratpoints_bit_array));
    /* the row shifts of the numerator classes, now that the primes are
     * known */
    class_offsets(&cls[0], sieve_list, args->sp3_max, &offsets[0]);
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

#ifdef DEBUG
        printf("\n  using squares\n");
        fflush(NULL);
#endif

        /* from the first square in the range; b*b <= b_high, written so
         * that the square cannot overflow */
        for(b = ceil_sqrt(args->b_low); b <= args->b_high/b; b++)
        { const rp_num_class *cl;

          bb = b*b;
          cl = &cls[bb & 0x3f];
          if(EXT0(cl->bits))
          { fill_bp_list(bb, cl->k, bp_list, args, sieve_list);
            total += sift(bb, survivors, args, cl,
                          sieve_list, &bp_list[0],
                          &quit, process, info);
            if(quit) { break; }
          }

#ifdef DEBUG
          else
          { printf("\nb = %ld: excluded mod 64\n", bb);
            fflush(NULL);
          }
#endif
        }
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
        {
#ifdef DEBUG
          printf("\n  divisor = %ld\n", *div);
          fflush(NULL);
#endif

          /* from the first multiple of the divisor by a square in the
           * range; d*b*b <= b_high, written so that the product cannot
           * overflow (the divisors are at most b_high, see setup_us1) */
          for(b = ceil_sqrt((args->b_low - 1)/(*div) + 1);
              b <= (args->b_high/(*div))/b; b++)
          { int flag = 1;
            const rp_num_class *cl;

            bb = (*div)*b*b;
            cl = &cls[bb & 0x3f];
            if(EXT0(cl->bits))
            { long i;

              for(i = 0; den_info[i].p; i++)
              { int v = valuation1(bb, den_info[i].p);
                if((v >= den_info[i].slope)
                     && ((v + (den_info[i].val)) & 1))
                { flag = 0; break; }
              }
              if(flag)
              { fill_bp_list(bb, cl->k, bp_list, args, sieve_list);
                total += sift(bb, survivors, args, cl,
                              sieve_list, &bp_list[0],
                              &quit, process, info);
                if(quit) { break; }
              }
            }

#ifdef DEBUG
            else
            { printf("\nb = %ld: excluded mod 64\n", bb);
              fflush(NULL);
            }
#endif
          }
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
        long w, w_low = args->b_low >> LONG_SHIFT;
        long w_high = args->b_high >> LONG_SHIFT;
        /* the Jacobi symbol test as a product of Legendre symbols, when the
         * leading coefficient allows it; see jacobi_setup */
        jacobi_info ji;
        unsigned char jtab[RP_JACOBI_TABLE];
        int fast_jacobi = (args->flags & RATPOINTS_USE_JACOBI)
                            && jacobi_setup(&ji, jtab, RP_JACOBI_TABLE,
                                            c[degree], work[0], args->b_high);

#ifdef DEBUG
        printf("\n  taking account of forbidden divisors of the denominator\n");
        if(args->flags & RATPOINTS_USE_JACOBI)
        { printf("  Jacobi symbol test %s\n",
                 fast_jacobi ? "by Legendre symbols" : "by jacobi1/jacobi");
        }
        fflush(NULL);
#endif

        /* The 2-adic test on a denominator depends on b mod 64 alone -- bit
         * b mod 64 of den_bits says whether its class has a numerator
         * pattern -- and the forbidden-divisor arrays are words indexed by b
         * mod 64 as well.  So the denominators are taken a word of 64 at a
         * time: the word of those that pass both tests is one AND per array,
         * and the loop below visits only the bits that are set, which on a
         * random curve are a third of the denominators.  Bit j of the word
         * for w stands for b = 64*w + j. */
        { forbidden_entry *fba = &forb_ba[0];

          while(fba->p)
          { fba->curr = fba->start + mod(w_low, fba->p);
            fba++;
          }
        }

#ifdef DEBUG
        printf("\n  den_bits = %*.*lx\n", WIDTH, WIDTH, den_bits);
        fflush(NULL);
#endif

        for(w = w_low; w <= w_high; w++)
        { unsigned long b_bits = den_bits;
          long base = w << LONG_SHIFT;

          { forbidden_entry *fba = &forb_ba[0];

            while(fba->p)
            { b_bits &= *(fba->curr);
              fba->curr++;
              if(fba->curr == fba->end) { fba->curr = fba->start; }
              fba++;
            }
          }
          /* the first and the last word may be entered part way */
          if(w == w_low) { b_bits &= ~0UL << (args->b_low & LONG_MASK); }
          if(w == w_high)
          { b_bits &= ~0UL >> (LONG_MASK - (args->b_high & LONG_MASK)); }

#ifdef DEBUG
          printf("\n  w = %ld: b_bits = %*.*lx\n", w, WIDTH, WIDTH, b_bits);
          fflush(NULL);
#endif

          while(b_bits)
          { const rp_num_class *cl;

            b = base + RP_CTZL(b_bits);
            b_bits &= b_bits - 1UL;
            cl = &cls[b & 0x3f];

            /* the Jacobi symbol test comes first when it is the cheap one:
             * a few multiplications against the divisions of the valuation
             * test, and it rejects half of what gets here */
            if(fast_jacobi && !jacobi_test(b, &ji))
            {
#ifdef DEBUG
              printf("\nb = %ld: excluded by Jacobi symbol\n", b);
              fflush(NULL);
#endif
              continue;
            }

            /* check if denominator is excluded: is v_p(b) one of the
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
                && (fast_jacobi || !(args->flags & RATPOINTS_USE_JACOBI)
                      || (use_c_long
                           ? jacobi1(b, c_long[degree])
                           : jacobi(b, work[0], c[degree])) == 1))
            { fill_bp_list(b, cl->k, bp_list, args, sieve_list);
              total += sift(b, survivors, args, cl,
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
          if(quit) { break; }
        }
      } /* if(args->flags & RATPOINTS_CHECK_DENOM) */
      else
      { long b;
        long bp_list[args->sp3_max > 0 ? args->sp3_max : 1];
          /* sp3_max, not sp3: adapt_primes may reach for a
           * further prime as the run goes on */

        for(b = args->b_low; b <= args->b_high; b++)
        { const rp_num_class *cl = &cls[b & 0x3f];

          if(EXT0(cl->bits))
          { fill_bp_list(b, cl->k, bp_list, args, sieve_list);
            total += sift(b, survivors, args, cl,
                          sieve_list, &bp_list[0],
                          &quit, process, info);
            if(quit) { break; }
          }

#ifdef DEBUG
          else
          { printf("\nb = %ld: excluded mod 64\n", b);
            fflush(NULL);
          }
#endif

          if(b == LONG_MAX) { break; } /* b++ would overflow */
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
    static unsigned long long last_cyc1 = 0, last_cyc2 = 0, last_cyc3 = 0;
    static unsigned long long last_and2 = 0, last_rows = 0, last_calls = 0;

    fprintf(stderr, "[runshape] Upred=%.6g Uact=%.6g Dpred=%.6g Dact=%.6g"
            " words=%lu arrays=%lu bits=%lu coprime=%lu checks=%lu kodd=%ld"
            /* and this curve's share of the phase counters (TODO item 30:
             * the cost of an AND against the footprint of the rows; rdtsc,
             * so scale by a pinned cycle count over the process before
             * comparing runs); and2 is 0 without RP_PHASE_COUNTS */
            " sp1=%ld sp2=%ld calls=%llu cyc1=%llu cyc2=%llu cyc3=%llu"
            " and2=%llu rows=%llu\n",
            args->run_words,
            (double)(_rp_arrays_swept - last_arrays)*(double)RBA_PACK,
            args->run_denoms, (double)(_rp_bp_dens - last_dens),
            args->n_words, args->n_arrays, args->n_bits,
            args->n_coprime, args->n_checks,
            EXT0(cls[1].bits) ? cls[1].k : -1L,
            args->sp1, args->sp2, _rp_sift0_calls - last_calls,
            _rp_phase1_cycles - last_cyc1, _rp_phase2_cycles - last_cyc2,
            _rp_check_cycles - last_cyc3,
#ifdef RP_PHASE_COUNTS
            _rp_and2 - last_and2,
#else
            0ULL,
#endif
            _rp_init_rows - last_rows);
    last_arrays = _rp_arrays_swept; last_dens = _rp_bp_dens;
    last_cyc1 = _rp_phase1_cycles; last_cyc2 = _rp_phase2_cycles;
    last_cyc3 = _rp_check_cycles; last_rows = _rp_init_rows;
    last_calls = _rp_sift0_calls;
#ifdef RP_PHASE_COUNTS
    last_and2 = _rp_and2;
#endif
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
