/***********************************************************************
 * ratpoints-3.2.0                                                     *
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
 * bench_check.c                                                       *
 *                                                                     *
 * Benchmark for the exact check, and the way to remeasure the four    *
 * constants RATPOINTS_CHECK_STEP, _LIMB, _CALL and _ROOT, and the two *
 * RATPOINTS_TEST_STEP and _CALL of the modular test before it.  Not   *
 * part of the library; build it with                                  *
 *   make bench_check                                                  *
 * and run it as                                                       *
 *   ./bench_check [reps]                                              *
 *                                                                     *
 * It performs exactly the gmp calls that _ratpoints_check_point makes *
 * for one surviving numerator -- the Horner loop over the             *
 * bc[k] = c[k]*b^(degree-k), the further multiplication by b for an   *
 * odd degree, the sign test and mpz_sqrtrem -- over a range of        *
 * degrees, coefficient sizes and height bounds, and prints beside     *
 * each measurement what the formula in find_points.c predicts.        *
 *                                                                     *
 * The per-denominator part, forming the bc[k] themselves, is set up   *
 * outside the timed loop, because the rule the constants feed charges *
 * only the marginal cost of one more check.                           *
 *                                                                     *
 * What to look at, and what not to.  Only the RATIOS between the      *
 * rows are used by the program, so a machine on which every check is  *
 * dearer needs no change to the constants; the absolute figures will  *
 * not match the formula and are not meant to.  What would need a      *
 * change is a gmp whose square root has its thresholds elsewhere:     *
 * that shows up as a jump in the measured column that does not line   *
 * up with a change in the limbs column, and the constant to move is   *
 * RATPOINTS_CHECK_ROOT.                                               *
 *                                                                     *
 * The error column is not expected to be zero.  The formula carries   *
 * two things a hot loop does not.  In the sieve a check is a cold     *
 * excursion out of the inner loop, and the fixed part of that is what *
 * RATPOINTS_CHECK_CALL holds, so the formula's ratios are flatter     *
 * than the ones here; and an odd degree restricts the denominators to *
 * squares, which leaves fewer survivors per denominator and so a      *
 * colder check still, which the odd-degree term absorbs.  On the      *
 * machine the shipped constants were measured on, the error column    *
 * runs from +37% on the odd-degree rows with small coefficients to    *
 * -13% at the top of the degree range.  Against the sieve's own       *
 * per-check figures, which is the comparison that matters, the same   *
 * formula is within 6% (there a check costs 44 cycles plus 1.44 times *
 * what this benchmark reports).                                       *
 *                                                                     *
 * The second table is the modular test (_ratpoints_check_point): one  *
 * modulus q of the sixteen, that is the reductions of a and b modulo  *
 * q, the loop over the residues of the coefficients and jacobi1(q, v) *
 * -- the same lines, copied from find_points.c --, over degrees 2 to  *
 * 100000, on random residues, so that each modulus rejects half of    *
 * the survivors, as on a real search; the cost is what the whole test *
 * costs a survivor over the moduli it tried.  "constant" reduces as   *
 * the check is compiled, with the loop over the moduli unrolled, so   *
 * that each modulus is a constant and the reductions multiply,        *
 * "divide" with the modulus read from a variable, which divides.      *
 * Beside them the formula, degree*TEST_STEP + TEST_CALL, and both as  *
 * ratios to the exact check of the reference row above, measured in   *
 * the same run: those are the ratios the switch between the two kinds *
 * of check compares.  Then jacobi1 alone, for each q and random v,    *
 * and the line degree*STEP + CALL fitted to the "constant" column.    *
 * The two constants are fitted to the sieve's own figures, where the  *
 * test of a survivor costs about 1.37 times what this benchmark       *
 * reports and the exact check of the reference curve 44 cycles plus   *
 * 1.44 times its figure here, and where its data are colder, so the   *
 * error column is not zero either: on the machine they were measured  *
 * on it runs from about -25% at low degree, where the test is mostly  *
 * jacobi1, to -8% at the top, where it is the loop.                   *
 ***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <gmp.h>
#include "ratpoints.h"

#ifdef __x86_64__
# include <x86intrin.h>
# define TICK() __rdtsc()
# define TICKNAME "rdtsc cycles"
#else
# include <time.h>
static unsigned long long TICK(void)
{ struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return((unsigned long long)t.tv_sec*1000000000ULL + (unsigned long long)t.tv_nsec);
}
# define TICKNAME "nanoseconds"
#endif

#define MAXDEG 24
#define NA 4096
#define LONG_LENGTH 64   /* as in rp-private.h */

static mpz_t c[MAXDEG+1], bc[MAXDEG+1], work[3];

/* the same formula as check_cost() in find_points.c */
static double predict(long degree, double cbits, double hbits)
{ double fbits = cbits + (double)degree*hbits;
  double limbs = fbits/64.0;
  double mid = 0.5*(cbits + fbits)/64.0;
  double root = ceil(0.5*limbs);
  double cost;

  if(root < 1.0) { root = 1.0; }
  cost = (double)degree*(RATPOINTS_CHECK_STEP + RATPOINTS_CHECK_LIMB*mid)
          + RATPOINTS_CHECK_CALL + RATPOINTS_CHECK_ROOT*(root - 1.0);
  if(degree & 1)
  { cost += RATPOINTS_CHECK_STEP + RATPOINTS_CHECK_LIMB*limbs; }
  return(cost);
}

/* the moduli of the modular test and the two functions it calls, as in
 * find_points.c */
#define RATPOINTS_TEST_MODULI 16
static const unsigned long test_moduli[RATPOINTS_TEST_MODULI] =
  { 2147483647, 2147483629, 2147483587, 2147483579, 2147483563, 2147483549,
    2147483543, 2147483497, 2147483489, 2147483477, 2147483423, 2147483399,
    2147483353, 2147483323, 2147483269, 2147483249 };
#if (defined(__GNUC__) && __GNUC__ >= 8) || defined(__clang__)
# define RP_UNROLL_TEST _Pragma("GCC unroll 16")
#else
# define RP_UNROLL_TEST
#endif

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

static inline int jacobi1(long b, const long lcf)
{
  long f;
  int neg = 0;

  while((b & 1) == 0) b >>= 1;
  f = lcf;
  if(f < 0)
  { if(f == LONG_MIN) { f /= 4; }
    f = -f; neg = 1;
  }
  if(b < 1UL<<(LONG_LENGTH - 5)) f = mod(f, b);
  if(f == 0) return(1);

  while(1)
  { long s = (neg && (b & 2)) ? -1 : 1;
    long n = f;
    long m = b;

    while(!(n & 1))
    { if(m & 2) s = -s;
      if(m & 4) s = -s;
      n >>= 1;
    }
    while(1)
    { if(n & m & 2) s = -s;
      while(m > n)
      { m -= n;
        do
        { if(n & 2) s = -s;
          if(n & 4) s = -s;
          m >>= 1;
        }
        while(!(m & 1));
      }
      if(m == n)
      { if(m == 1) return(s);
        b /= m;
        if(f == 0) return(1);
        break;
      }
      if(n & m & 2) s = -s;
      while(n > m)
      { n -= m;
        do
        { if(m & 2) s = -s;
          if(m & 4) s = -s;
          n >>= 1;
        }
        while(!(n & 1));
      }
      if(m == n)
      { if(m == 1) return(s);
        b /= m;
        if(f == 0) return(1);
        break;
      }
    }
  }
}

/* The test of one survivor (a, b): the moduli it tries, up to the first
 * modulo which F(a, b) is not a square, with the moduli constant as in the
 * check, and the same with them read from qvar, which makes the compiler
 * divide. */
static unsigned long qvar[RATPOINTS_TEST_MODULI];

static long __attribute__((noinline))
  test_constant(const unsigned long *cq0, long degree, long a, long b)
{ long j, k;

  RP_UNROLL_TEST
  for(j = 0; j < RATPOINTS_TEST_MODULI; j++)
  { const unsigned long q = test_moduli[j];
    const unsigned long *cq = cq0 + j*(degree + 1);
    long r = a % (long)q;
    unsigned long aq, bq = (unsigned long)b % q, v = cq[degree], bp = 1;

    if(r < 0) { r += (long)q; }
    aq = (unsigned long)r;
    for(k = degree - 1; k >= 0; k--)
    { bp = bp*bq % q; v = (v*aq + cq[k]*bp) % q; }
    if(degree & 1) { v = v*bq % q; }
    if(jacobi1((long)q, (long)v) < 0) { return(j + 1); }
  }
  return(RATPOINTS_TEST_MODULI);
}

static long __attribute__((noinline))
  test_divide(const unsigned long *cq0, long degree, long a, long b)
{ long j, k;

  for(j = 0; j < RATPOINTS_TEST_MODULI; j++)
  { const unsigned long q = qvar[j];
    const unsigned long *cq = cq0 + j*(degree + 1);
    long r = a % (long)q;
    unsigned long aq, bq = (unsigned long)b % q, v = cq[degree], bp = 1;

    if(r < 0) { r += (long)q; }
    aq = (unsigned long)r;
    for(k = degree - 1; k >= 0; k--)
    { bp = bp*bq % q; v = (v*aq + cq[k]*bp) % q; }
    if(degree & 1) { v = v*bq % q; }
    if(jacobi1((long)q, (long)v) < 0) { return(j + 1); }
  }
  return(RATPOINTS_TEST_MODULI);
}

static int __attribute__((noinline)) jacobi_alone(long q, long v)
{ return(jacobi1(q, v)); }
static volatile long sink; /* what is computed, so that it is computed */

/* What the test costs per modulus tried at this degree, by each of the two
 * reductions (into *divided the second), and the mean number of moduli a
 * survivor tried. */
static double test_row(long degree, long reps, double *divided,
                       double *tried, gmp_randstate_t rs)
{ long n = (degree < 10000) ? NA/4 : NA/64; /* survivors */
  long a[NA/4], b = 196613, i, j, k, r, nmod = 0;
  unsigned long *cq = malloc(RATPOINTS_TEST_MODULI*(degree + 1)
                             *sizeof(unsigned long));
  unsigned long long t, best = ~0ULL, best_d = ~0ULL;

  if(cq == NULL) { fprintf(stderr, "bench_check: no memory\n"); exit(1); }
  gmp_randseed_ui(rs, 1000*(unsigned long)degree + 1);
  for(j = 0; j < RATPOINTS_TEST_MODULI; j++)
  { for(k = 0; k <= degree; k++)
    { cq[j*(degree + 1) + k] = gmp_urandomm_ui(rs, test_moduli[j]); }
  }
  for(i = 0; i < n; i++)
  { a[i] = (long)gmp_urandomb_ui(rs, 18) - (1L << 17); }
  for(i = 0; i < n; i++) { nmod += test_constant(cq, degree, a[i], b); }
  if(degree >= 1000) { reps = reps/10 + 1; }
  for(r = 0; r < reps; r++)
  { long s = 0;

    t = TICK();
    for(i = 0; i < n; i++) { s += test_constant(cq, degree, a[i], b); }
    t = TICK() - t;
    if(s != nmod) { fprintf(stderr, "bench_check: the test is not the same\n"); exit(1); }
    if(t < best) { best = t; }
    s = 0;
    t = TICK();
    for(i = 0; i < n; i++) { s += test_divide(cq, degree, a[i], b); }
    t = TICK() - t;
    if(s != nmod) { fprintf(stderr, "bench_check: the two tests differ\n"); exit(1); }
    if(t < best_d) { best_d = t; }
  }
  free(cq);
  *divided = (double)best_d/(double)nmod;
  *tried = (double)nmod/(double)n;
  return((double)best/(double)nmod);
}

static double one_row(long degree, long cbits, long hbits, long reps,
                      gmp_randstate_t rs)
{ long k, i, r, a[NA], b;
  unsigned long long t, best = ~0ULL;

  /* the same coefficients every time this configuration is asked for, so
   * that the reference row below is the row it is compared against */
  gmp_randseed_ui(rs, 1000000*(unsigned long)degree
                       + 1000*(unsigned long)cbits + (unsigned long)hbits);
  for(k = 0; k <= degree; k++)
  { mpz_urandomb(c[k], rs, cbits);
    mpz_add_ui(c[k], c[k], 1);   /* positive, so F(a,b) >= 0 and the square
                                  * root always runs, as it nearly always
                                  * does in the sieve: the search is confined
                                  * to the intervals where f is positive */
  }
  b = (1L << (hbits-1)) + 1;     /* odd, like a real denominator */
  for(i = 0; i < NA; i++)
  { a[i] = (long)(gmp_urandomb_ui(rs, hbits)) + 1; }

  mpz_set_si(work[0], 1);
  for(k = degree-1; k >= 0; k--)
  { mpz_mul_ui(work[0], work[0], b);
    mpz_mul(bc[k], c[k], work[0]);
  }

  for(r = 0; r < reps; r++)
  { t = TICK();
    for(i = 0; i < NA; i++)
    { mpz_set(work[2], c[degree]);
      for(k = degree-1; k >= 0; k--)
      { mpz_mul_si(work[2], work[2], a[i]);
        mpz_add(work[2], work[2], bc[k]);
      }
      if(degree & 1) { mpz_mul_ui(work[2], work[2], b); }
      if(mpz_cmp_si(work[2], 0) >= 0)
      { mpz_sqrtrem(work[0], work[1], work[2]);
        if(mpz_cmp_si(work[1], 0) == 0) { mpz_add_ui(work[0], work[0], 0); }
      }
    }
    t = TICK() - t;
    if(t < best) { best = t; }
  }
  return((double)best/NA);
}

int main(int argc, char **argv)
{ static const long degrees[] = {3, 4, 5, 6, 7, 8, 10, 12, 16, 20, 0};
  static const long sizes[] = {4, 20, 60, 120, 250, 0};
  long reps = (argc > 1) ? atol(argv[1]) : 150;
  long k, hbits;
  double ref_meas = 0.0, ref_pred = 0.0;
  gmp_randstate_t rs;

  gmp_randinit_default(rs);
  gmp_randseed_ui(rs, 20260911);
  for(k = 0; k <= MAXDEG; k++) { mpz_init(c[k]); mpz_init(bc[k]); }
  for(k = 0; k < 3; k++) { mpz_init(work[k]); }

  printf("What one exact check costs, in %s, hot cache, best of %ld.\n",
         TICKNAME, reps);
  printf("\"measured\" is this program, \"formula\" is check_cost() in"
         " find_points.c;\nonly the ratios matter, and the reference row is"
         " degree 6 with 4-bit\ncoefficients at each height bound.\n\n");
  printf("%5s %6s %6s %5s %10s %8s %9s %8s %7s\n",
         "hbits", "degree", "cbits", "limbs", "measured", "rel", "formula",
         "rel", "error");

  for(hbits = 14; hbits <= 18; hbits += 4)
  { long di, si, nd = 0, ns = 0;
    double meas[8][12], pred[8][12];

    while(degrees[nd]) { nd++; }
    while(sizes[ns]) { ns++; }
    for(si = 0; si < ns; si++)
    { for(di = 0; di < nd; di++)
      { meas[si][di] = one_row(degrees[di], sizes[si], hbits, reps, rs);
        pred[si][di] = predict(degrees[di], (double)sizes[si], (double)hbits);
    } }
    /* every row as a ratio to degree 6 with 4-bit coefficients, which is
     * the curve RATPOINTS_CHECK_REFERENCE stands for */
    ref_meas = meas[0][3]; ref_pred = pred[0][3];
    for(si = 0; si < ns; si++)
    { for(di = 0; di < nd; di++)
      { double fbits = (double)sizes[si] + (double)degrees[di]*(double)hbits;

        printf("%5ld %6ld %6ld %5ld %10.1f %8.3f %9.1f %8.3f %+6.1f%%\n",
               hbits, degrees[di], sizes[si], (long)ceil(fbits/64.0),
               meas[si][di], meas[si][di]/ref_meas,
               pred[si][di], pred[si][di]/ref_pred,
               100.0*((pred[si][di]/ref_pred)/(meas[si][di]/ref_meas) - 1.0));
      }
      printf("\n");
    }
  }

  /* the modular test, against the reference row of the exact check at the
   * smaller height bound, measured again here so that both are taken at
   * the same moment */
  { static const long tdeg[] = {2, 3, 4, 6, 8, 12, 16, 20, 30, 50, 100, 200,
                                500, 1000, 3000, 10000, 30000, 100000, 0};
    double meas[32], sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0, sw = 0.0;
    double ref = one_row(6, 4, 14, reps, rs), ref_p = predict(6, 4.0, 14.0);
    long di, i;

    for(i = 0; i < RATPOINTS_TEST_MODULI; i++) { qvar[i] = test_moduli[i]; }
    printf("\nThe modular test: what one modulus costs, the cycles of the"
           " test of a survivor\nover the moduli it tried, \"constant\" as"
           " the check reduces, \"divide\" with\na division; \"formula\" is"
           " degree*RATPOINTS_TEST_STEP + RATPOINTS_TEST_CALL;\nthe ratios"
           " are to the exact check of degree 6 with 4-bit coefficients at\n"
           "hbits 14 (%.1f cycles here, %.1f by its formula).\n\n", ref, ref_p);
    printf("%6s %10s %10s %6s %9s %8s %9s %8s %7s\n", "degree", "constant",
           "divide", "tried", "formula", "rel", "formula", "rel", "error");
    for(di = 0; tdeg[di]; di++)
    { double div, tried, f;
      long d = tdeg[di];

      meas[di] = test_row(d, reps, &div, &tried, rs);
      f = (double)d*RATPOINTS_TEST_STEP + RATPOINTS_TEST_CALL;
      printf("%6ld %10.1f %10.1f %6.3f %9.1f %8.3f %9.1f %8.3f %+6.1f%%\n",
             d, meas[di], div, tried, f, meas[di]/ref, f, f/ref_p,
             100.0*((f/ref_p)/(meas[di]/ref) - 1.0));
      /* least squares of the relative error */
      sw += 1.0/(meas[di]*meas[di]); sx += d/(meas[di]*meas[di]);
      sxx += (double)d*d/(meas[di]*meas[di]); sy += 1.0/meas[di];
      sxy += d/meas[di];
    }
    { double det = sxx*sw - sx*sx;

      printf("\nfitted to the \"constant\" column: one modulus costs %.1f"
             " + %.2f*degree cycles here\n", (sxx*sy - sx*sxy)/det,
             (sxy*sw - sx*sy)/det);
    }
    { long q;
      double sum = 0.0;

      gmp_randseed_ui(rs, 20261122);
      printf("jacobi1(q, v) for random v below q, by q:");
      for(q = 0; q < RATPOINTS_TEST_MODULI; q++)
      { long v[NA], k, r;
        unsigned long long t, best = ~0ULL;
        int s = 0;

        for(k = 0; k < NA; k++)
        { v[k] = (long)gmp_urandomm_ui(rs, test_moduli[q]); }
        for(r = 0; r < reps; r++)
        { t = TICK();
          for(k = 0; k < NA; k++) { s += jacobi_alone((long)test_moduli[q], v[k]); }
          t = TICK() - t;
          if(t < best) { best = t; }
        }
        sink += s;
        printf("%s %.0f", (q % 8 == 0) ? "\n " : "", (double)best/NA);
        sum += (double)best/NA;
      }
      printf("\n mean %.1f cycles\n", sum/RATPOINTS_TEST_MODULI);
    }
  }
  return(0);
}
