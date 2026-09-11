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
 * bench_check.c                                                       *
 *                                                                     *
 * Benchmark for the exact check, and the way to remeasure the four    *
 * constants RATPOINTS_CHECK_STEP, _LIMB, _CALL and _ROOT.  Not part   *
 * of the library; build it with                                       *
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
 * only the marginal cost of one more check.                          *
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
 * formula is within 6% (see PARAM-NOTES.md).                          *
 ***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include "ratpoints.h"

#if defined(__x86_64__) || defined(__i386__)
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
  return(0);
}
