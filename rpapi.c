/***********************************************************************
 * ratpoints-2.2                                                       *
 *  - A program to find rational points on hyperelliptic curves        *
 * Copyright (C) 2026  Michael Stoll                                   *
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
 * rpapi.c                                                             *
 *                                                                     *
 * Tests of the library interface that the command-line program cannot *
 * reach: the checks on the arguments of find_points_work() and the    *
 * error codes they return, a callback that stops the search or        *
 * declines a point, the flags set through the API, and a sequence of  *
 * searches on one initialised structure.  "make testapi" compares the *
 * output with testbase-api.  Part of the test suite of TODO item 20;  *
 * the searches themselves are covered by rptest and test4.sh.         *
 *                                                                     *
 * Michael Stoll, September 19, 2026                                   *
 ***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "ratpoints.h"

static mpz_t c[RATPOINTS_MAX_DEGREE+1];
static ratpoints_interval domain[2*RATPOINTS_MAX_DEGREE];

/* the callbacks: count, count and print, decline, stop after the first */
static long seen = 0;

static int count(long a, long b, const mpz_t y, void *info, int *quit)
{ seen++; return(1); }

static int decline(long a, long b, const mpz_t y, void *info, int *quit)
{ seen++; return(0); }

static int stop(long a, long b, const mpz_t y, void *info, int *quit)
{ seen++; *quit = 1; return(1); }

static int show(long a, long b, const mpz_t y, void *info, int *quit)
{ seen++;
  printf("    (%ld : ", a); mpz_out_str(NULL, 10, y); printf(" : %ld)\n", b);
  return(1);
}

/* set the coefficients from a list of longs and the input fields to what
 * main.c starts from */
static void curve(ratpoints_args *args, long degree, const long *cof,
                  long height)
{ long k;

  for(k = 0; k <= degree; k++) { mpz_set_si(c[k], cof[k]); }
  args->cof = c; args->degree = degree; args->height = height;
  args->domain = domain; args->num_inter = 0;
  args->b_low = 1; args->b_high = height;
  args->sp1 = -1; args->sp2 = -1;
  args->survivors_per_word = -1.0; args->sp2_extra = -1; args->sp2_u0 = -1.0;
  args->cost_table = -1.0; args->adapt = -1; args->sp3_extra = -1;
  args->sp3_per_denom = -1.0; args->check_cost = -1.0;
  args->array_size = 0; args->sturm = RATPOINTS_DEFAULT_STURM;
  args->num_primes = -1; args->max_forbidden = -1;
  args->flags = 0;
}

static void report(const char *what, long ret)
{ printf("%-56s %4ld  (%ld seen)\n", what, ret, seen); seen = 0; }

int main(void)
{ ratpoints_args args;
  long n, ret;
  /* y^2 = x^2 + 1: the points at infinity, (0 : +-1 : 1) and the
   * primitive Pythagorean triples; 12 points up to height 10 */
  static const long pyth[] = {1, 0, 1};
  /* y^2 = x^6 + x^4 + x^2 + 1 = (x^2 + 1)(x^4 + 1) */
  static const long sextic[] = {1, 0, 1, 0, 1, 0, 1};
  static const long quintic[] = {1, 0, 0, 0, 0, 1};
  static const long dropped[] = {1, 2, 0, 0};     /* x^0 .. x^3, degree 1 */
  static const long constant[] = {5};

  for(n = 0; n <= RATPOINTS_MAX_DEGREE; n++) { mpz_init(c[n]); }

  /* initialise for degree 3 -- the work space then holds a polynomial of
   * degree up to 3 -- and search a conic */
  curve(&args, 2, pyth, 10); args.degree = 3;
  find_points_init(&args);
  args.degree = 2;
  report("find_points_work: y^2 = x^2 + 1, height 10", find_points_work(&args, count, NULL));
  printf("  moduli used: %ld in the first stage, %ld in both, %ld primes in the third\n",
         args.sp1_used, args.sp2_used, args.sp3_used - args.sp2_used);
  printf("  the region searched: %ld interval(s), [%g, %g]\n", args.num_inter,
         args.domain[0].low, args.domain[0].up);

  /* the checks on the arguments */
  curve(&args, 2, pyth, 10); args.cof = NULL;
  report("cof = NULL", find_points_work(&args, count, NULL));
  curve(&args, 2, pyth, 10); args.height = 0;
  report("height 0", find_points_work(&args, count, NULL));
  curve(&args, 2, pyth, 10); args.height = -7;
  report("height -7", find_points_work(&args, count, NULL));
  curve(&args, 2, pyth, 10); args.domain = NULL;
  report("domain = NULL", find_points_work(&args, count, NULL));
  curve(&args, 0, constant, 10);
  report("degree 0", find_points_work(&args, count, NULL));
  curve(&args, 2, pyth, 10); args.degree = -1;
  report("degree -1", find_points_work(&args, count, NULL));
  curve(&args, 2, pyth, 10); mpz_set_si(c[1], 0); mpz_set_si(c[2], 0);
  report("degree 2 with c[1] = c[2] = 0 (degree 0 after stripping)",
         find_points_work(&args, count, NULL));
  /* the work space was sized for degree 2 */
  curve(&args, 6, sextic, 10);
  report("degree 6 after initialising for degree 3",
         find_points_work(&args, count, NULL));
  curve(&args, 3, dropped, 20);
  report("degree 3 with two leading zeros (genus drops)",
         find_points_work(&args, count, NULL));
  /* the fields that are normalised */
  curve(&args, 2, pyth, 10); args.num_inter = -3;
  report("num_inter -3", find_points_work(&args, count, NULL));
  curve(&args, 2, pyth, 10); args.b_low = -2; args.b_high = 0;
  report("b_low -2, b_high 0", find_points_work(&args, count, NULL));
  curve(&args, 2, pyth, 10); args.b_low = 5; args.b_high = 3;
  report("b_low 5, b_high 3 (only the points at infinity)",
         find_points_work(&args, count, NULL));
  curve(&args, 2, pyth, 10); args.array_size = 1;
  report("array_size 1", find_points_work(&args, count, NULL));
  curve(&args, 2, pyth, 10); args.sturm = 1000;
  report("sturm 1000", find_points_work(&args, count, NULL));
  curve(&args, 2, pyth, 10); args.num_primes = 1000; args.sp2 = 500; args.sp1 = 600;
  report("num_primes 1000, sp2 500, sp1 600",
         find_points_work(&args, count, NULL));
  printf("  moduli used: %ld in the first stage, %ld in both\n", args.sp1_used, args.sp2_used);
  curve(&args, 2, pyth, 10); args.max_forbidden = 1000;
  report("max_forbidden 1000", find_points_work(&args, count, NULL));
  /* the input fields come back as they went in */
  curve(&args, 2, pyth, 10); args.sp1 = -1; args.sp2 = -1; args.num_primes = -1;
  args.b_low = 0; args.b_high = -1; args.sturm = 100; args.num_inter = -1;
  ret = find_points_work(&args, count, NULL);
  report("a search with the choices left to the library", ret);
  printf("  sp1 %ld sp2 %ld num_primes %ld max_forbidden %ld b_low %ld b_high %ld"
         " array_size %ld sturm %ld flags %#x, domain %s, %ld interval(s) searched\n",
         args.sp1, args.sp2, args.num_primes, args.max_forbidden, args.b_low,
         args.b_high, args.array_size, args.sturm,
         args.flags & RATPOINTS_FLAGS_INPUT_MASK,
         (args.domain == domain) ? "unchanged" : "CHANGED", args.num_inter);
  /* the flags */
  curve(&args, 2, pyth, 10); args.flags = RATPOINTS_NO_Y;
  report("RATPOINTS_NO_Y (one point of each pair)",
         find_points_work(&args, count, NULL));
  curve(&args, 2, pyth, 10); args.flags = RATPOINTS_NO_CHECK;
  args.sp1 = 15; args.sp2 = 30; args.num_primes = 30; /* what the sieve
                                   leaves depends on the moduli: pin them */
  report("RATPOINTS_NO_CHECK (the survivors, y = 0)",
         find_points_work(&args, count, NULL));
  curve(&args, 2, pyth, 10); args.flags = RATPOINTS_NO_REVERSE | RATPOINTS_NO_JACOBI;
  report("RATPOINTS_NO_REVERSE | RATPOINTS_NO_JACOBI",
         find_points_work(&args, count, NULL));
  /* the callback */
  curve(&args, 2, pyth, 10);
  report("a callback that declines every point", find_points_work(&args, decline, NULL));
  curve(&args, 2, pyth, 10);
  report("a callback that stops at the first point", find_points_work(&args, stop, NULL));
  curve(&args, 2, pyth, 10); args.num_inter = 1;
  domain[0].low = 0.5; domain[0].up = 2.0;
  report("the same, on the interval [0.5, 2]", find_points_work(&args, stop, NULL));
  /* a search within the intervals given */
  curve(&args, 2, pyth, 10); args.num_inter = 2;
  domain[0].low = -2.0; domain[0].up = -1.0; domain[1].low = 0.0; domain[1].up = 1.0;
  ret = find_points_work(&args, show, NULL);
  report("y^2 = x^2 + 1 on [-2, -1] U [0, 1]", ret);
  printf("  the region searched: %ld interval(s)\n", args.num_inter);
  find_points_clear(&args);

  /* the wrapper, initialising for the degree of each curve */
  curve(&args, 5, quintic, 30);
  report("find_points: y^2 = x^5 + 1, height 30", find_points(&args, count, NULL));
  curve(&args, 6, sextic, 30);
  report("find_points: y^2 = (x^2 + 1)(x^4 + 1), height 30",
         find_points(&args, count, NULL));
  curve(&args, 6, sextic, 30); args.num_inter = 1;
  domain[0].low = -0.5; domain[0].up = 0.5;
  report("the same on [-0.5, 0.5]", find_points(&args, show, NULL));

  for(n = 0; n <= RATPOINTS_MAX_DEGREE; n++) { mpz_clear(c[n]); }
  return(0);
}
