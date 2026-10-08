/***********************************************************************
 * ratpoints-3.1.0                                                     *
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
 * output with testbase-api.  The searches themselves are covered by   *
 * rptest and test4.sh.                                                *
 *                                                                     *
 * Michael Stoll, September 19, 2026                                   *
 ***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <limits.h>
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

/* keep the points of a run, to compare two runs point by point: what the
 * library delivers must not depend on the number of threads */
#define KEPT_MAX 8192
typedef struct { long a; long b; mpz_t y; } kept_point;
static kept_point kept[2][KEPT_MAX];
static long nkept[2] = {0, 0};
static int keep_to = 0;

static int keep(long a, long b, const mpz_t y, void *info, int *quit)
{ long n = nkept[keep_to];

  if(n < KEPT_MAX)
  { kept[keep_to][n].a = a; kept[keep_to][n].b = b;
    mpz_init_set(kept[keep_to][n].y, y);
    nkept[keep_to] = n + 1;
  }
  seen++;
  return(1);
}

/* the same, stopping at the first point */
static int keep_one(long a, long b, const mpz_t y, void *info, int *quit)
{ keep(a, b, y, info, quit); *quit = 1; return(1); }

/* the same, stopping at the first point with a denominator of at least
 * stop_b */
static long stop_b = 0;

static int keep_until(long a, long b, const mpz_t y, void *info, int *quit)
{ keep(a, b, y, info, quit); if(b >= stop_b) { *quit = 1; } return(1); }

/* whether the two kept runs agree, point by point and in order */
static void same_points(void)
{ long n, same = (nkept[0] == nkept[1]);

  for(n = 0; same && n < nkept[0]; n++)
  { same = (kept[0][n].a == kept[1][n].a && kept[0][n].b == kept[1][n].b
              && mpz_cmp(kept[0][n].y, kept[1][n].y) == 0);
  }
  printf("  the two runs: %ld and %ld points, %s\n", nkept[0], nkept[1],
         same ? "the same in the same order" : "DIFFERENT");
  for(n = 0; n < nkept[0]; n++) { mpz_clear(kept[0][n].y); }
  for(n = 0; n < nkept[1]; n++) { mpz_clear(kept[1][n].y); }
  nkept[0] = 0; nkept[1] = 0; keep_to = 0;
}

/* whether the run kept second found what the one kept first did, point by
 * point and in order; the second is emptied for the next run, the first
 * stays */
static int same_as_first(void)
{ long n, same = (nkept[0] == nkept[1]);

  for(n = 0; same && n < nkept[0]; n++)
  { same = (kept[0][n].a == kept[1][n].a && kept[0][n].b == kept[1][n].b
              && mpz_cmp(kept[0][n].y, kept[1][n].y) == 0);
  }
  for(n = 0; n < nkept[1]; n++) { mpz_clear(kept[1][n].y); }
  nkept[1] = 0; seen = 0;
  return(same);
}

/* the tuning fields, by number: those of type double and those of type
 * long (sturm and num_threads have tests of their own) */
static const char *const dfields[] = {"survivors_per_word", "sp2_u0",
                                      "cost_table", "sp3_per_denom",
                                      "check_cost"};
#define NDFIELDS 5
static void set_dfield(ratpoints_args *args, long f, double v)
{ switch(f)
  { case 0: args->survivors_per_word = v; break;
    case 1: args->sp2_u0 = v; break;
    case 2: args->cost_table = v; break;
    case 3: args->sp3_per_denom = v; break;
    case 4: args->check_cost = v; break;
  }
}
static const char *const lfields[] = {"sp1", "sp2", "sp2_extra", "sp3_extra",
                                      "adapt", "array_size", "num_primes",
                                      "max_forbidden"};
#define NLFIELDS 8
static void set_lfield(ratpoints_args *args, long f, long v)
{ switch(f)
  { case 0: args->sp1 = v; break;
    case 1: args->sp2 = v; break;
    case 2: args->sp2_extra = v; break;
    case 3: args->sp3_extra = v; break;
    case 4: args->adapt = v; break;
    case 5: args->array_size = v; break;
    case 6: args->num_primes = v; break;
    case 7: args->max_forbidden = v; break;
  }
}

/* nearly every one at once: those of type double at v, those of type long
 * at LONG_MAX but for sp2 and num_primes, which stay with the library */
static void all_extreme(ratpoints_args *args, double v)
{ long f;

  for(f = 0; f < NDFIELDS; f++) { set_dfield(args, f, v); }
  for(f = 0; f < NLFIELDS; f++) { set_lfield(args, f, LONG_MAX); }
  args->sp2 = -1; args->num_primes = -1;
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
  args->num_threads = 1;
  args->flags = 0;
}

/* y^2 = x - 9223372036854769900 at the height bound 2^63 - 1, over the
 * window at the top of the range that test4.sh searches as well: the
 * denominator 1 and the numerators from 2^63 - 1024 on, where fourteen
 * points lie, (9223372036854769900 + y^2 : +-y : 1) for y = 70, ..., 76,
 * besides the one at infinity */
static void top_window(ratpoints_args *args)
{ static const long top[] = {-9223372036854769900L, 1};

  curve(args, 1, top, LONG_MAX);
  args->b_high = 1; args->num_inter = 1;
  domain[0].low = 9223372036854774784.0; domain[0].up = 1e30;
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
  /* the five machine constants pinned, so that the report of what was
   * used depends neither on the values compiled in nor on a tuning.mk */
  args.survivors_per_word = 0.003; args.sp2_extra = 11; args.sp2_u0 = 1.6e6;
  args.cost_table = 38.0; args.sp3_per_denom = 0.013;
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
  /* beyond RATPOINTS_MAX_DEGREE, refused before anything is sized by it
   * (and before the coefficients are read, which c[] does not hold) */
  curve(&args, 2, pyth, 10); args.degree = RATPOINTS_MAX_DEGREE + 1;
  report("degree RATPOINTS_MAX_DEGREE + 1", find_points_work(&args, count, NULL));
  curve(&args, 2, pyth, 10); args.degree = LONG_MAX;
  report("degree LONG_MAX", find_points_work(&args, count, NULL));
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
  /* The intervals must be closed, l <= u, increasing and disjoint, with no
   * end NaN, as the command line has them; anything else is refused, and
   * num_inter and domain[] are left as they were.  At height 20 the curve
   * has two points at infinity and fourteen with x in [0, 3], six with x in
   * [-3, -1] and none with x in [5, 30]. */
  curve(&args, 2, pyth, 20); args.num_inter = 1;
  domain[0].low = 0.0; domain[0].up = 3.0;
  report("height 20 on [0, 3]", find_points_work(&args, count, NULL));
  curve(&args, 2, pyth, 20); args.num_inter = 2;
  domain[0].low = 0.0; domain[0].up = 2.0; domain[1].low = 1.0; domain[1].up = 3.0;
  report("height 20 on [0, 2], [1, 3] (overlapping)",
         find_points_work(&args, count, NULL));
  printf("  num_inter %ld, domain [%g, %g], [%g, %g]\n", args.num_inter,
         domain[0].low, domain[0].up, domain[1].low, domain[1].up);
  curve(&args, 2, pyth, 20); args.num_inter = 2;
  domain[0].low = 0.0; domain[0].up = 1.5; domain[1].low = 1.5; domain[1].up = 3.0;
  report("height 20 on [0, 1.5], [1.5, 3] (touching)",
         find_points_work(&args, count, NULL));
  curve(&args, 2, pyth, 20); args.num_inter = 2;
  domain[0].low = -3.0; domain[0].up = -1.0; domain[1].low = 5.0; domain[1].up = 30.0;
  report("height 20 on [-3, -1], [5, 30]", find_points_work(&args, count, NULL));
  curve(&args, 2, pyth, 20); args.num_inter = 2;
  domain[0].low = 5.0; domain[0].up = 30.0; domain[1].low = -3.0; domain[1].up = -1.0;
  report("height 20 on [5, 30], [-3, -1] (decreasing)",
         find_points_work(&args, count, NULL));
  curve(&args, 2, pyth, 20); args.num_inter = 1;
  domain[0].low = 2.0; domain[0].up = 1.0;
  report("height 20 on [2, 1] (empty)", find_points_work(&args, count, NULL));
  curve(&args, 2, pyth, 20); args.num_inter = 1;
  domain[0].low = 0.0; domain[0].up = NAN;
  report("height 20 on [0, NaN]", find_points_work(&args, count, NULL));
  curve(&args, 2, pyth, 20); args.num_inter = 2; args.sturm = -1;
  domain[0].low = NAN; domain[0].up = 1.0; domain[1].low = 2.0; domain[1].up = 3.0;
  report("height 20 on [NaN, 1], [2, 3], no Sturm step",
         find_points_work(&args, count, NULL));
  curve(&args, 2, pyth, 20); args.num_inter = 1;
  domain[0].low = -HUGE_VAL; domain[0].up = HUGE_VAL;
  report("height 20 on [-infinity, infinity]", find_points_work(&args, count, NULL));
  /* The default interval at a height bound that is not a double has its
   * ends rounded outward, so that it holds x = height: y^2 = x - 2^53 at
   * height 2^53 + 1 with the denominator 1 has the points x = 2^53 and
   * x = 2^53 + 1 and the one at infinity (the Sturm step to depth 55 cuts
   * the search down to the numerators from 2^53 on). */
  { static const long shifted[] = {-9007199254740992L, 1};

    curve(&args, 1, shifted, 9007199254740993L); args.b_high = 1; args.sturm = 55;
    report("y^2 = x - 2^53, height 2^53 + 1, denominator 1",
           find_points_work(&args, show, NULL));
  }
  /* Any number of intervals: 600000 of them, [n, n + 1/2] for n from
   * -300000 on, without the Sturm step and with it, which keeps a copy of
   * them.  Of the points of height up to 20, those with x = 0, 5/12, 4/3,
   * 12/5, -3/4, -8/15 and -15/8 lie in the intervals, and those with x = 3/4,
   * 8/15, 15/8, -5/12, -4/3 and -12/5 between them: fourteen points and the
   * two at infinity. */
  { long nbig = 600000, k;
    ratpoints_interval *big = malloc((nbig + 2)*sizeof(ratpoints_interval));

    if(big == NULL) { printf("no memory for the intervals\n"); return(1); }
    for(k = 0; k < nbig; k++)
    { big[k].low = (double)(k - nbig/2); big[k].up = big[k].low + 0.5; }
    curve(&args, 2, pyth, 20); args.sturm = -1;
    args.domain = big; args.num_inter = nbig;
    report("height 20 on [n, n + 1/2], 600000 of them, no Sturm step",
           find_points_work(&args, count, NULL));
    for(k = 0; k < nbig; k++)
    { big[k].low = (double)(k - nbig/2); big[k].up = big[k].low + 0.5; }
    curve(&args, 2, pyth, 20);
    args.domain = big; args.num_inter = nbig;
    report("height 20 on [n, n + 1/2], 600000 of them",
           find_points_work(&args, count, NULL));
    free(big);
  }
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

  /* The largest degree, RATPOINTS_MAX_DEGREE (100, divisible by 4):
   * y^2 = x^100 + 1 has the two points at infinity and (0 : +-1 : 1), and
   * no other, since a^100 + b^100 = y^2 with ab != 0 would be a solution of
   * u^4 + v^4 = w^2 in nonzero integers.  A larger degree is refused by
   * find_points; and given to find_points_init, a larger or a negative one
   * makes it allocate nothing and find_points_work refuse the search,
   * whatever degree that is then given (find_points_clear is called as
   * usual). */
  { static long big[RATPOINTS_MAX_DEGREE+1];
    static const long bad[] = {RATPOINTS_MAX_DEGREE + 1, LONG_MAX, -1, LONG_MIN};
    static const char *const bad_name[] = {"RATPOINTS_MAX_DEGREE + 1", "LONG_MAX",
                                           "-1", "LONG_MIN"};
    char what[80];
    long k;

    big[0] = 1; big[RATPOINTS_MAX_DEGREE] = 1;
    curve(&args, RATPOINTS_MAX_DEGREE, big, 10);
    sprintf(what, "find_points: y^2 = x^%d + 1, height 10", RATPOINTS_MAX_DEGREE);
    report(what, find_points(&args, show, NULL));
    curve(&args, 2, pyth, 10); args.degree = RATPOINTS_MAX_DEGREE + 1;
    report("find_points: degree RATPOINTS_MAX_DEGREE + 1", find_points(&args, count, NULL));
    curve(&args, 2, pyth, 10); args.degree = LONG_MAX;
    report("find_points: degree LONG_MAX", find_points(&args, count, NULL));
    for(k = 0; k < 4; k++)
    { curve(&args, 2, pyth, 10); args.degree = bad[k];
      find_points_init(&args);
      args.degree = 2;
      sprintf(what, "init for degree %s, then degree 2", bad_name[k]);
      report(what, find_points_work(&args, count, NULL));
      find_points_clear(&args);
    }
  }

  /* Threads.  Each search is run on one thread and on three, and the
   * points of the two runs are compared one by one: the library must
   * deliver the same points in the same order whatever the number of
   * threads, the same total, and under a callback that stops the search
   * the same first point.  (A library built without threads sieves on one
   * thread whatever num_threads says, and passes the same way.) */
  curve(&args, 2, pyth, 2000); args.degree = 6;
  find_points_init(&args);
  args.degree = 2;
  keep_to = 0; report("y^2 = x^2 + 1, height 2000, one thread",
                      find_points_work(&args, keep, NULL));
  curve(&args, 2, pyth, 2000); args.num_threads = 3;
  keep_to = 1; report("the same on three threads", find_points_work(&args, keep, NULL));
  printf("  num_threads after the search: %ld\n", args.num_threads);
  same_points();
  curve(&args, 6, sextic, 3000);
  keep_to = 0; report("(x^2 + 1)(x^4 + 1), height 3000, one thread",
                      find_points_work(&args, keep, NULL));
  curve(&args, 6, sextic, 3000); args.num_threads = 3;
  keep_to = 1; report("the same on three threads", find_points_work(&args, keep, NULL));
  same_points();
  curve(&args, 5, quintic, 500);
  keep_to = 0; report("x^5 + 1, height 500, stop at the first point, one thread",
                      find_points_work(&args, keep_one, NULL));
  curve(&args, 5, quintic, 500); args.num_threads = 3;
  keep_to = 1; report("the same on three threads", find_points_work(&args, keep_one, NULL));
  same_points();
  curve(&args, 2, pyth, 300); args.flags = RATPOINTS_NO_CHECK;
  keep_to = 0; report("the survivors (RATPOINTS_NO_CHECK), height 300, one thread",
                      find_points_work(&args, keep, NULL));
  curve(&args, 2, pyth, 300); args.flags = RATPOINTS_NO_CHECK; args.num_threads = 3;
  keep_to = 1; report("the same on three threads", find_points_work(&args, keep, NULL));
  same_points();
  curve(&args, 2, pyth, 2000); args.num_threads = -1;
  report("as many threads as there are processors", find_points_work(&args, count, NULL));
  curve(&args, 2, pyth, 2000); args.num_threads = 2;
  keep_to = 0; report("two threads after the pool of three", find_points_work(&args, keep, NULL));
  curve(&args, 2, pyth, 2000); args.num_threads = 1;
  keep_to = 1; report("and one thread after that", find_points_work(&args, keep, NULL));
  same_points();
  /* A search that the callback stops in a block sieved with fewer primes
   * than the block before it: on this curve of many points the correction
   * during the run drops a prime of the third stage, and the first point
   * with b >= 2700, (846 : +-y : 2755), lies in the first block sieved
   * without it.  The moduli reported must be what that block was sieved
   * with, on one thread and on three; the constants are pinned, as above.
   * Which block a correction first applies to depends on the width of a
   * bit array and on RP_ADAPT_LAG, and the stop lies in it as described
   * in the default build only, so what is printed is whether the two runs
   * report the same moduli, not the moduli. */
  { static const long rich[] = {16, 20, -131, -109, 152, 49, 4};
    long used[2][3], k;

    for(k = 0; k < 2; k++)
    { curve(&args, 6, rich, 50000); args.num_threads = 1 + 2*k;
      args.survivors_per_word = 0.003; args.sp2_extra = 11; args.sp2_u0 = 1.6e6;
      args.cost_table = 38.0; args.sp3_per_denom = 0.013;
      stop_b = 2700; keep_to = k;
      report(k ? "the same on three threads"
               : "a curve of many points, height 50000, stopped at b >= 2700",
             find_points_work(&args, keep_until, NULL));
      used[k][0] = args.sp1_used; used[k][1] = args.sp2_used;
      used[k][2] = args.sp3_used;
    }
    same_points();
    printf("  the moduli reported: %s\n",
           (used[0][0] == used[1][0] && used[0][1] == used[1][1]
              && used[0][2] == used[1][2])
             ? "the same on one and three threads" : "DIFFERENT");
  }
  find_points_clear(&args);
  curve(&args, 6, sextic, 3000); args.num_threads = 4;
  report("find_points on four threads (the pool made and freed inside)",
         find_points(&args, count, NULL));

  /* The tuning fields at the extremes of their range: each of type double
   * at NaN of either sign, at the infinities, at the largest double and at
   * the least positive one, each of type long at LONG_MAX and LONG_MIN, and
   * then several at once.  They decide how the search goes about it, never
   * what it finds, so every search has to deliver the points of the search
   * with the defaults, in the same order.  The leading coefficient of the
   * curve is not a square, so that the tests on the denominators are in
   * play.  Then the window at the top of the largest height bound that
   * test4.sh searches, where the height alone would let array_size be 2^58
   * bit arrays. */
  { static const long quad[] = {5, 3, 7};       /* y^2 = 7x^2 + 3x + 5 */
    static const double dval[] = {NAN, -NAN, INFINITY, -INFINITY, DBL_MAX,
                                  4.9406564584124654e-324};
    static const long lval[] = {LONG_MAX, LONG_MIN};
    long f, i;

    curve(&args, 2, quad, 1000);
    keep_to = 0; report("y^2 = 7x^2 + 3x + 5, height 1000, the defaults",
                        find_points(&args, keep, NULL));
    keep_to = 1;
    printf("  the fields of type double at NaN, -NaN, inf, -inf, the largest"
           " and the least\n  positive double, those of type long at"
           " LONG_MAX and LONG_MIN:\n");
    for(f = 0; f < NDFIELDS; f++)
    { long bad = -1;

      for(i = 0; i < 6; i++)
      { curve(&args, 2, quad, 1000); set_dfield(&args, f, dval[i]);
        find_points(&args, keep, NULL);
        if(!same_as_first() && bad < 0) { bad = i; }
      }
      if(bad < 0) { printf("    %-40s the same points\n", dfields[f]); }
      else { printf("    %-40s DIFFERENT at %g\n", dfields[f], dval[bad]); }
    }
    for(f = 0; f < NLFIELDS; f++)
    { long bad = -1;

      for(i = 0; i < 2; i++)
      { curve(&args, 2, quad, 1000); set_lfield(&args, f, lval[i]);
        find_points(&args, keep, NULL);
        if(!same_as_first() && bad < 0) { bad = i; }
      }
      if(bad < 0) { printf("    %-40s the same points\n", lfields[f]); }
      else { printf("    %-40s DIFFERENT at %ld\n", lfields[f], lval[bad]); }
    }
    /* the offset of the second stage at its largest, flat and scaled by a
     * run length so short that it stays near 2^63, and the number of
     * first-stage moduli at its largest, which the rule that looks at
     * further primes adds that offset to */
    curve(&args, 2, quad, 1000); args.sp2_extra = LONG_MAX; args.sp2_u0 = 0.0;
    find_points(&args, keep, NULL);
    printf("    %-40s %s\n", "sp2_extra LONG_MAX, sp2_u0 0",
           same_as_first() ? "the same points" : "DIFFERENT");
    curve(&args, 2, quad, 1000); args.sp2_extra = LONG_MAX; args.sp2_u0 = 1e-300;
    find_points(&args, keep, NULL);
    printf("    %-40s %s\n", "sp2_extra LONG_MAX, sp2_u0 1e-300",
           same_as_first() ? "the same points" : "DIFFERENT");
    curve(&args, 2, quad, 1000); args.sp1 = LONG_MAX; args.sp2_u0 = 0.0;
    find_points(&args, keep, NULL);
    printf("    %-40s %s\n", "sp1 LONG_MAX, sp2_u0 0",
           same_as_first() ? "the same points" : "DIFFERENT");
    /* nearly every field at once: the doubles at NaN, then at inf, and the
     * longs at LONG_MAX but for sp2 and num_primes, which are left to the
     * library, so that the second stage is chosen and the rule that looks
     * at further primes is in play */
    for(i = 0; i < 4; i += 2)
    { curve(&args, 2, quad, 1000); all_extreme(&args, dval[i]);
      find_points(&args, keep, NULL);
      printf("    %-40s %s\n", (i == 0) ? "all but sp2, num_primes; doubles NaN"
                                        : "all but sp2, num_primes; doubles inf",
             same_as_first() ? "the same points" : "DIFFERENT");
    }
    for(n = 0; n < nkept[0]; n++) { mpz_clear(kept[0][n].y); }

    top_window(&args);
    keep_to = 0; nkept[0] = 0;
    report("y^2 = x - 9223372036854769900, height 2^63 - 1",
           find_points(&args, keep, NULL));
    keep_to = 1;
    top_window(&args); args.array_size = LONG_MAX;
    find_points(&args, keep, NULL);
    printf("    %-40s %s\n", "array_size LONG_MAX",
           same_as_first() ? "the same points" : "DIFFERENT");
    top_window(&args); all_extreme(&args, INFINITY);
    find_points(&args, keep, NULL);
    printf("    %-40s %s\n", "all but sp2, num_primes; doubles inf",
           same_as_first() ? "the same points" : "DIFFERENT");
    for(n = 0; n < nkept[0]; n++) { mpz_clear(kept[0][n].y); }
    nkept[0] = 0; keep_to = 0;
  }

  for(n = 0; n <= RATPOINTS_MAX_DEGREE; n++) { mpz_clear(c[n]); }
  return(0);
}
