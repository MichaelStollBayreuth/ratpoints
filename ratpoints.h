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
 * ratpoints.h                                                         *
 *                                                                     *
 * Header file for the ratpoints program and library                   *
 *                                                                     *
 * Michael Stoll, September 21, 2009; January 7, 2022                  *
 * with changes by Bill Allombert, December 29, 2021                   *
 ***********************************************************************/

/* Use the GNU multiprecision library. */
#include <gmp.h>

#define RATPOINTS_MAX_DEGREE 100           /* max. degree of f(x) */

/* These were the fixed defaults up to version 2.2.3.  They are kept for
 * source compatibility and are no longer used by the library, which now
 * chooses sp1 and sp2 from the curve (see below).  They were tuned on
 * random genus 2 curves for Intel(R) processors with AVX2 and are a
 * reasonable fixed choice for that regime -- but only for that regime: on
 * curves with many rational points they cost about 50%. */
#define RATPOINTS_DEFAULT_SP1 11           /* Former fixed value for sp1 */
#define RATPOINTS_DEFAULT_SP2 19           /* Former fixed value for sp2 */

/* Normally sp1 and sp2 are not taken from the two values above, but are
 * chosen from the curve: setting them to a negative value (which is what
 * the command line program does when -n / -N are absent) asks for that.
 *
 * The first phase of the sieve costs one load and one AND per bit-array
 * per prime, whatever survives; the second phase costs nothing on an empty
 * bit-array but a good deal on a non-empty one.  So primes should be moved
 * into the first phase until few enough bit-arrays are left alive, and
 * measurements over curves spanning a factor of 600 in density show that
 * the best sp1 is where the expected number of surviving numerators per
 * 64-bit word falls below the constant below -- almost independently of the
 * curve.  Per word rather than per bit-array: measured across the register
 * widths, the best per-word figure is the same (0.0075, 0.0075, 0.010 at
 * 128, 256 and 512 bits) while the per-bit-array one doubles with each
 * doubling of the width.  The measurements behind this are in
 * PARAMETER-MODEL.md on the phases-by-register-width branch of the git
 * repository.
 *
 * The constant is a property of the machine, not of the curve: it is where
 * one more first-phase prime stops paying for itself against the cost of a
 * non-empty bit-array entering the second phase.  Anything within a factor
 * of about two of the value below costs less than 3%.
 *
 * To retune them for another machine, run "make tune", which measures both
 * over the two test sets and writes what it finds to tuning.mk; it takes
 * several minutes and wants an idle machine.  By hand, no rebuild is needed
 * either, since the two can be set on the command line:
 *   ./rptest -r <x> -R <n> -z          (random curves)
 *   ./rptest-many -r <x> -R <n> -z     (curves with many rational points)
 * minimising the sum of the two times.  Use both tests: they cover the two
 * regimes that matter, and a value that suits one can be poor for the other.
 * The offset for sp2 matters much less than the threshold, and at large
 * height bounds hardly at all. */
/* Both constants are compiled-in defaults only: they can be set per call
 * through the survivors_per_word and sp2_extra fields of ratpoints_args
 * (a negative value there means "use the compiled-in one"), and on the
 * command line with -r and -R. */
#ifndef RATPOINTS_SURVIVORS_PER_WORD
# define RATPOINTS_SURVIVORS_PER_WORD 0.0075 /* when to stop the first phase */
#endif
#ifndef RATPOINTS_SP2_EXTRA
# define RATPOINTS_SP2_EXTRA 9              /* sp2 = sp1 + this, capped */
#endif

/* ...but only for a run long enough that the fixed costs of a prime no
 * longer matter.  A prime is set up once and then used for the whole run:
 * its sieve table is built at most p times however many numerators there
 * are, and its entry in bp_list is stepped once per denominator.  Per word
 * of numerators sieved those fall like 1/U, where U is the number of such
 * words the run will sweep, so the marginal condition for one more prime
 * has the shape
 *   worth adding  <=>  a + (fixed)/U  <  what it saves,
 * and the number of primes worth using is therefore
 *   sp2_extra(U) = RATPOINTS_SP2_EXTRA / (1 + RATPOINTS_SP2_U0/U) .
 * RATPOINTS_SP2_U0 is the number of numerator words at which setting a
 * prime up costs as much as sieving with it; it is a property of the
 * machine, like the other two constants here.
 *
 * This is what makes one tuning serve every height bound.  The measured
 * best offset is 3 to 5 at a height bound of 16383 and 9 at 200000, and a
 * single fixed value costs 8% of the pair at the larger bound (17.6% of its
 * point-rich half), which is why it is worth a constant of its own.  The
 * first phase needs no such correction: a phase-1 prime costs one AND per
 * word unconditionally, so the fixed part is a far smaller share of it, and
 * the best threshold indeed hardly moves with the height bound.
 *
 * Setting RATPOINTS_SP2_U0 (or the sp2_u0 field, or -U) to zero switches the
 * correction off and restores a flat offset. */
#ifndef RATPOINTS_SP2_U0
# define RATPOINTS_SP2_U0 1.2e6
#endif

/* The third sieving stage tests one surviving numerator at a time against
 * further primes, computing (a * b^-1) mod p instead of reading a table.
 * It needs no table, so it can use primes that would never be worth
 * building one for; what it costs is a multiplication and a reduction per
 * survivor and prime, and one subtraction per denominator and prime to
 * carry b along.  How many primes it should use therefore depends on how
 * many survivors a denominator has, which is what the two constants below
 * express, both as fractions of what one exact check costs: the first is
 * the per-survivor cost of testing one prime, the second the per-denominator
 * cost of carrying it.  A prime is worth adding while
 *   S * (1 - q_s - r) > q_d ,
 * where S is the expected number of survivors per denominator still in play
 * and r is the prime's density.  The second of the two decides whether the
 * stage runs at all, so it is the one to tune; it can be set per call through
 * the sp3_per_denom field of ratpoints_args, and with -Q on the command line.
 * Since S falls by a factor of r with every prime added, this stops of its
 * own accord; with survivors thin on the ground it stops at once, which is
 * what should happen, because then the per-denominator cost is all there is.
 * The number of primes can also be fixed outright, through the sp3_extra
 * field of ratpoints_args or with -P on the command line. */
#ifndef RATPOINTS_SP3_PER_SURVIVOR
# define RATPOINTS_SP3_PER_SURVIVOR 0.055  /* one prime tested, per survivor */
#endif
#ifndef RATPOINTS_SP3_PER_DENOM
# define RATPOINTS_SP3_PER_DENOM 0.013     /* one prime carried, per denom. */
#endif
/* What fraction of the survivors of the first two phases the test for common
 * factors lets through.  It is about 6/pi^2 for the numerators that survive
 * by chance, a little less because the denominators that keep the most
 * survivors are those divisible by several of the sieving primes, which are
 * also the ones with the fewest coprime numerators.  It also carries the
 * halving that an even denominator brings, since only odd numerators are
 * then considered; neither factor is worth estimating separately. */
#ifndef RATPOINTS_SP3_COPRIME
# define RATPOINTS_SP3_COPRIME 0.7
#endif

#define RATPOINTS_DEFAULT_NUM_PRIMES 30    /* Default value for num_primes.
     Unless num_primes is set explicitly, this is where the search starts:
     sieving_info() looks at further primes when a curve does not leave
     enough of them informative for sp2 to be sp1 + RATPOINTS_SP2_EXTRA. */
#define RATPOINTS_DEFAULT_STURM 10         /* Default value for sturm_iter */

#define RATPOINTS_DEFAULT_MAX_FORBIDDEN 30 /* Default value for max_forbidden */

#define RATPOINTS_ARRAY_SIZE 256           /* Array size in bit-arrays */

/* data structure for intervals, used in finding the positivity region */
typedef struct {double low; double up;} ratpoints_interval;

/* main data structure for arguments and local data */
typedef struct { mpz_t *cof; long degree; long height;
                 ratpoints_interval *domain; long num_inter;
                 long b_low; long b_high; long sp1; long sp2; long sp3;
                 double survivors_per_word; long sp2_extra; double sp2_u0;
                 long sp3_extra; double sp3_per_denom;
                 long array_size;
                 long sturm; long num_primes; long max_forbidden;
                 unsigned int flags;
        /* from here: private data */
                 mpz_t *work; long work_length;
                 void *se_buffer; void *se_next;
                 void *ba_buffer; void *ba_next;
                 int *int_buffer; int *int_next;
                 void *sieve_list; void *stage3_list; void *magics;
                 void *den_info; void *divisors;
                 void *forb_ba; void *forbidden;
                 void *ba_buffer_na; long ba_buffer_primes;
                 double run_words; double run_denoms;
               }
        ratpoints_args;

/* Define the flag bits for the flags component: */
#define RATPOINTS_NO_CHECK        (unsigned int)0x0001
#define RATPOINTS_NO_Y            (unsigned int)0x0002
#define RATPOINTS_NO_REVERSE      (unsigned int)0x0004
#define RATPOINTS_NO_JACOBI       (unsigned int)0x0008
#define RATPOINTS_VERBOSE         (unsigned int)0x0010

#define RATPOINTS_FLAGS_INPUT_MASK \
 (RATPOINTS_NO_CHECK | RATPOINTS_NO_Y | RATPOINTS_NO_REVERSE | \
  RATPOINTS_NO_JACOBI | RATPOINTS_VERBOSE)

/* Flags bits for internal purposes */
#define RATPOINTS_REVERSED        (unsigned int)0x0100
#define RATPOINTS_CHECK_DENOM     (unsigned int)0x0200
#define RATPOINTS_USE_SQUARES     (unsigned int)0x0400
#define RATPOINTS_USE_SQUARES1    (unsigned int)0x0800
#define RATPOINTS_COMPUTE_BC      (unsigned int)0x2000

/* Return values of find_points() */
#define RATPOINTS_NON_SQUAREFREE (-1)
#define RATPOINTS_BAD_ARGS (-2)
#define RATPOINTS_WORK_LENGTH_TOO_SMALL (-3)

/* Function prototypes */
long find_points(ratpoints_args*,
                 int proc(long, long, const mpz_t, void*, int*), void*);

void find_points_init(ratpoints_args*);

long find_points_work(ratpoints_args*,
                      int proc(long, long, const mpz_t, void*, int*), void*);

void find_points_clear(ratpoints_args*);
