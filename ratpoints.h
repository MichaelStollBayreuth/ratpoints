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

/* The following seem to work well for Intel(R) processors with AVX2.
 * For just AVX without AVX2, 10/18/27 seems to be slightly faster.
 * This is for random genus 2 curves and height 10^5, as tested via
 * ./rptest -h 100000 > /dev/null .
 * In any case, the optimal values will depend on the CPU etc. */
#define RATPOINTS_DEFAULT_SP1 11           /* Default value for sp1 */
#define RATPOINTS_DEFAULT_SP2 19           /* Default value for sp2 */
#define RATPOINTS_DEFAULT_NUM_PRIMES 30    /* Default value for num_primes */
#define RATPOINTS_DEFAULT_STURM 10         /* Default value for sturm_iter */

#define RATPOINTS_DEFAULT_MAX_FORBIDDEN 30 /* Default value for max_forbidden */

#define RATPOINTS_ARRAY_SIZE 256           /* Array size in bit-arrays */

/* data structure for intervals, used in finding the positivity region */
typedef struct {double low; double up;} ratpoints_interval;

/* main data structure for arguments and local data */
typedef struct { mpz_t *cof; long degree; long height;
                 ratpoints_interval *domain; long num_inter;
                 long b_low; long b_high; long sp1; long sp2;
                 long array_size;
                 long sturm; long num_primes; long max_forbidden;
                 unsigned int flags;
        /* from here: private data */
                 mpz_t *work; long work_length;
                 void *se_buffer; void *se_next;
                 void *ba_buffer; void *ba_next;
                 int *int_buffer; int *int_next;
                 void *sieve_list;
                 void *den_info; void *divisors;
                 void *forb_ba; void *forbidden;
                 void *ba_buffer_na;
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
