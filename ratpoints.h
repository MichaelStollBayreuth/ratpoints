/***********************************************************************
 * ratpoints-3.1.0                                                     *
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
 * ratpoints.h                                                         *
 *                                                                     *
 * Header file for the ratpoints program and library                   *
 *                                                                     *
 * Michael Stoll, Sep 21, 2009; Jan 7, 2022; Sep 21, 2026              *
 * with changes by Bill Allombert, Dec 29, 2021                        *
 ***********************************************************************/

/* Use the GNU multiprecision library. */
#include <gmp.h>

#define RATPOINTS_MAX_DEGREE 100           /* max. degree of f(x) */

/* A fixed choice of the numbers of moduli in the first stage and in the
 * first two stages, as versions 2.x used it.  The library does not use
 * these values: when sp1 and sp2 are negative -- which is what the command
 * line program passes unless -n / -N are given -- it chooses both from the
 * curve, as described below.  They remain for source compatibility, and as
 * a reasonable fixed choice for random genus 2 curves at a height bound of
 * about 10^4; on curves with many rational points they cost about 50%. */
#define RATPOINTS_DEFAULT_SP1 11           /* a fixed value for sp1 */
#define RATPOINTS_DEFAULT_SP2 19           /* a fixed value for sp2 */

/* How the numbers of sieving moduli are chosen when sp1 and sp2 are
 * negative.
 *
 * The first stage of the sieve costs one load and one AND per bit array and
 * modulus, whatever survives; the second stage costs nothing on an empty
 * bit array and a good deal on a non-empty one.  So moduli are moved into
 * the first stage until few enough bit arrays are left alive.  The
 * candidates are ranked by cost per information (prime_cost and phase_1_key
 * in find_points.c), and the first stage takes a modulus of density r while
 *
 *   (survivors per word reaching it) * (1 - r) * (downstream cost factor)
 *       > RATPOINTS_SURVIVORS_PER_WORD * (its cost per word) .
 *
 * Its cost per word is one AND plus its fixed costs -- tables, entries per
 * denominator, overhead per call, row fetch -- spread over the run.  The
 * downstream factor is what a survivor costs in the later stages relative
 * to a long run: there the second stage kills it cheaply, while a run of a
 * few thousand words has no second stage and every survivor reaches the
 * extraction, which costs several times more.  Pricing the fixed costs is
 * what keeps a short run from taking moduli whose tables it cannot pay for:
 * at a height bound of 200 the stage stops after seven moduli instead of
 * thirteen whose tables would be a fifth of the run.
 *
 * The threshold is per 64-bit word rather than per bit array because the
 * best per-word value is the same at every register width, while the best
 * per-array value doubles with the width.  It is a property of the machine
 * -- where one more first-stage modulus stops paying for itself against the
 * cost of a non-empty bit array entering the second stage -- and nearly
 * independent of the curve: one value serves curves whose densities span a
 * factor of 600, and height bounds from 200 to 200000.  The optimum is
 * flat: a value within a factor of two of the best costs less than 3%.
 *
 * The second stage then gets RATPOINTS_SP2_EXTRA further moduli, fewer in a
 * short run (RATPOINTS_SP2_U0 below).
 *
 * Both constants are compiled-in defaults: they can be set per call through
 * the survivors_per_word and sp2_extra fields of ratpoints_args (a negative
 * value there means "use the compiled-in one") and on the command line with
 * -r and -R.  "make tune" measures them for the machine at hand, together
 * with the run length, the table cost and the third stage's cost per
 * denominator below, and writes what it finds to tuning.mk; it takes
 * several minutes and wants an idle machine.  By hand no rebuild is needed:
 *   ./rptest -r <x> -R <n> -U <u> -C <c> -z        (random curves)
 *   ./rptest-many -r <x> -R <n> -U <u> -C <c> -z   (curves with many points)
 * minimising the sum of the two times.  Use both tests: they cover the two
 * regimes that matter, and a value that suits one can be poor for the
 * other.  The offset matters much less than the threshold. */
#ifndef RATPOINTS_SURVIVORS_PER_WORD
# define RATPOINTS_SURVIVORS_PER_WORD 0.003 /* when to stop the first stage */
#endif
#ifndef RATPOINTS_SP2_EXTRA
# define RATPOINTS_SP2_EXTRA 11             /* sp2 = sp1 + this, capped */
#endif

/* The offset above is for a run long enough that the fixed costs of a
 * modulus no longer matter.  A modulus is set up once and then used for the
 * whole run: its sieving table is built at most m times however many
 * numerators there are, and its entry in the list of reduced denominators
 * (bp_list) is computed once per denominator.  Per word of numerators
 * sieved these costs fall like 1/U, where U is the number of words the run
 * will sweep, so the number of further moduli worth using is
 *   sp2_extra(U) = RATPOINTS_SP2_EXTRA / (1 + RATPOINTS_SP2_U0/U) .
 * RATPOINTS_SP2_U0 is the run length at which setting a modulus up costs as
 * much as sieving with it; it is a property of the machine, like the other
 * two constants.
 *
 * This is what makes one tuning serve every height bound.  The measured
 * best offset is 3 on random curves at a height bound of 16383, 5 on
 * point-rich ones there, 9 on random curves at 200000 and 12 on point-rich
 * ones, and the pair of constants gives 3, 8, 11 and 11.  A single fixed
 * value costs up to 12% at one of the two bounds.  Only tunings at two very
 * different run lengths ("make tune" and "make tunehigh") can separate the
 * two constants.
 *
 * Zero (also as the sp2_u0 field, or -U 0) makes the offset flat. */
#ifndef RATPOINTS_SP2_U0
# define RATPOINTS_SP2_U0 1.6e6
#endif

/* The third sieving stage tests one surviving numerator at a time against
 * further primes, computing (a * b^-1) mod p instead of reading a table.
 * It needs no table, so it can use primes for which building one would
 * never pay; what it costs is a multiplication and a reduction per survivor
 * and prime, and a set-up per prime for each denominator that brings a
 * survivor to the stage.  How many primes it should use therefore depends
 * on how many survivors a denominator has.  The two constants below are
 * these costs as fractions of what one exact check costs.  A prime of
 * density r is worth adding while
 *   S * (1 - RATPOINTS_SP3_PER_SURVIVOR - r) > RATPOINTS_SP3_PER_DENOM ,
 * where S is the expected number of survivors per denominator still in
 * play.  Since S falls by the factor r with every prime added, this stops
 * of its own accord, and with survivors thin on the ground it stops at
 * once, which is right, because then the cost per denominator is all there
 * is.  When no prime in hand pays but a prime of the mean density seen so
 * far on the curve would, the stage looks at a further prime beyond
 * num_primes, so the table buffer may grow past what the first two stages
 * needed.
 *
 * The second constant decides whether the stage runs at all, so it is the
 * one to tune: the sp3_per_denom field of ratpoints_args, -Q on the command
 * line, and the fifth constant "make tune" measures.  The running time is
 * flat in it from 0.003 to 0.05.  The number of primes can also be fixed
 * outright, through the sp3_extra field or with -P; then the stage takes
 * that many without the test above, and looks further only when the primes
 * in hand run out. */
#ifndef RATPOINTS_SP3_PER_SURVIVOR
# define RATPOINTS_SP3_PER_SURVIVOR 0.055  /* one prime tested, per survivor */
#endif
#ifndef RATPOINTS_SP3_PER_DENOM
# define RATPOINTS_SP3_PER_DENOM 0.013     /* one prime, per denominator */
#endif
/* The fraction of the survivors of the first two stages that the test for
 * common factors lets through, which enters the estimate of S.  It is about
 * 6/pi^2 for numerators that survive by chance; the denominators divisible
 * by several sieving primes keep the most survivors and have the fewest
 * coprime numerators, which lowers it, and the rows of the sieve that
 * already exclude numerators sharing a prime with the denominator raise it.
 * The value is fitted rather than derived: the separate effects are not
 * worth estimating. */
#ifndef RATPOINTS_SP3_COPRIME
# define RATPOINTS_SP3_COPRIME 0.7
#endif

/* What the operations of the sieve cost relative to each other.  The cost
 * constants from here to RATPOINTS_COST_CHECK (RATPOINTS_COMPOSITE_MAX
 * between them is a bound, not a cost) are per numerator word and in units
 * of what one first-stage modulus costs there, which is one AND per word --
 * a quarter of one AND on a 256-bit array, about 0.26 core cycles or 0.2
 * rdtsc cycles on the machine they were measured on (an Intel i7-1355U).
 * They are properties of the machine, measured rather than fitted: build
 * with -DRP_PHASE_TIMING and divide the cycles of each part by the number
 * of times it ran and by that unit.  At another register width one AND
 * covers another number of words and the unit changes with it (a 64-bit
 * build has these constants 1.6 times too small, a 128-bit build 1.2
 * times); that is within the flat part of every optimum, and "make tune"
 * covers the constants that matter.
 *
 * They are needed because moduli are not equally expensive.  The sieving
 * table of a modulus m has m rows and is built once for each class of
 * denominators modulo m that turns up, so at most m times: over a run of U
 * numerator words and D denominators that is COST_TABLE*m*min(D,m)/U per
 * word, which grows with the square of the modulus and falls as the run
 * gets longer.  At equal information a smaller modulus is therefore
 * better, and at a small height bound much better.  Setting COST_TABLE to
 * zero (the cost_table field, or -C 0) drops the term and ranks by
 * information alone.  Measured at 18 to 30; the running time is flat from
 * 20 to 70.  "make tune" sweeps it as its fourth constant. */
#ifndef RATPOINTS_COST_TABLE
# define RATPOINTS_COST_TABLE 38.0   /* building one row of a sieve table */
#endif
/* The largest composite modulus offered to the ranking.  A composite
 * modulus -- a prime power, or a product of primes and prime powers --
 * carries the information of all its factors for one AND per word; but its
 * row is m bit arrays and its table m rows of them, and the cost model
 * charges an AND the same whatever the modulus.  That holds only while the
 * rows of the first-stage moduli stay in the first-level cache: with rows
 * of some 200 bit arrays (7 KB) the cycles per AND double.  Products of two
 * mid-sized primes then save up to a fifth of the instructions and 3.7% of
 * the time on random curves at a height bound of 200000, but lose 4.5 to 7%
 * on point-rich curves, and telling the two cases apart would take a model
 * of the cache with several more constants.  Of the bounds 32, 64, 128 and
 * 255, 64 is the best on three of the four test suites and within 1% on
 * the fourth.  It admits the prime powers 9, 25, 27, 49 and the products
 * 15, 21, 33, 35, 39, 45, 51, 55, 57 and 63.  The cost of fetching the row
 * of a large modulus is charged by RATPOINTS_COST_LINE below. */
#ifndef RATPOINTS_COMPOSITE_MAX
# define RATPOINTS_COMPOSITE_MAX 64
#endif
#ifndef RATPOINTS_COST_BP
# define RATPOINTS_COST_BP 8.0       /* one entry of bp_list, per denominator:
     the denominator reduced modulo the modulus by a multiplication.
     Measured at 5.5 to 11 on the four test suites; the running time does
     not react to the value (8 against 24 within 0.3%). */
#endif
/* and, for a modulus that has a sieving table, filling in its sieve_spec
 * once for every denominator.  Measured at 17 to 32 per modulus on the four
 * test suites, three times the bp_list entry; between them the two are
 * about 7% of "make test1" -- a fixed cost per denominator and modulus,
 * which is what makes a short run want fewer moduli than a long one. */
#ifndef RATPOINTS_COST_SETUP
# define RATPOINTS_COST_SETUP 30.0
#endif
/* and, for a modulus of the first stage, what one call of the sieve costs
 * it beyond its ANDs: the reduction that sets its row pointer at the head
 * of the call and the narrower legs that sieve the bit arrays past the last
 * whole chunk (sift.c).  A call handles at most array_size bit arrays of
 * one interval of one denominator, so below a height bound of some 30000 it
 * is one call per interval and denominator, and the cost is of COST_SETUP's
 * kind.  Measured from the core cycles per first-stage AND on a 256-bit
 * array against the bit arrays per call: 1.03 at 200 of them, 1.40 at 56,
 * 3.0 to 3.4 at 10, 5 to 40 at 1 to 3.  The excess is 20 to 22 cycles per
 * call and modulus, of which the row fetch (COST_LINE below) accounts for 7
 * to 8 at the sizes of the primes the random curves use, leaving 12 to 13
 * core cycles; in the units of this block that is 50.  At a height bound of
 * 16383 (512 words per call) it is a tenth of an AND per word, at 200000
 * (three calls per denominator of 2000 words) 7%. */
#ifndef RATPOINTS_COST_CALL
# define RATPOINTS_COST_CALL 50.0
#endif
/* and the fetch of a first-stage modulus's table row for each denominator.
 * The row is m bit arrays, a denominator walks min(m, A) of them (A its
 * number of bit arrays), and they come from beyond the first-level cache,
 * the previous denominator's row having been another one; per word that is
 * COST_LINE*min(m,A)*(bytes per bit array / 64)*D/U.  Measured at 1.2 core
 * cycles per 64-byte cache line at a height bound of 16383, on random and
 * on point-rich curves alike (the latter use primes of 60 to 127 with 128
 * bit arrays per denominator and run at 1.85 cycles per AND against 1.40);
 * at 200000 the tables sit in the third-level cache and a line costs 2 to
 * 4 cycles, but there the term is a few per cent of an AND whatever the
 * value.  In the units of this block 1.2 cycles per line is 4.5.  At 16383
 * (A = 128 at full width) a modulus above 128 costs half again as much as
 * its ANDs (4.5/8 per word); at 200000 (A = 1500) every modulus costs
 * within 5% of them. */
#ifndef RATPOINTS_COST_LINE
# define RATPOINTS_COST_LINE 4.5
#endif
#ifndef RATPOINTS_COST_PHASE2
# define RATPOINTS_COST_PHASE2 110.0 /* one AND on a surviving bit-array */
#endif
/* what one survivor of the second stage costs from there on: the
 * extraction, the test for common factors, the third stage and, for the few
 * that get that far, the exact check.  Measured between 340 and 1900, the
 * small figures at a height bound of 16383 and the large ones at 200000,
 * where the numbers are bigger and the memory colder; the value here is for
 * the large bound, where the choice of the moduli reacts to it most. */
#ifndef RATPOINTS_COST_SURVIVOR
# define RATPOINTS_COST_SURVIVOR 1400.0
#endif
/* and what the exact check contributes to that, for the reference curve of
 * the check-cost estimate below.  It is the one term of COST_SURVIVOR that
 * the degree moves, and it is a small one: three per cent of a survivor at
 * a height bound of 200000 and eighteen per cent at 16383, because only one
 * survivor in twenty-five ever reaches the check.  In the same units as the
 * line above, so 306 rdtsc cycles divided by what one first-stage AND per
 * word costs, which is between 0.16 and 0.21 at the large bound. */
#ifndef RATPOINTS_COST_CHECK
# define RATPOINTS_COST_CHECK 1600.0
#endif

/* What one exact check costs, in rdtsc cycles, for the curve in hand.  The
 * two third-stage constants above are fractions of it, and what they are
 * fractions of is not the same for every curve: the check evaluates the
 * binary form F(a,b) at one numerator by Horner from the coefficients
 * bc[k] = c[k]*b^(degree-k), which are computed once per denominator, and
 * then takes an integer square root.  So it costs
 *   - one multiplication by a single word and one addition per Horner step,
 *     that is per degree, on numbers whose size runs from one limb up to the
 *     size of F and so averages half of it;
 *   - one more multiplication when the degree is odd, to make the form even;
 *   - a square root, which is nearly free while the root fits in one limb
 *     and costs about a fixed amount per further limb of the root after that.
 * F(a,b) has about  cbits + degree*hbits  bits, where cbits is the size of
 * the largest coefficient and hbits that of the height bound, so the degree
 * and the coefficients between them fix every quantity in the formula, and
 * all of it is known before the first prime is looked at.
 *
 * The four constants were measured by timing exactly the gmp calls the check
 * makes (bench_check.c), over degrees 2 to 20, coefficients of 4 to 250 bits
 * and height bounds of 2^14 and 2^18, with CHECK_CALL then set so that the
 * ratios agree with what the sieve itself shows (there a check costs 44
 * cycles plus 1.44 times the benchmark's figure); the estimate is within 6%
 * of the sieve's own figures over degrees 3 to 12.  Only the ratio of one
 * curve's check to another's is used, so a machine on which every check is
 * dearer needs no change here.  The estimate can be overridden per call
 * through the check_cost field of ratpoints_args, and with -W on the command
 * line; -W 306 assumes for every curve what a check costs on a degree-6
 * curve with small coefficients. */
#ifndef RATPOINTS_CHECK_STEP
# define RATPOINTS_CHECK_STEP 29.0   /* one Horner step */
#endif
#ifndef RATPOINTS_CHECK_LIMB
# define RATPOINTS_CHECK_LIMB 8.0    /* and per limb it carries */
#endif
#ifndef RATPOINTS_CHECK_CALL
# define RATPOINTS_CHECK_CALL 94.0   /* the call, and a one-limb square root */
#endif
#ifndef RATPOINTS_CHECK_ROOT
# define RATPOINTS_CHECK_ROOT 170.0  /* each further limb of the root */
#endif
/* What the formula gives for the curves the third-stage constants were
 * tuned on: degree 6, coefficients of a few bits, height bound between 2^14
 * and 2^18, where it ranges from 303 to 311.  The two fractions are divided
 * by the estimate over this, so they are unchanged on such a curve and
 * smaller on one whose check costs more. */
#ifndef RATPOINTS_CHECK_REFERENCE
# define RATPOINTS_CHECK_REFERENCE 306.0
#endif

#ifndef RATPOINTS_DEFAULT_NUM_PRIMES
#define RATPOINTS_DEFAULT_NUM_PRIMES 30    /* Default value for num_primes.
     Unless num_primes is set explicitly, this is where the search starts:
     sieving_info() looks at further primes when a curve does not leave
     enough of them informative for sp2 to be sp1 + RATPOINTS_SP2_EXTRA, and
     again when the third stage finds none in hand worth taking while a
     prime of the density the curve has been offering would be. */
#endif
#ifndef RATPOINTS_DEFAULT_STURM
#define RATPOINTS_DEFAULT_STURM 10         /* Default value for sturm_iter */
#endif

#ifndef RATPOINTS_DEFAULT_MAX_FORBIDDEN
#define RATPOINTS_DEFAULT_MAX_FORBIDDEN 64 /* Default value for max_forbidden:
     how many primes the denominators are tested against with bit arrays
     or by their valuation (option -F).  Since the denominators are taken a
     word of 64 at a time, one such prime costs a few instructions per word,
     so the cap is a generous one: the primes up to the square root of the
     height bound reach it at a height bound of about half a million. */
#endif

#ifndef RATPOINTS_ARRAY_SIZE
#define RATPOINTS_ARRAY_SIZE 256           /* Array size in bit-arrays */
#endif

/* data structure for intervals, used in finding the positivity region */
typedef struct {double low; double up;} ratpoints_interval;

/* main data structure for arguments and local data */
typedef struct { mpz_t *cof; long degree; long height;
                 ratpoints_interval *domain; long num_inter;
                 long b_low; long b_high; long sp1; long sp2; long sp3;
                 double survivors_per_word; long sp2_extra; double sp2_u0;
                 double cost_table; long adapt;
                 long sp3_extra; double sp3_per_denom; double check_cost;
                 long array_size;
                 long sturm; long num_primes; long max_forbidden;
                 long num_threads;
                   /* how many threads sieve: 1 (or less), the calling one;
                      n > 1: n threads besides the calling one, which hands
                      out the work and delivers the points, in the order a
                      single thread would find them; negative: as many as
                      there are online processors.  A library built without
                      threads (RATPOINTS_NO_THREADS) ignores it. */
                 unsigned int flags;
                 long sp1_used; long sp2_used; long sp3_used;
                   /* output: the number of sieving moduli (primes and
                      composite moduli) the last search used in the first
                      stage and in the first two, and the number of moduli
                      and primes in all three stages; the input fields
                      above come back as they went in */
        /* from here: private data */
                 mpz_t *work; long work_length;
                 void *se_buffer; void *se_next;
                 void *ba_buffer; void *ba_next;
                 int *int_buffer; int *int_next;
                 void *sieve_list; void *stage3_list; void *magics;
                 void *den_info; void *divisors;
                 void *forb_ba; void *forbidden;
                 void *forb_words; long forb_words_len;
                 void *ba_buffer_na; long ba_buffer_primes;
                 long ba_buffer_arrays; void *pw_buffer;
                 double run_words; double run_denoms;
                 unsigned long adapt_at; long sp3_max;
                 double check_rel;
                 void *pool;
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
#define RATPOINTS_USE_JACOBI      (unsigned int)0x1000
  /* the Jacobi symbol test on the denominators applies: even degree, the
     leading coefficient is not a square, and RATPOINTS_NO_JACOBI is not set */
#define RATPOINTS_NO_REVERSE_AUTO (unsigned int)0x4000
  /* the program itself decided not to reverse (intervals given, or bounds
     on the denominator); kept apart from RATPOINTS_NO_REVERSE so that the
     caller's input is not changed */

/* Return values of find_points() */
#define RATPOINTS_NON_SQUAREFREE (-1)
#define RATPOINTS_BAD_ARGS (-2)
#define RATPOINTS_WORK_LENGTH_TOO_SMALL (-3)
#define RATPOINTS_NO_MEMORY (-4)  /* a sieving thread ran out of memory for
                                     the points of a block */

/* Function prototypes */
long find_points(ratpoints_args*,
                 int proc(long, long, const mpz_t, void*, int*), void*);

void find_points_init(ratpoints_args*);

long find_points_work(ratpoints_args*,
                      int proc(long, long, const mpz_t, void*, int*), void*);

void find_points_clear(ratpoints_args*);
