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
 * repository.  (Those figures were 0.0075 for the rule as it then was; the
 * tuning session of 2.3 made the rule count the survivors a modulus
 * removes rather than the ones it meets, a factor 1 - r of about a half,
 * and charge the modulus its per-call cost and its row fetch, so the same
 * measurements put the constant at 0.003 -- a sweep of the threshold on
 * the corrected model at 200 to 200000 has a flat basin from 0.002 to
 * 0.004 at 16383, and "make tune" confirmed it.)
 *
 * The constant is a property of the machine, not of the curve: it is where
 * one more first-phase prime stops paying for itself against the cost of a
 * non-empty bit-array entering the second phase.  Anything within a factor
 * of about two of the value below costs less than 3%.
 *
 * To retune them for another machine, run "make tune", which measures these
 * two, the run length below and the table cost further down, over the two
 * test sets and writes what it finds to tuning.mk; it takes several minutes
 * and wants an idle machine.  By hand, no rebuild is needed either, since
 * all four can be set on the command line:
 *   ./rptest -r <x> -R <n> -U <u> -C <c> -z        (random curves)
 *   ./rptest-many -r <x> -R <n> -U <u> -C <c> -z   (curves with many points)
 * minimising the sum of the two times.  Use both tests: they cover the two
 * regimes that matter, and a value that suits one can be poor for the other.
 * The offset for sp2 matters much less than the threshold, and at large
 * height bounds hardly at all. */
/* Both constants are compiled-in defaults only: they can be set per call
 * through the survivors_per_word and sp2_extra fields of ratpoints_args
 * (a negative value there means "use the compiled-in one"), and on the
 * command line with -r and -R. */
#ifndef RATPOINTS_SURVIVORS_PER_WORD
# define RATPOINTS_SURVIVORS_PER_WORD 0.003 /* when to stop the first phase */
#endif
#ifndef RATPOINTS_SP2_EXTRA
# define RATPOINTS_SP2_EXTRA 11             /* sp2 = sp1 + this, capped */
#endif

/* ...but only for a run long enough that the fixed costs of a prime no
 * longer matter.  A prime is set up once and then used for the whole run:
 * its sieve table is built at most p times however many numerators there
 * are, and its entry in bp_list is computed once per denominator.  Per word
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
 * This is what makes one tuning serve every height bound.  The measured best
 * offset is 3 on random curves at a height bound of 16383, 5 on point-rich
 * ones there, 9 on random curves at 200000 and 12 on point-rich ones, and
 * the pair above lands on 3, 8, 11 and 11.  A single fixed value costs 12%
 * of the point-rich half at the larger bound.
 *
 * Setting RATPOINTS_SP2_U0 (or the sp2_u0 field, or -U) to zero switches the
 * correction off and restores a flat offset.
 *
 * The first phase has a correction of its own since 2.3 (TODO item 29).
 * A modulus is taken while the survivors per word it removes -- the
 * survivors it meets times 1 - r, r its density -- exceed the threshold
 * times what the modulus costs per word -- one AND, plus its tables, its
 * per-denominator and per-call entries and the fetch of its row spread
 * over the run (prime_cost and phase_1_key in find_points.c) -- with the
 * survivors weighted by what one costs
 * downstream relative to a long run: there the second phase's eleven
 * moduli kill it cheaply, while in a run of a few thousand words there is
 * no second phase and every survivor reaches the extraction, which costs
 * several times more (RATPOINTS_COST_PHASE2 and RATPOINTS_COST_SURVIVOR
 * below say how much).  At 200000 both factors are one.  At 16383, where
 * "make tune" fits the threshold, both are about three at the modulus
 * where the phase stops and nearly cancel; the phase takes 0.4 moduli
 * more than before there, which the next fit will weigh.  At a height
 * bound of 200 the cost factor is a hundred for the small primes and
 * several hundred at the stopping modulus, the downstream factor six, the
 * phase stops after seven moduli where it took thirteen whose tables were
 * a fifth of the run, and the run is a quarter shorter.  (The factor
 * 1 - r, the per-call cost and the row fetch came with the tuning
 * session's estimate corrections, refitted as one group with the other
 * corrections of that session; see the manual.) */
#ifndef RATPOINTS_SP2_U0
# define RATPOINTS_SP2_U0 1.6e6
#endif

/* The third sieving stage tests one surviving numerator at a time against
 * further primes, computing (a * b^-1) mod p instead of reading a table.
 * It needs no table, so it can use primes that would never be worth
 * building one for; what it costs is a multiplication and a reduction per
 * survivor and prime, and a set-up for each denominator that brings a
 * survivor to the stage at all (up to 2.3's item 24 it was set up for every
 * denominator, and b carried along modulo each of its primes, which is what
 * the second constant below was tuned for).  How many primes it should use
 * therefore depends on how many survivors a denominator has, which is what
 * the two constants below express, both as fractions of what one exact
 * check costs: the first is the per-survivor cost of testing one prime, the
 * second the per-denominator cost of carrying it.  A prime is worth adding
 * while
 *   S * (1 - q_s - r) > q_d ,
 * where S is the expected number of survivors per denominator still in play
 * and r is the prime's density.  The second of the two decides whether the
 * stage runs at all, so it is the one to tune; it can be set per call through
 * the sp3_per_denom field of ratpoints_args, and with -Q on the command line.
 * Since S falls by a factor of r with every prime added, this stops of its
 * own accord; with survivors thin on the ground it stops at once, which is
 * what should happen, because then the per-denominator cost is all there is.
 * When no prime in hand pays but a prime of the mean density seen so far on
 * the curve would, the stage looks at a further prime beyond num_primes
 * (item 27), so the table buffer may grow past what the first two stages
 * needed.  The number of primes can also be fixed outright, through the
 * sp3_extra field of ratpoints_args or with -P on the command line; then
 * the stage takes that many without the test above, and looks further only
 * when the primes in hand run out. */
#ifndef RATPOINTS_SP3_PER_SURVIVOR
# define RATPOINTS_SP3_PER_SURVIVOR 0.055  /* one prime tested, per survivor */
#endif
#ifndef RATPOINTS_SP3_PER_DENOM
# define RATPOINTS_SP3_PER_DENOM 0.013     /* one prime, per denominator */
#endif
/* What fraction of the survivors of the first two phases the test for common
 * factors lets through.  It is about 6/pi^2 for the numerators that survive
 * by chance, a little less because the denominators that keep the most
 * survivors are those divisible by several of the sieving primes, which are
 * also the ones with the fewest coprime numerators.  It was fitted when the
 * even denominators alone were packed, at half width, and carries what that
 * left over; since 2.3 every class of the denominator is packed with its
 * own stride (rp_num_class in rp-private.h) and run_shape counts that, so
 * this is a fitted fudge; neither factor is worth estimating separately. */
#ifndef RATPOINTS_SP3_COPRIME
# define RATPOINTS_SP3_COPRIME 0.7
#endif

/* What the sieve's operations cost, relative to each other.  The cost
 * constants from here to RATPOINTS_COST_CHECK (RATPOINTS_COMPOSITE_MAX
 * between them is a bound, not a cost) are per numerator word and in units
 * of what one first-phase prime costs there, which is one AND per word --
 * a quarter of one AND on a 256-bit array, about 0.26 core cycles or 0.2
 * rdtsc cycles on the machine they were measured on; they are properties of
 * the machine, measured rather than fitted, by building with
 * -DRP_PHASE_TIMING and dividing the cycles of each part by the number of
 * times it ran and by that unit.  (At another register width one AND
 * covers another number of words and the unit changes with it, so these
 * constants want measuring again there; see tune.sh on the table cost.)
 *
 * They are here because the primes are not equally expensive and the rule
 * that picks them used to act as though they were.  The sieve table for p
 * has p rows and is built once for each denominator class that turns up, so
 * at most p times: over a run of U numerator words that is
 * COST_TABLE*p*min(D,p)/U per word, which grows with the square of the prime
 * and falls as the run gets longer.  At equal information a smaller prime is
 * therefore strictly better, and at a small height bound it is much better.
 * Setting COST_TABLE to zero (the cost_table field, or -C 0) drops the term
 * and restores ranking by information alone.  "make tune" sweeps it as its
 * fourth stage, after the three constants above. */
#ifndef RATPOINTS_COST_TABLE
# define RATPOINTS_COST_TABLE 38.0   /* building one row of a sieve table */
#endif
/* The largest composite modulus offered to the ranking (2.3, TODO item
 * 21).  A composite modulus -- a prime power, or a product of primes and
 * prime powers -- carries the information of all its factors for one AND
 * per word; but its row is m bit arrays and its table m rows of them, and
 * the cost model charges an AND the same whatever the modulus.  That holds
 * while the rows are small: measured here, the first phase's cycles per
 * AND double when rows of some 200 bit arrays (7 KB) take the place of the
 * small primes' rows, which live in the first-level cache and whose tables
 * fit the second, and the products of two mid-size primes then save a
 * fifth of the instructions and a tenth of the time on random curves at
 * height 200000 while losing 7% on the point-rich ones.  With the moduli
 * bounded by 64 -- the powers 9, 25, 27, 49 and the products of the primes
 * up to 21 -- the cycles follow the instructions on the random curves (32,
 * 128 and 255 were measured too; 64 wins on test1, testhigh and
 * testhighmany and loses a point to 32 on test1many).  Item 30 measured
 * what a per-AND cost rising with the modulus would buy: the count rule
 * is at the fixed-count optimum already, and the larger products pay only
 * on random curves at 200000 (3.7% of testhigh) while costing 4.5% of the
 * point-rich suite, which a set-aware model with several cache constants
 * would be needed to tell apart; so the bound stays, and the row fetch a
 * modulus above the bit arrays of a denominator pays for is charged by
 * RATPOINTS_COST_LINE below. */
#ifndef RATPOINTS_COMPOSITE_MAX
# define RATPOINTS_COMPOSITE_MAX 64
#endif
#ifndef RATPOINTS_COST_BP
# define RATPOINTS_COST_BP 8.0       /* one entry of bp_list, per denominator:
     the denominator reduced modulo the prime by a multiplication.  Measured
     at 5.5 to 11 on the four test suites in September 2026, when the entry
     stopped being stepped from the previous denominator (24 before, measured
     10 to 38); the runs do not react to the value within the noise. */
#endif
/* and, for a prime that has a sieving table, filling in its sieve_spec once
 * for every denominator.  Measured at 17 to 32 per prime on the four test
 * suites in September 2026, three times the bp_list entry; between them the
 * two are about 7% of "make test1" (13% before item 24) -- a fixed cost per
 * denominator and per prime, which is what makes a short run want fewer
 * primes than a long one. */
#ifndef RATPOINTS_COST_SETUP
# define RATPOINTS_COST_SETUP 30.0
#endif
/* and, for a modulus of the first phase, what one call of the sieve costs
 * it beyond its ANDs: the reduction that sets its row pointer at the head
 * of the call and the narrower legs that sieve the bit arrays past the last
 * whole chunk (sift.c).  A call handles at most array_size bit arrays of
 * one denominator's interval, so below a height bound of some 30000 it is
 * one call per interval and denominator, and the cost is of COST_SETUP's
 * kind.  Measured in September 2026 (item 30) from the counter build's
 * core cycles per first-phase AND on a 256-bit array against the bit
 * arrays per call: 1.03 at 200 of them, 1.40 at 56, 3.0-3.4 at 10, 5 to 40
 * at 1 to 3.  The excess over 1.03 times the arrays per call is 20-22
 * cycles per call and modulus, of which the row fetch (COST_LINE below)
 * accounts for 7-8 at the sizes of the primes the random curves use,
 * leaving some 12-13 core cycles; in the units of this block that is 50.
 * At 16383 (512 words per call) it is a tenth of an AND per word, at
 * 200000 (three calls per denominator of 2000 words) 7%. */
#ifndef RATPOINTS_COST_CALL
# define RATPOINTS_COST_CALL 50.0
#endif
/* and the fetch of a first-phase modulus's table row for each denominator.
 * The row is p bit arrays, a denominator walks min(p, A) of them (A its
 * bit arrays), and they come from beyond the first-level cache, the
 * previous denominator's row having been another one; per word that is
 * COST_LINE*min(p,A)*(bytes per bit array / 64)*D/U.  Measured in
 * September 2026 (item 30): the fit of the counter build's cycles per
 * first-phase AND at 16383 gives 1.2 core cycles per 64-byte line on the
 * random curves, and the point-rich curves there (primes of 60-127 with
 * 128 bit arrays per denominator) run at 1.85 cycles per AND against
 * 1.40, the same 1.2 per line; at 200000 the tables sit in the
 * third-level cache and the fit gives 2-4, but there the term is a few per
 * cent of an AND whatever the value.  In the units of this block 1.2
 * cycles per line is 4.5.  At 16383 (A = 128 at full width) a modulus
 * above 128 costs half again as much as its ANDs (4.5/8 per word), at
 * 200000 (A = 1500) every modulus costs within 5% of them: this is the
 * cost that item 21's bound on the composite moduli stood in for at the
 * smaller bound. */
#ifndef RATPOINTS_COST_LINE
# define RATPOINTS_COST_LINE 4.5
#endif
#ifndef RATPOINTS_COST_PHASE2
# define RATPOINTS_COST_PHASE2 110.0 /* one AND on a surviving bit-array */
#endif
/* what one survivor of the second phase costs from there on: the extraction,
 * the test for common factors, the third stage and, for the few that get
 * that far, the exact check.  Measured between 340 and 1900, the small
 * figures at a height bound of 16383 and the large ones at 200000, where the
 * numbers are bigger and the memory colder; the value here is for the large
 * bound, since that is where the correction of item 14 ever fires. */
#ifndef RATPOINTS_COST_SURVIVOR
# define RATPOINTS_COST_SURVIVOR 1400.0
#endif
/* and what the exact check contributes to that, for the reference curve of
 * check_cost() below.  It is the one term of COST_SURVIVOR that the degree
 * moves, and it is a small one: three per cent of a survivor at a height
 * bound of 200000 and eighteen per cent at 16383, because only one survivor
 * in twenty-five ever reaches the check.  In the same units as the line
 * above, so 306 rdtsc cycles divided by what one first-phase AND per word
 * costs, which is between 0.16 and 0.21 at the large bound. */
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
 * makes, over degrees 2 to 20, coefficients of 4 to 250 bits and height
 * bounds of 2^14 and 2^18, and then setting CHECK_CALL so that the ratios
 * agree with what the sieve itself shows; see PARAM-NOTES.md on the
 * adaptive-parameters branch of the git repository.  Only the ratio of one
 * curve's check to another's is used, so a machine on which every check is
 * dearer needs no change here.  The estimate can be overridden per call
 * through the check_cost field of ratpoints_args, and with -W on the command
 * line; -W 306 is what versions before 2.3 did, which is to assume that every
 * curve costs what a degree-6 one with small coefficients costs. */
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
     so the cap is a generous one; it was 30 up to 2.2.4, which the primes
     of the compiled table alone reach at a height bound of 200000. */
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
                 unsigned int flags;
                 long sp1_used; long sp2_used; long sp3_used;
                   /* output: the number of sieving moduli (primes, and
                      since 2.3 composite moduli) the last search used in
                      the first stage and in the first two, and the number
                      of moduli and primes in all three stages; the input
                      fields above come back as they went in */
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
                 unsigned long n_words; unsigned long n_arrays;
                 unsigned long n_bits; unsigned long n_coprime;
                 unsigned long n_checks; unsigned long n_sifts;
                 unsigned long n_words_2;
                 unsigned long adapt_at; long sp3_max; int stage3_filled;
                 double check_rel;
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
#define RATPOINTS_COMPUTE_BC      (unsigned int)0x2000
#define RATPOINTS_NO_REVERSE_AUTO (unsigned int)0x4000
  /* the program itself decided not to reverse (intervals given, or bounds
     on the denominator); kept apart from RATPOINTS_NO_REVERSE so that the
     caller's input is not changed */

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
