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
 * sift.c                                                              *
 *                                                                     *
 * The sieving procedure for ratpoints                                 *
 *                                                                     *
 * Michael Stoll, Apr 14; 2009, January 7, 2022                        *
 * with changes by Bill Allombert, Dec 29, 2021                        *
 ***********************************************************************/

#include "rp-private.h"


/* ---------------------------------------------------------------------
 * Development instrumentation: split the work of _ratpoints_sift0 into
 * its three stages and record how much survives each of them.  Build
 * with -DRP_PHASE_TIMING for the timings, and additionally with
 * -DRP_PHASE_COUNTS for the survivor counts, which need a popcount pass
 * over the survivors and two more per surviving unit in phase 2.  The
 * two flags should not be combined when the timings are what is wanted.
 * A machine-readable line is written to stderr when the program exits.
 *
 * The three stages are
 *   1  sieving with the first sp1 primes,
 *   2  scanning the survivors and sieving them with the next sp2-sp1,
 *   3  the exact check with gmp, which runs inside stage 2 and is
 *      subtracted from it.
 *
 * The timer is rdtsc, which counts reference cycles at a constant rate
 * and is therefore a clock, not a core-cycle counter.  Frequency drift
 * affects the stages alike, so the ratios are sound; for core cycles run
 * the whole thing under "perf stat -e cpu_core/cycles/" and split that.
 * --------------------------------------------------------------------- */

#ifdef RP_PHASE_TIMING

#include <x86intrin.h>
#include <string.h>

unsigned long long _rp_phase1_cycles = 0, _rp_phase2_cycles = 0;
unsigned long long _rp_check_cycles = 0, _rp_check_calls = 0;
/* The part of the exact check that is per denominator rather than per
 * survivor: the powers of b multiplying the coefficients, recomputed
 * on the first check for a new b.  A third stage saves that part only
 * when it removes every survivor of a denominator, so the marginal
 * cost of one check is (cyc3 - cycbc)/checks, not cyc3/checks.
 * Written by _ratpoints_check_point in find_points.c . */
unsigned long long _rp_bc_cycles = 0, _rp_bc_calls = 0;
/* the per-denominator loop that steps b modulo each sieving prime */
unsigned long long _rp_bp_cycles = 0, _rp_bp_dens = 0, _rp_bp_steps = 0;
/* building a sieve table: once for each pair (prime, denominator class), so
 * at most p times for the prime p however long the run is.  rows counts the
 * bit-arrays written, which is what the cost should be proportional to. */
unsigned long long _rp_init_cycles = 0, _rp_init_calls = 0, _rp_init_rows = 0;
/* filling in sieve_spec and check_spec, which is done once per denominator
 * for each prime of the first two phases and each prime of the third */
unsigned long long _rp_setup_cycles = 0, _rp_setup_dens = 0;
/* clearing the two boundary words and the bit arrays that only exist to
 * make the count a multiple of RATPOINTS_CHUNK.  Until 2.3 this was a pass
 * over the whole array, writing the 2-adic pattern into every bit array
 * before the first phase ANDed anything into it; the first phase's first
 * prime does that now, and this is what is left.  So little is left that the
 * two rdtsc reads around it are a good part of what it reports: take the
 * figure as an upper bound on the cost, not as a measurement of it. */
unsigned long long _rp_fill_cycles = 0, _rp_fill_arrays = 0;
/* and the whole of sift(), so that what is left of it once the set-up and
 * the two phases are taken out can be seen */
unsigned long long _rp_sift_cycles = 0, _rp_sift_calls = 0;
unsigned long long _rp_sift0_calls = 0, _rp_arrays_swept = 0;
long _rp_sp1 = -1, _rp_sp2 = -1, _rp_sp3 = -1;
/* the whole run, so that core cycles from "perf stat" can be apportioned
 * to the phases: cycles_i = perf_cycles * cyc_i / cyctot */
unsigned long long _rp_t_start = 0;
static void _rp_t_begin(void) __attribute__((constructor));
static void _rp_t_begin(void) { _rp_t_start = __rdtsc(); }
#ifdef RP_PHASE_COUNTS
/* bits set on entry to phase 1, after phase 1, and after phase 2;
 * and the number of units (bit-arrays, or words under
 * USE_LONG_IN_PHASE_2) that are non-zero after phase 1 */
unsigned long long _rp_bits_in = 0, _rp_bits_1 = 0, _rp_bits_2 = 0;
unsigned long long _rp_units_surviving = 0;
/* work actually done in phase 2: AND steps in the sp2-sp1 loop (which
 * stops early once nums is empty) and iterations of the bit-extraction
 * loops.  The latter is now one per set bit: the loops used to run up to the
 * highest set bit, and the gap between ext2 and bits_2 was what said that was
 * worth changing. */
unsigned long long _rp_and2 = 0, _rp_ext2 = 0;
/* and in phase 1, where every prime is applied to every word: the product of
 * the arrays swept and sp1, which is what the phase-1 cost is proportional
 * to and hence what divides into it to give the cost of one AND. */
unsigned long long _rp_and1 = 0;

static inline unsigned _rp_popcnt(const ratpoints_bit_array *a)
{ unsigned long w[RBA_PACK]; unsigned i, c = 0;

  memcpy(w, a, sizeof(*a));
  for(i = 0; i < RBA_PACK; i++) { c += __builtin_popcountl(w[i]); }
  return c;
}
#endif

# define RP_TIC(t) unsigned long long t = __rdtsc()
# define RP_TOC(t, acc) (acc) += __rdtsc() - (t)
/* wrap one call to _ratpoints_check_point, so that the exact verification
 * can be subtracted from the cost of the second phase */
# define RP_CHECK(call) \
    ({ unsigned long long t_ = __rdtsc(); long r_ = (call); \
       _rp_check_cycles += __rdtsc() - t_; _rp_check_calls++; r_; })

static void _rp_phase_report(void) __attribute__((destructor));

static void _rp_phase_report(void)
{ unsigned long long c2;

  if(_rp_phase1_cycles + _rp_phase2_cycles == 0) { return; }
  c2 = _rp_phase2_cycles - _rp_check_cycles;
  /* one parseable line; "bits_*" are only meaningful with RP_PHASE_COUNTS */
  fprintf(stderr,
          "[phasedata] width=%d chunk=%d long2=%d sp1=%ld sp2=%ld sp3=%ld"
          " calls=%llu arrays=%llu bits_in=%llu bits_1=%llu bits_2=%llu"
          " units_1=%llu and2=%llu ext2=%llu checks=%llu bc=%llu"
          " dens=%llu bpsteps=%llu tabs=%llu rows=%llu cyctab=%llu"
          " and1=%llu cycsetup=%llu setupdens=%llu"
          " cycfill=%llu fillarrays=%llu cycsift=%llu siftcalls=%llu"
          " cyc1=%llu cyc2=%llu cyc3=%llu cycbc=%llu cycbp=%llu cyctot=%llu\n",
          (int)(8*(int)sizeof(ratpoints_bit_array)), (int)RATPOINTS_CHUNK,
#ifdef USE_LONG_IN_PHASE_2
          1,
#else
          0,
#endif
          _rp_sp1, _rp_sp2, _rp_sp3, _rp_sift0_calls, _rp_arrays_swept,
#ifdef RP_PHASE_COUNTS
          _rp_bits_in, _rp_bits_1, _rp_bits_2, _rp_units_surviving,
          _rp_and2, _rp_ext2,
#else
          0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL,
#endif
          _rp_check_calls, _rp_bc_calls, _rp_bp_dens, _rp_bp_steps,
          _rp_init_calls, _rp_init_rows, _rp_init_cycles,
#ifdef RP_PHASE_COUNTS
          _rp_and1,
#else
          0ULL,
#endif
          _rp_setup_cycles, _rp_setup_dens,
          _rp_fill_cycles, _rp_fill_arrays, _rp_sift_cycles, _rp_sift_calls,
          _rp_phase1_cycles, c2, _rp_check_cycles, _rp_bc_cycles,
          _rp_bp_cycles, __rdtsc() - _rp_t_start);
}

#else

# define RP_TIC(t)
# define RP_TOC(t, acc)
# define RP_CHECK(call) (call)

#endif /* RP_PHASE_TIMING */

/* ---------------------------------------------------------------------
 * Development instrumentation: cut the pipeline short, so that the cost
 * of one stage can be had as a difference of two whole-program cycle
 * counts and nothing has to be attributed by a timer inside the loop.
 * Build with -DRP_STOP_AFTER=<n>:
 *   1  stop after the first phase (no second phase at all)
 *   2  ... after locating the survivors (no sieving, no extraction)
 *   3  ... after sieving them with the next sp2-sp1 primes
 *   4  ... after extracting the bits, but without the exact check
 *   0 or unset: the whole thing, i.e. the real program.
 * Differences between consecutive levels give the stages, and because
 * every level sieves the first phase identically, level 1 cancels out of
 * all of them.  The truncated levels accumulate what they would have
 * used into _rp_sink, which is printed at exit so that the compiler
 * cannot drop the work whose cost is being measured.
 * --------------------------------------------------------------------- */
#ifndef RP_STOP_AFTER
# define RP_STOP_AFTER 0
#endif

#if RP_STOP_AFTER == 4
/* level 4 runs the extraction loop and its relprime test, but not the exact
 * check that they lead to */
# define RP_CHECK_POINT(a, b) (_rp_sink += 1 + (unsigned long)(a), 0L)
#else
# define RP_CHECK_POINT(a, b) \
    RP_CHECK(_ratpoints_check_point((a), (b), args, quit, process, info))
#endif

#if RP_STOP_AFTER
#include <stdio.h>
unsigned long _rp_sink = 0;
static void _rp_sink_report(void) __attribute__((destructor));
static void _rp_sink_report(void)
{ fprintf(stderr, "[stopafter] level=%d sink=%lu\n", RP_STOP_AFTER, _rp_sink); }
#endif

/* The reductions modulo a prime by a multiplication, RP_MULMOD and
 * RP_MULDIV, and the limit RP_MULMOD_LIMIT below which they are exact, are
 * in rp-private.h; find_points.c uses them too. */

/* Development switch: with -DRP_MOD_CHOICE the two ways of reducing a word
 * number modulo a prime live in the same binary, chosen by the environment
 * variable RP_MOD_MUL, so that they can be timed against each other without
 * the code-alignment difference that two builds would bring (see the note in
 * the Makefile).  The test is per prime and per call and is perfectly
 * predicted, and it is present in both arms, so it does not favour either. */
#ifdef RP_MOD_CHOICE
#include <stdlib.h>
int _rp_use_mod_mul = 1;
static void _rp_mod_choose(void) __attribute__((constructor));
static void _rp_mod_choose(void)
{ const char *e = getenv("RP_MOD_MUL");

  if(e) { _rp_use_mod_mul = atoi(e); }
  fprintf(stderr, "[modchoice] mod_mul=%d\n", _rp_use_mod_mul);
}
# define RP_USE_MOD_MUL _rp_use_mod_mul
#else
# define RP_USE_MOD_MUL 1
#endif

/* (a mod p) in [0, p), for |a| <= RP_MULMOD_LIMIT.  The sign is put back
 * afterwards rather than removed beforehand by a shift, because the callers
 * cannot bound a tightly enough for a shift to stay inside 32 bits. */
static inline long mod_mul(long a, long p, unsigned long m)
{ long r = RP_MULMOD((a < 0) ? -a : a, p, m);

  return((a < 0 && r) ? p - r : r);
}

/* Walking the set bits of a word of survivors.
 *
 * The bits still set after the second sieving stage are the numerators that
 * have to be checked exactly; bit t of a word stands for the numerator
 * a0 + d*t.  Stepping through every bit position up to the highest one set
 * costs about thirty iterations for each bit that is actually there, at the
 * survival rate the parameters aim at, so go straight to the set bits
 * instead: the lowest is at RP_CTZL(w) (rp-private.h), and w &= w-1 clears
 * it.
 */

/* Body runs once per set bit of w, with a set to that bit's numerator and t
 * to its position.  w is consumed. */
#define RP_EACH_SET_BIT(w, first, step, a, t) \
  for(; (w) && (((t) = RP_CTZL(w)), ((a) = (first) + (step)*(t)), 1); \
      (w) &= (w) - 1UL)

/* ---------------------------------------------------------------------
 * The third stage.
 *
 * What is left after the second phase is a handful of numerators per
 * denominator, and each of them would go straight to the exact check, which
 * homogenises f, evaluates it in multi-precision arithmetic and takes an
 * integer square root.  That is two orders of magnitude dearer than one
 * more sieving step, so it pays to test a few more primes first -- but not
 * by building a sieve table for them, which costs O(p) per denominator and
 * is only worth it while a whole bit array is still in play.  Here the
 * condition is evaluated one numerator at a time instead:
 *
 *   f(a/b) is a square mod p  <==>  is_f_square[(a * b^-1) mod p]
 *
 * which is what the tables encode as well, so this rejects exactly what a
 * table for that prime would have rejected, and never a genuine point --
 * the table of squares counts zero as a square.
 *
 * The stage runs after the test for common factors, not before it.  A
 * numerator sharing a factor with the denominator stands for a fraction
 * that a smaller denominator has already dealt with, so f takes the same
 * value there and every prime accepts it; those survivors can only be
 * removed by the gcd, and there is no point in testing them here first.
 * --------------------------------------------------------------------- */

static inline int stage3(long a, const check_spec *csp, long n)
{ long i;

  for(i = 0; i < n; i++)
  { long p = csp[i].p;
    long binv = csp[i].binv;
    long bias = csp[i].bias;
    long am;

    if(!binv) { continue; } /* p divides the denominator: no information */
    if(bias)
    { am = RP_MULMOD(a*binv + bias, p, csp[i].magic); }
    else
    { am = a % p;
      if(am < 0) { am += p; }
      am = (am*binv) % p;
    }
    if(!csp[i].is_f_square[am]) { return(0); }
  }
  return(1);
}

/**************************************************************************
 * check if m and n are relatively prime                                  *
 **************************************************************************/

static inline int relprime(long m, long n)
{
  /* n (the denominator) is always positive here */
  if(m == 0) { return(n == 1); }
  if(m < 0) { m = -m; }
  if(!(m & 1)) /* m is even */
  { if(!(n & 1)) { return(0); } /* n is also even */
    m >>= 1; while(!(m & 1)) { m >>= 1; } /* n odd: replace m by odd part */
  }
  while(!(n & 1)) { n >>= 1; } /* replace n by odd part */
  /* successively subtract the smaller from the larger
   * and replace the result by its odd part,
   * until both are equal (to their gcd) */
  while(n != m)
  { if(n > m)
    { n -= m; n >>= 1; while(!(n & 1)) { n >>= 1; } }
    else
    { m -= n; m >>= 1; while(!(m & 1)) { m >>= 1; } }
  }
  return(m == 1);
}

/**************************************************************************
 * Try to avoid divisions                                                 *
 **************************************************************************/

#ifdef RP_MOD_COUNTS
/* Development instrumentation: how often the helper below actually divides,
 * and how far outside [-16b, 16b) its argument is when it does. */
unsigned long long _rp_mod_calls = 0, _rp_mod_divs = 0, _rp_mod_quot = 0;
static void _rp_mod_report(void) __attribute__((destructor));
static void _rp_mod_report(void)
{ fprintf(stderr, "[moddata] sift calls=%llu divs=%llu (%.2f%%)"
                  " mean|a/b|=%.1f\n",
          _rp_mod_calls, _rp_mod_divs,
          _rp_mod_calls ? 100.0*(double)_rp_mod_divs/(double)_rp_mod_calls : 0.0,
          _rp_mod_divs ? (double)_rp_mod_quot/(double)_rp_mod_divs : 0.0);
}
# define RP_MOD_TICK(a, b) do { _rp_mod_calls++; } while(0)
# define RP_MOD_DIV(a, b) do { _rp_mod_divs++; \
      _rp_mod_quot += (unsigned long long)(((a) < 0 ? -(a) : (a))/(b)); } while(0)
#else
# define RP_MOD_TICK(a, b)
# define RP_MOD_DIV(a, b)
#endif

/* returns a mod b (for b positive) in [0,b) */
static inline long mod(long a, long b)
{
  long b1 = b << 4; /* b1 = 16*b */

  RP_MOD_TICK(a, b);
  /* if a is outside [-16*b, 16*b), then use divison */
  if(a < -b1) { RP_MOD_DIV(a, b); a %= b; if(a < 0) { a += b; } return(a); }
  if(a < 0) { a += b1; }
  else { if(a >= b1) { RP_MOD_DIV(a, b); return(a % b); } }
  /* otherwise subtract 2-power multiples of b if necessary
   * to obtain the remainder. */
  b1 >>= 1; /* b1 = 8*b */
  if(a >= b1) { a -= b1; }
  b1 >>= 1; /* b1 = 4*b */
  if(a >= b1) { a -= b1; }
  b1 >>= 1; /* b1 = 2*b */
  if(a >= b1) { a -= b1; }
  if(a >= b) { a -= b; }
  return(a);
}

/* What happens to one surviving numerator: the test for common factors,
 * then the third stage.  The two counters say how many numerators got that
 * far, which is what lets the number of primes be corrected during the run
 * instead of predicted before it; they cost one increment each on paths
 * taken a few times in a million. */
static inline int accepted(long a, long b, const check_spec *csp, long n,
                           ratpoints_args *args)
{ if(!relprime(a, b)) { return(0); }
  args->n_coprime++;
  if(!stage3(a, csp, n)) { return(0); }
  args->n_checks++;
  return(1);
}

/**************************************************************************
 * The inner loop of the sieving procedure                                *
 **************************************************************************/

/* b is the denominator;
 * the bit-arrays to be dealt with are indexed w_low..w_high-1,
 * where index 0 is the array whose zeroth bit corresponds to 0
 * (or to 1 when using only odd numerators, as specified by which_bits).
 * survivors points to space to be used for the sieving.
 * sieves points to the sieving information.
 * quit will be set when the search is stopped (because a point was found).
 * process is the function used to deal with a point that was found;
 * it is passed the pointer info, which can be used for data that
 * should persist between calls. */
long _ratpoints_sift0(long b, long w_low, long w_high,
           ratpoints_args *args, bit_selection which_bits,
           ratpoints_bit_array *survivors, ratpoints_bit_array bits16,
           long mask_low, long mask_high, long n_pad, sieve_spec *sieves,
           check_spec *checks, int *quit,
           int process(long, long, const mpz_t, void*, int*), void *info)
{
  long total = 0;
  long sp1 = args->sp1; /* number of primes in first stage */
  long sp2 = args->sp2; /* number of primes in first and second stage combined */
  long nchecks = args->sp3 - sp2; /* further primes, for the third stage */
  const unsigned long *magics = (const unsigned long *)args->magics;
  /* Whether the reductions modulo the primes below can multiply by the
   * reciprocal, which is exact for values below 2^32.  What is reduced is
   * a word number plus a sieve_spec offset, and the offsets carry a
   * multiple of the prime between 2^31 and 2^31 + p (RP_ROW_BIAS), so the
   * values lie in [0, 2^32) for every word number in this range -- which
   * reaches a height bound above 10^11.  Beyond it, mod() divides. */
  int small = (w_low >= -RP_ROW_BIAS
                && w_high <= RP_ROW_BIAS - 2*RATPOINTS_MAX_PRIME);

#ifdef DEBUG
  /* There is nothing in the survivors array to print: since 2.3 the first
   * phase's first prime writes it, so on entry it holds either nothing at
   * all (the first call) or the previous denominator's leavings.  What goes
   * into it is this pattern, with the two ends cleared after the phase. */
  { printf("\nsift0(b = %ld) @ start: %ld bit arrays from ",
           b, w_high - w_low);
    PRINT_RBA(bits16);
    printf("\n  mask_low = %ld, mask_high = %ld, padding = %ld\n",
           mask_low, mask_high, n_pad);
    fflush(NULL);
  }
#endif

  /* now do the sieving (fast!) */

  args->n_words += (unsigned long)(w_high - w_low)*RBA_PACK;

#ifdef RP_PHASE_TIMING
  _rp_sift0_calls++;
  _rp_arrays_swept += w_high - w_low;
  _rp_sp1 = sp1; _rp_sp2 = sp2; _rp_sp3 = args->sp3;
#endif
#ifdef RP_PHASE_COUNTS
  /* Every bit array starts from the same 2-adic pattern, so this is a
   * multiplication rather than a pass over the array.  It used to be a
   * population count on each of them, which was 18% of "make testhigh" and
   * sat outside every timed region -- so the instrumentation was itself the
   * larger part of the "quarter of the run in no phase" that prompted this
   * change.  It is not quite the same quantity: the two boundary words are
   * counted unmasked, so the count is high by whatever those two ends would
   * have lost.  That is negligible over a long numerator interval and is not
   * negligible over a short one, which is the case at a small height bound
   * -- so read bits_in as an upper bound, and compare it with bits_1 rather
   * than trusting it absolutely. */
  _rp_bits_in += (unsigned long long)(w_high - w_low - n_pad)
                   *_rp_popcnt(&bits16);
  _rp_and1 += (unsigned long long)(w_high - w_low)*sp1;
#endif
  RP_TIC(_rp_t1);

#ifdef DEBUG
  printf("\nsift0: sp1 = %ld, sp2 = %ld\n\n", sp1, sp2);
  fflush(NULL);
#endif

#if (defined(RATPOINTS_CHUNK) && (RATPOINTS_CHUNK > 1) && (RATPOINTS_CHUNK <= 16))
  /* Use separate variables for the individual bit-arrays;
   * they should be mapped to CPU registers.
   * This saves load/store instructions.
   * The more registers can be used, the better!
   * The code here is for up to 16 registers.
   * It will need to be extended in the obvious way to allow more,
   * e.g., 32 registers when using 512-bit vector operations. */

  /* First set the start fields for the first phase of sieving.  (The second
   * phase finds its rows directly; see there.)
   *
   * This is the busiest reduction in the program: sp1 of them for every call,
   * and a call handles at most RATPOINTS_ARRAY_SIZE bit arrays, so at a large
   * height bound there are of the order of a billion.  mod() is the wrong
   * shape for it.  Its conditional-subtraction chain only avoids the division
   * while the word number is within sixteen primes of zero, and the word
   * number grows with the height bound while the primes do not: at height
   * 2*10^5 it divides on 13% of the calls, and the rest walk a chain of
   * data-dependent branches.  Multiplying by the reciprocal does the whole
   * job in a few cycles and branchlessly. */
  { long n;

    for(n = 0; n < sp1; n++)
    { long a = w_low + sieves[n].offset;

      sieves[n].start = sieves[n].ptr
                          + ((small && RP_USE_MOD_MUL)
                               ? mod_mul(a, sieves[n].p, magics[n])
                               : mod(a, sieves[n].p));
    }
  }

  if(sp1 == 0)
  { /* No first phase at all, which -n 0 asks for.  Then nothing has written
     * the bit arrays, since the loop below folds that into the first prime
     * and there is no first prime; write the pattern here instead. */
    long i;

    for(i = w_high - w_low; i; i--) { survivors[i-1] = bits16; }
  }
  else
  { ratpoints_bit_array *surv = survivors;
    long w_low_new;

    /* Take RATPOINTS_CHUNK bit-arrays and apply phase 1 to them,
     * then repeat with the next RATPOINTS_CHUNK bit-arrays. */
    for(w_low_new = w_low; w_low_new < w_high; surv += RATPOINTS_CHUNK, w_low_new += RATPOINTS_CHUNK)
    { long n;
      /* The first prime writes the registers instead of reading them back
       * from memory, ANDing in the 2-adic pattern every bit array starts
       * from as it goes.  That is what makes the pass that used to fill the
       * array before any of this unnecessary: one store and one load per bit
       * array, on 1.7e10 of them in "make testhigh".  The boundary words and
       * the padding are dealt with after the phase instead, which comes to
       * the same thing because AND is commutative. */
      ratpoints_bit_array *siv0 = sieves[0].start;
#if (RATPOINTS_CHUNK >= 1)
      ratpoints_bit_array reg0 = bits16 & *siv0++;
#endif
#if (RATPOINTS_CHUNK >= 2)
      ratpoints_bit_array reg1 = bits16 & *siv0++;
#endif
#if (RATPOINTS_CHUNK >= 3)
      ratpoints_bit_array reg2 = bits16 & *siv0++;
#endif
#if (RATPOINTS_CHUNK >= 4)
      ratpoints_bit_array reg3 = bits16 & *siv0++;
#endif
#if (RATPOINTS_CHUNK >= 5)
      ratpoints_bit_array reg4 = bits16 & *siv0++;
#endif
#if (RATPOINTS_CHUNK >= 6)
      ratpoints_bit_array reg5 = bits16 & *siv0++;
#endif
#if (RATPOINTS_CHUNK >= 7)
      ratpoints_bit_array reg6 = bits16 & *siv0++;
#endif
#if (RATPOINTS_CHUNK >= 8)
      ratpoints_bit_array reg7 = bits16 & *siv0++;
#endif
#if (RATPOINTS_CHUNK >= 9)
      ratpoints_bit_array reg8 = bits16 & *siv0++;
#endif
#if (RATPOINTS_CHUNK >= 10)
      ratpoints_bit_array reg9 = bits16 & *siv0++;
#endif
#if (RATPOINTS_CHUNK >= 11)
      ratpoints_bit_array reg10 = bits16 & *siv0++;
#endif
#if (RATPOINTS_CHUNK >= 12)
      ratpoints_bit_array reg11 = bits16 & *siv0++;
#endif
#if (RATPOINTS_CHUNK >= 13)
      ratpoints_bit_array reg12 = bits16 & *siv0++;
#endif
#if (RATPOINTS_CHUNK >= 14)
      ratpoints_bit_array reg13 = bits16 & *siv0++;
#endif
#if (RATPOINTS_CHUNK >= 15)
      ratpoints_bit_array reg14 = bits16 & *siv0++;
#endif
#if (RATPOINTS_CHUNK >= 16)
      ratpoints_bit_array reg15 = bits16 & *siv0++;
#endif

      while(siv0 >= sieves[0].end) { siv0 -= sieves[0].p; }
      sieves[0].start = siv0;

#ifdef DEBUG
      /* the same trace the loop below prints for every other prime */
      { printf("\nsift0 after prime p = %ld, w_low_new = %ld"
               " [high numerators to the left]:\n\n",
               sieves[0].p, w_low_new);
#if (RATPOINTS_CHUNK >= 16)
        PRINT_RBA(reg15); printf("\n");
#endif
#if (RATPOINTS_CHUNK >= 15)
        PRINT_RBA(reg14); printf("\n");
#endif
#if (RATPOINTS_CHUNK >= 14)
        PRINT_RBA(reg13); printf("\n");
#endif
#if (RATPOINTS_CHUNK >= 13)
        PRINT_RBA(reg12); printf("\n");
#endif
#if (RATPOINTS_CHUNK >= 12)
        PRINT_RBA(reg11); printf("\n");
#endif
#if (RATPOINTS_CHUNK >= 11)
        PRINT_RBA(reg10); printf("\n");
#endif
#if (RATPOINTS_CHUNK >= 10)
        PRINT_RBA(reg9); printf("\n");
#endif
#if (RATPOINTS_CHUNK >= 9)
        PRINT_RBA(reg8); printf("\n");
#endif
#if (RATPOINTS_CHUNK >= 8)
        PRINT_RBA(reg7); printf("\n");
#endif
#if (RATPOINTS_CHUNK >= 7)
        PRINT_RBA(reg6); printf("\n");
#endif
#if (RATPOINTS_CHUNK >= 6)
        PRINT_RBA(reg5); printf("\n");
#endif
#if (RATPOINTS_CHUNK >= 5)
        PRINT_RBA(reg4); printf("\n");
#endif
#if (RATPOINTS_CHUNK >= 4)
        PRINT_RBA(reg3); printf("\n");
#endif
#if (RATPOINTS_CHUNK >= 3)
        PRINT_RBA(reg2); printf("\n");
#endif
#if (RATPOINTS_CHUNK >= 2)
        PRINT_RBA(reg1); printf("\n");
#endif
#if (RATPOINTS_CHUNK >= 1)
        PRINT_RBA(reg0); printf("\n");
#endif
        fflush(NULL);
      }
#endif


      for(n = 1; n < sp1; n++)
      { /* retrieve the pointer to the beginning of the relevant bits */
        ratpoints_bit_array *siv1 = sieves[n].start;
        /* This points to >= RATPOINTS_CHUNK consecutive bit-arrays
         * of information (see init.c, gen_find_point_h.c),
         * so we can safely step siv1 that many times. */

        /* perform the sieving on RATPOINTS_CHUNK registers */
#if (RATPOINTS_CHUNK >= 1)
        AND(reg0, *siv1++);
#endif
#if (RATPOINTS_CHUNK >= 2)
        AND(reg1, *siv1++);
#endif
#if (RATPOINTS_CHUNK >= 3)
        AND(reg2, *siv1++);
#endif
#if (RATPOINTS_CHUNK >= 4)
        AND(reg3, *siv1++);
#endif
#if (RATPOINTS_CHUNK >= 5)
        AND(reg4, *siv1++);
#endif
#if (RATPOINTS_CHUNK >= 6)
        AND(reg5, *siv1++);
#endif
#if (RATPOINTS_CHUNK >= 7)
        AND(reg6, *siv1++);
#endif
#if (RATPOINTS_CHUNK >= 8)
        AND(reg7, *siv1++);
#endif
#if (RATPOINTS_CHUNK >= 9)
        AND(reg8, *siv1++);
#endif
#if (RATPOINTS_CHUNK >= 10)
        AND(reg9, *siv1++);
#endif
#if (RATPOINTS_CHUNK >= 11)
        AND(reg10, *siv1++);
#endif
#if (RATPOINTS_CHUNK >= 12)
        AND(reg11, *siv1++);
#endif
#if (RATPOINTS_CHUNK >= 13)
        AND(reg12, *siv1++);
#endif
#if (RATPOINTS_CHUNK >= 14)
        AND(reg13, *siv1++);
#endif
#if (RATPOINTS_CHUNK >= 15)
        AND(reg14, *siv1++);
#endif
#if (RATPOINTS_CHUNK >= 16)
        AND(reg15, *siv1++);
#endif

        /* update the pointer for the next round
         * (RATPOINTS_CHUNK-1 bit-arrays after sieves[n].end) */
        while(siv1 >= sieves[n].end) { siv1 -= sieves[n].p; }
        sieves[n].start = siv1;

#ifdef DEBUG
        { printf("\nsift0 after prime p = %ld, w_low_new = %ld [high numerators to the left]:\n\n",
                 sieves[n].p, w_low_new);
#if (RATPOINTS_CHUNK >= 16)
          PRINT_RBA(reg15); printf("\n");
#endif
#if (RATPOINTS_CHUNK >= 15)
          PRINT_RBA(reg14); printf("\n");
#endif
#if (RATPOINTS_CHUNK >= 14)
          PRINT_RBA(reg13); printf("\n");
#endif
#if (RATPOINTS_CHUNK >= 13)
          PRINT_RBA(reg12); printf("\n");
#endif
#if (RATPOINTS_CHUNK >= 12)
          PRINT_RBA(reg11); printf("\n");
#endif
#if (RATPOINTS_CHUNK >= 11)
          PRINT_RBA(reg10); printf("\n");
#endif
#if (RATPOINTS_CHUNK >= 10)
          PRINT_RBA(reg9); printf("\n");
#endif
#if (RATPOINTS_CHUNK >= 9)
          PRINT_RBA(reg8); printf("\n");
#endif
#if (RATPOINTS_CHUNK >= 8)
          PRINT_RBA(reg7); printf("\n");
#endif
#if (RATPOINTS_CHUNK >= 7)
          PRINT_RBA(reg6); printf("\n");
#endif
#if (RATPOINTS_CHUNK >= 6)
          PRINT_RBA(reg5); printf("\n");
#endif
#if (RATPOINTS_CHUNK >= 5)
          PRINT_RBA(reg4); printf("\n");
#endif
#if (RATPOINTS_CHUNK >= 4)
          PRINT_RBA(reg3); printf("\n");
#endif
#if (RATPOINTS_CHUNK >= 3)
          PRINT_RBA(reg2); printf("\n");
#endif
#if (RATPOINTS_CHUNK >= 2)
          PRINT_RBA(reg1); printf("\n");
#endif
#if (RATPOINTS_CHUNK >= 1)
          PRINT_RBA(reg0); printf("\n");
#endif
          fflush(NULL);
        }
#endif

      }
      /* store the contents of the registers back into memory */
#if (RATPOINTS_CHUNK >= 1)
      surv[0] = reg0;
#endif
#if (RATPOINTS_CHUNK >= 2)
      surv[1] = reg1;
#endif
#if (RATPOINTS_CHUNK >= 3)
      surv[2] = reg2;
#endif
#if (RATPOINTS_CHUNK >= 4)
      surv[3] = reg3;
#endif
#if (RATPOINTS_CHUNK >= 5)
      surv[4] = reg4;
#endif
#if (RATPOINTS_CHUNK >= 6)
      surv[5] = reg5;
#endif
#if (RATPOINTS_CHUNK >= 7)
      surv[6] = reg6;
#endif
#if (RATPOINTS_CHUNK >= 8)
      surv[7] = reg7;
#endif
#if (RATPOINTS_CHUNK >= 9)
      surv[8] = reg8;
#endif
#if (RATPOINTS_CHUNK >= 10)
      surv[9] = reg9;
#endif
#if (RATPOINTS_CHUNK >= 11)
      surv[10] = reg10;
#endif
#if (RATPOINTS_CHUNK >= 12)
      surv[11] = reg11;
#endif
#if (RATPOINTS_CHUNK >= 13)
      surv[12] = reg12;
#endif
#if (RATPOINTS_CHUNK >= 14)
      surv[13] = reg13;
#endif
#if (RATPOINTS_CHUNK >= 15)
      surv[14] = reg14;
#endif
#if (RATPOINTS_CHUNK >= 16)
      surv[15] = reg15;
#endif
    }
  }

#else /* RATPOINTS_CHUNK not between 2 and 16 */

  { long n;
    long range = w_high - w_low;

    /* Write the 2-adic pattern into the bit arrays.  The chunked arm above
     * folds this into its first prime and so does not need the pass; this
     * one would have to duplicate three loops to do the same, and it is not
     * the arm anything is built with, so it pays for the pass instead. */
    { long i;

      for(i = range; i; i--) { survivors[i-1] = bits16; }
    }

    for(n = 0; n < sp1; n++)
    { ratpoints_bit_array *sieve_n = sieves[n].ptr;
        /* points to the bit-array with sieve information */
      long p = sieves[n].p; /* the prime */
      long a = -w_low - sieves[n].offset;
      long r = small ? mod_mul(a, p, magics[n]) : mod(a, p);
        /* r is such that the relevant information starts
         * at sieve_n[p-r] */
      ratpoints_bit_array *surv = survivors;
        /* pointer stepping through the survivors array */

      if(w_high < w_low + r)
      { /* If we get here, r > 0, since w_high >= w_low always.
         * In this case, we can just process range bit-arrays
         * in a row, without overshooting the end of sieve_n. */
        ratpoints_bit_array *siv1 = &sieve_n[p-r];
        ratpoints_bit_array *siv0 = siv1 + range;

        while(siv1 != siv0)
        { AND(*surv, *siv1++); surv++; }
      }
      else
      { /* Otherwise, we first have to do r steps,
         * then move the pointer to the sieve information
         * back by p, and continue in packets of p;
         * finally, there may be another partial run. */
        ratpoints_bit_array *siv1 = &sieve_n[p-r];
        ratpoints_bit_array *surv_end = &survivors[range - p];

        { long i;

          for(i = r; i; i--)
          { AND(*surv, *siv1++); surv++; }
        }
        siv1 -= p;
        while(surv <= surv_end)
        { long i;

          for(i = p; i; i--)
          { AND(*surv, *siv1++); surv++; }
          siv1 -= p;
        }
        surv_end += p;
        while(surv < surv_end)
        { AND(*surv, *siv1++); surv++; }
      }

#ifdef DEBUG
      { long k, c = 0;

        printf("\nsift0 after prime p = %ld [high numerators to the left]:", p);
        for(k = range - 1; k >= 0; k--, c++)
        { if((c & (0xff >> RBA_SHIFT)) == 0) { printf("\n"); }
          PRINT_RBA(survivors[k]);
        }
        printf("\n");
        fflush(NULL);
      }
#endif

    }

  }
#endif /* RATPOINTS_CHUNK */

  RP_TOC(_rp_t1, _rp_phase1_cycles);

  /* The two ends of the numerator interval, and the bit arrays that only
   * exist to make the count a multiple of RATPOINTS_CHUNK.  This used to be
   * done before the first phase, on the same pass that wrote the 2-adic
   * pattern into every bit array; since AND is commutative, clearing the
   * bits afterwards clears the same bits, and doing it here is two bit
   * arrays per call rather than all of them. */
  { RP_TIC(t_fill);

    if(mask_low) { MASKL(survivors, mask_low); }
    if(mask_high) { MASKU(&survivors[w_high - w_low - n_pad - 1], mask_high); }
    if(n_pad)
    { long i;

      for(i = w_high - w_low - n_pad; i < w_high - w_low; i++)
      { survivors[i] = zero; }
    }
    RP_TOC(t_fill, _rp_fill_cycles);
#ifdef RP_PHASE_TIMING
    /* what this region writes, which since 2.3 is at most the two boundary
     * words and the padding -- not the whole range, which is what
     * _rp_arrays_swept already counts */
    _rp_fill_arrays += (unsigned long long)((mask_low ? 1 : 0)
                                             + (mask_high ? 1 : 0) + n_pad);
#endif
  }

#ifdef DEBUG
  { long n, c = 0;

    printf("\nsift0(b = %ld) after phase 1 [high numerators to the left]:\n", b);
    for(n = w_high - w_low - 1; n >= 0; n--, c++)
    { if((c & (0xff >> RBA_SHIFT)) == 0) { printf("\n"); }
      PRINT_RBA(survivors[n]);
    }
    printf("\n\n");
    fflush(NULL);
  }
#endif

  RP_TIC(_rp_t2);

#if RP_STOP_AFTER != 1

  /* Second phase of the sieve: test each surviving bit array with more primes */
  { ratpoints_bit_array *surv0 = &survivors[0];
    ratpoints_bit_array *surv_end = &survivors[w_high - w_low];
    long i;

    /* Step through the survivors array.  sp1 is chosen so that only a few
     * per cent of the bit-arrays are non-empty here, so nearly all the work
     * is stepping over the empty ones; do that in a tight loop of its own
     * rather than re-entering the body below every time.  (Or-ing several
     * bit-arrays together and stepping over the whole group at once was
     * tried, and is slower at this survival rate: a group that contains a
     * survivor has wasted its whole "or", and that happens often enough to
     * cost more than the tests it saves.)
     * The loop has no bound test either.  It steps until it meets a
     * non-zero bit array, and the one written here, just past the range,
     * is what it meets when nothing survived (find_points_work leaves room
     * for it).  With the bound test the scan was nine instructions per bit
     * array, two of them branches; now it is four -- the load, the test,
     * the branch and the pointer step -- and the position is recovered from
     * the pointer once per survivor instead. */
    *surv_end = ~zero;
    for(;;)
    { ratpoints_bit_array nums;
      long n;
#ifndef USE_LONG_IN_PHASE_2
      sieve_spec *ssp = &sieves[sp1];
#endif

      while(TESTZ(*surv0)) { surv0++; }
      if(surv0 >= surv_end) { break; }
      i = w_low + (surv0 - survivors);
      nums = *surv0++;
      args->n_arrays++;

#if RP_STOP_AFTER == 2
      if(TEST(nums)) { _rp_sink += EXT0(nums); }
      continue;
#endif

#ifdef RP_PHASE_COUNTS
      if(TEST(nums))
      { _rp_units_surviving++; _rp_bits_1 += _rp_popcnt(&nums); }
#endif

#ifdef DEBUG
      if(TEST(nums))
      { printf("\nsurviving word ");
        PRINT_RBA(nums);
        printf(" @ i = %ld\n", i);
        fflush(NULL);
      }
#endif

#ifdef USE_LONG_IN_PHASE_2
      /* Keep the first phase and the scan at the full register width, but do
       * the rest one 64-bit word at a time: of a surviving bit-array, only
       * the words that are themselves non-zero are sieved with the remaining
       * primes and then extracted.  No different table layout is needed --
       * the tables hold RBA_PACK copies of the pattern, so the word wanted is
       * word k of the table row for the bit array, found as in the other
       * arm.
       * This is slower than sieving the whole bit-array at once, by a few per
       * cent at 128 and 256 bits, and it is instructive that it is: the
       * number of AND steps is the same either way, because the other words
       * of the array were already zero, and a narrow and a wide read of the
       * same table come from one cache line.  See PHASE-NOTES.md on the
       * phases-by-register-width branch. */
      { long a0, da, d, k;

        if(which_bits == num_all)
        { d = 1; a0 = i * RBA_LENGTH; da = LONG_LENGTH; }
        else
        { d = 2; a0 = i * (2*RBA_LENGTH); da = 2*LONG_LENGTH;
          if(which_bits == num_odd) { a0++; }
        }

        for(k = 0; k < RBA_PACK; k++)
        { unsigned long numsk = EXT(nums, k);
          sieve_spec *sspk = &sieves[sp1];
          const unsigned long *mg = &magics[sp1];
          long a, t, a0k = a0 + k*da;  /* first numerator of word no. k */

          if(!numsk) { continue; }

          for(n = sp2-sp1; n && numsk; n--, sspk++, mg++)
          { /* word k of the table row for this bit array; the row is found
             * as in the other arm below */
            long p = sspk->p;
            long v = i + sspk->offset;
            long row = small ? RP_MULMOD(v, p, *mg) : mod(v, p);

#ifdef RP_PHASE_COUNTS
            _rp_and2++;
#endif
            numsk &= ((const unsigned long *)sspk->ptr)[RBA_PACK*row + k];
          }

#if RP_STOP_AFTER == 3
          _rp_sink += numsk; continue;
#endif
#ifdef RP_PHASE_COUNTS
          _rp_bits_2 += __builtin_popcountl(numsk);
#endif

          RP_EACH_SET_BIT(numsk, a0k, d, a, t)
          {
#ifdef RP_PHASE_COUNTS
            _rp_ext2++;
#endif
            args->n_bits++;
            if(accepted(a, b, checks, nchecks, args))
            { total += RP_CHECK_POINT(a, b);
              if(*quit) return(total);
            }
          }
        }
      }
#else
      /* Sieve with the next sp2-sp1 primes while some bits are set.  The
       * table row for word number i is the one at index (i + offset) mod p
       * (see sieve_spec in rp-private.h), and the reduction multiplies by
       * the reciprocal wherever small says that is exact.  Until 2.3 the
       * row was reached from a pointer set up at the head of the call, by
       * subtracting p until it pointed back into the table: one to three
       * data-dependent branches per AND, most of them mispredicted, and one
       * reduction per prime and call whether or not a single bit array had
       * survived. */
      { const unsigned long *mg = &magics[sp1];

        for(n = sp2-sp1; n && TEST(nums); n--, ssp++, mg++)
        { long p = ssp->p;
          long v = i + ssp->offset;

#ifdef RP_PHASE_COUNTS
          _rp_and2++;
#endif
          AND(nums, ssp->ptr[small ? RP_MULMOD(v, p, *mg) : mod(v, p)]);

#ifdef DEBUG
          printf("after prime p = %ld:\n ", p);
          PRINT_RBA(nums);
          printf("\n");
          fflush(NULL);
#endif
        }
      }

#if RP_STOP_AFTER == 3
      _rp_sink += EXT0(nums); continue;
#endif

#ifdef RP_PHASE_COUNTS
      _rp_bits_2 += _rp_popcnt(&nums);
#endif

      /* Check the survivors of the sieve if they really give points. */
      if(TEST(nums))
      { long a0, a, da, d, t;
        /* a  := the numerator corresponding to the selected bit
         * a0 := numerator corresponding to the lowest bit
         * d  := step size in numerators from one bit to the next
         * da := step size in numerators from one word to the next */

        /* Set d, a0, da according to which_bits. */
        if(which_bits == num_all)
        { d = 1; a0 = i * RBA_LENGTH; da = LONG_LENGTH; }
        else
        { d = 2; a0 = i * (2*RBA_LENGTH); da = 2*LONG_LENGTH;
          if(which_bits == num_odd) { a0++; }
        }

        { /* extract the first word */
          unsigned long nums0 = EXT0(nums);

          RP_EACH_SET_BIT(nums0, a0, d, a, t)
          { /* one bit that is set */

#ifdef RP_PHASE_COUNTS
            _rp_ext2++;
#endif
            args->n_bits++;

#ifdef DEBUG
            printf("\nsurviving bit no. %ld --> a = %ld. ", t, a);
            if(accepted(a, b, checks, nchecks, args))
            { printf("Check point...\n");
              fflush(NULL);
              total += RP_CHECK_POINT(a, b);
              if(*quit) return(total); /* if quit was set, stop */
            }
            else
            { printf("Not in lowest terms, or rejected by the third stage"
                     " --> skip.\n"); fflush(NULL); }
#else
            if(accepted(a, b, checks, nchecks, args))
            /* the fraction a/b is in lowest terms and survives the third
             * stage: check if we really get a point, and if so, process it. */
            { total += RP_CHECK_POINT(a, b);
              if(*quit) return(total); /* if quit was set, stop */
            }
#endif

          }

          { /* process the remaining words */
            long k;

            for (k = 1; k < RBA_PACK; k++)
            {
              unsigned long nums1 = EXT(nums,k);

              a0 += da; /* numerator corresponding to first bit of word no. k */
              RP_EACH_SET_BIT(nums1, a0, d, a, t)
              { /* one bit that is set */

#ifdef RP_PHASE_COUNTS
                _rp_ext2++;
#endif
                args->n_bits++;

#ifdef DEBUG
                printf("\nsurviving bit no. %ld --> a = %ld. ",
                       LONG_LENGTH*k + t, a);
                if(accepted(a, b, checks, nchecks, args))
                { printf("Check point...\n");
                  fflush(NULL);
                  total += RP_CHECK_POINT(a, b);
                  if(*quit) return(total); /* if quit was set, stop */
                }
                else
                { printf("Not in lowest terms, or rejected by the third stage"
                         " --> skip.\n"); fflush(NULL); }
#else
                if(accepted(a, b, checks, nchecks, args))
                { total += RP_CHECK_POINT(a, b);
                  if(*quit) return(total);
                }
#endif

              }
            }
          }
        }
      }
#endif /* USE_LONG_IN_PHASE_2 */
    }
  }
#endif /* RP_STOP_AFTER != 1 */

  RP_TOC(_rp_t2, _rp_phase2_cycles);

  return(total);
}
