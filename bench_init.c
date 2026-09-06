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
 * bench_init.c                                                        *
 *                                                                     *
 * Correctness check and benchmark for the sieve_init_<p> functions    *
 * in init.c .  Not part of the library; build it with                 *
 *   make bench_init                                                   *
 * and run it as                                                       *
 *   ./bench_init [reps [pmin [pmax]]]                                 *
 *                                                                     *
 * It first rebuilds every sieve table (all primes, all denominator    *
 * classes) and compares it word for word against a reference computed *
 * directly from the definition                                        *
 *   si[x] bit i  <==>  is_f_square[((x*LONG_LENGTH + i) * b^-1) mod p] *
 * so that a faster implementation of sieve_init can be checked        *
 * without running the whole test suite; then it times the same sweep. *
 * Restricting the primes with pmin/pmax is useful because the primes  *
 * above LONG_LENGTH (CODE_INIT_SIEVE2) dominate the total cost.       *
 *                                                                     *
 * Michael Stoll, Sep 6, 2026                                          *
 ***********************************************************************/
#include "rp-private.h"
#include "primes.h"
#include <time.h>

extern ratpoints_init_fun sieve_init[RATPOINTS_NUM_PRIMES];
extern void *pointer_align(void *xx, long m);

static double now(void)
{ struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + 1e-9*ts.tv_nsec; }

static long inv_mod_p(long p, long b)
{ long i = 1, n = b; while(1) { if(n%p == 1) return i; i++; n += b; } }

/* Independent reference for one table, straight from the definition:
 *   si[x] bit i  <==>  is_f_square[ ((x*LONG_LENGTH + i) * d) mod p ],  d = b^-1
 * then replicated with period p over RBA_PACK*p words, then CHUNK-1 wrap copies. */
static void reference(long p, long d, const int *isfs, unsigned long *ref)
{ long x, i, k;
  for(x = 0; x < p; x++)
  { unsigned long w = 0;
    for(i = 0; i < LONG_LENGTH; i++)
    { long j = ((x*LONG_LENGTH + i) % p) * d % p;
      if(isfs[j]) { w |= 1UL << i; }
    }
    ref[x] = w;
  }
  for(k = 1; k < RBA_PACK; k++) { for(x = 0; x < p; x++) ref[x + k*p] = ref[x]; }
  for(k = 0; k < (RATPOINTS_CHUNK-1)*RBA_PACK; k++) ref[p*RBA_PACK + k] = ref[k];
}

int main(int argc, char *argv[])
{
  long reps = (argc > 1) ? atol(argv[1]) : 40;
  long pmin = (argc > 2) ? atol(argv[2]) : 0;
  long pmax = (argc > 3) ? atol(argv[3]) : 1000;
  long need = 0, pn, b, seed = 12345;
  int  *isfs_buf = malloc(RATPOINTS_NUM_PRIMES*(RATPOINTS_MAX_PRIME+1)*sizeof(int));
  long *inv_buf  = malloc(RATPOINTS_NUM_PRIMES*RATPOINTS_MAX_PRIME*sizeof(long));
  unsigned long *ref = malloc((RATPOINTS_MAX_PRIME*RBA_PACK
                               + RATPOINTS_CHUNK*RBA_PACK)*sizeof(unsigned long));
  ratpoints_sieve_entry *se = malloc(RATPOINTS_NUM_PRIMES*sizeof(ratpoints_sieve_entry));
  ratpoints_args args;
  void *ba_na;
  long bad = 0, tables = 0;
  double t;

  for(pn = 0; pn < RATPOINTS_NUM_PRIMES; pn++)
    need += prime[pn]*(prime[pn] + RATPOINTS_CHUNK-1);
  ba_na = malloc((need+1)*sizeof(ratpoints_bit_array));
  args.ba_buffer = pointer_align(ba_na, sizeof(ratpoints_bit_array));

  /* pseudo-random but fixed is_f_square data, and the inverse tables */
  for(pn = 0; pn < RATPOINTS_NUM_PRIMES; pn++)
  { long p = prime[pn], i;
    int  *isfs = isfs_buf + pn*(RATPOINTS_MAX_PRIME+1);
    long *inv  = inv_buf  + pn*RATPOINTS_MAX_PRIME;
    for(i = 0; i <= p; i++)
    { seed = seed*6364136223846793005L + 1442695040888963407L;
      isfs[i] = (int)((seed >> 40) & 1); }
    inv[0] = 0;
    for(i = 1; i < p; i++) inv[i] = inv_mod_p(p, i);
    se[pn].init = sieve_init[pn]; se[pn].p = p;
    se[pn].is_f_square = isfs; se[pn].inverses = inv; se[pn].offset = 0;
  }

  /* ---- correctness: every prime, every denominator class ---- */
  for(pn = 0; pn < RATPOINTS_NUM_PRIMES; pn++)
  { long p = prime[pn];
    for(b = 1; b < p; b++)
    { unsigned long *si;
      long w = p*RBA_PACK + (RATPOINTS_CHUNK-1)*RBA_PACK, j;
      args.ba_next = args.ba_buffer;
      si = (unsigned long *)(*(se[pn].init))(&se[pn], b, &args);
      reference(p, inv_buf[pn*RATPOINTS_MAX_PRIME + b], se[pn].is_f_square, ref);
      for(j = 0; j < w; j++)
        if(si[j] != ref[j])
        { if(bad < 5) printf("MISMATCH p=%ld b=%ld word %ld: got %016lx want %016lx\n",
                             p, b, j, si[j], ref[j]);
          bad++; break; }
      tables++;
  } }
  printf("%s: %ld tables checked, %ld mismatching\n",
         bad ? "*** FAIL ***" : "OK", tables, bad);

  /* ---- benchmark: sweep all (p,b) pairs, reps times ---- */
  t = now();
  { long r;
    for(r = 0; r < reps; r++)
    { args.ba_next = args.ba_buffer;
      for(pn = 0; pn < RATPOINTS_NUM_PRIMES; pn++)
      { long p = prime[pn];
        if(p < pmin || p > pmax) continue;
        for(b = 1; b < p; b++) (*(se[pn].init))(&se[pn], b, &args);
      }
    }
  }
  t = now() - t;
  { long nt = 0; for(pn = 0; pn < RATPOINTS_NUM_PRIMES; pn++)
    { if(prime[pn] >= pmin && prime[pn] <= pmax) nt += prime[pn]-1; }
    printf("  %ld reps x %ld tables : %7.4f s  -> %7.1f ns/table\n",
         reps, nt, t, 1e9*t/(reps*(double)nt)); }
  return bad ? 1 : 0;
}
