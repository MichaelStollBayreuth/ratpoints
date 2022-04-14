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

/* returns a mod b (for b positive) in [0,b) */
static inline long mod(long a, long b)
{
  long b1 = b << 4; /* b1 = 16*b */

  /* if a is outside [-16*b, 16*b), then use divison */
  if(a < -b1) { a %= b; if(a < 0) { a += b; } return(a); }
  if(a < 0) { a += b1; }
  else { if(a >= b1) { return(a % b); } }
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
           ratpoints_bit_array *survivors, sieve_spec *sieves, int *quit,
           int process(long, long, const mpz_t, void*, int*), void *info)
{
  long total = 0;
  long sp1 = args->sp1; /* number of primes in first stage */
  long sp2 = args->sp2; /* number of primes in first and second stage combined */

#ifdef DEBUG
  { long n, c = 0;
    printf("\nsift0(b = %ld) @ start [high numerators to the left]:\n", b);
    for(n = w_high - w_low - 1; n >= 0; n--, c++)
    { if((c & (0xff >> RBA_SHIFT)) == 0) { printf("\n"); }
      PRINT_RBA(survivors[n]);
    }
    printf("\n");
    fflush(NULL);
  }
#endif

  /* now do the sieving (fast!) */

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

  /* first set the start fields for the first and second phases of sieving */
  { long n;

    for(n = 0; n < sp2; n++)
    { sieves[n].start = sieves[n].ptr + mod(w_low + sieves[n].offset, sieves[n].p); }
  }

  { ratpoints_bit_array *surv = survivors;
    long w_low_new;

    /* Take RATPOINTS_CHUNK bit-arrays and apply phase 1 to them,
     * then repeat with the next RATPOINTS_CHUNK bit-arrays. */
    for(w_low_new = w_low; w_low_new < w_high; surv += RATPOINTS_CHUNK, w_low_new += RATPOINTS_CHUNK)
    { long n;
      /* read data from memory into registers */
#if (RATPOINTS_CHUNK >= 1)
      ratpoints_bit_array reg0 = surv[0];
#endif
#if (RATPOINTS_CHUNK >= 2)
      ratpoints_bit_array reg1 = surv[1];
#endif
#if (RATPOINTS_CHUNK >= 3)
      ratpoints_bit_array reg2 = surv[2];
#endif
#if (RATPOINTS_CHUNK >= 4)
      ratpoints_bit_array reg3 = surv[3];
#endif
#if (RATPOINTS_CHUNK >= 5)
      ratpoints_bit_array reg4 = surv[4];
#endif
#if (RATPOINTS_CHUNK >= 6)
      ratpoints_bit_array reg5 = surv[5];
#endif
#if (RATPOINTS_CHUNK >= 7)
      ratpoints_bit_array reg6 = surv[6];
#endif
#if (RATPOINTS_CHUNK >= 8)
      ratpoints_bit_array reg7 = surv[7];
#endif
#if (RATPOINTS_CHUNK >= 9)
      ratpoints_bit_array reg8 = surv[8];
#endif
#if (RATPOINTS_CHUNK >= 10)
      ratpoints_bit_array reg9 = surv[9];
#endif
#if (RATPOINTS_CHUNK >= 11)
      ratpoints_bit_array reg10 = surv[10];
#endif
#if (RATPOINTS_CHUNK >= 12)
      ratpoints_bit_array reg11 = surv[11];
#endif
#if (RATPOINTS_CHUNK >= 13)
      ratpoints_bit_array reg12 = surv[12];
#endif
#if (RATPOINTS_CHUNK >= 14)
      ratpoints_bit_array reg13 = surv[13];
#endif
#if (RATPOINTS_CHUNK >= 15)
      ratpoints_bit_array reg14 = surv[14];
#endif
#if (RATPOINTS_CHUNK >= 16)
      ratpoints_bit_array reg15 = surv[15];
#endif

      for(n = 0; n < sp1; n++)
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

    for(n = 0; n < sp1; n++)
    { ratpoints_bit_array *sieve_n = sieves[n].ptr;
        /* points to the bit-array with sieve information */
      long p = sieves[n].p; /* the prime */
      long r = mod(-w_low - sieves[n].offset, p);
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

    /* initialize pointers in sieve for the second phase */
    for(n = sp1; n < sp2; n++)
    { sieves[n].start = sieves[n].ptr + mod(w_low + sieves[n].offset, sieves[n].p); }
  }
#endif /* RATPOINTS_CHUNK */

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

#ifdef USE_LONG_IN_PHASE_2
  /* Work with unsigned longs instead of bit-arrays */
  /* Second phase of the sieve: test each surviving word in the bit-array
   * with more primes */
  { unsigned long *surv0 = (unsigned long *)survivors;
    long i, base = 0;

    for(i = RBA_PACK*w_low; i < RBA_PACK*w_high; i++, base++)
    { unsigned long nums = *surv0++;
      sieve_spec *ssp = &sieves[sp1];
      long n;

#ifdef DEBUG
      if(nums)
      { printf("\nsurviving word %*.*lx @ i = %ld\n", WIDTH, WIDTH, nums, i);
        fflush(NULL);
      }
#endif

      for(n = sp2-sp1; n && nums; n--)
      { unsigned long *ptr = (unsigned long *)ssp->start;
        long pp = RBA_PACK*ssp->p;

        ptr += base;
        while(ptr >= (unsigned long *)ssp->end) { ptr -= pp; }
        nums &= *ptr;

#ifdef DEBUG
        printf("after prime p = %ld: %*.*lx\n ", ssp->p, WIDTH, WIDTH, nums);
        fflush(NULL);
#endif

        ssp++;
      }

      /* Check the survivors of the sieve if they really give points */
      if(nums)
      { long a0, a, d;
             /* a will be the numerator corresponding to the selected bit */
#ifdef DEBUG
        long bit = 0;
#endif

        /* a0 := numerator corresponding to lowest bit,
         * d  := step size in numerators from one bit to the next */
        if(which_bits == num_all)
        { d = 1; a0 = i << LONG_SHIFT; }
        else
        { d = 2; a0 = i << (LONG_SHIFT+1);
          if(which_bits == num_odd) { a0++; }
        }

        for(a = a0; nums; a += d, nums >>= 1)
        { /* test one bit */

#ifdef DEBUG
          if(nums & 1)
          { printf("\nsurviving bit no. %ld --> a = %ld. ", bit, a);
            if(relprime(a, b))
            { printf("Check point...\n");
              fflush(NULL);
              total += _ratpoints_check_point(a, b, args, quit, process, info);
              if(*quit) return(total); /* if quit was set, stop */
            }
            else
            { printf("Not in lowest terms --> skip.\n"); fflush(NULL); }
          }
          bit++;
#else
          if((nums & 1) && relprime(a, b))
          { total += _ratpoints_check_point(a, b, args, quit, process, info);
            if(*quit) return(total); /* if quit was set, stop */
          }
#endif

        }
      }
    }
  }

#else /* not defined(USE_LONG_IN_PHASE_2) */

  /* Second phase of the sieve: test each surviving bit array with more primes */
  { ratpoints_bit_array *surv0 = &survivors[0];
    long i, base = 0;

    /* Step through the survivors array */
    for(i = w_low; i < w_high; i++, base++)
    { ratpoints_bit_array nums = *surv0++;
      sieve_spec *ssp = &sieves[sp1];
      long n;

#ifdef DEBUG
      if(TEST(nums))
      { printf("\nsurviving word ");
        PRINT_RBA(nums);
        printf(" @ i = %ld\n", i);
        fflush(NULL);
      }
#endif

      /* Sieve with the next sp2-sp1 primes while some bits are set. */
      for(n = sp2-sp1; n && TEST(nums); n--)
      { ratpoints_bit_array *ptr = (ssp->start) + base;
        long p = ssp->p;

        while(ptr >= ssp->end) { ptr -= p; }
        AND(nums, *ptr);

#ifdef DEBUG
        printf("after prime p = %ld:\n ", p);
        PRINT_RBA(nums);
        printf("\n");
        fflush(NULL);
#endif

        ssp++;
      }

      /* Check the survivors of the sieve if they really give points. */
      if(TEST(nums))
      { long a0, a, da, d;
        /* a  := the numerator corresponding to the selected bit
         * a0 := numerator corresponding to the lowest bit
         * d  := step size in numerators from one bit to the next
         * da := step size in numerators from one word to the next */

#ifdef DEBUG
        long bit = 0; /* counter for the bits, used for output */
#endif

        /* Set d, a0, da according to which_bits. */
        if(which_bits == num_all)
        { d = 1; a0 = i << RBA_SHIFT; da = LONG_LENGTH; }
        else
        { d = 2; a0 = i << (RBA_SHIFT+1); da = 2*LONG_LENGTH;
          if(which_bits == num_odd) { a0++; }
        }

        { /* extract the first word */
          unsigned long nums0 = EXT0(nums);

          for(a = a0; nums0; a += d, nums0 >>= 1)
          { /* test one bit */

#ifdef DEBUG
            if(nums0 & 1)
            { printf("\nsurviving bit no. %ld --> a = %ld. ", bit, a);
              if(relprime(a, b))
              { printf("Check point...\n");
                fflush(NULL);
                total += _ratpoints_check_point(a, b, args, quit, process, info);
                if(*quit) return(total); /* if quit was set, stop */
              }
              else
              { printf("Not in lowest terms --> skip.\n"); fflush(NULL); }
            }
            bit++;
#else
            if((nums0 & 1) && relprime(a, b))
            /* bit is set and fraction a/b is in lowest terms:
             * check if we really get a point, and if so, process it. */
            { total += _ratpoints_check_point(a, b, args, quit, process, info);
              if(*quit) return(total); /* if quit was set, stop */
            }
#endif

          }

          { /* process the remaining words */
            long k;

            for (k = 1; k < RBA_PACK; k++)
            {
              unsigned long nums1 = EXT(nums,k);

#ifdef DEBUG
              bit = LONG_LENGTH * k;
#endif

              a0 += da; /* numerator corresponding to first bit of word no. k */
              for (a = a0; nums1; a += d, nums1 >>= 1)
              { /* test one bit */

#ifdef DEBUG
                if(nums1 & 1)
                { printf("\nsurviving bit no. %ld --> a = %ld. ", bit, a);
                  if(relprime(a, b))
                  { printf("Check point...\n");
                    fflush(NULL);
                    total += _ratpoints_check_point(a, b, args, quit, process, info);
                    if(*quit) return(total); /* if quit was set, stop */
                  }
                  else
                  { printf("Not in lowest terms --> skip.\n"); fflush(NULL); }
                }
                bit++;
#else
                if((nums1 & 1) && relprime(a, b))
                { total += _ratpoints_check_point(a, b, args, quit, process, info);
                  if(*quit) return(total);
                }
#endif

              }
            }
          }
        }
      }
      /* Attempt to save some subtractions, but no improvement... */
      /*
      if(base == BASE_REPEAT)
      { for(n = sp1; n < sp2; n++)
        { ratpoints_bit_array *start = sieves[n].start + BASE_REPEAT;

          while(start >= sieves[n].end) { start -= sieves[n].p; }
          sieves[n].start = start;
        }
        base = 0;
      }
      */
    }
  }
#endif /* USE_LONG_IN_PHASE_2 */

  return(total);
}
