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
 * gen_find_points_h.c                                                 *
 *                                                                     *
 * This program writes the file find_points.h                          *
 *                                                                     *
 * Michael Stoll, Mar 8, 2009                                          *
 * with changes by Bill Allombert, Dec 29, 2021                        *
 ***********************************************************************/

#include "rp-private.h"
#include "primes.h"

long inv_mod_p(long p, long b)
{ /* invert b mod p */
  /* doesn't have to be fast...
     -- actually this will be faster in our range than the usual XGCD. */
  long i = 1, n = b;

  while(1)
  { if(n%p == 1) { return(i); }
    i++;
    n += b;
  }
}

int main(int argc, char *argv[])
{
  long n;

  { int work[RATPOINTS_MAX_PRIME];

    printf("static const int "
           "squares[RATPOINTS_NUM_PRIMES+1][RATPOINTS_MAX_PRIME] =\n{\n");
    for(n = 0; n < RATPOINTS_NUM_PRIMES; n++)
    { long p = prime[n];
      long i;

      work[0] = 1;
      for(i = 1; i < p; i++) work[i] = 0;
      /* record non-zero squares mod p, p odd */
      for(i = 1; i < p; i += 2) work[(i*i) % p] = 1;
      printf("{");
      for(i = 0; i < p; i++)
      { printf("%d", work[i]);
        if(i < p-1) printf(",");
      }
      printf((n < RATPOINTS_NUM_PRIMES - 1) ? "},\n " : "}\n};\n");
    }
  }

  printf("static const long offsets[RATPOINTS_NUM_PRIMES] =\n{");
  for(n = 0; n < RATPOINTS_NUM_PRIMES; n++)
  { long p = prime[n];

    { printf("%ld", inv_mod_p(p, (2*RBA_LENGTH)%p)); }
    printf((n < RATPOINTS_NUM_PRIMES - 1) ? "," : "};\n\n");
  }

  printf("static const long "
         "inverses[RATPOINTS_NUM_PRIMES][RATPOINTS_MAX_PRIME] =\n{");
  for(n = 0; n < RATPOINTS_NUM_PRIMES; n++)
  { long p = prime[n];
    long i;

    printf("{0");
    for(i = 1; i < p; i++)
    { printf(",%ld", inv_mod_p(p, i)); }
    printf((n < RATPOINTS_NUM_PRIMES - 1) ? "},\n " : "}\n};\n");
  }

  { unsigned long work[RATPOINTS_MAX_PRIME];

    printf("unsigned long "
           "sieves0[RATPOINTS_NUM_PRIMES][RBA_PACK*(RATPOINTS_MAX_PRIME_EVEN +  RATPOINTS_CHUNK-1)] =\n{\n");
    for(n = 0; n < RATPOINTS_NUM_PRIMES; n++)
    { long p = prime[n];
      long i, j;

      for(i = 0; i < p; i++) { work[i] = ~0UL; }
      for(i = 0; i < LONG_LENGTH; i++)
      { work[(p*i)>>LONG_SHIFT] &= ~(1UL<<((p*i) & LONG_MASK)); }

      printf("{");
      for(j = 0; j < RBA_PACK; j++)
      { for(i = 0; i < p; i++)
        { printf("0x%*.*lx, ", 16, 16, work[i]); }
      }
      for(i = 0; i < RBA_PACK*(RATPOINTS_CHUNK-1); i++)
      { printf("0x%*.*lx", 16, 16, work[i%p]);
        if(i < RBA_PACK*(RATPOINTS_CHUNK-1) - 1) { printf(", "); }
      }
      printf((n < RATPOINTS_NUM_PRIMES - 1) ? "},\n" : "}\n");
    }
    printf("};\n\n");
  }

  return(0);
}
