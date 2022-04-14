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
 * gen_init_sieve_h.c                                                  *
 *                                                                     *
 * This program writes the file init_sieve.h                           *
 *                                                                     *
 * Michael Stoll, Jan 9, 2008                                          *
 ***********************************************************************/

#include "rp-private.h"
#include "primes.h"

ratpoints_bit_array work[RATPOINTS_MAX_PRIME];

int main(int argc, char *argv[])
{
  long n;

  for(n = 0; n < RATPOINTS_NUM_PRIMES; n++)
  { long p = prime[n];
    if(p < LONG_LENGTH) { printf("CODE_INIT_SIEVE1(%ld)\n", p); }
    else { printf("CODE_INIT_SIEVE2(%ld)\n", p); }
  }
  printf("\n");

  printf("ratpoints_init_fun sieve_init[RATPOINTS_NUM_PRIMES] = \n{");
  for(n = 0; n < RATPOINTS_NUM_PRIMES; n++)
  { long p = prime[n];

    { printf("&sieve_init_%ld", p); }
    if(n < RATPOINTS_NUM_PRIMES - 1) printf(",\n ");
  }
  printf("};\n\n");

  return(0);
}
