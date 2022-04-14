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
 * rptest.c                                                            *
 *                                                                     *
 * Test program for ratpoints                                          *
 *                                                                     *
 * Michael Stoll, May 27, 2009; January 2, 2022                        *
 ***********************************************************************/

#include <stdlib.h>
#include <stdio.h>

#include "ratpoints.h"

#include "testdata.h"

mpz_t c[RATPOINTS_MAX_DEGREE+1];  /* The coefficients of f */

ratpoints_interval domain[2*RATPOINTS_MAX_DEGREE];

/**************************************************************************
 * function that processes the points                                     *
 **************************************************************************/

/* data structure to be passed to process() */
typedef struct {int print_between; int no_output; int one_point; int no_y;}
        data;

/* process a point (a:y:b) by printing "[a,y,b]" (unless no_output is set) */
int process(long a, long b, const mpz_t y, void *info0, int *quit)
{ data *info = (data *)info0;

  *quit = 0;
  if(info->no_output) return(1);
  if(info->print_between) { printf(","); }
  else { info->print_between = 1; }
  printf("[%ld,", a);
  mpz_out_str((FILE *)NULL, 10, y);
  printf(",%ld]", b);
  return(1);
}

int main(int argc, char *argv[])
{
  long n;
  ratpoints_args args;

  /* parameters; see documentation */
  long degree        = 6;
  long height        = 16383;
  long sieve_primes1 = RATPOINTS_DEFAULT_SP1;
  long sieve_primes2 = RATPOINTS_DEFAULT_SP2;
  long num_primes    = RATPOINTS_DEFAULT_NUM_PRIMES;
  long max_forbidden = RATPOINTS_DEFAULT_MAX_FORBIDDEN;
  long b_low         = 1;
  long b_high        = height;
  long sturm_iter    = RATPOINTS_DEFAULT_STURM;
  long array_size    = 0; /* will be set to default value in find_points_work */
  int no_check       = 0;
  int no_y           = 0;
  int no_reverse     = 0;
  int no_jacobi      = 0;
  int no_output      = 0;
  int iterations     = 1; /* number of repetitions for each curve, set by -m */

  int b_high_set     = 0;

  unsigned int flags = 0;

  data *info = malloc(sizeof(data)); /* allocate an instance of data */

  /* initialize multi-precision integer variables */
  for(n = 0; n <= degree; n++) { mpz_init(c[n]); }

/**************************************************************************
 * get at the input                                                       *
 **************************************************************************/

  /* recognise optional args */
  { long i = 1;
    while(i < argc)
    {
      if(*(argv[i]) != '-') return(-6);
      switch(argv[i][1])
      { case 'h': /* height bound */
          if(argc == i) return(-6);
          i++;
          if(sscanf(argv[i], " %ld", &height) != 1) return(-6);
          i++;
          break;
        case 'p': /* max number of primes used */
          if(argc == i) return(-6);
          i++;
          if(sscanf(argv[i], " %ld", &num_primes) != 1) return(-6);
          i++;
          break;
        case 'F': /* max number of "forbidden divisors of denominator" */
          if(argc == i) return(-6);
          i++;
          if(sscanf(argv[i], " %ld", &max_forbidden) != 1) return(-6);
          i++;
          break;
        case 'n': /* number of primes used for first stage of sieving */
          if(argc == i) return(-6);
          i++;
          if(sscanf(argv[i], " %ld", &sieve_primes1) != 1) return(-6);
          i++;
          break;
        case 'N': /* number of primes used for sieving altogether */
          if(argc == i) return(-6);
          i++;
          if(sscanf(argv[i], " %ld", &sieve_primes2) != 1) return(-6);
          i++;
          break;
        case 'j': /* do not use Jacobi sum test */
          no_jacobi = 1;
          i++;
          break;
        case 'J': /* do use Jacobi sum test */
          no_jacobi = 0;
          i++;
          break;
        case 'k': /* keep: do not reverse polynomial */
          no_reverse = 1;
          i++;
          break;
        case 'K': /* allow reversal of polynomial */
          no_reverse = 0;
          i++;
          break;
        case 'x': /* no check */
          no_check = 1;
          i++;
          break;
        case 'X': /* do check points */
          no_check = 0;
          i++;
          break;
        case 'y': /* no y */
          no_y = 1;
          i++;
          break;
        case 'Y': /* print complete points */
          no_y = 0;
          i++;
          break;
        case 'z': /* no output */
          no_output = 1;
          i++;
          break;
        case 'Z': /* do print points */
          no_output = 0;
          i++;
          break;
        case 's': /* no Sturm sequence computation */
          sturm_iter = -1;
          i++;
          break;
        case 'S': /* Sturm sequence */
          i++;
          if(i <= argc && argv[i][0] != '-')
          { if(sscanf(argv[i], " %ld", &sturm_iter) != 1) return(-6);
            i++;
          }
          else sturm_iter = RATPOINTS_DEFAULT_STURM;
          break;
        case 'd': /* Bounds for denom */
          switch(argv[i][2])
          { case 'l': /* lower bound */
              if(argc == i) return(-6);
              i++;
              if(sscanf(argv[i], " %ld", &b_low) != 1) return(-6);
              i++;
              break;
            case 'u': /* upper bound */
              if(argc == i) return(-6);
              i++;
              if(sscanf(argv[i], " %ld", &b_high) != 1) return(-6);
              i++;
              b_high_set = 1; /* note that upper bound was set */
              break;
            default: return(-6);
          }
          break;
        case 'm': /* Number of iterations */
          if(argc == i) return(-6);
          i++;
          if(sscanf(argv[i], " %d", &iterations) != 1) return(-6);
          i++;
          break;
        default: return(-6);
  } } }

  /* In case height was changed: */
  if(!b_high_set) { b_high = height; }

  /* initialize */
  args.degree = degree; /* this information is needed for the initialization */
  find_points_init(&args);

  if(no_check)   { flags |= RATPOINTS_NO_CHECK;   }
  if(no_y)       { flags |= RATPOINTS_NO_Y;       }
  if(no_reverse) { flags |= RATPOINTS_NO_REVERSE; }
  if(no_jacobi)  { flags |= RATPOINTS_NO_JACOBI;  }

  /* Repeat computation iterations times */
  { long count;

    for(count = iterations; count; count--)
    { for(n = 0; n < NUM_TEST; n++)
      { /* set up polynomial */
        long k;

        for(k = 0; k < 7; k++) { mpz_set_si(c[k], testdata[n][k]); }

        /* Fill args structure. Recall:
          typedef struct {mpz_t *cof; long degree; long height;
                          ratpoints_interval *domain; long num_inter;
                          long b_low; long b_high; long sp1; long sp2;
                          double ratio1; double ratio2; long array_size;
                          long sturm; long num_primes; long max_forbidden;
                          unsigned int flags; ...}
            ratpoints_args; */

        args.cof           = &c[0];
        args.degree        = degree;
        args.height        = height;
        args.domain        = &domain[0];
        args.num_inter     = 0;
        args.b_low         = b_low;
        args.b_high        = b_high;
        args.sp1           = sieve_primes1;
        args.sp2           = sieve_primes2;
        args.array_size    = array_size;
        args.sturm         = sturm_iter;
        args.num_primes    = num_primes;
        args.max_forbidden = max_forbidden;
        args.flags         = flags;

        /* Fill info. Recall:
          typedef struct {int print_between; int no_output; int one_point;
                          int no_y;} data; */

        info->print_between = 0;
        info->no_output     = no_output;
        info->no_y          = 0;

        if(no_output == 0) { printf("{"); }
        find_points_work(&args, process, (void *)info);
        if(no_output == 0) { printf("}\n"); }
        /* fflush(NULL); */
      }
    }
  }

  /* clean up */
  find_points_clear(&args);
  free(info);
  for(n = 0; n <= degree; n++) { mpz_clear(c[n]); }

  return(0);
}
