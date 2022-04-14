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
 * sturm.c                                                             *
 *                                                                     *
 * Sturm sequence and positivity intervals                             *
 *                                                                     *
 * Michael Stoll, Jan 9, 2008                                          *
 ***********************************************************************/

#include "ratpoints.h"

/**************************************************************************
 * Arguments of _ratpoints_compute_sturm() : (from the args argument)     *
 *                                                                        *
 * + cofs   - points to an array of mpz_t's holding the coefficients of   *
 *            the polynomial                                              *
 * + degree - the degree of the polynomial                                *
 * + iter   - the number of iteration steps in the refinement of the      *
 *            intervals                                                   *
 * + ivlist - points to an array of intervals giving the current search   *
 *            domain                                                      *
 * + num_iv - the number of intervals in ivlist, must be > 0              *
 *                                                                        *
 * NOTE: ivlist must be able to store  1 + floor(degree/2)  additional    *
 *       intervals.                                                       *
 *                                                                        *
 * Return values :                                                        *
 * + >0 - normal operation, ivlist modified (intersection with positivity *
 *        domain), return value is number of intervals                    *
 * +  0 - polynomial is everywhere negative, ivlist may be changed        *
 * + -1 - polynomial is not squarefree, ivlist unchanged                  *
 **************************************************************************/

/* A helper function: evaluate the polynomial in cofs[] of given degree
   at num/2^denexp and return the sign. */

static long eval_sign(ratpoints_args *args, mpz_t *cofs, long degree,
                      long num, long denexp)
{
  long n, e, s;
  mpz_t *work = args->work;

  /* Horner scheme... */
  mpz_set(work[0], cofs[degree]);
  for(n = degree-1, e = denexp; n >= 0; n--, e += denexp)
  { mpz_mul_si(work[0], work[0], num);
    mpz_mul_2exp(work[1], cofs[n], e);
    mpz_add(work[0], work[0], work[1]);
  }
  s = mpz_cmp_si(work[0], 0);
  return(s);
}

long _ratpoints_compute_sturm(ratpoints_args *args)
{
  mpz_t *cofs = args->cof;
  long degree = args->degree;
  long iter = args->sturm;
  ratpoints_interval *ivlist = args->domain;
  long num_iv = args->num_inter;
  long n, m, k, new_num;
  mpz_t sturm[degree+1][degree+1]; /* Array to hold the polynomials */
  long sturm_degs[degree+1]; /* The degrees of the polynomials */
  mpz_t *work = args->work;
  long count1 = 0, count2 = 0;

  /* first initialize */
  for(n = 0; n <= degree; n++)
  { for(m = 0; m <= degree; m++)
    { mpz_init(sturm[n][m]); }
  }
  /* copy polynomial f into first entry */
  for(n = 0; n <= degree; n++) { mpz_set(sturm[0][n], cofs[n]); }
  sturm_degs[0] = degree;
  /* compute derivative in second entry */
  for(n = 0; n < degree; n++)
  { mpz_set(work[2], cofs[n+1]);
    mpz_mul_si(sturm[1][n], work[2], n+1);
  }
  sturm_degs[1] = degree - 1;

  /* now do polynomial divisions ... */
  for(k = 2; k <= degree; k++)
  { long d1 = sturm_degs[k-1], d2 = sturm_degs[k-2];
    /* first copy sturm[k-2] into sturm[k] */
    for(n = 0; n <= degree - (k-2); n++) { mpz_set(sturm[k][n], sturm[k-2][n]); }
    /* now build linear combination that reduces the degree */
    while(d2 >= d1)
    { mpz_gcd(work[2], sturm[k-1][d1], sturm[k][d2]);
      mpz_fdiv_q(work[0], sturm[k-1][d1], work[2]);
      mpz_fdiv_q(work[1], sturm[k][d2], work[2]);
      if(mpz_cmp_si(work[0], 0) < 0)
      { mpz_neg(work[0], work[0]); mpz_neg(work[1], work[1]); }
      /* sturm[k] = work[0] * sturm[k] - work[1] * x^(d2-d1) * sturm[k-1] */
      for(n = 0; n <= d1; n++)
      { mpz_mul(sturm[k][n+d2-d1], sturm[k][n+d2-d1], work[0]);
        mpz_submul(sturm[k][n+d2-d1], work[1], sturm[k-1][n]);
      }
      for(n = 0; n < d2-d1; n++)
      { mpz_mul(sturm[k][n], sturm[k][n], work[0]); }
      d2--;
      while(mpz_cmp_si(sturm[k][d2], 0) == 0 && d2 >= 0) { d2--; }
      if(d2 < 0)  /* not squarefree */
      { for(n = 0; n <= degree; n++)
        { for(m = 0; m <= degree; m++)
            { mpz_clear(sturm[n][m]); }
        }
        return(-1);
      }
    }
    /* change sign */
    for(n = 0; n <= d2; n++) { mpz_neg(sturm[k][n], sturm[k][n]); }
    /* normalize */
    mpz_set_ui(work[2], 0);
    for(n = 0; n <= d2; n++)
    { mpz_gcd(work[2], work[2], sturm[k][n]);
      if(mpz_cmp_ui(work[2], 1) == 0) { break; }
    }
    if(mpz_cmp_ui(work[2], 1) != 0)
    { for(n = 0; n <= d2; n++)
      { mpz_fdiv_q(sturm[k][n], sturm[k][n], work[2]); }
    }
    sturm_degs[k] = d2;
    if(d2 == 0) { break; } /* sturm[k] is constant */
  }

  /* compute number of real zeros */
  for(n = 0; n < k; n++)
  { long d1 = sturm_degs[n], d2 = sturm_degs[n+1];
    int s1 = mpz_cmp_si(sturm[n][d1], 0),
        s2 = mpz_cmp_si(sturm[n+1][d2], 0);

    if(s1 != s2) { count1++; }
    if(d1 & 1) { s1 = -s1; }
    if(d2 & 1) { s2 = -s2; }
    if(s1 != s2) { count2++; }
  }

  if(count2 == count1 && mpz_cmp_si(cofs[0], 0) < 0)
  { /* no real roots, negative constant term ==> no points */
    for(n = 0; n <= degree; n++)
    { for(m = 0; m <= degree; m++)
        { mpz_clear(sturm[n][m]); }
    }
    args->num_inter = 0;
    return(0);
  }

  /* Find list of intervals that may contain points */
  /* recall: typedef struct {double low; double up;} ratpoints_interval; */
  { ratpoints_interval ivlocal[1 + (degree>>1)];
    ratpoints_interval *iptr = &ivlocal[0];
    long max = (long)(((unsigned long)(-1))>>1);
    long min = -max;
    long num_intervals;
    long slcf = mpz_cmp_si(cofs[degree], 0);

    /* recursive helper function */
    void iterate(long nl, long nr, long del, long der, long cleft, long cright,
                 long sl, long sr, long depth)
    { /* nl/2^del, nr/2^der : interval left/right endpoints,
         cleft, cright: sign change counts at endpoints,
         sl, sr: signs at endpoints,
         depth: iteration depth */
      if(cleft == cright && sl < 0) { return; }
      /* here we know the polynomial is negative on the interval */
      if((cleft == cright && sl > 0) || depth >= iter)
      /* we have to add/extend an interval if we either know that
         the polynomial is positive on the interval (first condition)
         or the maximal iteration depth has been reached (second condition) */
      { double l = ((double)nl)/((double)(1<<del));
        double u = ((double)nr)/((double)(1<<der));

        if(iptr == &ivlocal[0])
        { iptr->low = l; iptr->up  = u; iptr++; }
        else
        { if((iptr-1)->up == l) /* extend interval */
          { (iptr-1)->up = u; }
          else /* new interval */
          { iptr->low = l; iptr->up  = u; iptr++; }
        }
        return;
      }
      /* now we must split the interval and evaluate the sturm sequence
         at the midpoint */
      { long nm, dem, s0, s1, s2, s, cmid = 0, n;

        if(nl == min)
        { if(nr == max) { nm = 0; dem = 0; }
          else { nm = (nr == 0) ? -1 : 2*nr; dem = 0; }
        }
        else
        { if(nr == max) { nm = (nl == 0) ? 1 : 2*nl; dem = 0; }
          else /* "normal" case */
          { if(del == der) /* then both are zero */
            { if(((nl+nr) & 1) == 0) { nm = (nl+nr)>>1; dem = 0; }
              else { nm = nl+nr; dem = 1; }
            }
            else /* here one de* is greater */
            { if(del > der) { nm = nl + (nr<<(del-der)); dem = del+1; }
              else { nm = (nl<<(der-del)) + nr; dem = der+1; }
            }
          }
        }
        s0 = eval_sign(args, sturm[0], sturm_degs[0], nm, dem);
        s1 = eval_sign(args, sturm[1], sturm_degs[1], nm, dem);
        if(s0*s1 == -1) { cmid++; }
        s = (s1 == 0) ? s0 : s1;
        for(n = 2; n <= k; n++)
        { s2 = eval_sign(args, sturm[n], sturm_degs[n], nm, dem);
          if(s2 == -s) { cmid++; s = s2; }
          else if(s2 != 0) { s = s2; }
        }
        /* now recurse */
        iterate(nl, nm, del, dem, cleft, (s0==0) ? (cmid+1) : cmid,
                sl, (s0==0) ? -s1 : s0, depth+1);
        iterate(nm, nr, dem, der, cmid, cright,
                (s0==0) ? s1 : s0, sr, depth+1);
      }
    } /* end iterate() */

    iterate(min, max, 0, 0, count2, count1,
            (degree & 1) ? -slcf : slcf, slcf, 0);
    num_intervals = iptr - &ivlocal[0];
    /* intersect with given intervals */
    { ratpoints_interval local_copy[num_iv];
      long n, n1, n2;

      /* make a copy of the given list */
      for(n = 0; n < num_iv; n++) { local_copy[n] = ivlist[n]; }
      n1 = 0; n2 = 0; n = 0;
      while(n1 < num_intervals && n2 < num_iv)
      { if(ivlocal[n1].low <= local_copy[n2].low)
        { if(ivlocal[n1].up < local_copy[n2].low)
          { n1++; } /* can forget this interval */
          else
          { if(ivlocal[n1].up <= local_copy[n2].up)
            { /* note intersection */
              ivlist[n].low = local_copy[n2].low;
              ivlist[n].up = ivlocal[n1].up;
              n++;
              n1++;
            }
            else
            { /* note intersection */
              ivlist[n] = local_copy[n2];
              n++;
              n2++;
        } } }
        else /* here, ivlocal[n1].low > local_copy[n2].low */
        { if(local_copy[n2].up < ivlocal[n1].low)
          { n2++; } /* can forget this interval */
          else
          { if(local_copy[n2].up <= ivlocal[n1].up)
            { /* note intersection */
              ivlist[n].low = ivlocal[n1].low;
              ivlist[n].up = local_copy[n2].up;
              n++;
              n2++;
            }
            else
            { /* note intersection */
              ivlist[n] = ivlocal[n1];
              n++;
              n1++;
      } } } }
      args->num_inter = new_num = n;
    }
  }
  for(n = 0; n <= degree; n++)
  { for(m = 0; m <= degree; m++)
      { mpz_clear(sturm[n][m]); }
  }
  return(new_num);
}

