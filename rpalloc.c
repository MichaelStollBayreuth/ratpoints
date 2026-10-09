/***********************************************************************
 * ratpoints-3.1.0                                                     *
 *  - A program to find rational points on hyperelliptic curves        *
 * Copyright (C) 2026  Michael Stoll                                   *
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
 * rpalloc.c                                                           *
 *                                                                     *
 * What the library and the program do when memory runs out.  Linked   *
 * with -Wl,--wrap=malloc,--wrap=calloc,--wrap=realloc (an option of   *
 * GNU ld), every call of malloc, calloc and realloc made by the       *
 * library, which is linked statically, or by main.c goes to the       *
 * functions below, which make the n-th of them fail, or every one     *
 * from the n-th on.  The memory gmp takes for its integers is         *
 * allocated inside the gmp library and is not affected; by default    *
 * gmp ends the program when that runs out.                            *
 *                                                                     *
 * Built as rpalloc, this runs a sequence of searches through the      *
 * library, for n = 1, 2, 3, ... until the sequence runs through       *
 * without an allocation failing, each n in a process of its own, and  *
 * judges every search against the same search without failures: it    *
 * must deliver the same points in the same order and return the same  *
 * value, or return RATPOINTS_NO_MEMORY having delivered the first so  *
 * many of those points -- only when an allocation failed during it,   *
 * and not once the callback has stopped it.  After a failed           *
 * find_points_init, every search on that structure must return        *
 * RATPOINTS_NO_MEMORY without delivering a point.  A crash, or a      *
 * search that does not end within a minute, fails too.                *
 *                                                                     *
 * Compiled with -DRPALLOC_PROGRAM and linked with main.c, as          *
 * ratpoints-alloc, it is the program with the allocation that         *
 * RPALLOC_FAIL=n names failing, or every one from the one that        *
 * RPALLOC_FROM=n names; it says on standard error, when it exits, how *
 * many allocations it made and how many of them failed.               *
 *                                                                     *
 * "make testalloc" runs both (testalloc.sh).                          *
 *                                                                     *
 * Michael Stoll, October 9, 2026                                      *
 ***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <signal.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include "ratpoints.h"

void *__real_malloc(size_t size);
void *__real_calloc(size_t n, size_t size);
void *__real_realloc(void *p, size_t size);
void *__wrap_malloc(size_t size);
void *__wrap_calloc(size_t n, size_t size);
void *__wrap_realloc(void *p, size_t size);

/* The allocation that fails is the fail_at-th, counted from 1, and with
 * fail_on set every later one as well; none when fail_at is 0.  The
 * sieving threads allocate too, so the counts are kept atomically. */
static long fail_at = 0;
static int fail_on = 0;
static long n_alloc = 0, n_failed = 0;

static int fail_this(void)
{ long k = __atomic_add_fetch(&n_alloc, 1, __ATOMIC_SEQ_CST);

  if(fail_at > 0 && (k == fail_at || (fail_on && k > fail_at)))
  { __atomic_add_fetch(&n_failed, 1, __ATOMIC_SEQ_CST);
    return(1);
  }
  return(0);
}

void *__wrap_malloc(size_t size)
{ return(fail_this() ? NULL : __real_malloc(size)); }

void *__wrap_calloc(size_t n, size_t size)
{ return(fail_this() ? NULL : __real_calloc(n, size)); }

void *__wrap_realloc(void *p, size_t size)
{ return(fail_this() ? NULL : __real_realloc(p, size)); }

#ifdef RPALLOC_PROGRAM

/* The program: which allocation fails comes from the environment, and the
 * counts go to standard error when the program exits, through error() or
 * from main, so that the caller can tell when n has passed the last
 * allocation.  A run that has not ended after a minute is ended by
 * SIGALRM, so that a search that hangs shows as a failure. */
static void rpalloc_report(void)
{ fprintf(stderr, "rpalloc: %ld allocations, %ld failed\n",
          __atomic_load_n(&n_alloc, __ATOMIC_SEQ_CST),
          __atomic_load_n(&n_failed, __ATOMIC_SEQ_CST));
}

__attribute__((constructor)) static void rpalloc_start(void)
{ const char *s;

  if((s = getenv("RPALLOC_FAIL")) != NULL && *s) { fail_at = atol(s); }
  if((s = getenv("RPALLOC_FROM")) != NULL && *s)
  { fail_at = atol(s); fail_on = 1; }
  atexit(rpalloc_report);
  alarm(60);
}

#else /* the driver of the library */

static mpz_t c[RATPOINTS_MAX_DEGREE+1];
static ratpoints_interval domain[2*RATPOINTS_MAX_DEGREE];

/* A point as the search delivers it: a, b, and y by |y| modulo the prime
 * 2^61 - 1 and the sign of y (no |y| here is that large, so the two are y
 * itself). */
typedef struct { long a, b; unsigned long y; int s; } rp_pt;

#define MAXPTS 1500
static rp_pt got[MAXPTS];      /* the points of the search under way */
static long ngot = 0;          /* how many it has delivered */
static long stop_at = 0;       /* stop it at this point; 0: never */

static int record(long a, long b, const mpz_t y, void *info, int *quit)
{ if(ngot < MAXPTS)
  { got[ngot].a = a; got[ngot].b = b;
    got[ngot].y = mpz_tdiv_ui(y, 2305843009213693951UL);
    got[ngot].s = mpz_sgn(y);
  }
  ngot++;
  if(stop_at > 0 && ngot >= stop_at) { *quit = 1; }
  return(1);
}

static long failures(void)
{ return(__atomic_load_n(&n_failed, __ATOMIC_SEQ_CST)); }

/* The searches.  The first NSHARED run one after the other on one
 * structure, initialised for degree 6, and with find_points_init they make
 * every allocation there is in the library: the arrays of every search; on
 * several threads the pool of threads and the arrays in which the threads
 * keep the points of their blocks, which the first search grows from
 * nothing and stops, so that a block past the stop may be the one that
 * runs out; the block of sieve tables, enlarged for a curve of many points
 * that looks at more primes than the default (the second curve of
 * testdata-high-many.h), for the tables of composite moduli and for every
 * prime there is; the patterns of the forbidden divisors beyond the table
 * of primes (denominators from 251^2 on).  After a failure the structure
 * must serve the next search; the eighth repeats the second.  The ninth
 * goes through find_points, which initialises and clears a structure of
 * its own.  The numbers of points of the searches that are not stopped
 * were counted by brute force, over every coprime pair (a, b) in the
 * range. */
#define NSEARCH 9
#define NSHARED 8
static const char *const sname[NSEARCH] =
{ "y^2 = x^2 + 1, height 1000, stopped at the 700th point",
  "y^2 = x^2 + 1, height 1000",
  "a curve of many points, height 16383",
  "y^2 = 441x^4 + 126x^2 + 1, height 500, composite moduli",
  "y^2 = (x^2 + 1)(x^4 + 1), height 1000, every prime",
  "y^2 = 7x^2 + 3x + 5, denominators 66000 to 66100",
  "y^2 = x^5 + 1, height 500, stopped at the first point",
  "y^2 = x^2 + 1, height 1000, again",
  "y^2 = (x^2 + 1)(x^4 + 1), height 500, find_points" };
static const long expect[NSEARCH] = {700, 1436, 330, 16, 8, 2, 1, 1436, 8};

/* the reference: what each search delivers when nothing fails */
static long ref_ret[NSEARCH], ref_n[NSEARCH];
static rp_pt ref[NSEARCH][MAXPTS];

/* set the coefficients from a list of longs and the input fields to what
 * main.c starts from */
static void curve(ratpoints_args *args, long degree, const long *cof,
                  long height)
{ long k;

  for(k = 0; k <= degree; k++) { mpz_set_si(c[k], cof[k]); }
  args->cof = c; args->degree = degree; args->height = height;
  args->domain = domain; args->num_inter = 0;
  args->b_low = 1; args->b_high = height;
  args->sp1 = -1; args->sp2 = -1;
  args->survivors_per_word = -1.0; args->sp2_extra = -1; args->sp2_u0 = -1.0;
  args->cost_table = -1.0; args->adapt = -1; args->sp3_extra = -1;
  args->sp3_per_denom = -1.0; args->check_cost = -1.0;
  args->array_size = RATPOINTS_ARRAY_SIZE; args->sturm = RATPOINTS_DEFAULT_STURM;
  args->num_primes = -1; args->max_forbidden = -1;
  args->num_threads = 1;
  args->flags = 0;
}

static void setup(ratpoints_args *args, long s, long nt)
{ static const long pyth[] = {1, 0, 1};
  static const long rich[] = {456976, -448032, -255200, 208380, 61033, -12834, 81};
  static const long quartic[] = {1, 0, 126, 0, 441};
  static const long sextic[] = {1, 0, 1, 0, 1, 0, 1};
  static const long quad[] = {5, 3, 7};
  static const long quintic[] = {1, 0, 0, 0, 0, 1};

  stop_at = 0;
  switch(s)
  { case 0: curve(args, 2, pyth, 1000); stop_at = 700; break;
    case 1: case 7: curve(args, 2, pyth, 1000); break;
    case 2: curve(args, 6, rich, 16383); break;
    case 3: curve(args, 4, quartic, 500);
            args->num_primes = 40; args->sp1 = 25; args->sp2 = 35;
            args->cost_table = 0.0;
            break;
    case 4: curve(args, 6, sextic, 1000); args->num_primes = 1000; break;
    case 5: curve(args, 2, quad, 70000);
            args->b_low = 66000; args->b_high = 66100; args->num_inter = 1;
            domain[0].low = 0.0; domain[0].up = 0.01;
            break;
    case 6: curve(args, 5, quintic, 500); stop_at = 1; break;
    case 8: curve(args, 6, sextic, 500); break;
  }
  args->num_threads = nt;
}

static int same_points(const rp_pt *p, const rp_pt *q, long n)
{ long k;

  for(k = 0; k < n; k++)
  { if(p[k].a != q[k].a || p[k].b != q[k].b || p[k].y != q[k].y
       || p[k].s != q[k].s)
    { return(0); }
  }
  return(1);
}

/* Whether search s did what it must; says what is wrong if not.  A search
 * that the callback has stopped must not fail: the point that could not be
 * kept, if any, comes after the stop. */
static int judge(long s, long ret, int failed, int init_failed)
{ long n = ngot;
  int prefix = (n <= ref_n[s] && same_points(got, ref[s], n));
  int ok;

  if(init_failed) { ok = (ret == RATPOINTS_NO_MEMORY && n == 0); }
  else if(ret == ref_ret[s] && n == ref_n[s] && prefix) { ok = 1; }
  else
  { ok = (failed && ret == RATPOINTS_NO_MEMORY && prefix
            && (stop_at == 0 || n < stop_at));
  }
  if(!ok)
  { printf("  n = %ld, %s: returned %ld after %ld point(s)%s, %s; without"
           " failures it returns %ld after %ld\n", fail_at, sname[s], ret, n,
           prefix ? "" : " that are not the first of those without failures",
           init_failed ? "find_points_init having failed"
                       : (failed ? "an allocation having failed in it"
                                 : "no allocation having failed in it"),
           ref_ret[s], ref_n[s]);
  }
  return(ok);
}

/* Run the sequence on nt threads.  With keep set, the points are kept as
 * the reference; otherwise each search is judged against it.  Returns the
 * number of searches judged wrong. */
static long sequence(long nt, int keep)
{ static const long sextic[] = {1, 0, 1, 0, 1, 0, 1};
  ratpoints_args args;
  long s, f, wrong = 0;
  int init_failed;

  curve(&args, 6, sextic, 1000);
  f = failures();
  find_points_init(&args);
  init_failed = (failures() > f);
  for(s = 0; s < NSEARCH; s++)
  { long ret;

    setup(&args, s, nt);
    ngot = 0;
    f = failures();
    ret = (s < NSHARED) ? find_points_work(&args, record, NULL)
                        : find_points(&args, record, NULL);
    if(keep)
    { ref_ret[s] = ret; ref_n[s] = ngot;
      memcpy(ref[s], got, ((ngot < MAXPTS) ? ngot : MAXPTS)*sizeof(rp_pt));
    }
    else if(!judge(s, ret, failures() > f, init_failed && s < NSHARED))
    { wrong++; }
    if(s == NSHARED - 1) { find_points_clear(&args); }
  }
  return(wrong);
}

int main(int argc, char *argv[])
{ long nt, n, first = 1, last = 0, total, bad = 0, k;
  int from;

  if(argc < 3 || (strcmp(argv[2], "one") != 0 && strcmp(argv[2], "from") != 0))
  { printf("usage: rpalloc threads one|from [n]\n"
           "  fails the n-th allocation (one), or every one from the n-th on"
           " (from),\n  for n = 1, 2, 3, ... or for the n given\n");
    return(2);
  }
  nt = atol(argv[1]);
  from = (strcmp(argv[2], "from") == 0);
  if(argc > 3) { first = last = atol(argv[3]); }
  for(k = 0; k <= RATPOINTS_MAX_DEGREE; k++) { mpz_init(c[k]); }

  /* the reference, and its numbers of points against the brute-force ones */
  sequence(nt, 1);
  total = n_alloc;
  for(k = 0; k < NSEARCH; k++)
  { if(ref_n[k] != expect[k] || ref_ret[k] != expect[k] || ref_n[k] > MAXPTS)
    { printf("without failures, %s: returned %ld after %ld point(s), expected %ld\n",
             sname[k], ref_ret[k], ref_n[k], expect[k]);
      bad++;
    }
  }
  if(bad) { return(1); }

  /* each n in a process of its own, which inherits the reference and has
   * no thread besides its own: the pools of the reference are gone.  It
   * exits with bit 0 set when a search was wrong and bit 1 set when no
   * allocation failed, which means that n has passed the last one. */
  for(n = first; last == 0 || n <= last; n++)
  { pid_t pid;
    int st;

    fflush(stdout);
    pid = fork();
    if(pid < 0) { perror("rpalloc: fork"); return(2); }
    if(pid == 0)
    { long wrong;

      alarm(60);
      fail_at = n; fail_on = from; n_alloc = 0; n_failed = 0;
      wrong = sequence(nt, 0);
      fail_at = 0; /* the sequence is over: what exit does gets its memory */
      exit((wrong ? 1 : 0) | (n_failed > 0 ? 0 : 2));
    }
    while(waitpid(pid, &st, 0) < 0)
    { if(errno != EINTR) { perror("rpalloc: waitpid"); return(2); } }
    if(WIFSIGNALED(st))
    { printf("  n = %ld: ended by signal %d%s\n", n, WTERMSIG(st),
             (WTERMSIG(st) == SIGALRM) ? ", not having ended within a minute" : "");
      bad++;
    }
    else
    { if(WEXITSTATUS(st) & 1) { bad++; }
      if(WEXITSTATUS(st) & 2) { break; }
    }
    if(n > 10*total + 1000)
    { printf("  n = %ld: the allocations do not end\n", n); bad++; break; }
  }
  printf("%ld thread%s, %s: n = %ld to %ld, %s\n", nt, (nt == 1) ? "" : "s",
         from ? "every allocation from the n-th on failing"
              : "the n-th allocation failing",
         first, (last > 0) ? last : n - 1,
         bad ? "WRONG" : "every search right");

  for(k = 0; k <= RATPOINTS_MAX_DEGREE; k++) { mpz_clear(c[k]); }
  return(bad ? 1 : 0);
}

#endif
