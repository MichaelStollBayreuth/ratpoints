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
 * main.c                                                              *
 *                                                                     *
 * Main program file for the ratpoints executable                      *
 *                                                                     *
 * Michael Stoll, May 27, 2009; Jan 7-18, 2022; Sep 6-21, 2026         *
 ***********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>

#include "ratpoints.h"

/**************************************************************************
 * define                                                                 *
 **************************************************************************/

#define RATPOINTS_VERSION \
  "This is ratpoints-3.1.0 Copyright (C) 2008,2009,2022,2023,2026 by Michael Stoll.\n\n" \
  "This program comes with ABSOLUTELY NO WARRANTY.\n" \
  "This is free software, and you are welcome to redistribute it under the\n" \
  "terms of the GNU General Public License version 2 or later.\n\n" \
  "Please acknowledge use of this program in published work.\n\n"

/**************************************************************************
 * global variables                                                       *
 **************************************************************************/

mpz_t c[RATPOINTS_MAX_DEGREE+1];  /* The coefficients of f */
mpz_t tmp, tmp2;      /* for intermediate results with gmp integers */

int quiet;            /* A flag saying whether to suppress messages */
int one_point;        /* A flag saying if one point is enough */
int points_at_infty;  /* A flag saying if we should look for points at infty */
int no_output;        /* A flag that indicates that no points should be printed */
char *print_format;   /* The printf format for printing points */
char *string_before;  /* String to be printed before the points */
char *string_between; /* String to be printed between the points */
char *string_after;   /* String to be printed after the points */

ratpoints_interval domain[2*RATPOINTS_MAX_DEGREE];
                      /* This contains the intervals representing the search region */

char *usage_str =
    "Usage: ratpoints 'a_0 a_1 ... a_n' max_height\n"
    "                 [-dl low_den] [-du up_den]\n"
    "                 [-f format] [-fs str] [-fm str] [-fe str] [-y] [-Y]\n"
    "                 [[-l low1] -u up1 ... -l lown [-u upn]]\n"
    "                 [-n num_primes1] [-N num_primes2] [-p max_primes]\n"
    "                 [-r survivors_per_word] [-R extra_primes] [-U words]\n"
    "                 [-C table_cost] [-A adapt]\n"
    "                 [-P stage3_primes] [-Q stage3_cost] [-W check_cost]\n"
    "                 [-F max_forbidden] [-t threads] [-s] [-S [iter]]\n"
    "                 [-q] [-v] [-z] [-Z] [-1] [-i] [-I]\n"
    "                 [-k] [-K] [-j] [-J] [-x] [-X]\n\n";

/**************************************************************************
 * prototypes                                                             *
 **************************************************************************/

int read_input(long, char *argv[], ratpoints_args*);
char *scan_mpz(char*, mpz_t);
int scan_long(char*, long*);
int scan_double(char*, double*);
void print_poly(mpz_t*, long);
void print_string(char*);
void message(long n, long total, ratpoints_args *args);
void error(long);

/**************************************************************************
 * function that processes the points                                     *
 **************************************************************************/

/* structure containing data to be passed to the point-processing function */
typedef struct {int print_between; int no_output; int one_point;
                int points_at_infty; int no_y; char *pf;
                char* string_before; char *string_between;} data;

int process(long a, long b, const mpz_t y, void *info0, int *quit)
{ data *info = (data *)info0;
  char *fmt = info->pf;
  if(b == 0 && !(info->points_at_infty)) return(0);
  *quit = info->one_point;
  if(info->no_output) return(1);
  if(info->print_between) { print_string(info->string_between); }
  else { print_string(info->string_before); info->print_between = 1; }
  while(*fmt)
  { char c = *fmt++;
    switch(c)
    { case '%':
        c = *fmt++;
        switch(c)
        { case 0: putchar('%'); return(1);
          case 'x': printf("%ld", a); break;
          case 'y': if(info->no_y)
                    { putchar('%'); putchar('y'); }
                    else
                    { mpz_out_str((FILE *)NULL, 10, y); }
                    break;
          case 'z': printf("%ld", b); break;
          default: putchar('%'); putchar(c); break;
        }
        break;
      case '\\':
        c = *fmt++;
        switch(c)
        { case 0: putchar('\\'); return(1);
          case 't': printf("\t"); break;
          case 'n': printf("\n"); break;
          case '\\': putchar('\\'); break;
          case '%': putchar('%'); break;
          default: putchar('\\'); putchar(c); break;
        }
        break;
      default: putchar(c); break;
  } }
  fflush(stdout);
  return(1);
}

/**************************************************************************
 * main                                                                   *
 **************************************************************************/

int main(int argc, char *argv[])
{
  long total;
  int before_printed; /* whether process has printed string_before */
  ratpoints_args args;

  { long n;
    /* initialize multi-precision integer variables */
    for(n = 0; n <= RATPOINTS_MAX_DEGREE; n++) { mpz_init(c[n]); }
  }
  mpz_init(tmp); mpz_init(tmp2);;

  /* read input */
  if(argc < 3) { error(2); }
  if(read_input(argc-1, &argv[0], &args))
  { message(1, 0, &args);
    return(0);
  }
  if(!quiet)
  { message(0, 0, &args);
    message(5, args.degree, &args);
    if(!(args.flags & RATPOINTS_VERBOSE))
    { message(6, args.height, &args); printf("\n"); }
  }
  /* use default array size */
  args.array_size = RATPOINTS_ARRAY_SIZE;

  if(!quiet && args.num_inter > 0) message(3, 0, &args);

  /* typedef struct {int print_between; int no_output; int one_point;
                     int points_at_infty; int no_y; char *pf;
                     char* string_before; char *string_between;} data; */

  /* initialize data structure to be passed to process... */
  { data *info = malloc(sizeof(data));

    if(info == NULL) { error(10); }
    info->print_between   = 0;
    info->no_output       = no_output;
    info->pf              = print_format;
    info->one_point       = one_point;
    info->points_at_infty = points_at_infty;
    info->string_before   = string_before;
    info->string_between  = string_between;
    info->no_y            = args.flags & (RATPOINTS_NO_CHECK | RATPOINTS_NO_Y);

    /* ...and call find_points */
    total = find_points(&args, process, (void *)info);
    before_printed = info->print_between;
    /* clean up */
    free(info);
  }
  if(total < 0) /* error signaled by find_points? */
  { mpz_clear(tmp); mpz_clear(tmp2);
    if(total == RATPOINTS_NON_SQUAREFREE) error(8);
    if(total == RATPOINTS_NO_MEMORY) error(10);
    error(9);
  }
  /* process prints string_before before the first point it prints;
   * if it printed none (there was no point, or -z suppressed them),
   * print string_before now, in place of the points. */
  if(!before_printed) { print_string(string_before); }
  /* terminate the output of the points */
  print_string(string_after);
  if(!quiet)
  { if(!(args.flags & RATPOINTS_VERBOSE))
    { printf("\n\n");
      message(4, 0, &args);
      message(7, 0, &args);
      message(8, 0, &args);
    }
    message(2, total, &args);
  }

  /* clear multi-precision integer variables */
  mpz_clear(tmp); mpz_clear(tmp2);
  { long n;
    for(n = 0; n <= RATPOINTS_MAX_DEGREE; n++) { mpz_clear(c[n]); }
  }
  return(0);
}

/**************************************************************************
 * procedures                                                             *
 **************************************************************************/

/**************************************************************************
 * get at the input                                                       *
 **************************************************************************/

/* The end that stands for a missing first -l (negated) or last -u: the
 * height bound as a double, rounded up when it is not one, so that the
 * interval holds every x = a/b with |a| <= height (the library takes the
 * ends of the intervals exactly). */
static double height_end(long height)
{ double h = (double)height;

  if(h < 9223372036854775808.0 && (long)h < height)
  { h = nextafter(h, HUGE_VAL); }
  return(h);
}

int read_input(long argc, char *argv[], ratpoints_args *args)
{
  int l_seen = 0;     /* flag for dealing with -l -u */
  mpz_t fff;          /* temporary variable for gmp integers */
  long num_inter = 0; /* number of positivity intervals */

  /* Initialize args */
  /* typedef struct {mpz_t *cof; long degree; long height;
                    ratpoints_interval *domain; long num_inter;
                    long b_low; long b_high; long sp1; long sp2;
                    double ratio1; double ratio2; long array_size;
                    long sturm; long num_primes; long max_forbidden;
                    unsigned int flags; ...}
            ratpoints_args;
  */
  args->cof           = &c[0];
  args->degree        = 0;
  args->height        = 0;
  args->domain        = &domain[0];
  args->num_inter     = 0;  /* No interval up to now */
  args->b_low         = 1;  /* denominators go from 1 to h */
  args->b_high        = -1;
  args->sp1           = -1; /* negative: choose from the curve */
  args->sp2           = -1; /* negative: choose from the curve */
  args->survivors_per_word = -1.0; /* negative: compiled-in default */
  args->sp2_extra     = -1; /* negative: compiled-in default */
  args->sp2_u0        = -1.0; /* negative: compiled-in default */
  args->cost_table    = -1.0; /* negative: compiled-in default */
  args->adapt         = -1; /* negative: adapt when free to */
  args->sp3_extra     = -1; /* negative: choose from the curve */
  args->sp3_per_denom = -1.0; /* negative: compiled-in default */
  args->check_cost    = -1.0; /* negative: estimate it from the curve */
  args->array_size    = RATPOINTS_ARRAY_SIZE;    /* default */
  args->sturm         = RATPOINTS_DEFAULT_STURM; /* default */
  args->num_primes    = -1; /* gives default value */
  args->max_forbidden = -1; /* gives default value */
  args->num_threads   = 1;  /* sieve on this thread; see -t */
  args->flags         = 0;  /* do the check by default */
                            /* list y-coordinates by default */
                            /* allow reversal of polynomial */
                            /* use Jacobi symbol test */

  mpz_init(fff); /* initialize gmp integer */
  /* read coefficients: integers separated by white space (see scan_mpz),
   * at most RATPOINTS_MAX_DEGREE + 1 of them.  Those beyond go into fff,
   * so that a list with something in it that is not an integer is refused
   * as such, however long it is. */
  { char *s = argv[1];
    long n = 0; /* the number of integers read */

    for(;;)
    { while(isspace((unsigned char)*s)) { s++; }
      if(*s == 0) { break; }
      s = scan_mpz(s, (n <= RATPOINTS_MAX_DEGREE) ? c[n] : fff);
      if(s == NULL) { error(11); }
      n++;
    }
    if(n > RATPOINTS_MAX_DEGREE + 1) { error(3); }
    if(n <= 1) { error(5); }
    args->degree = n - 1;
  }

  /* read height bound */
  if(!scan_long(argv[2], &(args->height)) || args->height < 1)
  { error(4); }
  args->b_high = args->height;

  /* Set global variables to their default values */
  no_output       = 0;    /* print points by default */
  quiet           = 0;    /* don't be quiet */
  one_point       = 0;    /* look for all points */
  points_at_infty = 1;    /* also find points at infinity */
  print_format    = NULL;
  string_before   = NULL;
  string_between  = NULL;
  string_after    = NULL;

  /* recognize optional args */
  { long i = 3;
    while(i <= argc)
    {
      if(*(argv[i]) != '-') { error(6); }
      switch(argv[i][1])
      { case 'l': /* lower intevral endpoint */
          if(argc == i) error(7);
          i++;
          if(l_seen) { error(7); } /* -l -l */
          if(num_inter == RATPOINTS_MAX_DEGREE) { error(7); }
          if(!scan_double(argv[i], &domain[num_inter].low)) { error(7); }
          if(num_inter > 0 && domain[num_inter-1].up >= domain[num_inter].low)
          { error(7); }
          i++;
          l_seen = 1;
          break;
        case 'u': /* upper interval endpoint */
          if(argc == i) { error(7); }
          i++;
          if(!l_seen)
          { if(num_inter == 0) { domain[0].low = -height_end(args->height); }
            else { error(7); } /* -u -u */
          }
          if(!scan_double(argv[i], &domain[num_inter].up)) { error(7); }
          /* a first interval with no -l that ends below -height holds no
           * x = a/b of the search, and it is made [-infinity, u], which is
           * what a missing first -l stands for, rather than turned around,
           * which the library would refuse */
          if(!l_seen && domain[num_inter].up < domain[num_inter].low)
          { domain[num_inter].low = -HUGE_VAL; }
          if(domain[num_inter].low > domain[num_inter].up) { error(7); }
          i++;
          l_seen = 0;
          num_inter++;
          break;
        case 'p': /* max number of primes used */
          if(argc == i) { error(6); }
          i++;
          if(!scan_long(argv[i], &(args->num_primes))) { error(6); }
          i++;
          break;
        case 'F': /* max number of "forbidden divisors of denominator" */
          if(argc == i) { error(6); }
          i++;
          if(!scan_long(argv[i], &(args->max_forbidden))) { error(6); }
          i++;
          break;
        case 't': /* sieving threads; 0: as many as there are processors */
          if(argc == i) { error(6); }
          i++;
          if(!scan_long(argv[i], &(args->num_threads))) { error(6); }
          if(args->num_threads == 0) { args->num_threads = -1; }
          i++;
          break;
        case 'n': /* number of primes used for first stage of sieving */
          if(argc == i) { error(6); }
          i++;
          if(!scan_long(argv[i], &(args->sp1))) { error(6); }
          i++;
          break;
        case 'N': /* number of primes used for sieving altogether */
          if(argc == i) { error(6); }
          i++;
          if(!scan_long(argv[i], &(args->sp2))) { error(6); }
          i++;
          break;
        case 'r': /* the survivors of the first stage per 64-bit word swept
                   * that one more modulus has to remove, for a modulus that
                   * costs one AND per word; decides sp1 when -n is absent
                   * (ratpoints.h, RATPOINTS_SURVIVORS_PER_WORD) */
          if(argc == i) { error(6); }
          i++;
          if(!scan_double(argv[i], &(args->survivors_per_word)))
          { error(6); }
          i++;
          break;
        case 'R': /* how many primes to add to sp1 to get sp2, when -N is
                   * absent */
          if(argc == i) { error(6); }
          i++;
          if(!scan_long(argv[i], &(args->sp2_extra))) { error(6); }
          i++;
          break;
        case 'U': /* the number of numerator words at which a second-stage
                   * prime pays for setting itself up; scales -R down for
                   * short runs.  Zero leaves -R flat */
          if(argc == i) { error(6); }
          i++;
          if(!scan_double(argv[i], &(args->sp2_u0))) { error(6); }
          i++;
          break;
        case 'A': /* how much of the choice to correct during the run from
                   * what the sieve is doing: 0 none, 1 (the default) the
                   * third stage, 2 the second stage as well */
          if(argc == i) { error(6); }
          i++;
          if(!scan_long(argv[i], &(args->adapt))) { error(6); }
          i++;
          break;
        case 'C': /* what building one row of a sieve table costs, in units
                   * of one first-stage AND per word; zero ranks the primes
                   * by what they say alone */
          if(argc == i) { error(6); }
          i++;
          if(!scan_double(argv[i], &(args->cost_table))) { error(6); }
          i++;
          break;
        case 'P': /* how many primes the third stage adds to sp2; negative
                   * means choose it from the curve */
          if(argc == i) { error(6); }
          i++;
          if(!scan_long(argv[i], &(args->sp3_extra))) { error(6); }
          i++;
          break;
        case 'Q': /* what one third-stage prime costs per denominator, as a
                   * fraction of one exact check; decides -P when it is
                   * absent */
          if(argc == i) { error(6); }
          i++;
          if(!scan_double(argv[i], &(args->sp3_per_denom))) { error(6); }
          i++;
          break;
        case 'W': /* what one exact check costs, in rdtsc cycles; negative
                   * means estimate it from the degree and the coefficients,
                   * and 306 is the value for a degree-6 curve with small
                   * coefficients */
          if(argc == i) { error(6); }
          i++;
          if(!scan_double(argv[i], &(args->check_cost))) { error(6); }
          i++;
          break;
        case 'f': /* printing format */
          switch(argv[i][2])
          { case 0: /* just -f */
              if(argc == i) { error(6); }
              i++;
              { long l = strlen(argv[i]);
                print_format = malloc((l+1)*sizeof(char));
                if(print_format == NULL) { error(10); }
                strcpy(print_format, argv[i]);
              }
              i++;
              break;
            case 's': /* starting string */
              if(argc == i) { error(6); }
              i++;
              { long l = strlen(argv[i]);
                string_before = malloc((l+1)*sizeof(char));
                if(string_before == NULL) { error(10); }
                strcpy(string_before, argv[i]);
              }
              i++;
              break;
            case 'm': /* in-between string */
              if(argc == i) { error(6); }
              i++;
              { long l = strlen(argv[i]);
                string_between = malloc((l+1)*sizeof(char));
                if(string_between == NULL) { error(10); }
                strcpy(string_between, argv[i]);
              }
              i++;
              break;
            case 'e': /* ending string */
              if(argc == i) { error(6); }
              i++;
              { long l = strlen(argv[i]);
                string_after = malloc((l+1)*sizeof(char));
                if(string_after == NULL) { error(10); }
                strcpy(string_after, argv[i]);
              }
              i++;
              break;
            default: { error(6); }
          }
          break;
        case 'q': /* quiet */
          quiet = 1;
          i++;
          break;
        case 'v': /* verbose */
          args->flags |= RATPOINTS_VERBOSE;
          i++;
          break;
        case 'j': /* do not use Jacobi symbol */
          args->flags |= RATPOINTS_NO_JACOBI;
          i++;
          break;
        case 'J': /* do use Jacobi symbol */
          args->flags &= ~RATPOINTS_NO_JACOBI;
          i++;
          break;
        case 'k': /* keep: do not reverse polynomial */
          args->flags |= RATPOINTS_NO_REVERSE;
          i++;
          break;
        case 'K': /* allow reversal of polynomial */
          args->flags &= ~RATPOINTS_NO_REVERSE;
          i++;
          break;
        case 'x': /* no check */
          args->flags |= RATPOINTS_NO_CHECK;
          i++;
          break;
        case 'X': /* do check points */
          args->flags &= ~RATPOINTS_NO_CHECK;
          i++;
          break;
        case 'y': /* print only x-coordinates */
          args->flags |= RATPOINTS_NO_Y;
          i++;
          break;
        case 'Y': /* print all points */
          args->flags &= ~RATPOINTS_NO_Y;
          i++;
          break;
        case 'z': /* no output */
          no_output = 1;
          i++;
          break;
        case 'Z': /* output the points */
          no_output = 0;
          i++;
          break;
        case '1': /* only one point */
          one_point = 1;
          i++;
          break;
        case 'i': /* no points at infty */
          points_at_infty = 0;
          i++;
          break;
        case 'I': /* print points at infty */
          points_at_infty = 1;
          i++;
          break;
        case 's': /* don't use Sturm sequence computation */
          args->sturm = -1;
          i++;
          break;
        case 'S': /* Sturm sequence */
          args->sturm = RATPOINTS_DEFAULT_STURM;
          i++;
          if(i <= argc && argv[i][0] != '-')
          { if(!scan_long(argv[i], &(args->sturm))) { error(6); }
            i++;
          }
          break;
        case 'd': /* Bounds for denom */
          switch(argv[i][2])
          { case 'l': /* lower bound */
              if(argc == i) { error(6); }
              i++;
              if(!scan_long(argv[i], &(args->b_low))) { error(6); }
              i++;
              break;
            case 'u': /* upper bound */
              if(argc == i) { error(6); }
              i++;
              if(!scan_long(argv[i], &(args->b_high))) { error(6); }
              i++;
              break;
            default: { error(6); }
          }
          break;
        default: { error(6); }
  } } }
  if(l_seen)
  /* complete last interval; one that starts above the height bound holds no
   * x = a/b of the search, and it is made [l, +infinity] rather than turned
   * around, which the library would refuse */
  { domain[num_inter].up = height_end(args->height);
    if(domain[num_inter].up < domain[num_inter].low)
    { domain[num_inter].up = HUGE_VAL; }
    num_inter++;
  }
  args->num_inter = num_inter;

  if(!print_format)
  /* default print format */
  { print_format = (args->flags & (RATPOINTS_NO_CHECK | RATPOINTS_NO_Y))
                     ? "(%x : %z)\n" : "(%x : %y : %z)\n";
  }

  /* quiet implies not verbose */
  if(quiet) { args->flags &= ~RATPOINTS_VERBOSE; }

  mpz_clear(fff); /* clean up */
  return(0);
}

/* Read the integer that s starts with -- an optional sign and at least one
 * digit, ended by white space or by the end of s -- into x, and return a
 * pointer to what follows it; return NULL if s does not start with one.
 * White space is what isspace() counts as such in the C locale, in which
 * the program runs: spaces, tabs, line feeds, carriage returns, vertical
 * tabs and form feeds. */
char *scan_mpz(char *s, mpz_t x)
{
  int neg = 0; /* flag for sign */

  if(*s == '-') { neg = 1; s++; }
  else if(*s == '+') { s++; }
  if(!('0' <= *s && *s <= '9')) { return NULL; }
  mpz_set_ui(tmp2, 0);
  while('0' <= *s && *s <= '9')
  { mpz_mul_ui(tmp2, tmp2, 10);
    mpz_add_ui(tmp2, tmp2, (unsigned long)(*s - '0'));
    s++;
  }
  if(*s != 0 && !isspace((unsigned char)*s)) { return NULL; }
  if(neg) { mpz_neg(tmp2, tmp2); }
  mpz_set(x, tmp2);
  return s;
}

/* Whether nothing but white space follows the number that ends at end:
 * a number on the command line may have white space before it (strtol and
 * strtod skip that) and after it, so that the carriage return that a
 * script with CRLF line ends leaves on the last argument of a line does
 * no harm. */
static int only_space_after(const char *end)
{
  while(isspace((unsigned char)*end)) { end++; }
  return(*end == 0);
}

/* Read the integer that is the whole of s, but for white space before and
 * after it, into x.  Returns 0 if s is not an integer or if the integer
 * does not fit a long. */
int scan_long(char *s, long *x)
{
  char *end;

  errno = 0;
  *x = strtol(s, &end, 10);
  return(end != s && errno == 0 && only_space_after(end));
}

/* Read the real number that is the whole of s, but for white space before
 * and after it, into x.  Returns 0 if s is not a number or if it is a NaN.
 * An infinity is accepted, since -l -inf and -u inf are interval ends that
 * bound nothing, and a number beyond the range of a double is read as one. */
int scan_double(char *s, double *x)
{
  char *end;

  *x = strtod(s, &end);
  return(end != s && !isnan(*x) && only_space_after(end));
}


/**************************************************************************
 * output routines                                                        *
 **************************************************************************/

/* Print a polynomial with mpz_t coefficients. */
void print_poly(mpz_t *coeffs, long degree)
{
  int flag = 0; /* 0 if at beginning */
  int i;
  char *s;

  for(i = degree; i >= 0; i--)
  { mpz_set(tmp, coeffs[i]);
    if(mpz_cmp_si(tmp, 0) != 0)
    { if(mpz_cmp_si(tmp, 0) > 0)
      { printf(flag ? " + " : ""); }
      else
      { printf(flag ? " - " : "- ");
        mpz_neg(tmp, tmp);
      }
      flag = 1;
      switch(i)
      { case 0: s = mpz_get_str((char *) 0, 10, tmp);
                printf("%s", s); free(s); break;
        case 1: if(mpz_cmp_si(tmp, 1) == 0)
                { printf("x"); }
                else
                { s = mpz_get_str((char *) 0, 10, tmp);
                  printf("%s x", s); free(s);
                }
                break;
        default: if(mpz_cmp_si(tmp, 1) == 0)
                 { printf("x^%d", i); }
                 else
                 { s = mpz_get_str((char *) 0, 10, tmp);
                   printf("%s x^%d", s, i); free(s);
                 }
                 break;
  } } }
  printf("\n");
  fflush(stdout);
}

/* Print a string, recognizing \t, \n, \0, \%, \\. */
void print_string(char *str)
{
  if(str)
  { while(*str)
    { char c = *str++;
      if(c == '\\')
      { c = *str++;
        switch(c)
        { case 0: putchar('\\'); return;
          case 't': printf("\t"); break;
          case 'n': printf("\n"); break;
          case '\\': putchar('\\'); break;
          case '%': putchar('%'); break;
          default: putchar('\\'); putchar(c); break;
      } }
      else putchar(c);
} } }


/* Print message number n. */
void message(long n, long total, ratpoints_args *args)
{
  switch(n)
  { case 0: printf("\n%s\n", RATPOINTS_VERSION);
            break;
    case 1: printf("\n%s\n%s", RATPOINTS_VERSION, usage_str);
            break;
    case 2: if(args->flags & (RATPOINTS_NO_CHECK | RATPOINTS_NO_Y))
            { printf("\n%ld rational point pairs found.\n\n", total); }
            else
            { printf("\n%ld rational points found.\n\n", total); }
            break;
    case 3: printf("Search region:\n  ");
            { long i;
              for(i = 0; i < args->num_inter; i++)
              { if(i) { printf(" U "); }
                printf("[%f, %f]", args->domain[i].low, args->domain[i].up);
            } }
            printf("\n");
            break;
    case 4: printf("%ld moduli used for first stage of sieving,\n",
                   args->sp1_used);
            printf("%ld moduli used for both stages of sieving together,\n",
                   args->sp2_used);
            printf("%ld further primes used in the third stage.\n",
                   args->sp3_used - args->sp2_used);
            break;
    case 5: printf("\nCurve equation is  y^2 = ");
            print_poly(args->cof, args->degree);
            printf("\n");
            break;
    case 6: printf("max. Height = %ld\n", args->height);
            break;
    case 7: if(args->flags & RATPOINTS_REVERSED)
            { printf("Polynomial was reversed for better performance.\n"); }
            if(args->flags & RATPOINTS_USE_SQUARES)
            { printf("(Reversed) polynomial is +-monic of odd degree:\n");
              printf("  could restrict denominators to squares.\n");
            }
            if(args->flags & RATPOINTS_USE_SQUARES1)
            { printf("(Reversed) polynomial has odd degree:\n");
              printf("  could restrict denominators essentially to squares.\n");
            }
            if(args->flags & RATPOINTS_NO_JACOBI)
            { printf("Jacobi symbol test was not used.\n"); }
            if(args->flags & RATPOINTS_NO_CHECK)
            { printf("Points were not checked exactly:\n");
              printf("   some of the printed x-coordinates may not give points.\n");
            }
            break;
    case 8: printf("Search intervals:\n");
            { long i;
              for(i = 0; i < args->num_inter; i++)
              { printf("[%lf, %lf]", args->domain[i].low, args->domain[i].up);
                if(i < args->num_inter -1) { printf(" U "); }
              }
              printf("\n");
            }
            break;
  }
  fflush(stdout);
}

/* Print error message number n and exit. */
void error(long n)
{
  switch(n)
  { case 3: printf("\nToo many coefficients.\n\n"); break;
    case 4: printf("\nIncorrect height argument.\n\n"); break;
    case 5: printf("\nThe polynomial must have degree at least 1.\n\n"); break;
    case 6: printf("\nWrong syntax for optional arguments:\n\n");
    case 2: printf("\n%s\n%s", RATPOINTS_VERSION, usage_str); break;
    case 7: printf("\nIncorrect interval arguments (not numbers, ");
            printf("not alternating, too many, or not ordered).\n\n");
            break;
    case 8: printf("\nPolynomial is not square-free.\n\n"); break;
    case 9: printf("\nBug no. 1 - please report!\n\n"); break;
    case 10: printf("\nNot enough memory.\n\n"); break;
    case 11: printf("\nIncorrect coefficients (not integers separated by ");
             printf("blanks).\n\n");
             break;
  }
  fflush(stdout);
  exit(n);
}
