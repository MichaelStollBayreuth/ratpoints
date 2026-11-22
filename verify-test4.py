#!/usr/bin/env python3
"""Check the reference output of test4.sh independently of the sieve.

usage: verify-test4.py [testbase4]

testbase4 (or the file named) is read invocation by invocation -- each
begins with the line test4.sh prints, "# <arg> <arg> ...".  For every
invocation whose output is a list of points, the points are compared with
a brute-force search: every pair (a, b) with gcd(a, b) = 1 in the
denominator range and the search intervals of the command line, within
the height bound, is tested for F(a, b) = sum c_i a^i b^(D-i) being a
square, in arbitrary precision and without any sieve.  The ends of the
intervals are the numbers given, taken exactly as rational numbers (1.6
is 8/5), as the program takes them (see c_end).  The same holds
for the x-coordinates -x and -y print, for the single point of -1 and for
the count -z reports.  A search too large for that (the two curves at
heights 30000 to 100000 in the third part) is checked against PARI/GP's
hyperellratpoints when gp is on the PATH, which is another implementation
of the search, and otherwise reported as not checked -- except when the
degree is odd and the denominator range is what makes the search large:
then only the denominators a point can have are tried, e m^2 with e a
divisor of the leading coefficient (see brute_points), which is how the
runs at the top of the range of a long are checked; and a run of -1 whose
point has the denominator 1, the first of every search, is checked by
trying the numerators of that denominator in the order of the search up
to the point printed, which has to be the first point found (see
first_point_b1).  The pairs -x prints beyond the points are survivors of
the sieve, which need not be points, but they must lie in the search:
coprime, in the height box, in the denominator range and in the search
intervals.

Whether the program refuses its arguments, and with which exit status, is
decided by the rules of main.c, which parse_args follows: the coefficients
are integers separated by white space, at most RATPOINTS_MAX_DEGREE + 1 of
them, in the first argument or, when that is @name, in the file name,
which test4.sh writes and FILES records (a name that FILES does not have
cannot be read, and a zero byte in a file is refused like a word that is
no integer); a number is the whole argument but for white space around it, an
integer within the range of a long or a real number that is not NaN; and
the search intervals alternate and are ordered by their ends as the
numbers they are, however many there are.  An
invocation that reports an exit status must report the one these rules
give, with its message, and one they refuse must report it; one they
accept may end with the message for a polynomial that is not squarefree,
which is checked, or else is checked like any other.  Invocations that
test an output format or the report of -v are not checked here beyond
what they print, except that the points a -v report contains are still
compared, and that wherever the equation of the curve (without -q) or a
list of the coefficients (the report of -v, at its start and after a
reversal) is printed, it must be the polynomial in full below the degree
RATPOINTS_PRINT_DEGREE and its short form from there on (see
curve_line), and that wherever the report gives the numbers of the
moduli of the first stage and of the first two stages, they must keep
within -n and -N, which are hard limits also where the first stage is
chosen by its rule, and that wherever the report says what one exact check
is put at, it says it in the form and with the numbers that the constants
of ratpoints.h give, with the modular test before the exact check where it
costs less (see check_report).  A run whose standard output refuses every write
(announced with ">/dev/full" or ">&-" after the arguments) must report
the failure on standard error and exit with status 13 when the program
writes anything there -- the message of a refusal, the banner, a point
(see writes) -- and must print nothing and exit with 0 when it writes
nothing.

Prints one line per invocation that fails and a summary; exits 1 on any
failure.  Takes a minute or two.
"""
import sys, os, re, math, subprocess, shutil
from fractions import Fraction

MAX_BRUTE = 40_000_000   # candidate pairs a brute-force search is allowed
rule_used = 0            # searches over the denominators e m^2 only

# the list of the coefficients that test4.sh reads from a file has two of a
# million digits, beyond the 4300 to which Python limits the conversion of
# a string to an int by default
if hasattr(sys, 'set_int_max_str_digits'): sys.set_int_max_str_digits(0)

# The files that test4.sh writes for the runs that read the list of the
# coefficients from a file, @name, by name, with what each holds; '-' is
# what the run of @- is given on standard input.  A name that is not here
# (no file, a directory, the empty name of @ alone) cannot be read.
BIG_ZEROS = '0'*999_999
FILES = {'lf-tab': ' 1\t0\n\n\t1\r\n', 'empty': '', 'nul': '1 0 1\0 2\n',
         'word': '1 0 x\n', '-': '1\n0\n1\n',
         'big': '-1' + BIG_ZEROS + '6 0 1' + BIG_ZEROS + '7\n',
         'deg250': '1 -1 2' + ' 0'*122 + ' 2' + ' 0'*124 + ' 1\n',
         'deg251': '1 2 1' + ' 0'*247 + ' -1 2\n',
         'deg999': '0 1 2' + ' 0'*995 + ' -1 2\n',
         'deg999z': '0 1 2' + ' 0'*995 + ' -1 2 0\n',
         'deg1000': '1 -1 2' + ' 0'*497 + ' 2' + ' 0'*499 + ' 1\n',
         'deg1001': '0 1 2' + ' 0'*997 + ' -1 2\n',
         'deg3000': '%d -1 2' % (10**100 + 7)**2 + ' 0'*1497 + ' %d' % (2*(10**100 + 7))
                    + ' 0'*1499 + ' 1\n'}

def macro(name):
    """the macro name of the ratpoints.h next to this script, the default
    that its #ifndef gives; without it the rules of the command line or of
    the output are not known"""
    fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ratpoints.h')
    try:
        return int(re.search(r'^#\s*define\s+%s\s+(\d+)' % name, open(fn).read(), re.M)[1])
    except (OSError, TypeError):
        sys.exit('verify-test4.py: no %s in %s' % (name, fn))
MAX_DEGREE = macro('RATPOINTS_MAX_DEGREE')
# from this degree on, the curve and the lists of its coefficients are
# printed in their short forms
PRINT_DEGREE = macro('RATPOINTS_PRINT_DEGREE')

def macro_double(name):
    """the constant of type double that the macro name of ratpoints.h is,
    the default of its #ifndef"""
    fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ratpoints.h')
    try:
        return float(re.search(r'^#\s*define\s+%s\s+([0-9.eE+-]+)' % name, open(fn).read(), re.M)[1])
    except (OSError, TypeError, ValueError):
        sys.exit('verify-test4.py: no %s in %s' % (name, fn))
# what one exact check and the modular test before it cost (check_report)
CHECK = {n: macro_double('RATPOINTS_CHECK_' + n) for n in ('STEP', 'LIMB', 'CALL', 'ROOT', 'REFERENCE')}
TEST = {n: macro_double('RATPOINTS_TEST_' + n) for n in ('STEP', 'CALL', 'TRIED')}
TEST_MODULI = 16   # the primes of the modular test (find_points.c)

# White space as isspace() has it in the C locale, in which ratpoints runs;
# a number may have it before and after it (strtol and strtod skip it
# before, main.c after), and it separates the coefficients
WS = ' \t\n\v\f\r'
INT = re.compile(r'[ \t\n\v\f\r]*[+-]?[0-9]+[ \t\n\v\f\r]*')
REAL = re.compile(r'[ \t\n\v\f\r]*[+-]?(?:'
                  r'(?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[eE][+-]?[0-9]+)?'
                  r'|0[xX](?:[0-9a-fA-F]+\.?[0-9a-fA-F]*|\.[0-9a-fA-F]+)(?:[pP][+-]?[0-9]+)?'
                  r'|[iI][nN][fF](?:[iI][nN][iI][tT][yY])?'
                  r')[ \t\n\v\f\r]*')

def c_long(s):
    """main.c's scan_long: the integer that s is, None if s is not one as a
    whole or if it does not fit a long"""
    if not INT.fullmatch(s): return None
    v = int(s.strip(WS))
    return v if -2**63 <= v < 2**63 else None

def c_double(s):
    """main.c's scan_double: the double that strtod makes of s (rounded to
    nearest; beyond the range of a double, an infinity), None if s is not a
    number as a whole or if it is a NaN"""
    if not REAL.fullmatch(s): return None
    t = s.strip(WS)
    if t.lstrip('+-')[:2] in ('0x', '0X'):
        try: return float.fromhex(t)
        except OverflowError: return -math.inf if t.startswith('-') else math.inf
    return float(t)

HALF_TINY = Fraction(1, 2**1075)   # half the least positive double

def c_end(s):
    """main.c's scan_end: the number s denotes as the end of a search
    interval, exactly -- a decimal number m 10^k as the rational number it
    is (1.6 is 8/5), a hexadecimal one m 2^k as the dyadic one (which
    float.fromhex would round to a double when it has more digits than a
    double holds) -- or an infinity, which a number that strtod reads as one
    (beyond the range of a double) also is; a number other than 0 that
    strtod reads as 0 is taken as plus or minus half the least positive
    double.  None if s is not a number as a whole or if it is a NaN"""
    x = c_double(s)
    if x is None or math.isinf(x): return x
    t = s.strip(WS)
    sign = -1 if t.startswith('-') else 1
    t = t.lstrip('+-')
    if t[:2] in ('0x', '0X'):
        base, m = 16, re.fullmatch(r'0[xX]([0-9a-fA-F]*)(?:\.([0-9a-fA-F]*))?(?:[pP]([+-]?[0-9]+))?', t)
    else:
        base, m = 10, re.fullmatch(r'([0-9]*)(?:\.([0-9]*))?(?:[eE]([+-]?[0-9]+))?', t)
    frac = m[2] or ''
    digits = int(m[1] + frac, base)
    # the exponent is looked at only when strtod has found a double other
    # than 0 for the number, so that it is not too large to compute with
    if digits == 0: return Fraction(0)
    if x == 0: return sign*HALF_TINY
    if base == 16: return sign*digits*Fraction(2)**(int(m[3] or 0) - 4*len(frac))
    return sign*digits*Fraction(10)**(int(m[3] or 0) - len(frac))

def coefficients(s):
    """main.c's reading of the coefficient list: integers, each an optional
    sign and at least one digit, separated by white space; None if the list
    is not of that form"""
    toks = [t for t in re.split(r'[ \t\n\v\f\r]+', s) if t]
    if not all(re.fullmatch(r'[+-]?[0-9]+', t) for t in toks): return None
    return [int(t) for t in toks]

def parse_args(args):
    """The arguments as main.c reads them: what is searched and printed, and
    in 'exit' the status the program exits with when it refuses them (None
    when it does not)"""
    o = dict(q=False, v=False, z=False, one=False, inf=True, x=False, y=False,
             s=False, S=None, dl=1, du=None, ivs=[], fmt=None, other={},
             exit=None, strs=False, keep=False)
    if len(args) < 2: o['exit'] = 2; return o
    lst = args[0]
    if lst.startswith('@'):
        lst = FILES.get(lst[1:])
        if lst is None: o['exit'] = 12; return o
        if '\0' in lst: o['exit'] = 11; return o
    c = coefficients(lst)
    if c is None: o['exit'] = 11; return o
    if len(c) > MAX_DEGREE + 1: o['exit'] = 3; return o
    if len(c) <= 1: o['exit'] = 5; return o
    o['c'] = c
    H = c_long(args[1])
    if H is None or H < 1: o['exit'] = 4; return o
    o['H'] = H
    n = len(args); i = 2; l_seen = False; low = None
    while i < n:
        a = args[i]
        if not a.startswith('-'): o['exit'] = 6; return o
        k = a[1:2]; k2 = a[2:3]   # main.c looks at these two characters only
        if k in ('l', 'u'):
            if i + 1 >= n: o['exit'] = 7; return o
            i += 1
            if k == 'l':
                if l_seen: o['exit'] = 7; return o
                low = c_end(args[i])
                if low is None: o['exit'] = 7; return o
                if o['ivs'] and o['ivs'][-1][1] >= low: o['exit'] = 7; return o
                l_seen = True
            else:
                if not l_seen:
                    if o['ivs']: o['exit'] = 7; return o
                    low = None   # no bound below: the height bound or -infinity
                up = c_end(args[i])
                if up is None or (low is not None and low > up): o['exit'] = 7; return o
                o['ivs'].append((low, up)); l_seen = False
        elif k in ('p', 'F', 't', 'n', 'N', 'R', 'A', 'P') or (k == 'd' and k2 in ('l', 'u')):
            if i + 1 >= n: o['exit'] = 6; return o
            i += 1
            v = c_long(args[i])
            if v is None: o['exit'] = 6; return o
            if k == 'd': o['d' + k2] = v
            else: o['other'][k] = v
        elif k in ('r', 'U', 'C', 'Q', 'W'):
            if i + 1 >= n: o['exit'] = 6; return o
            i += 1
            if c_double(args[i]) is None: o['exit'] = 6; return o
            o['other'][k] = args[i]
        elif k == 'f' and k2 in ('', 's', 'm', 'e'):
            if i + 1 >= n: o['exit'] = 6; return o
            i += 1
            if k2 == '': o['fmt'] = args[i]
            else: o['strs'] = True
        elif k == 'S':
            o['s'] = False; o['S'] = 10
            if i + 1 < n and not args[i+1].startswith('-'):
                i += 1
                o['S'] = c_long(args[i])
                if o['S'] is None: o['exit'] = 6; return o
        elif k in ('q', 'v', 'z', 'Z', '1', 'i', 'I', 'x', 'X', 'y', 'Y', 's',
                   'k', 'K', 'j', 'J'):
            if k in ('q', 'v', 'z', 'x', 'y', 's'): o[k] = True
            elif k in ('Z', 'X', 'Y'): o[k.lower()] = False
            elif k == '1': o['one'] = True
            elif k in ('i', 'I'): o['inf'] = (k == 'I')
            elif k in ('k', 'K'): o['keep'] = (k == 'k')
        else: o['exit'] = 6; return o
        i += 1
    if l_seen: o['ivs'].append((low, None))  # no bound above
    return o

# the messages of the refusals, by exit status (2 prints the usage only)
MESSAGE = {3: 'Too many coefficients.',
           4: 'Incorrect height argument.',
           5: 'The polynomial must have degree at least 1.',
           6: 'Wrong syntax for optional arguments:',
           7: 'Incorrect interval arguments (not numbers, not alternating, or not ordered).',
           11: 'Incorrect coefficients (not integers separated by blanks).',
           12: 'Cannot read the file of the coefficients.'}

def denominators(H, dl, du):
    """the denominator range the library searches: b_low and b_high below 1
    stand for 1 and the height bound, and b_high is at most that"""
    return max(dl, 1), (H if (du is None or du < 1) else min(du, H))

def numerators(b, lo, up, H):
    """the a with lo <= a/b <= up and |a| <= H: ceil(b lo) .. floor(b up);
    an end None is no bound"""
    if lo is None or lo == -math.inf: alo = -H
    elif lo == math.inf: alo = H + 1
    else: alo = max(-H, math.ceil(Fraction(lo)*b))
    if up is None or up == math.inf: ahi = H
    elif up == -math.inf: ahi = -H - 1
    else: ahi = min(H, math.floor(Fraction(up)*b))
    return alo, ahi

def in_search(a, b, H, dl, du, ivs):
    """whether (a : b) is a pair the search looks at: coprime, b in the
    denominator range, |a| <= H, and a/b in one of the intervals"""
    b_low, b_high = denominators(H, dl, du)
    if not (b_low <= b <= b_high and abs(a) <= H and math.gcd(a, b) == 1): return False
    return any(alo <= a <= ahi for alo, ahi in
               (numerators(b, lo, up, H) for lo, up in (ivs or [(None, None)])))

def brute_points(c, H, dl, du, ivs, want_inf):
    """The points the program should print, as a dict (a, b) -> set of y
    (both signs), plus the points at infinity.  None if the search is too
    large."""
    # strip leading zeros as the library does; the genus must not drop
    d = len(c) - 1
    while d > 0 and c[d] == 0: d -= 1
    D = d + (d & 1)
    pts = {}
    if want_inf:
        if d & 1: pts[(1, 0)] = {0}
        else:
            r = math.isqrt(c[d]) if c[d] >= 0 else -1
            if r >= 0 and r*r == c[d]: pts[(1, 0)] = {r, -r}
    b_low, b_high = denominators(H, dl, du)
    # the search region: the intervals of the command line, closed, their
    # ends the numbers given, taken exactly (see c_end); a missing first -l
    # or last -u, or no interval at all, is no bound (the program gives the
    # library the height bound rounded outward, which holds every a/b with
    # |a| <= H, or an infinite end where that would turn the interval
    # around)
    if not ivs: ivs = [(None, None)]
    # the denominators to try: every one in the range -- or, when the range
    # is too long for that and the degree is odd, only those a point can
    # have.  For odd d, F(a, b) = b (c_d a^d + b (...)), so with gcd(a, b) = 1
    # a prime not dividing c_d divides F(a, b) exactly as often as b, which
    # for a square must be an even number of times: b = e m^2 with e a
    # divisor of |c_d| (every divisor is tried, not only the squarefree ones)
    bs = range(b_low, b_high + 1)
    if b_high - b_low >= 10**7:
        if not (d & 1) or abs(c[d]) > 10**6: return None
        divs = [e for e in range(1, abs(c[d]) + 1) if c[d] % e == 0]
        ms = [(e, math.isqrt((b_low - 1)//e) + 1, math.isqrt(b_high//e)) for e in divs]
        if sum(max(0, hi - lo + 1) for e, lo, hi in ms) > 10**6: return None
        bs = sorted(e*m*m for e, lo, hi in ms for m in range(lo, hi + 1))
        global rule_used
        rule_used += 1
    # the size of the search, estimated in floating point
    total = 0
    for b in bs:
        for lo, up in ivs:
            flo = -float(H) if lo is None else max(float(b)*lo, -float(H))
            fup = float(H) if up is None else min(float(b)*up, float(H))
            total += max(0, fup - flo) + 1
        if total > MAX_BRUTE: return None
    for b in bs:
        pw = [b**(D - i) for i in range(D + 1)]   # b^(D-i)
        cb = [c[i]*pw[i] if i <= d else 0 for i in range(D + 1)]
        for lo, up in ivs:
            alo, ahi = numerators(b, lo, up, H)
            for a in range(alo, ahi + 1):
                if math.gcd(a, b) != 1: continue
                v = 0
                for i in range(D, -1, -1): v = v*a + cb[i]
                if v < 0: continue
                r = math.isqrt(v)
                if r*r == v: pts[(a, b)] = {r, -r}
    return pts

def gp_points(c, H):
    """PARI/GP's search, for the two searches too large to brute-force;
    None if gp is not available or fails"""
    if shutil.which('gp') is None: return None
    poly = '+'.join('(%d)*x^%d' % (ci, i) for i, ci in enumerate(c))
    script = 'v=hyperellratpoints(%s,%d); for(i=1,#v,print(v[i][1]," ",v[i][2])); quit' % (poly, H)
    try:
        out = subprocess.run(['gp', '-q'], input=script, capture_output=True, text=True, timeout=600).stdout
    except Exception:
        return None
    d = len(c) - 1
    while d > 0 and c[d] == 0: d -= 1
    pts = {}
    for line in out.split('\n'):
        f = line.split()
        if len(f) != 2: continue
        x = Fraction(f[0]); a, b = x.numerator, x.denominator
        if abs(a) > H or b > H: continue
        D = d + (d & 1)
        v = sum(c[i]*a**i*b**(D-i) for i in range(d + 1))
        r = math.isqrt(v)
        if r*r != v: return None
        pts[(a, b)] = {r, -r}
    if d & 1: pts[(1, 0)] = {0}
    else:
        r = math.isqrt(c[d]) if c[d] >= 0 else -1
        if r >= 0 and r*r == c[d]: pts[(1, 0)] = {r, -r}
    return pts

def first_point_b1(c, H, ivs, a_end):
    """For -1 on a search too large for brute_points, when the point printed
    has the denominator 1: the first point of the search with that
    denominator, as a dict like brute_points', found by trying its
    numerators in the order of the search -- the intervals in turn, each
    upward -- up to a_end; empty if there is none, None if they are too
    many.  The denominator 1 is the first of every search that includes it,
    whatever the shape of the loop over the denominators (the squares, the
    squares times the divisors of the leading coefficient, the divisor 1
    first, or every one), so this is the point -1 has to print."""
    d = len(c) - 1
    while d > 0 and c[d] == 0: d -= 1
    if not ivs: ivs = [(None, None)]
    tried = 0
    for lo, up in ivs:
        alo, ahi = numerators(1, lo, up, H)
        for a in range(alo, min(ahi, a_end) + 1):
            tried += 1
            if tried > 10**6: return None
            v = sum(c[i]*a**i for i in range(d + 1))   # F(a, 1)
            if v < 0: continue
            r = math.isqrt(v)
            if r*r == v: return {(a, 1): {r, -r}}
    return {}

def squarefree_mod(c, p):
    """whether f = sum c_i x^i, of degree len(c) - 1 modulo the prime p,
    is squarefree modulo p: the gcd of f and f' modulo p is a constant"""
    def strip(a):
        while a and a[-1] == 0: a.pop()
        return a
    f = strip([x % p for x in c])
    g = strip([(i*c[i]) % p for i in range(1, len(c))])
    while g:
        inv = pow(g[-1], p - 2, p)
        while len(f) >= len(g):
            q = f[-1]*inv % p; k = len(f) - len(g)
            for i in range(len(g)): f[k+i] = (f[k+i] - q*g[i]) % p
            strip(f)
        f, g = g, f
    return len(f) == 1

def squarefree(c):
    """whether f = sum c_i x^i, leading zeros dropped, is squarefree as a
    binary form of even degree: no repeated root, and no double root at
    infinity (the degree dropped from odd to even or by two or more)"""
    d0 = len(c) - 1; d = d0
    while d > 0 and c[d] == 0: d -= 1
    if d <= 0: return False
    if (d + 1) >> 1 < (d0 + 1) >> 1: return False
    # f is squarefree when it is so modulo a prime that does not divide its
    # leading coefficient, its discriminant then not being 0 modulo the
    # prime; this decides it fast at a large degree, and the gcd over Q
    # below, which can take minutes there, decides the rest
    for p in (2**61 - 1, 2**31 - 1):
        if c[d] % p and squarefree_mod(c[:d+1], p): return True
    f = [Fraction(x) for x in c[:d+1]]
    g = [Fraction(i*c[i]) for i in range(1, d+1)]
    # polynomial gcd over Q; squarefree iff it is a constant
    while g and any(g):
        while g and g[-1] == 0: g.pop()
        if not g: break
        while len(f) >= len(g) and any(f):
            q = f[-1]/g[-1]; k = len(f) - len(g)
            for i in range(len(g)): f[k+i] -= q*g[i]
            while f and f[-1] == 0: f.pop()
        f, g = g, f
    return len(f) == 1

def poly_text(c):
    """the polynomial sum c_i x^i as main.c's print_poly prints it"""
    s = ''
    for i in range(len(c) - 1, -1, -1):
        a = c[i]
        if a == 0: continue
        if a > 0: s += ' + ' if s else ''
        else: s += ' - ' if s else '- '; a = -a
        if i == 0: s += '%d' % a
        elif i == 1: s += 'x' if a == 1 else '%d x' % a
        else: s += ('x^%d' % i) if a == 1 else '%d x^%d' % (a, i)
    return s

def curve_line(line, c):
    """What a line that prints the curve or the list of its coefficients
    must be, on the coefficients c given; None for any other line.  The
    equation of the curve (message 5 of main.c) is the polynomial in full
    below the degree RATPOINTS_PRINT_DEGREE, its leading zero coefficients
    not counted, and from there on that degree and the number of bits of
    its largest coefficient (mpz_sizeinbase, which counts 1 for a 0); so is
    the list of the coefficients of the report of -v, with its leading
    zeros dropped, and so is the list after the reversal of an odd degree
    with a zero constant term, which swaps c_n and c_(d+1-n) for n = 1 to
    d.  The equation of the curve is followed by a blank line in both forms
    (see the check in main)"""
    d = len(c) - 1
    while d > 0 and c[d] == 0: d -= 1
    bits = max(abs(a).bit_length() or 1 for a in c[:d+1])
    if line.startswith('Curve equation is  y^2 = '):
        if d < PRINT_DEGREE: return 'Curve equation is  y^2 = ' + poly_text(c)
        return ('Curve equation is  y^2 = f(x), where f has degree %d and coefficients'
                ' of at most %d bits' % (d, bits))
    for head, cs in (('  coefficients: ', c[:d+1]),
                     ('  coefficients now: ', c[:1] + [c[d+1-n] for n in range(1, d+1)])):
        if line.startswith(head):
            if d < PRINT_DEGREE: return head[:-1] + ''.join(' %d' % a for a in cs)
            return head + 'of at most %d bits (not printed)' % bits
    return None

def check_report(o):
    """The lines in which the report of -v says what one exact check is put
    at, for the run that parse_args read into o: find_points.c's
    check_cost() -- the -W given when it is positive, and otherwise the
    estimate from the degree, the largest coefficient and the height bound
    --, compared in sieving_info with what the modular test costs a
    survivor, TEST_TRIED*(degree*TEST_STEP + TEST_CALL): where the test
    costs less and there is a check at all (no -x), the test is on and the
    third stage is told what it costs, in three lines, and otherwise what
    the check costs, in two.  The same operations in doubles, in the same
    order; the degree is that of the search, the leading zeros dropped and
    one less after the reversal of an even degree with a zero constant term
    (which leaves the coefficients that are not zero as they were)"""
    c, H = o['c'], o['H']
    d = len(c) - 1
    while d > 0 and c[d] == 0: d -= 1
    du = o['du'] if o['du'] is not None and o['du'] >= 1 else H
    if not o['keep'] and not o['ivs'] and o['dl'] <= 1 and du >= H \
       and d % 2 == 0 and c[0] == 0:
        d -= 1
    w = c_double(o['other']['W']) if 'W' in o['other'] else -1.0
    if w > 0.0: exact = w
    else:
        hbits = math.log(float(H) + 1.0)/math.log(2.0)
        cbits = 1.0
        for a in c:
            if a != 0 and float(abs(a).bit_length()) > cbits: cbits = float(abs(a).bit_length())
        fbits = cbits + float(d)*hbits
        limbs = fbits/64.0
        mid = 0.5*(cbits + fbits)/64.0
        root = math.ceil(0.5*limbs)
        if root < 1.0: root = 1.0
        exact = float(d)*(CHECK['STEP'] + CHECK['LIMB']*mid) + CHECK['CALL'] \
                + CHECK['ROOT']*(root - 1.0)
        if d & 1: exact += CHECK['STEP'] + CHECK['LIMB']*limbs
    test = TEST['TRIED']*(float(d)*TEST['STEP'] + TEST['CALL'])
    modular = not o['x'] and test < exact
    rel = (test if modular else exact)/CHECK['REFERENCE']
    if rel <= 0.0: rel = 1.0
    if modular:
        return ['  one exact check is put at %.0f cycles; survivors are tested modulo %d'
                % (exact, TEST_MODULI),
                '    primes first, at %.0f cycles, %.2f times what a check costs'
                % (rel*CHECK['REFERENCE'], rel),
                '    on the curves the third stage was tuned on']
    return ['  one exact check is put at %.0f cycles, %.2f times what it costs'
            % (rel*CHECK['REFERENCE'], rel),
            '    on the curves the third stage was tuned on']

POINT = re.compile(r'^\((-?\d+) : (-?\d+) : (\d+)\)$')
XONLY = re.compile(r'^\((-?\d+) : (\d+)\)$')
FOUND = re.compile(r'^(\d+) rational points? (?:pairs )?found\.$')
# the numbers of the moduli of the first stage and of the first two stages
# together in the report of -v
MODULI = re.compile(r'==> use (\d+) moduli in the first stage, (\d+) altogether,$')

# the runs of w() in test4.sh, whose standard output refuses every write
WRITE_FAILS = ('>/dev/full', '>&-')

def writes(o):
    """Whether the program writes anything to its standard output on the
    command line that parse_args read into o: the message of a refusal (for
    status 2 the usage), the banner and the report without -q, the message
    for a polynomial that is not squarefree, a point the search finds unless
    -z suppresses it.  None if that is not decided here: with a string of
    -fs, -fm or -fe, under -x (which prints survivors of the sieve as well),
    or when the search is too large for brute_points"""
    if o['exit'] is not None or not o['q']: return True
    if o['strs'] or o['x']: return None
    if not squarefree(o['c']): return True
    if o['z']: return False
    pts = brute_points(o['c'], o['H'], o['dl'], o['du'], o['ivs'], o['inf'])
    return None if pts is None else bool(pts)

def main():
    fn = sys.argv[1] if len(sys.argv) > 1 else 'testbase4'
    # newline='': a carriage return in an argument stays in its line
    lines = open(fn, newline='').read().split('\n')
    # split into invocations
    invs = []
    cur = None
    for line in lines:
        # an announce line "# <arg> <arg> ..."; a format test whose output
        # ends without a newline glues the next announce onto its last line
        k = line.find('# <')
        if k > 0 and cur is not None:
            cur[1].append(line[:k]); line = line[k:]
        if line.startswith('# <') or line == '#':
            # a filtered run announces its filter after " | "
            head, _, filt = line.partition(' | ')
            cur = (re.findall(r'<([^>]*)>', head), [], filt)
            invs.append(cur)
        elif line.startswith('====') or line.startswith('----'):
            cur = None
        elif cur is not None:
            cur[1].append(line)
    checked = failed = skipped = gp_used = refusals = write_fails = curves = 0
    moduli = reports = 0
    notes = []
    for args, out, filt in invs:
        o = parse_args(args)
        def fail(msg):
            nonlocal failed
            failed += 1; print('FAIL:', ' '.join('<%s>' % a for a in args)); print('   ', msg)
        exit_line = [l for l in out if l.startswith('exit ')]
        body = [l for l in out if not l.startswith('exit ')]
        # the curve and the lists of its coefficients, wherever they are
        # printed, in full or in their short forms
        if 'c' in o:
            for i, l in enumerate(body):
                want = curve_line(l, o['c'])
                if want is None: continue
                if l != want: fail('printed "%s", where it must be "%s"' % (l[:100], want[:100])); break
                # message 5 ends in a blank line, in both forms, after which
                # comes the height bound (message 6) or, under -v, the search
                # region (message 3) when intervals are given and otherwise
                # the report, which begins with a blank line of its own
                # (unless a filter of the run has taken the lines out)
                after = body[i+1:i+3]
                blank2 = o['v'] and not o['ivs']
                if not filt and l.startswith('Curve equation is') and \
                   (after[:1] != [''] or (len(after) > 1 and (after[1] == '') != blank2)):
                    fail('"%s" is followed by %s, not by a blank line and %s'
                         % (l[:100], after, 'another' if blank2 else 'a line of text')); break
                curves += 1
        # -n and -N are hard limits on the moduli of the first stage and of
        # the first two stages, also where the rule chooses the first stage
        # (a filter of the run may have cut the line before the "==>")
        for l in body:
            m = MODULI.search(l if '==>' in l else '==> ' + l)
            if not m or 'c' not in o: continue
            s1, s2 = int(m[1]), int(m[2])
            lim1, lim2 = o['other'].get('n', -1), o['other'].get('N', -1)
            if s1 > s2 or (lim2 >= 0 and s2 > lim2) or (lim1 >= 0 and s1 > lim1):
                fail('%d moduli in the first stage and %d altogether, beyond -n %d or -N %d'
                     % (s1, s2, lim1, lim2))
            moduli += 1
        # what one exact check is put at, and whether the modular test is on
        for i, l in enumerate(body):
            if not l.startswith('  one exact check is put at') or 'c' not in o: continue
            want = check_report(o)
            got = body[i:i+len(want)]
            if got != want: fail('the report of the check: %s, where it must be %s' % (got, want))
            reports += 1
        if filt in WRITE_FAILS:
            # what the program prints on standard error, where it reports
            # that a write failed, and its exit status
            w = writes(o)
            code = int(exit_line[-1].split()[1]) if exit_line else None
            text = [l for l in body if l]
            if w is None:
                skipped += 1; notes.append('not checked (output %s): %s' % (filt, ' '.join(args[:2])))
            elif w and (code != 13 or text != ['Error writing to standard output.']):
                fail('output %s: exit %s and %s, where the program writes and must fail with 13'
                     % (filt, code, text))
            elif not w and (code != 0 or text):
                fail('output %s: exit %s and %s, where the program writes nothing and must exit with 0'
                     % (filt, code, text))
            else: write_fails += 1
            continue
        if exit_line:
            code = int(exit_line[-1].split()[1])
            # a test of the arguments: refused with the status and the
            # message the rules give -- or accepted, and then checked like
            # any other run, one that ended in the message for a polynomial
            # that is not squarefree included
            if o['exit'] is not None or code not in (0, 8):
                if code != o['exit']:
                    fail('exit %d, where the rules of the command line give %s'
                         % (code, o['exit'] or 'none'))
                elif code in MESSAGE and MESSAGE[code] not in body:
                    fail('exit %d without its message "%s"' % (code, MESSAGE[code]))
                else: refusals += 1
                continue
        elif o['exit'] is not None:
            fail('the rules of the command line refuse these arguments (exit %d)' % o['exit'])
            continue
        if filt and 'found' not in filt:
            skipped += 1; continue   # a filtered report: no points to check
        c = o['c']; H = o['H']
        # a run that ended in the message for a non-squarefree polynomial:
        # right when the polynomial is one, wrong otherwise -- and a
        # squarefree polynomial must not be refused
        refused = any(l == 'Polynomial is not square-free.' for l in body)
        sqf = squarefree(c)
        if refused != (not sqf):
            fail('refused as not squarefree: %s, squarefree: %s' % (refused, sqf))
            continue
        if refused: checked += 1; continue
        want_inf = o['inf']
        printed = {}
        xlist = []; count = None
        for l in body:
            m = POINT.match(l)
            if m:
                a, y, b = int(m[1]), int(m[2]), int(m[3])
                printed.setdefault((a, b), []).append(y); continue
            m = XONLY.match(l)
            if m: xlist.append((int(m[1]), int(m[2]))); continue
            m = FOUND.match(l)
            if m: count = int(m[1])
        if o['fmt'] or o['strs']:
            skipped += 1; continue   # formats: the points are checked elsewhere
        if o['z'] and count is None:
            skipped += 1; continue
        # the expected set
        big = H > 5000 and len(str(H)) < 15
        pts = None
        if not big or o['du'] is not None or o['ivs']:
            pts = brute_points(c, H, o['dl'], o['du'], o['ivs'], want_inf)
        if pts is None and o['one'] and o['dl'] == 1 and len(printed) == 1 \
           and not o['x'] and not o['y'] and list(printed)[0][1] == 1:
            pts = first_point_b1(c, H, o['ivs'], list(printed)[0][0])
        if pts is None:
            if o['dl'] == 1 and o['du'] is None and not o['ivs']:
                pts = gp_points(c, H)
                if pts is not None:
                    gp_used += 1
                    if not want_inf: pts.pop((1, 0), None)
            if pts is None:
                skipped += 1; notes.append('not checked (too large, no gp): %s' % ' '.join(args[:2])); continue
        checked += 1
        if o['z']:
            # the count: one per point with y = 0, two otherwise; one per
            # point (x-coordinate, or survivor) with -y or -x
            exp = 1 if o['one'] and pts else 0
            if not o['one']:
                exp = sum(1 if (0 in ys and len(ys) == 1) or o['x'] or o['y'] else 2 for ys in pts.values())
            if count != exp: fail('count %d, expected %d' % (count, exp))
            continue
        xset = set(xlist)
        if len(xset) != len(xlist): fail('an x-coordinate printed twice'); continue
        if o['one']:
            allp = printed or {k: None for k in xset}
            if len(allp) > 1 or sum(len(v) for v in printed.values()) > 1:
                fail('more than one point with -1'); continue
            if pts and not allp: fail('-1 printed nothing, but there are points')
            elif allp:
                (a, b), = allp.keys()
                if (a, b) not in pts:
                    fail('-1 printed (%d : %d), not a point (or, where only the numerators '
                         'of the denominator 1 up to it are tried, not the first)' % (a, b))
                elif printed and printed[(a, b)][0] not in pts[(a, b)]: fail('-1: wrong y')
            continue
        if o['x'] or o['y']:
            exp = set(pts.keys())
            # what -x prints beyond the points are survivors of the sieve,
            # and the sieve looks only at the pairs of the search (the
            # points at infinity aside, which a reversed polynomial reaches
            # as x = 0)
            stray = [p for p in xset - exp if p != (1, 0) and not
                     in_search(p[0], p[1], H, o['dl'], o['du'], o['ivs'])] if o['x'] else []
            if xset != exp:
                extra = xset - exp; missing = exp - xset
                if missing: fail('missing x-coordinates: %s' % sorted(missing)[:5])
                elif stray: fail('-x printed pairs outside the search: %s' % sorted(stray)[:5])
                elif o['x']: notes.append('-x printed %d survivor(s) that are not points: %s' % (len(extra), ' '.join(args[:2])))
                else: fail('extra x-coordinates: %s' % sorted(extra)[:5])
            continue
        # the full list: for each point both y unless y = 0
        exp = {k: sorted(v) for k, v in pts.items()}
        got = {k: sorted(v) for k, v in printed.items()}
        if exp != got:
            missing = [k for k in exp if k not in got]; extra = [k for k in got if k not in exp]
            wrong = [k for k in exp if k in got and exp[k] != got[k]]
            fail('missing %s extra %s wrong y %s' % (sorted(missing)[:5], sorted(extra)[:5], wrong[:5]))
    print('%d invocations checked (%d by gp, %d over the denominators e m^2 only), '
          '%d refusals of the arguments checked, %d runs on a failing standard output checked, '
          '%d lines of the curve or its coefficients checked, '
          '%d numbers of the moduli of the first two stages checked, '
          '%d reports of the cost of a check checked, '
          '%d skipped (formats, reports), %d failed'
          % (checked, gp_used, rule_used, refusals, write_fails, curves, moduli, reports,
             skipped, failed))
    for n in notes: print(' ', n)
    sys.exit(1 if failed else 0)

if __name__ == '__main__':
    main()
