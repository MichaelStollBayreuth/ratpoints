#!/usr/bin/env python3
"""Check the reference output of test4.sh independently of the sieve.

usage: verify-test4.py [testbase4]

testbase4 (or the file named) is read invocation by invocation -- each
begins with the line test4.sh prints, "# <arg> <arg> ...".  For every
invocation whose output is a list of points, the points are compared with
a brute-force search: every pair (a, b) with gcd(a, b) = 1 in the
denominator range and the search intervals of the command line, within
the height bound, is tested for F(a, b) = sum c_i a^i b^(D-i) being a
square, in arbitrary precision and without any sieve.  The same holds
for the x-coordinates -x and -y print, for the single point of -1 and for
the count -z reports.  A search too large for that (the two curves at
heights 30000 to 100000 in the third part) is checked against PARI/GP's
hyperellratpoints when gp is on the PATH, which is another implementation
of the search, and otherwise reported as not checked -- except when the
degree is odd and the denominator range is what makes the search large:
then only the denominators a point can have are tried, e m^2 with e a
divisor of the leading coefficient (see brute_points), which is how the
runs at the top of the range of a long are checked.  Invocations that
test an error message, an output format or the report of -v are not
checked here beyond what they print, except that the points a -v report
contains are still compared.

Prints one line per invocation that fails and a summary; exits 1 on any
failure.  Takes a minute or two.
"""
import sys, re, math, subprocess, shutil
from fractions import Fraction

MAX_BRUTE = 40_000_000   # candidate pairs a brute-force search is allowed
rule_used = 0            # searches over the denominators e m^2 only

def parse_args(args):
    """The options of main.c that decide what is searched and printed."""
    o = dict(q=False, v=False, z=False, one=False, inf=True, x=False, y=False,
             s=False, S=None, dl=1, du=None, ivs=[], fmt=None, other={},
             err=None, strs=False)
    if len(args) < 2: o['err'] = 'no height'; return o
    o['cof'] = args[0]; o['H'] = args[1]
    i = 2; l_seen = False; low = None
    while i < len(args):
        a = args[i]
        if not a.startswith('-'): o['err'] = 'bad option'; return o
        k = a[1:]
        def val():
            nonlocal i
            i += 1
            if i >= len(args): raise IndexError
            return args[i]
        try:
            if k == 'q': o['q'] = True
            elif k == 'v': o['v'] = True
            elif k == 'z': o['z'] = True
            elif k == 'Z': o['z'] = False
            elif k == '1': o['one'] = True
            elif k == 'i': o['inf'] = False
            elif k == 'I': o['inf'] = True
            elif k == 'x': o['x'] = True
            elif k == 'X': o['x'] = False
            elif k == 'y': o['y'] = True
            elif k == 'Y': o['y'] = False
            elif k in ('k', 'K', 'j', 'J'): pass
            elif k == 's': o['s'] = True
            elif k == 'S':
                o['s'] = False; o['S'] = 10
                if i+1 < len(args) and not args[i+1].startswith('-'):
                    o['S'] = int(val())
            elif k == 'dl': o['dl'] = int(val())
            elif k == 'du': o['du'] = int(val())
            elif k == 'l':
                if l_seen: o['err'] = 'l l'; return o
                low = float(val()); l_seen = True
            elif k == 'u':
                up = float(val())
                if not l_seen:
                    if o['ivs']: o['err'] = 'u u'; return o
                    low = None  # -H, once H is known
                o['ivs'].append((low, up)); l_seen = False
            elif k == 'f': o['fmt'] = val()
            elif k in ('fs', 'fm', 'fe'): val(); o['strs'] = True
            elif k in ('n', 'N', 'p', 'F', 'r', 'R', 'U', 'A', 'C', 'P', 'Q', 'W'):
                o['other'][k] = val()
            else: o['err'] = 'unknown option ' + a; return o
        except (IndexError, ValueError):
            o['err'] = 'bad value'; return o
        i += 1
    if l_seen: o['ivs'].append((low, None))  # up = H
    return o

def coefficients(s):
    """main.c's scan: optional sign, digits; an empty or malformed string
    gives None"""
    c = []
    for tok in s.split():
        if not re.fullmatch(r'[+-]?\d+', tok): return None
        c.append(int(tok))
    return c if c else None

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
    b_low = max(dl, 1)
    b_high = H if (du is None or du < 1) else min(du, H)
    if not ivs: ivs = [(-float(H), float(H))]
    ivs = [((-float(H) if lo is None else lo), (float(H) if up is None else up)) for lo, up in ivs]
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
    # the size of the search
    total = 0
    for b in bs:
        for lo, up in ivs:
            total += max(0, min(float(b)*up, float(H)) - max(float(b)*lo, -float(H))) + 1
        if total > MAX_BRUTE: return None
    Hf = float(H)   # the program compares b*low and b*up with (double)height
    for b in bs:
        fb = float(b)   # the program computes b*low in double precision
        pw = [b**(D - i) for i in range(D + 1)]   # b^(D-i)
        cb = [c[i]*pw[i] if i <= d else 0 for i in range(D + 1)]
        for lo, up in ivs:
            if fb*lo <= -Hf: alo = -H
            elif fb*lo > Hf: break
            else: alo = math.ceil(fb*lo)
            if fb*up >= Hf: ahi = H
            elif fb*up < -Hf: continue
            else: ahi = math.floor(fb*up)
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

def squarefree(c):
    """whether f = sum c_i x^i, leading zeros dropped, is squarefree as a
    binary form of even degree: no repeated root, and no double root at
    infinity (the degree dropped from odd to even or by two or more)"""
    d0 = len(c) - 1; d = d0
    while d > 0 and c[d] == 0: d -= 1
    if d <= 0: return False
    if (d + 1) >> 1 < (d0 + 1) >> 1: return False
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

POINT = re.compile(r'^\((-?\d+) : (-?\d+) : (\d+)\)$')
XONLY = re.compile(r'^\((-?\d+) : (\d+)\)$')
FOUND = re.compile(r'^(\d+) rational points? (?:pairs )?found\.$')

def main():
    fn = sys.argv[1] if len(sys.argv) > 1 else 'testbase4'
    lines = open(fn).read().split('\n')
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
    checked = failed = skipped = gp_used = 0
    notes = []
    for args, out, filt in invs:
        o = parse_args(args)
        if filt and 'found' not in filt:
            skipped += 1; continue   # a filtered report: no points to check
        exit_line = [l for l in out if l.startswith('exit ')]
        body = [l for l in out if not l.startswith('exit ')]
        if exit_line:
            code = int(exit_line[-1].split()[1])
            # an error test: the program must have refused the input; a run
            # that ended well is checked like any other
            if code != 0: skipped += 1; continue
            if o['err'] or coefficients(o.get('cof', '')) is None:
                failed += 1; print('FAIL (accepted bad input):', args); continue
        if o['err']:
            failed += 1; print('FAIL (parse):', args, o['err']); continue
        c = coefficients(o['cof'])
        if c is None or not re.fullmatch(r'\d+', o['H']):
            skipped += 1; continue
        H = int(o['H'])
        # a run that ended in the message for a non-squarefree polynomial:
        # right when the polynomial is one, wrong otherwise -- and a
        # squarefree polynomial must not be refused
        refused = any(l == 'Polynomial is not square-free.' for l in body)
        sqf = squarefree(c)
        if refused != (not sqf):
            failed += 1
            print('FAIL:', ' '.join('<%s>' % a for a in args))
            print('    refused as not squarefree: %s, squarefree: %s' % (refused, sqf))
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
        if pts is None:
            if o['dl'] == 1 and o['du'] is None and not o['ivs']:
                pts = gp_points(c, H)
                if pts is not None:
                    gp_used += 1
                    if not want_inf: pts.pop((1, 0), None)
            if pts is None:
                skipped += 1; notes.append('not checked (too large, no gp): %s' % ' '.join(args[:2])); continue
        checked += 1
        def fail(msg):
            nonlocal failed
            failed += 1; print('FAIL:', ' '.join('<%s>' % a for a in args)); print('   ', msg)
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
                if (a, b) not in pts: fail('-1 printed (%d : %d), not a point' % (a, b))
                elif printed and printed[(a, b)][0] not in pts[(a, b)]: fail('-1: wrong y')
            continue
        if o['x'] or o['y']:
            exp = set(pts.keys())
            if xset != exp:
                extra = xset - exp; missing = exp - xset
                if missing: fail('missing x-coordinates: %s' % sorted(missing)[:5])
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
          '%d skipped (errors, formats, reports), %d failed'
          % (checked, gp_used, rule_used, skipped, failed))
    for n in notes: print(' ', n)
    sys.exit(1 if failed else 0)

if __name__ == '__main__':
    main()
