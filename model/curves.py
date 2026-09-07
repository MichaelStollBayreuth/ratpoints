"""Parse the example curve files and predict the sieve's survivor rates.

The prediction mirrors sieving_info() in find_points.c: for each odd prime p,
  np = #{ a mod p : f(a) is a square mod p (0 counts) }
  inf = (deg f odd) or (leading coefficient is a square mod p)
  r   = (np*(p-1) + p)/p^2   if inf, else np/p
primes with np == p carry no information and are dropped; the rest are sorted
by increasing r, and ratpoints uses the first sp1 of them in phase 1 and the
first sp2 altogether.  So the predicted survivor rate after k primes is the
product of the k smallest r.
"""
import re, os

PRIMES = [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,
          97,101,103,107,109,113,127]

def parse_poly(txt):
    """'9*x^6 + 72*x^5 - 306*x^3 + 9' -> [9,0,0,-306,0,72,9] (ascending)."""
    t = txt.replace(' ', '')
    if t.startswith('(') and t.endswith(')'): t = t[1:-1]
    if not re.fullmatch(r'[-+0-9x*^]+', t) or 'x' not in t: return None
    terms = re.findall(r'[-+]?[^-+]+', t)
    coeffs = {}
    for term in terms:
        m = re.fullmatch(r'([-+]?)(\d*)\*?(?:x(?:\^(\d+))?)?', term)
        if not m: return None
        sign, num, exp = m.groups()
        if 'x' not in term: k = 0
        else: k = int(exp) if exp else 1
        c = int(num) if num else 1
        if sign == '-': c = -c
        coeffs[k] = coeffs.get(k, 0) + c
    d = max(coeffs)
    if d < 3 or d > 6: return None
    return [coeffs.get(i, 0) for i in range(d+1)]

def read_curves(path, tag):
    # some polynomials are wrapped over two lines, the first ending in + or -
    raw, pending, start = [], '', 0
    for ln, line in enumerate(open(path), 1):
        t = line.split('//')[0].strip()
        if pending:
            raw.append((start, pending + t)); pending = ''; continue
        if t.endswith('+') or t.endswith('-'): pending, start = t, ln
        else: raw.append((ln, t))
    out = []
    for ln, line in raw:
        if not line or line.startswith('{') or ':' in line: continue
        if 'points' in line or 'ht' in line: continue
        c = parse_poly(line)
        if c: out.append(('%s:%d' % (tag, ln), c))
    return out

def densities(coeffs):
    """The sorted list of r values, as sieving_info computes them."""
    deg = len(coeffs) - 1
    rs = []
    for p in PRIMES:
        sq = set((x*x) % p for x in range(p))          # squares incl. 0
        np_ = sum(1 for a in range(p)
                  if sum(c * pow(a, i, p) for i, c in enumerate(coeffs)) % p in sq)
        if np_ == p and (deg % 2 == 1 or coeffs[deg] % p in sq):
            continue                                    # no information
        inf = (deg % 2 == 1) or (coeffs[deg] % p in sq)
        r = (np_*(p-1) + p)/float(p*p) if inf else np_/float(p)
        rs.append((r, p))
    rs.sort()
    return rs

def rate(rs, k):
    x = 1.0
    for r, _ in rs[:k]: x *= r
    return x

if __name__ == '__main__':
    base = '/home/mstoll/software/git/ratpoints/examples'
    cs = read_curves(os.path.join(base, 'examples'), 'ex')
    cs += read_curves(os.path.join(base, 'best-curves'), 'bc')
    print('parsed %d curves' % len(cs))
    rows = []
    for name, c in cs:
        rs = densities(c)
        if len(rs) < 20: continue
        rows.append((rate(rs, 11), name, c, rs))
    rows.sort()
    print('predicted survivor rate after 11 primes:')
    print('  min %.3e   median %.3e   max %.3e'
          % (rows[0][0], rows[len(rows)//2][0], rows[-1][0]))
    import pickle
    pickle.dump(rows, open('/tmp/claude-1000/-home-mstoll-software-git-ratpoints/83b9b35a-2986-4d7f-a258-6c32978b4583/scratchpad/curves.pkl','wb'))
