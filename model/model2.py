"""Cost model for the three stages, fitted and cross-validated.

Everything the model needs is known before sieving starts:
  calls    number of sift0 calls           (one per denominator interval)
  arrays   bit-arrays to sweep             (height and Sturm bounds)
  R(k)     product of the k smallest prime densities np/p, which
           sieving_info() already computes
  W        bits per bit-array
"""
import csv, sys, math, numpy as np
S='/tmp/claude-1000/-home-mstoll-software-git-ratpoints/83b9b35a-2986-4d7f-a258-6c32978b4583/scratchpad'
sys.path.insert(0,S); from curves import densities
W = 256

def load(path):
    rows = list(csv.DictReader(open(path)))
    for r in rows:
        for k, v in r.items():
            if k not in ('curve','f'): r[k] = float(v)
        for i in (1,2,3): r['c%d'%i] = r['perf_cycles']*r['cyc%d'%i]/r['cyctot']
        r['ctot'] = r['c1']+r['c2']+r['c3']
    return rows

rows = load(S+'/sweep2.csv')
# calls depends only on (curve, height); take it from sweep2 for sweep1's rows
callmap = {(r['curve'], r['height']): r['calls'] for r in rows}
try:
    old = [r for r in load(S+'/sweep.csv')
           if (r['curve'], r['height']) in callmap]
    for r in old: r['calls'] = callmap[(r['curve'], r['height'])]
    rows = rows + old
except IOError: pass
print('%d rows, %d curves, heights %s'
      % (len(rows), len(set(r['curve'] for r in rows)),
         sorted(set(int(r['height']) for r in rows))))

D = {}
for r in rows:
    if r['curve'] not in D:
        Rk = [1.0]
        for d,_ in densities([int(x) for x in r['f'].split()]): Rk.append(Rk[-1]*d)
        D[r['curve']] = Rk
def R(cur,k):
    Rk = D[cur]; return Rk[min(int(k), len(Rk)-1)]

# F: the floor of non-reduced representations (k*a, k*b) of genuine points,
# which no amount of sieving removes.  Estimated per (curve, height) from the
# most heavily sieved row available.
floor = {}
for r in sorted(rows, key=lambda r: r['sp2']):
    floor[(r['curve'],r['height'])] = max(0.0, r['bits_2'] -
                                          r['bits_in']*R(r['curve'], r['sp2']))
rho = float(np.median([r['checks']/r['bits_2'] for r in rows if r['bits_2']>100]))

NAMES = ['calls','arrays','sp1arr','units_1','and2','bits_2','checks']
def regs(r, sp1=None, sp2=None):
    cur = r['curve']; sp1 = r['sp1'] if sp1 is None else sp1
    sp2 = r['sp2'] if sp2 is None else sp2
    F = floor[(cur, r['height'])]
    R1 = R(cur, sp1)
    b2 = r['bits_in']*R(cur, sp2) + F
    return {'calls': r['calls'], 'arrays': r['arrays'],
            'sp1arr': r['arrays']*sp1,
            'units_1': r['arrays']*(1-(1-R1)**W),
            'and2': r['arrays']*sum(1-(1-R(cur, sp1+j))**W
                                    for j in range(int(sp2-sp1))),
            'bits_2': b2, 'checks': rho*b2}

def fit(y, names, data, verbose=None):
    A = np.array([[regs(r)[n] for n in names] for r in data])
    b = np.array([y(r) for r in data]); keep = list(range(len(names)))
    for _ in range(len(names)):
        x,*_ = np.linalg.lstsq(A[:,keep], b, rcond=None)
        if (x>=0).all(): break
        keep = [k for k,xi in zip(keep,x) if xi>0]
        if not keep: return {}
    if verbose:
        pred = A[:,keep]@x
        print('\n  %s:  R^2 = %.4f' % (verbose,
              1-((b-pred)**2).sum()/((b-b.mean())**2).sum()))
        for k,xi in zip(keep,x): print('      %-10s %12.2f cycles' % (names[k],xi))
    return dict(zip([names[k] for k in keep], x))

def fit_all(data, verbose=False):
    c = {}
    for lab, y, ns in [('stage 1', lambda r: r['c1'], ['calls','arrays','sp1arr']),
                       ('stage 2', lambda r: r['c2'],
                        ['calls','arrays','units_1','and2','bits_2']),
                       ('stage 3', lambda r: r['c3'], ['checks'])]:
        for k,v in fit(y, ns, data, lab if verbose else None).items():
            c[k] = c.get(k,0.0)+v
    return c

print('\nrho (fraction of candidates in lowest terms) = %.3f' % rho)
print('\n=== coefficients, fitted on all curves ===')
COEF = fit_all(rows, verbose=True)
pred = np.array([sum(COEF.get(n,0)*regs(r)[n] for n in NAMES) for r in rows])
meas = np.array([r['ctot'] for r in rows])
print('\n  total cost: R^2 = %.4f, median |rel err| %.1f%%, 90th %.1f%%'
      % (1-((meas-pred)**2).sum()/((meas-meas.mean())**2).sum(),
         100*np.median(np.abs(pred/meas-1)), 100*np.percentile(np.abs(pred/meas-1),90)))

print('\n=== leave-one-curve-out: does the model choose good parameters? ===')
groups = {}
for r in rows: groups.setdefault((r['curve'],r['height']), []).append(r)
print('  %-9s %7s %8s %8s %9s %9s' %
      ('curve','height','best','model','model cost','default'))
exc, dexc = [], []
for (cur,H), g in sorted(groups.items()):
    C = fit_all([r for r in rows if r['curve'] != cur])   # held out
    best = min(g, key=lambda r: r['ctot'])
    pick = min(g, key=lambda r: sum(C.get(n,0)*regs(r)[n] for n in NAMES))
    dflt = min(g, key=lambda r: abs(r['sp1']-11)+abs(r['sp2']-19))
    e = 100*(pick['ctot']/best['ctot']-1); exc.append(e)
    de = 100*(dflt['ctot']/best['ctot']-1); dexc.append(de)
    print('  %-9s %7d  %2d/%-2d    %2d/%-2d   %+8.1f%%  %+8.1f%%'
          % (cur, H, best['sp1'], best['sp2'], pick['sp1'], pick['sp2'], e, de))
print('\n  held-out model: median %+.1f%%, mean %+.1f%%, worst %+.1f%% above the best'
      % (np.median(exc), np.mean(exc), max(exc)))
print('  shipped default: median %+.1f%%, mean %+.1f%%, worst %+.1f%%'
      % (np.median(dexc), np.mean(dexc), max(dexc)))

# ---------------------------------------------------------------------------
print('\n=== the marginal rule (what sieving_info could actually do) ===')
B  = COEF.get('sp1arr',0); Cc = COEF.get('arrays',0); Dd = COEF.get('units_1',0)
E  = COEF.get('and2',0)
Fb = COEF.get('bits_2',0) + rho*COEF.get('checks',0)
print('  cost of one more phase-1 prime, per bit-array:      %.2f cycles' % B)
print('  cost of one AND step in phase 2:                    %.2f cycles' % E)
print('  cost of entering phase 2 with a non-empty unit:     %.2f cycles' % Dd)
print('  cost of one candidate surviving phase 2:            %.2f cycles' % Fb)
def u(cur,k): return 1-(1-R(cur,k))**W
def choose(cur, bits_in, arrays, F, kmax=30):
    """sp1: keep moving a prime into phase 1 while that is cheaper.
       sp2: keep adding a phase-2 prime while it removes more than it costs."""
    sp1 = 1
    while sp1 < kmax and B - Dd*(u(cur,sp1)-u(cur,sp1+1)) - E*u(cur,sp1) < 0:
        sp1 += 1
    sp2 = sp1
    while sp2 < kmax and E*u(cur,sp2)*arrays < Fb*bits_in*(R(cur,sp2)-R(cur,sp2+1)):
        sp2 += 1
    return sp1, sp2
print('\n  %-9s %7s %8s %9s %9s' % ('curve','height','best','rule','rule cost'))
rexc = []
for (cur,H), g in sorted(groups.items()):
    F = floor[(cur,H)]; r0 = g[0]
    s1, s2 = choose(cur, r0['bits_in'], r0['arrays'], F)
    near = min(g, key=lambda r: abs(r['sp1']-s1)+abs(r['sp2']-s2))
    best = min(g, key=lambda r: r['ctot'])
    e = 100*(near['ctot']/best['ctot']-1); rexc.append(e)
    print('  %-9s %7d  %2d/%-2d     %2d/%-2d    %+8.1f%%'
          % (cur, H, best['sp1'], best['sp2'], s1, s2, e))
print('\n  marginal rule: median %+.1f%%, mean %+.1f%%, worst %+.1f%% above the best'
      % (np.median(rexc), np.mean(rexc), max(rexc)))

print('\n=== what is the best single fixed choice? ===')
combos = sorted(set((int(r['sp1']),int(r['sp2'])) for r in rows))
score = []
for (s1,s2) in combos:
    ex = []
    for (cur,H), g in groups.items():
        m = [r for r in g if r['sp1']==s1 and r['sp2']==s2]
        b = min(g, key=lambda r: r['ctot'])['ctot']
        if m: ex.append(m[0]['ctot']/b)
    if len(ex) == len(groups): score.append((np.mean(ex), np.median(ex), s1, s2))
score.sort()
print('  %-8s %10s %10s' % ('sp1/sp2','mean','median'))
for mean, med, s1, s2 in score[:5]:
    print('  %2d/%-5d %9.2fx %9.2fx' % (s1, s2, mean, med))
d = [s for s in score if (s[2],s[3])==(11,18)]
if d: print('  %2d/%-5d %9.2fx %9.2fx   <- the shipped default' % (11,18,d[0][0],d[0][1]))
