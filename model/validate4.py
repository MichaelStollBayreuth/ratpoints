"""Validate the cost model against the random-curve sample in testdata.h.

The strong test: fit the cost coefficients on the point-rich curves from
examples/ only, then use them to choose (sp1, sp2) for random curves that
the fit has never seen -- a different population, not just a held-out curve.
"""
import csv, sys, numpy as np
S='/tmp/claude-1000/-home-mstoll-software-git-ratpoints/83b9b35a-2986-4d7f-a258-6c32978b4583/scratchpad'
sys.path.insert(0,S); from curves import densities
W = 256
def mbits(r):
    """mean number of set bits per bit-array on entry to phase 1; this is
    W for which_bits == num_all and about W/2 otherwise, and ratpoints
    knows which case it is before sieving."""
    return r['bits_in']/r['arrays']

def load(path, callmap=None):
    rows = list(csv.DictReader(open(path)))
    for r in rows:
        for k,v in r.items():
            if k not in ('curve','f'): r[k] = float(v)
        if 'calls' not in r and callmap: r['calls'] = callmap[(r['curve'],r['height'])]
        for i in (1,2,3): r['c%d'%i] = r['perf_cycles']*r['cyc%d'%i]/r['cyctot']
        r['ctot'] = r['c1']+r['c2']+r['c3']
    return rows

rich = load(S+'/sweep2.csv')
cm = {(r['curve'],r['height']): r['calls'] for r in rich}
rich += [r for r in load(S+'/sweep.csv', cm) if (r['curve'],r['height']) in cm]
rand = load(S+'/sweep4.csv')
print('point-rich sample: %d rows, %d curves' % (len(rich), len(set(r['curve'] for r in rich))))
print('random sample:     %d rows, %d curves' % (len(rand), len(set(r['curve'] for r in rand))))

D = {}
def Rk(r):
    if r['curve'] not in D:
        v=[1.0]
        for d,_ in densities([int(x) for x in r['f'].split()]): v.append(v[-1]*d)
        D[r['curve']] = v
    return D[r['curve']]
def R(r,k):
    v = Rk(r); return v[min(int(k), len(v)-1)]

def floors(rows):
    f = {}
    for r in sorted(rows, key=lambda r: r['sp2']):
        f[(r['curve'],r['height'])] = max(0.0, r['bits_2'] - r['bits_in']*R(r,r['sp2']))
    return f
NAMES = ['calls','arrays','sp1arr','units_1','and2','bits_2','checks']
def regs(r, F, rho, sp1=None, sp2=None):
    sp1 = r['sp1'] if sp1 is None else sp1; sp2 = r['sp2'] if sp2 is None else sp2
    b2 = r['bits_in']*R(r,sp2) + F
    return {'calls': r['calls'], 'arrays': r['arrays'], 'sp1arr': r['arrays']*sp1,
            'units_1': r['arrays']*(1-(1-R(r,sp1))**mbits(r)),
            'and2': r['arrays']*sum(1-(1-R(r,sp1+j))**mbits(r) for j in range(int(sp2-sp1))),
            'bits_2': b2, 'checks': rho*b2}
def fit_all(rows, F, rho):
    coef = {}
    for y, ns in [(lambda r: r['c1'], ['calls','arrays','sp1arr']),
                  (lambda r: r['c2'], ['calls','arrays','units_1','and2','bits_2']),
                  (lambda r: r['c3'], ['checks'])]:
        A = np.array([[regs(r,F[(r['curve'],r['height'])],rho)[n] for n in ns] for r in rows])
        b = np.array([y(r) for r in rows]); keep=list(range(len(ns)))
        for _ in range(len(ns)):
            x,*_ = np.linalg.lstsq(A[:,keep],b,rcond=None)
            if (x>=0).all(): break
            keep=[k for k,xi in zip(keep,x) if xi>0]
            if not keep: break
        for k,xi in zip(keep,x): coef[ns[k]] = coef.get(ns[k],0.0)+xi
    return coef

Fr, Fn = floors(rich), floors(rand)
rho_r = float(np.median([r['checks']/r['bits_2'] for r in rich if r['bits_2']>100]))
rho_n = float(np.median([r['checks']/r['bits_2'] for r in rand if r['bits_2']>100]))
print('\nfraction of candidates in lowest terms: %.3f (point-rich) vs %.3f (random)'
      % (rho_r, rho_n))

print('\n=== do the a-priori predictors hold on random curves? ===')
for lab, meas, pred in [
    ('bits after phase 1', lambda r: r['bits_1'], lambda r: r['bits_in']*R(r,r['sp1'])),
    ('units entering phase 2', lambda r: r['units_1'],
     lambda r: r['arrays']*(1-(1-R(r,r['sp1']))**mbits(r))),
    ('AND steps in phase 2', lambda r: r['and2'],
     lambda r: r['arrays']*sum(1-(1-R(r,r['sp1']+j))**mbits(r) for j in range(int(r['sp2']-r['sp1']))))]:
    a=np.array([meas(r) for r in rand]); b=np.array([pred(r) for r in rand]); m=a>0
    print('  %-24s predicted/measured: median %.3f (10th %.3f, 90th %.3f)'
          % (lab, np.median(b[m]/a[m]), np.percentile(b[m]/a[m],10), np.percentile(b[m]/a[m],90)))

COEF_RICH = fit_all(rich, Fr, rho_r)
COEF_RAND = fit_all(rand, Fn, rho_n)
print('\n=== fitted costs: does the machine behave the same on both samples? ===')
print('  %-10s %14s %14s' % ('', 'point-rich', 'random'))
for n in NAMES:
    print('  %-10s %14.2f %14.2f' % (n, COEF_RICH.get(n,0), COEF_RAND.get(n,0)))

groups = {}
for r in rand: groups.setdefault((r['curve'],r['height']), []).append(r)
def report(label, coef, rho):
    exc=[]
    for (cur,H),g in sorted(groups.items()):
        F = Fn[(cur,H)]
        best = min(g, key=lambda r: r['ctot'])
        pick = min(g, key=lambda r: sum(coef.get(n,0)*regs(r,F,rho)[n] for n in NAMES))
        exc.append(100*(pick['ctot']/best['ctot']-1))
    print('  %-42s median %+5.1f%%  mean %+5.1f%%  worst %+6.1f%%'
          % (label, np.median(exc), np.mean(exc), max(exc)))
    return exc
print('\n=== choosing (sp1, sp2) for random curves ===')
report('model, coefficients from the POINT-RICH curves', COEF_RICH, rho_r)
report('model, coefficients from the random curves', COEF_RAND, rho_n)
dexc=[]
for (cur,H),g in sorted(groups.items()):
    best=min(g,key=lambda r:r['ctot'])
    d=[r for r in g if r['sp1']==11 and r['sp2']==18]
    if d: dexc.append(100*(d[0]['ctot']/best['ctot']-1))
print('  %-42s median %+5.1f%%  mean %+5.1f%%  worst %+6.1f%%'
      % ('the shipped default 11/19 (grid 11/18)', np.median(dexc), np.mean(dexc), max(dexc)))

print('\n=== where do the optima sit for random curves? ===')
from collections import Counter
opt = Counter()
for (cur,H),g in sorted(groups.items()):
    b = min(g, key=lambda r: r['ctot']); opt[(int(b['sp1']),int(b['sp2']))] += 1
for k,v in sorted(opt.items()): print('  sp1/sp2 = %2d/%-2d : %d of %d' % (k[0],k[1],v,len(groups)))

print('\n=== which coordinate does the model get wrong? ===')
from collections import Counter
c1 = Counter(); c2 = Counter()
for (cur,H),g in sorted(groups.items()):
    F = Fn[(cur,H)]
    best = min(g, key=lambda r: r['ctot'])
    pick = min(g, key=lambda r: sum(COEF_RICH.get(n,0)*regs(r,F,rho_r)[n] for n in NAMES))
    c1[(int(best['sp1']), int(pick['sp1']))] += 1
    c2[(int(best['sp2']), int(pick['sp2']))] += 1
print('  sp1  best -> model : ' + ', '.join('%d->%d x%d' % (a,b,n) for (a,b),n in sorted(c1.items())))
print('  sp2  best -> model : ' + ', '.join('%d->%d x%d' % (a,b,n) for (a,b),n in sorted(c2.items())))

print('\n=== a hybrid: keep the default unless the model predicts a big win ===')
def hybrid(rows_, floors_, rho, thresh):
    gs = {}
    for r in rows_: gs.setdefault((r['curve'],r['height']), []).append(r)
    exc = []
    for (cur,H),g in sorted(gs.items()):
        F = floors_[(cur,H)]
        cost = lambda r: sum(COEF_RICH.get(n,0)*regs(r,F,rho)[n] for n in NAMES)
        best = min(g, key=lambda r: r['ctot'])
        pick = min(g, key=cost)
        dflt = [r for r in g if r['sp1']==11 and r['sp2']==18]
        if dflt and cost(pick) > (1-thresh)*cost(dflt[0]): pick = dflt[0]
        exc.append(100*(pick['ctot']/best['ctot']-1))
    return exc
for th in (0.0, 0.10, 0.20, 0.30):
    a = hybrid(rand, Fn, rho_n, th); b = hybrid(rich, Fr, rho_r, th)
    print('  keep default unless model predicts >%2d%% gain:  random %+5.1f%% / %+5.1f%%   point-rich %+5.1f%% / %+5.1f%%   (median/mean)'
          % (100*th, np.median(a), np.mean(a), np.median(b), np.mean(b)))
ex_r = []
for (cur,H),g in sorted(groups.items()):
    best=min(g,key=lambda r:r['ctot']); d=[r for r in g if r['sp1']==11 and r['sp2']==18]
    if d: ex_r.append(100*(d[0]['ctot']/best['ctot']-1))
gs={}
for r in rich: gs.setdefault((r['curve'],r['height']),[]).append(r)
ex_p=[]
for (cur,H),g in sorted(gs.items()):
    best=min(g,key=lambda r:r['ctot']); d=[r for r in g if r['sp1']==11 and r['sp2']==18]
    if d: ex_p.append(100*(d[0]['ctot']/best['ctot']-1))
print('  the plain default everywhere:                    random %+5.1f%% / %+5.1f%%   point-rich %+5.1f%% / %+5.1f%%'
      % (np.median(ex_r), np.mean(ex_r), np.median(ex_p), np.mean(ex_p)))

print('\n=== use the model for sp1 only, with a simple rule for sp2 ===')
def policy(rows_, floors_, rho, sp2rule):
    gs={}
    for r in rows_: gs.setdefault((r['curve'],r['height']),[]).append(r)
    exc=[]
    for (cur,H),g in sorted(gs.items()):
        F=floors_[(cur,H)]
        cost=lambda r: sum(COEF_RICH.get(n,0)*regs(r,F,rho)[n] for n in NAMES)
        # the model's sp1: minimise with sp2 held at the default
        cands=[r for r in g if r['sp2']==18] or g
        s1=int(min(cands,key=cost)['sp1'])
        s2=sp2rule(s1)
        pick=min(g,key=lambda r: abs(r['sp1']-s1)*100+abs(r['sp2']-s2))
        best=min(g,key=lambda r:r['ctot'])
        exc.append(100*(pick['ctot']/best['ctot']-1))
    return exc
for lab, rule in [('sp2 = 19 (the default)', lambda s: 19),
                  ('sp2 = sp1 + 8', lambda s: s+8),
                  ('sp2 = sp1 + 11', lambda s: s+11),
                  ('sp2 = max(19, sp1 + 8)', lambda s: max(19, s+8))]:
    a=policy(rand,Fn,rho_n,rule); b=policy(rich,Fr,rho_r,rule)
    print('  model sp1, %-24s random %+5.1f%% / %+5.1f%%   point-rich %+5.1f%% / %+5.1f%%'
          % (lab, np.median(a), np.mean(a), np.median(b), np.mean(b)))
