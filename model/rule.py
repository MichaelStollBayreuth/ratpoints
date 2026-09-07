"""Test the simplest possible rule: choose sp1 as the first k with R(k) <= T."""
import csv, sys, numpy as np
S='/tmp/claude-1000/-home-mstoll-software-git-ratpoints/83b9b35a-2986-4d7f-a258-6c32978b4583/scratchpad'
sys.path.insert(0,S); from curves import densities
def load(p, cm=None):
    rows=list(csv.DictReader(open(p)))
    for r in rows:
        for k,v in r.items():
            if k not in ('pop','curve','f'): r[k]=float(v)
        if 'calls' not in r and cm: r['calls']=cm[(r['curve'],r['height'])]
        r['ctot']=sum(r['perf_cycles']*r['cyc%d'%i]/r['cyctot'] for i in (1,2,3))
    return rows
rich=load(S+'/sweep2.csv')
cm={(r['curve'],r['height']):r['calls'] for r in rich}
rich+=[r for r in load(S+'/sweep.csv',cm) if (r['curve'],r['height']) in cm]
rand=load(S+'/sweep4.csv'); fine=load(S+'/sweep5.csv')

D={}
def R(r,k):
    if r['curve'] not in D:
        v=[1.0]
        for d,_ in densities([int(x) for x in r['f'].split()]): v.append(v[-1]*d)
        D[r['curve']]=v
    v=D[r['curve']]; return v[min(int(k),len(v)-1)]

def sp1_from(r, T, use_u):
    m = r['bits_in']/r['arrays']
    for k in range(1, 31):
        x = 1-(1-R(r,k))**m if use_u else R(r,k)
        if x <= T: return k
    return 30

def evaluate(rows, T, use_u, label):
    g={}
    for r in rows: g.setdefault((r['curve'],r['height']),[]).append(r)
    exc=[]
    for key,rs in sorted(g.items()):
        s1=sp1_from(rs[0],T,use_u)
        near=min(rs,key=lambda r:(abs(r['sp1']-s1), r['ctot']))
        cand=[r for r in rs if r['sp1']==near['sp1']]
        pick=min(cand,key=lambda r:r['ctot']) if len(cand)>1 else near
        # be honest: choose sp2 by the fixed rule, not by hindsight
        s2=max(19,int(near['sp1'])+8)
        pick=min(cand,key=lambda r:abs(r['sp2']-s2))
        best=min(rs,key=lambda r:r['ctot'])
        exc.append(100*(pick['ctot']/best['ctot']-1))
    return np.array(exc)

print('  rule: sp1 = first k with R(k) <= T,  sp2 = max(19, sp1+8)\n')
print('  %-12s %22s %22s' % ('T','point-rich (27 groups)','random (28 groups)'))
for T in (5e-5, 1e-4, 1.33e-4, 2e-4, 4e-4):
    a=evaluate(rich,T,False,''); b=evaluate(rand,T,False,'')
    print('  R <= %.2e %11.1f%% / %6.1f%% %14.1f%% / %6.1f%%'
          % (T, np.median(a), np.mean(a), np.median(b), np.mean(b)))
print('\n  for comparison, the u-based form (u = 1-(1-R)^m):')
for T in (0.013, 0.02, 0.026, 0.04):
    a=evaluate(rich,T,True,''); b=evaluate(rand,T,True,'')
    print('  u <= %.3f     %11.1f%% / %6.1f%% %14.1f%% / %6.1f%%'
          % (T, np.median(a), np.mean(a), np.median(b), np.mean(b)))
print('\n  (median / mean excess over the best setting in each grid)')
