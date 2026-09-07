"""With a fine sp1 scan: locate the optimum by a parabolic fit through the
three cheapest points, and report what has survived phase 1 there."""
import csv, numpy as np
S='/tmp/claude-1000/-home-mstoll-software-git-ratpoints/83b9b35a-2986-4d7f-a258-6c32978b4583/scratchpad'
rows=list(csv.DictReader(open(S+'/sweep5.csv')))
for r in rows:
    for k,v in r.items():
        if k not in ('pop','curve','f'): r[k]=float(v)
    r['ctot']=sum(r['perf_cycles']*r['cyc%d'%i]/r['cyctot'] for i in (1,2,3))
g={}
for r in rows: g.setdefault((r['pop'],r['curve']),[]).append(r)

print('  %-4s %-9s %5s %7s %7s   %10s %8s   %8s' %
      ('pop','curve','m','sp1*','fine','R(sp1*)','u(sp1*)','flatness'))
out={}
for (pop,cur),rs in sorted(g.items()):
    rs.sort(key=lambda r: r['sp1'])
    c=np.array([r['ctot'] for r in rs]); s=np.array([r['sp1'] for r in rs])
    i=int(np.argmin(c))
    # parabola through the minimum and its neighbours -> sub-integer optimum
    if 0 < i < len(c)-1:
        d=c[i-1]-2*c[i]+c[i+1]
        sub=s[i] + (0.5*(c[i-1]-c[i+1])/d if d>0 else 0.0)
    else: sub=s[i]
    r0=rs[i]
    R=r0['bits_1']/r0['bits_in']; u=r0['units_1']/r0['arrays']
    # how flat: how many sp1 values are within 2% of the minimum
    flat=int((c <= 1.02*c[i]).sum())
    print('  %-4s %-9s %5.0f %7d %7.1f   %10.2e %8.4f   %d values within 2%%' %
          (pop,cur,r0['bits_in']/r0['arrays'],s[i],sub,R,u,flat))
    out.setdefault(pop,[]).append((R,u,sub))
print()
for pop,lab in (('rich','point-rich'),('rand','random')):
    v=out[pop]; R=np.array([x[0] for x in v]); u=np.array([x[1] for x in v])
    print('  %-11s surviving bits R: median %.2e  (%.1e .. %.1e, factor %.1f)'
          % (lab, np.median(R), R.min(), R.max(), R.max()/R.min()))
    print('  %-11s non-empty arrays u: median %.4f  (%.4f .. %.4f, factor %.1f)'
          % ('', np.median(u), u.min(), u.max(), u.max()/u.min()))
allv=out['rich']+out['rand']
R=np.array([x[0] for x in allv]); u=np.array([x[1] for x in allv])
print('\n  ALL 12 curves, spanning a factor %.0f in R(11):' % 630)
print('    surviving BITS   at the optimum : median %.2e, spread factor %.1f' % (np.median(R), R.max()/R.min()))
print('    non-empty ARRAYS at the optimum : median %.4f, spread factor %.1f' % (np.median(u), u.max()/u.min()))
