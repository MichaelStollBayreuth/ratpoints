import re, sys, pickle, math
S='/tmp/claude-1000/-home-mstoll-software-git-ratpoints/83b9b35a-2986-4d7f-a258-6c32978b4583/scratchpad'
sys.path.insert(0,S); from curves import densities, rate
txt = open('/home/mstoll/software/git/ratpoints/testdata.h').read()
body = txt[txt.index('long testdata'):]
cs = [[int(x) for x in m.split(',')] for m in re.findall(r'\{(-?\d+(?:,-?\d+)+)\}', body)]
print('parsed %d curves, %d coefficients each' % (len(cs), len(cs[0])))
rows = []
for i, c in enumerate(cs):
    while len(c) > 1 and c[-1] == 0: c = c[:-1]        # drop a zero leading coeff
    if len(c) < 4: continue
    rs = densities(c)
    if len(rs) < 20: continue
    rows.append((rate(rs, 11), 'td:%d' % i, c, rs))
rows.sort()
print('predicted survivor rate after 11 primes, over %d random curves:' % len(rows))
import statistics
vals = [r[0] for r in rows]
for q in (0, 5, 25, 50, 75, 95, 100):
    print('   %3d%%  %.3e' % (q, vals[min(len(vals)-1, int(q*(len(vals)-1)/100))]))
pickle.dump(rows, open(S+'/testdata.pkl','wb'))
