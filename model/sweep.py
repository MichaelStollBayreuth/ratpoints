"""Grid over curve x height x (sp1, sp2), recording the survivor rates and
the cycles of each of the three stages.  Counts come from the deterministic
RP_PHASE_COUNTS build (one run); cycles from the RP_PHASE_TIMING build under
perf, pinned to CPU 0, minimum of two runs."""
import pickle, subprocess, sys, re, math, csv, time
S='/tmp/claude-1000/-home-mstoll-software-git-ratpoints/83b9b35a-2986-4d7f-a258-6c32978b4583/scratchpad'
sys.path.insert(0,S); from curves import rate
good = pickle.load(open(S+'/good.pkl','rb'))
# nine curves, log-spaced in the predicted survivor rate
idx = [int(round(i*(len(good)-1)/8.0)) for i in range(9)]
picks = [good[i] for i in idx]

HEIGHTS = [25000, 75000, 150000]
COMBOS  = [(sp1, sp2) for sp1 in (8, 11, 14, 17) for sp2 in (14, 18, 22, 26)
           if sp2 >= sp1]

def counts(arg, h, sp1, sp2):
    p = subprocess.run([S+'/bin/rp256-count', arg, str(h), '-n', str(sp1),
                        '-N', str(sp2), '-q'], capture_output=True, text=True,
                       timeout=300)
    m = dict(re.findall(r'(\w+)=(-?\d+)', p.stderr))
    return {k: int(v) for k, v in m.items()} if 'arrays' in m else None

def cycles(arg, h, sp1, sp2):
    """(perf core cycles, phasedata dict) for the faster of two runs."""
    best = None
    for _ in range(2):
        p = subprocess.run(['perf','stat','-x,','-e','cpu_core/cycles/',
                            'taskset','-c','0', S+'/bin/rp256-time', arg,
                            str(h), '-n', str(sp1), '-N', str(sp2), '-q'],
                           capture_output=True, text=True, timeout=300)
        c = re.search(r'^(\d+),,cpu_core/cycles', p.stderr, re.M)
        d = dict(re.findall(r'(\w+)=(-?\d+)', p.stderr))
        if not c or 'cyctot' not in d: continue
        d = {k:int(v) for k,v in d.items()}
        if best is None or int(c.group(1)) < best[0]: best = (int(c.group(1)), d)
    return best

out = open(S+'/sweep.csv','w', newline='')
w = csv.writer(out)
w.writerow(['curve','f','R1_pred','R2_pred','height','sp1','sp2','arrays',
            'bits_in','bits_1','bits_2','units_1','and2','ext2','checks',
            'cyc1','cyc2','cyc3','cyctot','perf_cycles'])
t0 = time.time(); n = 0
for r11, name, c, rs, _, _ in picks:
    arg = ' '.join(map(str, c))
    for h in HEIGHTS:
        for sp1, sp2 in COMBOS:
            n += 1
            try:
                cnt = counts(arg, h, sp1, sp2)
                cyc = cycles(arg, h, sp1, sp2)
            except subprocess.TimeoutExpired:
                print('timeout %s h=%d %d/%d' % (name,h,sp1,sp2)); continue
            if not cnt or not cyc: continue
            pc, d = cyc
            w.writerow([name, arg, rate(rs,sp1), rate(rs,sp2), h, sp1, sp2,
                        cnt['arrays'], cnt['bits_in'], cnt['bits_1'],
                        cnt['bits_2'], cnt['units_1'], cnt['and2'],
                        cnt['ext2'], cnt['checks'],
                        d['cyc1'], d['cyc2'], d['cyc3'], d['cyctot'], pc])
            out.flush()
    print('%s done (%d configs, %.0f s elapsed)' % (name, n, time.time()-t0))
out.close(); print('finished %d configs in %.0f s' % (n, time.time()-t0))
