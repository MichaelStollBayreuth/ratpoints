"""Fine scan in sp1 (step 1) along the ridge sp2 = min(30, sp1+8), to locate
the optimum precisely and see what has survived phase 1 there."""
import pickle, subprocess, sys, re, csv, time
S='/tmp/claude-1000/-home-mstoll-software-git-ratpoints/83b9b35a-2986-4d7f-a258-6c32978b4583/scratchpad'
sys.path.insert(0,S); from curves import rate
rich = pickle.load(open(S+'/good.pkl','rb'))          # point-rich, from examples/
rand = [r for r in pickle.load(open(S+'/testdata.pkl','rb')) if r[0] > 0]
lo, hi = int(0.05*len(rand)), int(0.95*len(rand))
picks  = [('rich', rich[i*(len(rich)-1)//5][1], rich[i*(len(rich)-1)//5][2],
           rich[i*(len(rich)-1)//5][3], 75000) for i in range(6)]
picks += [('rand', rand[lo+(hi-lo)*i//5][1], rand[lo+(hi-lo)*i//5][2],
           rand[lo+(hi-lo)*i//5][3], 150000) for i in range(6)]

def counts(arg,h,s1,s2):
    p=subprocess.run([S+'/bin/rp256-count',arg,str(h),'-n',str(s1),'-N',str(s2),'-q'],
                     capture_output=True,text=True,timeout=900)
    m=dict(re.findall(r'(\w+)=(-?\d+)',p.stderr))
    return {k:int(v) for k,v in m.items()} if 'arrays' in m else None
def cycles(arg,h,s1,s2,reps=3):
    best=None
    for _ in range(reps):
        p=subprocess.run(['perf','stat','-x,','-e','cpu_core/cycles/','taskset','-c','0',
                          S+'/bin/rp256-time',arg,str(h),'-n',str(s1),'-N',str(s2),'-q'],
                         capture_output=True,text=True,timeout=900)
        c=re.search(r'^(\d+),,cpu_core/cycles',p.stderr,re.M)
        d=dict(re.findall(r'(\w+)=(-?\d+)',p.stderr))
        if not c or 'cyctot' not in d: continue
        d={k:int(v) for k,v in d.items()}
        if best is None or int(c.group(1))<best[0]: best=(int(c.group(1)),d)
    return best

out=open(S+'/sweep5.csv','w',newline=''); w=csv.writer(out)
w.writerow(['pop','curve','f','R1_pred','height','sp1','sp2','arrays','calls',
            'bits_in','bits_1','bits_2','units_1','and2','ext2','checks',
            'cyc1','cyc2','cyc3','cyctot','perf_cycles'])
t0=time.time(); n=0
for pop,name,c,rs,h in picks:
    arg=' '.join(map(str,c))
    for s1 in range(6,31):
        s2=min(30,s1+8); n+=1
        try:
            cnt=counts(arg,h,s1,s2); cyc=cycles(arg,h,s1,s2)
        except subprocess.TimeoutExpired: continue
        if not cnt or not cyc: continue
        pc,d=cyc
        w.writerow([pop,name,arg,rate(rs,s1),h,s1,s2,cnt['arrays'],cnt['calls'],
                    cnt['bits_in'],cnt['bits_1'],cnt['bits_2'],cnt['units_1'],
                    cnt['and2'],cnt['ext2'],cnt['checks'],
                    d['cyc1'],d['cyc2'],d['cyc3'],d['cyctot'],pc]); out.flush()
    print('%s %s done (%d, %.0f s)'%(pop,name,n,time.time()-t0))
out.close(); print('finished %d in %.0f s'%(n,time.time()-t0))
