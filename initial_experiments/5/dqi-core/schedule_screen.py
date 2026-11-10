"""Exact local count screen and DP over degree-feasible block schedules."""
from collections import Counter
import json
from pathlib import Path
import time
from import_published import load
from blocked_monic_euclid import apply as block
from monomial_serial_euclid import apply as serial
from components.trimmed_middle import TrimmedArithmetic
from components.budget_inverse import bind

HERE=Path(__file__).resolve().parent

def component(k,size,kind,inverse=False):
    m=load();b=m.Builder();W=[b.qubits(12) for _ in range(70)];flags=b.qubits(k+size)
    class CountEmitter(m.Em):
        T=0
        def emit_mac(self,rt):self.T+=sum(t[0]=='P' for t in rt)
        def emit_free(self,rt):pass
    E=CountEmitter(b,*b.qubits(3),m.Alloc(b));m.inv_word=bind(m)
    math=TrimmedArithmetic(m,E);math.virtual_zero=True
    args=(m,E,math,W[:k-1],W[k-1:2*k-1],W[2*k-1:-1],W[-1])
    if kind=='serial':serial(*args,flags[k],inverse,linear_padding=True)
    else:block(*args,size,inverse,flags=flags[k:],monomial=True,suffix=True,linear_padding=True)
    return dict(k=k,size=size,kind=kind,inverse=inverse,T=E.T+b._cur_stats['toffoli'],Q=b.peak_qubits)

def main():
    start=time.time();rows=[]
    for k in range(1,35):
        for size,kind in [(1,'serial'),(1,'block'),(2,'block'),(3,'block')]:
            if k+size>35:continue
            if size==3 and k>19:continue
            f=component(k,size,kind);r=component(k,size,kind,True)
            row=dict(k=k,size=size,kind=kind,forward_T=f['T'],reverse_T=r['T'],Q=max(f['Q'],r['Q']))
            rows.append(row)
        print(json.dumps(dict(k=k,rows=rows[-4:])),flush=True)
    # A larger head's minimum peak already exceeds the selected limit.
    wide=component(1,4,'block');rows.append(dict(wide_screen=wide))
    costs={35:(0,[])}
    for k in reversed(range(1,35)):
        options=[]
        for row in rows:
            if row.get('k')!=k or row['Q']>961:continue
            tail,path=costs[k+row['size']]
            options.append((row['forward_T']+row['reverse_T']+tail,[row]+path))
        costs[k]=min(options,key=lambda x:x[0])
    selected=json.loads((HERE/'results/attribution_b3_s22_trimmed_retain.json').read_text())
    fixed=sum(c['T'] for c in selected['calls'] if c['name'].endswith('/first') or c['name'].endswith('/scale'))
    report=dict(status='COUNTED',scope='exact emitted-product local components and optimal schedule within tested serial/b1/b2/b3 family; no full decoder acceptance',
                limit_Q=961,optimal_core_T=costs[1][0]+fixed,fixed_first_and_scale_T=fixed,
                optimal_schedule=costs[1][1],rows=rows,seconds=time.time()-start)
    (HERE/'results/schedule_screen.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['rows','optimal_schedule']}),flush=True)

if __name__=='__main__':main()
