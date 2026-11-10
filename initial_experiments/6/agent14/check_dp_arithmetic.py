"""Independent exact polynomial checks of optimized split products."""
import json
from pathlib import Path
import numpy as np
from triangular_dp import DPArithmetic
from import_published import load
from native_tools import finish,check

def run(n,degree,q,start,holes):
 m=load();rng=np.random.default_rng(14553+n+q+start);shots=256
 av=rng.integers(0,4096,(degree+1,shots));av[0]=1
 for i in holes:av[i]=0
 tv=rng.integers(0,4096,(n,shots));want=tv.copy()
 for shot in range(shots):
  for i in range(start,n):
   for j in range(i):
    h=i-j
    if h%q==0 and h//q<=degree:want[i,shot]^=m.gmul(m.gpow(int(av[h//q,shot]),q),int(tv[j,shot]))
 reports=[]
 for mode in ['forward','inverse','both','inverse_both']:
  b=m.Builder();AA=[b.qubits(12) for _ in range(degree+1)];TT=[b.qubits(12) for _ in range(n)];E=m.Em(b,*b.qubits(3),m.Alloc(b));math=DPArithmetic(m,E)
  for i in holes:math.known_zero.add(tuple(AA[i]))
  for inv in ([False,True] if mode=='both' else [True,False] if mode=='inverse_both' else [mode=='inverse']):math.triangular_suffix(AA,q,TT,start,inv)
  c=finish(b,[(f'A{i}',w) for i,w in enumerate(AA)]+[(f'T{i}',w) for i,w in enumerate(TT)],f'dp_n{n}_d{degree}_q{q}_s{start}_{mode}')
  v0=want if mode in ['inverse','inverse_both'] else tv;v1=want if mode in ['forward','inverse_both'] else tv
  src={**{f'A{i}':v for i,v in enumerate(av)},**{f'T{i}':v for i,v in enumerate(v0)}};dst={**{f'A{i}':v for i,v in enumerate(av)},**{f'T{i}':v for i,v in enumerate(v1)}}
  r=check(c,src,dst,shots);r['mode']=mode;reports.append(r);print(json.dumps(r),flush=True)
 result=dict(n=n,degree=degree,q=q,start=start,holes=holes,shots=shots,circuits=reports,status='PASS' if all(r['status']=='PASS' for r in reports) else 'FAIL')
 assert result['status']=='PASS',result
 return result
rows=[run(*args) for args in [(13,4,1,0,[]),(17,7,1,5,[2,5]),(15,4,2,3,[1]),(11,2,4,0,[])]]
Path('results/dp_arithmetic_checks.json').write_text(json.dumps(dict(status='PASS',cases=rows),indent=2)+'\n')
