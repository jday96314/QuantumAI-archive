import json
from pathlib import Path
import numpy as np
from import_published import load
from gold_arithmetic import GoldArithmetic
from native_tools import finish,check
m=load();b=m.Builder();X=b.qubits(12);Y=b.qubits(12);E=m.Em(b,*b.qubits(3),m.Alloc(b));math=GoldArithmetic(m,E)
math.product(math.spec([(1,X,1)]),math.spec([(1,X,2)]),lambda tb:[Y[i] for i in tb])
c=finish(b,[('X',X),('Y',Y)],'cube28_kernel');rng=np.random.default_rng(142026);x=np.arange(4096);y=rng.integers(0,4096,4096);want=y^np.array([m.gpow(int(v),3) for v in x]);reports=[]
for seed in [142026,142027]:
 r=check(c,dict(X=x,Y=y),dict(X=x,Y=want),4096,seed);reports.append(r)
 assert r['status']=='PASS',r
assert c.operation_statistics().num_toffoli_like==28
res=dict(status='PASS',scope='all 4096 inputs, arbitrary target word, two measurement seeds, exact phase and scratch restoration',shots=8192,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,results=reports)
Path('results/cube28_native.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res),flush=True)
