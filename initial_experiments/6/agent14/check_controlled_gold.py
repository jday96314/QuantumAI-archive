import json
from pathlib import Path
import numpy as np
from import_published import load
from gold_arithmetic import GoldArithmetic
from quantum_control import selected
from native_tools import finish,check
m=load();b=m.Builder();X=b.qubits(12);Y=b.qubits(12);flag=b.qubits(1);E=m.Em(b,*b.qubits(3),m.Alloc(b))
with selected(m,E,flag[0]) as (cm,ce):
 math=GoldArithmetic(cm,ce);math.product(math.spec([(1,X,1)]),math.spec([(1,X,2)]),lambda tb:[Y[i] for i in tb])
c=finish(b,[('X',X),('Y',Y),('F',flag)],'controlled_cube40_kernel');rng=np.random.default_rng(142028);x=np.tile(np.arange(4096),2);f=np.repeat([0,1],4096);y=rng.integers(0,4096,8192);want=y^np.array([m.gpow(int(v),3)*int(t) for v,t in zip(x,f)])
reports=[]
for seed in [142028,142029]:
 r=check(c,dict(X=x,Y=y,F=f),dict(X=x,Y=want,F=f),8192,seed);reports.append(r);assert r['status']=='PASS',r
assert c.operation_statistics().num_toffoli_like==40
res=dict(status='PASS',scope='all 4096 field inputs times both quantum controls; arbitrary target, two measurement seeds, exact phase and scratch restoration',shots=16384,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,results=reports)
Path('results/cube40_controlled_native.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res),flush=True)
