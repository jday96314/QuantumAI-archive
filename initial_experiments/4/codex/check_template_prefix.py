"""Bind the classical streaming model to the actual compiled BM prefix values.

Only internal values are checked here, not the unfinished Wronskian decoder.
"""
import json,subprocess
from pathlib import Path
import numpy as np
import kickmix as km
from probe_stream_template import model
from probe_wronskian import p,mul
root=Path(__file__).resolve().parent;tag='stream10_offset3_select_noscale_compact'
hook=json.loads((root/f'parameters/{tag}_hook.json').read_text())
b=km.CircuitBuilder();syn=[b.create_quantum_register(12,f'S{i}') for i in range(70)]
c=b.finish_circuit();pre=root/'parameters/prefix.pre.kmx';pre.parent.mkdir(exist_ok=True);c.to_file(str(pre),format='kmx')
empty=root/'parameters/prefix.empty';empty.write_text('');phases=root/'parameters/prefix.phases';phases.write_text('')
path=root/'candidates/template_prefix.kmx'
subprocess.run([str(root/'translate_prefix'),str(root/f'../../../initial_experiments/3/parameters/{tag}_tape.npy'),str(hook['index']),'0',str(pre),str(empty),str(empty),str(path),str(hook['outcomes']),str(phases),str(hook['index'])],check=True)
c=km.Circuit.from_file(str(path),format='kmx');N=1024;rng=np.random.default_rng(812739)
S=np.zeros((70,N),dtype=np.int64)
for shot in range(N):
 L=shot%36;positions=rng.choice(4095,L,replace=False);values=rng.integers(1,4096,L)
 for r in range(70):
  for pos,v in zip(positions,values):S[r,shot]^=mul(int(v),p.EXP[(r+1)*int(pos)%4095])
sim=km.Simulator(batch_size=N,seed=1582);sim.use_same_registers_as(c)
bits=lambda a,n:((a[None,:]>>np.arange(n)[:,None])&1).astype(bool)
for reg,values in zip(syn,S):sim.write_across_shots(reg,bits(values,12))
sim.do(c)
def read_word(word):
 values=sim.read_across_shots(km.array([km.q(i) for i in word]))
 return (values.astype(np.int64)*(1<<np.arange(len(word)))[:,None]).sum(0)
actuals={name:np.stack([read_word(word) for word in hook[name]]) for name in ['sigma','template']}
actualS=np.stack([read_word(list(range(12*i,12*i+12))) for i in range(70)])
counter=read_word(hook['counter']);failures=[]
for shot in range(N):
 sig,bank,d,synvals,reason=model(S[:,shot].tolist())
 for name,expected in [('sigma',sig),('template',bank),('S',synvals)]:
  actual=actualS[:,shot] if name=='S' else actuals[name][:,shot]
  if not np.array_equal(actual,expected):failures.append(dict(shot=shot,part=name,different=np.flatnonzero(actual!=expected).tolist()))
 if counter[shot]!=(d-1)%256:failures.append(dict(shot=shot,part='counter',actual=int(counter[shot]),expected=(d-1)%256))
result=dict(status='PASS' if not failures else 'FAIL',scope='BM prefix quantum/classical value comparison',shots=N,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,failures=failures)
print(json.dumps(result,indent=2))
