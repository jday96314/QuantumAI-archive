"""Validate H=E V+O z U from the in-place quantum Sylvester prototype."""
import json,random
from pathlib import Path
import numpy as np
import kickmix as km
from bezout_algebra import product
from probe_euclid_linear import g
root=Path(__file__).resolve().parent
c=km.Circuit.from_file(str(root/'candidates/sylvester.kmb'),format='kmb');N=512;rng=random.Random(517229)
Es=[];Os=[];Us=[];Vs=[];Hs=[];rejected=0
while len(Es)<N:
 E=[1]+[rng.randrange(1,4096) for _ in range(17)];O=[rng.randrange(1,4096) for _ in range(18)]
 U=[rng.randrange(4096) for _ in range(17)];V=[rng.randrange(4096) for _ in range(18)]
 try:
  f,h,eps,idx,pre,log=g.forward(E,O)
  if eps:rejected+=1;continue
 except g.ClaimError:rejected+=1;continue
 H=[a^b for a,b in zip(product(E,V),product(O,[0]+U))]
 Es.append(E);Os.append(O);Us.append(U);Vs.append(V);Hs.append(H)
arrays={'E':np.array(Es).T,'O':np.array(Os).T,'X':np.array(Us).T,'V':np.array(Vs).T};H=np.array(Hs).T
bits=lambda a:((a[None,:]>>np.arange(12)[:,None])&1).astype(bool)
sim=km.Simulator(batch_size=N,seed=529);sim.use_same_registers_as(c)
for name,values in arrays.items():
 for i,row in enumerate(values):sim.write_across_shots(c.register_data[f'{name}{i}'],bits(row))
sim.do(c);failures=[]
for name,values in arrays.items():
 expected=H[:18] if name=='V' else H[18:] if name=='X' else values
 for i,row in enumerate(expected):
  out=sim.read_across_shots(c.register_data[f'{name}{i}']);bad=np.flatnonzero(np.any(out!=bits(row),axis=0))
  if len(bad):failures.append(dict(part=name,word=i,bad_shots=bad[:10].tolist(),count=len(bad)))
public={q for reg in c.register_data.values() for q in reg};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
junk=np.flatnonzero(sim.read_across_shots(scratch).any(axis=0));phase=np.flatnonzero(sim.read_phase_flipped_across_shots())
result=dict(status='PASS' if not failures and not len(junk) and not len(phase) else 'FAIL',scope='generic degree-17 polynomials, E[0]=1, no Euclidean degree gaps',shots=N,rejected_degenerate_pairs=rejected,Q=c.num_qubits,public_qubits=len(public),extra_qubits=c.num_qubits-len(public),T=c.operation_statistics().num_toffoli_like,failures=failures,scratch_errors=junk[:20].tolist(),phase_errors=phase[:20].tolist())
print(json.dumps(result,indent=2))
