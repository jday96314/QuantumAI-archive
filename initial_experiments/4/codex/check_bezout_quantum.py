"""Validate the in-place Bezout prototype on generic coprime polynomial pairs."""
import json,random
from pathlib import Path
import numpy as np
import kickmix as km
from bezout_algebra import bezout,matvec
from probe_euclid_linear import g
root=Path(__file__).resolve().parent;hook=json.loads((root/'parameters/bezout_hook.json').read_text())
c=km.Circuit.from_file(str(root/'candidates/bezout.kmb'),format='kmb');N=512;rng=random.Random(519120)
Es=[];Os=[];Xs=[];Ys=[];rejected=0
while len(Es)<N:
 E=[rng.randrange(1,4096) for _ in range(18)];O=[rng.randrange(1,4096) for _ in range(18)];x=[rng.randrange(4096) for _ in range(17)]
 try:
  f,h,eps,idx,pre,log=g.forward(E,O)
  if eps:rejected+=1;continue
 except g.ClaimError:rejected+=1;continue
 Es.append(E);Os.append(O);Xs.append(x);Ys.append(matvec(bezout(E,O),x))
arrays={'E':np.array(Es).T,'O':np.array(Os).T,'X':np.array(Xs).T};expected=np.array(Ys).T
bits=lambda a:((a[None,:]>>np.arange(12)[:,None])&1).astype(bool)
sim=km.Simulator(batch_size=N,seed=721);sim.use_same_registers_as(c)
for name,values in arrays.items():
 for i,row in enumerate(values):sim.write_across_shots(c.register_data[f'{name}{i}'],bits(row))
sim.do(c);failures=[]
for name,values in arrays.items():
 for i,row in enumerate(expected if name=='X' else values):
  out=sim.read_across_shots(c.register_data[f'{name}{i}']);bad=np.flatnonzero(np.any(out!=bits(row),axis=0))
  if len(bad):failures.append(dict(part=name,word=i,bad_shots=bad[:10].tolist(),count=len(bad)))
public={q for reg in c.register_data.values() for q in reg};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
junk=np.flatnonzero(sim.read_across_shots(scratch).any(axis=0));phase=np.flatnonzero(sim.read_phase_flipped_across_shots())
result=dict(status='PASS' if not failures and not len(junk) and not len(phase) else 'FAIL',scope='generic coprime polynomials, no Euclidean degree gaps',shots=N,rejected_degenerate_pairs=rejected,Q=c.num_qubits,public_qubits=len(public),extra_qubits=c.num_qubits-len(public),T=c.operation_statistics().num_toffoli_like,failures=failures,scratch_errors=junk[:20].tolist(),phase_errors=phase[:20].tolist())
print(json.dumps(result,indent=2))
