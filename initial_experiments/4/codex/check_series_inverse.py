"""Check normalized-series reciprocal, including all public and scratch bits."""
import argparse,json,random
from pathlib import Path
import numpy as np
import kickmix as km
from import_published import load

def check(n,shots,optimized=False):
 root=Path(__file__).resolve().parent;m=load();rng=random.Random(1745+n)
 tag=f'series_inverse_{n}'+('_optimized' if optimized else '')
 c=km.Circuit.from_file(str(root/f'candidates/{tag}.kmb'),format='kmb')
 values=[[1]+[rng.randrange(4096) for _ in range(n-1)] for _ in range(shots)]
 # Exercise sparse, zero-tail, and special all-one series as well.
 values[:min(3,shots)]=([[1]+[0]*(n-1),[1]*n,[1]+[4095]*(n-1)])[:min(3,shots)]
 expected=[]
 for F in values:
  G=[1]
  for k in range(1,n):
   v=0
   for i in range(1,k+1):v^=m.gmul(F[i],G[k-i])
   G.append(v)
  expected.append(G)
 values=np.asarray(values,dtype=np.uint16).T;expected=np.asarray(expected,dtype=np.uint16).T
 bits=lambda a:((a[None,:]>>np.arange(12)[:,None])&1).astype(bool)
 sim=km.Simulator(batch_size=shots,seed=721);sim.use_same_registers_as(c)
 for i,row in enumerate(values):sim.write_across_shots(c.register_data[f'F{i}'],bits(row))
 sim.do(c);failures=[]
 for i,row in enumerate(expected):
  bad=np.flatnonzero(np.any(sim.read_across_shots(c.register_data[f'F{i}'])!=bits(row),axis=0))
  if len(bad):failures.append(dict(word=i,bad_shots=bad[:10].tolist(),count=len(bad)))
 public={q for reg in c.register_data.values() for q in reg}
 scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
 junk=np.flatnonzero(sim.read_across_shots(scratch).any(axis=0));phase=np.flatnonzero(sim.read_phase_flipped_across_shots())
 out=dict(status='PASS' if not failures and not len(junk) and not len(phase) else 'FAIL',scope='normalized reciprocal modulo z^n, constant coefficient one',n=n,shots=shots,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,failures=failures,scratch_errors=junk[:20].tolist(),phase_errors=phase[:20].tolist())
 (root/f'results/{tag}_check.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
 return out

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=35);p.add_argument('--shots',type=int,default=512);p.add_argument('--optimized',action='store_true');a=p.parse_args();check(a.n,a.shots,a.optimized)
