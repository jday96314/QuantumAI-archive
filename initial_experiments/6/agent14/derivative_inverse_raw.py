"""Overwrite A by A/B' modulo a squarefree monic odd-degree B.

Native implementation on the generic Euclidean E/O chart. It restores B and
has a literal inverse even on reversibly armed unsupported degree gaps.
"""
import argparse,json
from pathlib import Path
import numpy as np
import kickmix as km
from extension_inverse import Adapter,setup
from cubic_low_workspace_inverse import bind
from tower_inverse import apply_map,map_ops
from triangular_dp import DPArithmetic
from bezout_payload_monic import apply as bezout_apply
from fast_dqi import FastBuilder
from lean_emitter import create as lean_create


def apply(m,E,A,B,one,inverse=False):
 n=len(A);assert n%2 and len(B)==n;t=n//2;b=E.b.raw;even=A[::2];odd=A[1::2];f=B[::2];g=B[1::2];math=DPArithmetic(m,E);rootops=map_ops([m.gpow(1<<i,2048) for i in range(12)])
 def root(back=False):
  for word in A+B:apply_map(b,km.array([km.q(q) for q in word]),rootops,back)
 def cross(first=True):
  for i in range(t):E.mac(m.prod_rt(even[i] if first else even[-1],g[i],even[-1] if first else even[i]))
 def bez(reverse=False):
  with m.opened(E,'derivative_gap_index',max(1,t.bit_length())) as index,m.opened(E,'derivative_gap_flag',1) as flag:bezout_apply(m,E,f,g,[even[:-1],odd],reverse,history=(index,flag[0]))
 def swap():
  for a,c in zip(even[:-1],odd):
   for q,r in zip(a,c):E.b.swap(q,r)
 def reverse():
  for i in range(n//2):
   for q,r in zip(A[i],A[-1-i]):E.b.swap(q,r)
 if inverse:
  reverse();math.triangular([one]+list(reversed(B[1:])),1,A,True);root();cross(False);swap();bez(True);cross();root(True)
 else:
  root();cross();bez();swap();cross(False);root(True);math.triangular([one]+list(reversed(B[1:])),1,A);reverse()


def build(n,inverse=False):
 m=setup()[0];m.inv_word=bind(m);b=FastBuilder(km.CircuitBuilder());aa=[b.create_quantum_register(12,f'A{i}') for i in range(n)];bb=[b.create_quantum_register(12,f'B{i}') for i in range(n)];ab=Adapter(b);sc,E=lean_create(m,ab);apply(m,E,[[q.id for q in w] for w in aa],[[q.id for q in w] for w in bb],[],inverse);ab.release_recycle('derivative_emitter',sc);c=b.finish_circuit();return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)))


def test(n=7,shots=256):
 from derivative_inverse_model import invmod,divide
 from quadratic_schur_model import F
 m,_,fw,bk,tr=setup();rng=np.random.default_rng(142086+n);cases=[];rejected=0
 while len(cases)<shots:
  roots=rng.choice(np.arange(1,4096),n,replace=False).tolist();B=[1]
  for a in roots:B=F.product(B,[a,1])
  f=[F.pow(x,2048) for x in B[::2]];g=[F.pow(x,2048) for x in B[1::2]]
  if len(cases)%2 and n>=7:
   f=rng.integers(0,4096,n//2+1).tolist();g=rng.integers(0,4096,n//2).tolist()+[1];f[-2]=F.mul(f[-1],g[-2]);B=[0]*(n+1)
   for i,v in enumerate(f):B[2*i]=F.mul(v,v)
   for i,v in enumerate(g):B[2*i+1]=F.mul(v,v)
  valid=True;used=False;_,r=divide(f,g);f,g=g,r;expected=n//2
  while expected:
   if len(g)!=expected:
    if not used and len(g)==expected-1 and f[0]:g=[0]+g;used=True
    else:valid=False;break
   if expected==1:break
   _,r=divide(f,g);f,g=g,r;expected-=1
  if not valid or not g or not g[0]:rejected+=1;continue
  A=rng.integers(0,4096,n).tolist();I=invmod(F.derivative(B),B);want=divide(F.product(A,I),B)[1];want+=[0]*(n-len(want));cases.append((A,B[:-1],want))
 cf=build(n);ci=build(n,True);bits=lambda v:((np.asarray(v)[None,:]>>np.arange(12)[:,None])&1).astype(bool);rows=[]
 for inverse,c in [(False,cf),(True,ci)]:
  sim=km.Simulator(batch_size=shots,seed=142086+inverse);sim.use_same_registers_as(c)
  for i in range(n):sim.write_across_shots(c.register_data[f'A{i}'],bits([tr(v[2 if inverse else 0][i],fw) for v in cases]));sim.write_across_shots(c.register_data[f'B{i}'],bits([tr(v[1][i],fw) for v in cases]))
  sim.do(c);errors={}
  for name,key in [('A',0 if inverse else 2),('B',1)]:
   for i in range(n):
    if not np.array_equal(sim.read_across_shots(c.register_data[f'{name}{i}']),bits([tr(v[key][i],fw) for v in cases])):errors[f'{name}{i}']=True
  public={q for reg in c.register_data.values() for q in reg};junk=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);r=dict(inverse=inverse,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,errors=errors,scratch=bool(sim.read_across_shots(junk).any()),phase=bool(sim.read_phase_flipped_across_shots().any()));rows.append(r);print(json.dumps(r),flush=True);assert not errors and not r['scratch'] and not r['phase'];sim.do(cf if inverse else ci);assert not sim.read_across_shots(junk).any();assert not sim.read_phase_flipped_across_shots().any()
 result=dict(status='PASS',n=n,shots=shots*2,rows=rows,rejected_generic_charts=rejected,scope='distinct nonzero split roots and forced one-gap squarefree nonsplit polynomials, arbitrary numerator, exact reciprocal derivative overwrites, fresh inverse and paired phase/scratch cleanup on generic or one repaired leading-zero charts');Path(f'results/derivative_inverse_raw_n{n}_native.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=7);p.add_argument('--shots',type=int,default=256);a=p.parse_args();test(a.n,a.shots)
