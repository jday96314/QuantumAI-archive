"""Contract the proper outer rational tail without storing its two zero tops."""
import argparse,json
from pathlib import Path
import numpy as np
import kickmix as km
from extension_inverse import Adapter,setup
from cubic_low_workspace_inverse_compact import bind
from fast_dqi import FastBuilder
from lean_emitter import create
from triangular_dp import DPArithmetic
from schur_division_lowtemp import create as quotient
from raw_cf_step import evaluate
from tower_inverse import map_ops,apply_map
from schur_fold import axpy


def build(n,kappa,inverse=False):
 m=setup()[0];m.inv_word=bind(m);b=FastBuilder(km.CircuitBuilder());N=[b.create_quantum_register(12,f'N{i}') for i in range(n+1)];B=[b.create_quantum_register(12,f'B{i}') for i in range(n)];ghost=b.create_classical_register(24,'outer_ghost');ab=Adapter(b);sc,E=create(m,ab);a=E.A.open('outer_a',12);h=E.A.open('outer_b',12);math=DPArithmetic(m,E);pts,I=m.vinv(n+2)
 def qscale(back=False):
  ops=map_ops([m.gmul(kappa,1<<j) for j in range(12)])
  for i in (range(n) if back else reversed(range(n))):
   if not back:apply_map(b,B[i],ops)
   for j in [i-1,i-2]:
    if j>=0:b.cx(B[j],B[i])
   if back:apply_map(b,B[i],ops,True)
 def update():
  for col,x in enumerate(pts):
   left=math.spec([(1,a,1),(x,h,1)]);right=math.spec([(m.gpow(x,i),[q.id for q in w],1) for i,w in enumerate(N)])
   def targets(tb):return [B[i][k].id for i in range(n) for k in m.bitlist(m.gmul(I[i][col],m.tmask(tb)),12)]
   math.product(left,right,targets)
 if inverse:
  with evaluate(b,B,N,kappa) as (nv,dv):quotient(b,E,nv,dv,a+h,kappa,[])
  for q,bit in zip(a+h,ghost):ab.z(q,bit)
  b.bit_store(ghost,False);update();qscale(True)
 else:
  lead=[q.id for q in N[-1]];nextlead=[q.id for q in N[-2]]
  with m.inv_word(E,lead) as ni:
   # The two head fields are still zero. Use a as the product workspace.
   m.mon_slot(E,nextlead,ni,lead,a)
   E.mac(m.prod_rt(h,[q.id for q in B[-1]],ni))
   for q,r in zip(B[-2],B[-1]):b.cx(q,r)
   E.mac(m.prod_rt(a,[q.id for q in B[-1]],ni))
   for q,r in zip(B[-2],B[-1]):b.cx(q,r)
   E.mac(m.prod_rt(a,h,nextlead))
   # Recycle h for restoring nextlead; repay its measurement on recreation.
   replay=b.alloc_bits(12);b.hmr(km.array([km.q(q) for q in h]),replay)
   m.mon_slot(E,nextlead,lead,ni,h)
   E.mac(m.prod_rt(h,[q.id for q in B[-1]],ni))
   for q,bit in zip(h,replay):E.b.z(q,bit)
   b.bit_store(replay,False);b.free(replay)
  qscale();update();b.hmr(km.array([km.q(q) for q in a+h]),ghost)
 if inverse:
  lead=[q.id for q in N[-1]];nextlead=[q.id for q in N[-2]]
  with m.inv_word(E,lead) as ni:
   E.free(m.prod_rt(h,[q.id for q in B[-1]],ni))
   m.mon_slot(E,nextlead,ni,lead,h)
   E.mac(m.prod_rt(h,[q.id for q in B[-1]],ni))
   for q,r in zip(B[-2],B[-1]):b.cx(q,r)
   E.free(m.prod_rt(a,[q.id for q in B[-1]],ni)+m.prod_rt(a,h,nextlead))
   for q,r in zip(B[-2],B[-1]):b.cx(q,r)
   m.mon_slot(E,nextlead,lead,ni,a)
   E.free(m.prod_rt(h,[q.id for q in B[-1]],ni))
 E.A.close('outer_b',h);E.A.close('outer_a',a);ab.release_recycle('outer_emitter',sc);c=b.finish_circuit();return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(24,c.num_bits)))


def test(N=6,shots=256):
 from raw_cf_model import endpoint,extract
 from quadratic_schur_model import F,node_maps
 n=N//2;m,k,fw,bk,tr=setup();_,kap,_,_=node_maps(n);rng=np.random.default_rng(142124+N);cases=[]
 for case in range(shots):
  d=case%(n+1);pos=rng.choice(np.arange(1,4096),d,replace=False).tolist();val=rng.integers(1,4096,d).tolist();heads,_=extract(F.syndrome(pos,val,N));nn,bb,exc,_=endpoint(heads)
  if exc:continue
  for i in reversed(range(1,n)):
   a,h=heads[i];nn,bb=F.add(F.product([a,h],nn),F.product([kap[i],1,1],bb)),nn
   while nn and not nn[-1]:nn.pop()
  if len(nn)-1<n:nn=[0]+nn;bb=[0]+bb
  a,h=heads[0];dd=F.add(F.product([a,h],nn),F.product([kap[0],1,1],bb));dd+= [0]*(n+2-len(dd));assert not any(dd[n:]);cases.append((nn+[0]*(n+1-len(nn)),bb+[0]*(n-len(bb)),dd[:n],a,h))
 shots=len(cases);cf=build(n,tr(kap[0],fw));ci=build(n,tr(kap[0],fw),True);bits=lambda v:((np.array(v)[None,:]>>np.arange(12)[:,None])&1).astype(bool);sim=km.Simulator(batch_size=shots,seed=142124);sim.use_same_registers_as(ci)
 for name,ix,count in [('N',0,n+1),('B',1,n)]:
  for i in range(count):sim.write_across_shots(cf.register_data[f'{name}{i}'],bits([tr(v[ix][i],fw) for v in cases]))
 sim.do(cf)
 for i in range(n):assert np.array_equal(sim.read_across_shots(cf.register_data[f'B{i}']),bits([tr(v[2][i],fw) for v in cases])),i
 masks=sim.read_across_shots(cf.register_data['outer_ghost']);head=np.concatenate([bits([tr(v[3],fw) for v in cases]),bits([tr(v[4],fw) for v in cases])]);expected=((head&masks).sum(0)%2).astype(bool);assert np.array_equal(sim.read_phase_flipped_across_shots(),expected)
 sim.do(ci)
 for name,ix,count in [('N',0,n+1),('B',1,n)]:
  for i in range(count):assert np.array_equal(sim.read_across_shots(cf.register_data[f'{name}{i}']),bits([tr(v[ix][i],fw) for v in cases])),(name,i)
 public={q for name,reg in cf.register_data.items() if name!='outer_ghost' for q in reg};junk=km.array([km.q(q) for q in range(max(cf.num_qubits,ci.num_qubits)) if km.q(q) not in public]);assert not sim.read_across_shots(junk).any();assert not sim.read_phase_flipped_across_shots().any();assert not sim.read_across_shots(cf.register_data['outer_ghost']).any();r=dict(status='PASS',N=N,shots=shots,forward=dict(Q=cf.num_qubits,T=cf.operation_statistics().num_toffoli_like),inverse=dict(Q=ci.num_qubits,T=ci.operation_statistics().num_toffoli_like));Path(f'results/raw_cf_outer_lowtemp_n{N}_native.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=6);a=p.parse_args();test(a.n)
