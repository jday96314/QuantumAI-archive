"""Affine-chain step: (u,t), equal-degree N/B -> T D (N/B).

Both polynomials have degree d and B is monic. One extra numerator
leading field stays explicit. Completion/boundary flags describe a known
zero tail; the omitted t component is repaid against 1 or a supplied
boundary factor. This is a native component, not a complete decoder.
"""
from contextlib import contextmanager
import argparse,json
from pathlib import Path
import numpy as np
import kickmix as km
from extension_inverse import Adapter,setup
from extension_division import create as divide_create,phase as divide_phase
from schur_dynamic_fold import dapply,powers
from schur_fold import axpy,xword
from schur_control import caxpy
from low_workspace_inverse import bind
from fast_dqi import FastBuilder


@contextmanager
def evaluate(b,N,B,kappa):
 m=setup()[0];deg=len(B);assert len(N)==deg+1 and deg>=1;extra=None;nv=N[:2]
 if deg==1:extra=b.alloc_qubits(12);bv=[B[0],extra]
 else:bv=B[:2]
 pp=powers(deg,kappa,m)
 def linear():
  for j in range(2,deg+1):axpy(b,N[j],nv[0],pp[j][0],m);axpy(b,N[j],nv[1],pp[j][1],m)
  for j in range(2,deg):axpy(b,B[j],bv[0],pp[j][0],m);axpy(b,B[j],bv[1],pp[j][1],m)
  xword(b,bv[0],pp[deg][0]);xword(b,bv[1],pp[deg][1])
 linear();yield [q.id for w in nv for q in w],[q.id for w in bv for q in w];linear()
 if extra is not None:b.free(extra)


def affine(b,E,u,t,N,B,inverse=False):
 m=setup()[0];m.inv_word=bind(m);U=[q.id for q in u];T=[q.id for q in t]
 with m.inv_word(E,T) as ti,m.opened(E,'affine_chain_product',12) as temp:
  def scale(back=False):
   for w in N:m.mon_slot(E,[q.id for q in w],ti if back else T,T if back else ti,temp)
  def shear():
   for n,v in zip(N[:-1],B):E.mac(m.prod_rt([q.id for q in n],U,[q.id for q in v]))
   b.cx(u,N[-1])
  if inverse:shear();scale(True)
  else:scale();shear()


def build(d,kappa,inverse=False):
 m=setup()[0];b=FastBuilder(km.CircuitBuilder());headers=[];rid=0
 def reg(n,name):
  nonlocal rid
  q=b.alloc_qubits(n)
  for x in q:headers.append(f'APPEND_TO_REGISTER q{x.id} r{rid}')
  headers.append(f'REGISTER r{rid} "{name}"');rid+=1;return q
 u=reg(12,'u');t=reg(12,'t');N=[reg(12,f'N{i}') for i in range(d+1)];B=[reg(12,f'B{i}') for i in range(d)];nl=reg(12,'Nlead');bt=reg(12,'Btop');completion=reg(1,'completed')[0];boundary=reg(1,'boundary')[0];factor=reg(12,'factor');ab=Adapter(b);sc=ab.qubits(3);E=m.Em(ab,*sc,m.Alloc(ab));NN=N+[nl];BB=B+[bt]
 def correction():
  for bit in range(12):
   if kappa>>bit&1:b.cx(completion,N[-1][bit]);b.cx(boundary,N[-1][bit])
  b.cx(completion,bt[0]);b.cx(boundary,bt[0])
 if inverse:
  with evaluate(b,NN,BB,kappa) as (nv,bv):divide_create(b,E,nv,bv,[q.id for w in [u,t] for q in w],kappa)
  b.cx(completion,t[0]);caxpy(m,E,boundary.id,[q.id for q in factor],[q.id for q in t],1);affine(b,E,u,t,NN,BB,True);correction();dapply(b,E,N,BB,kappa,[],True);b.x(bt[0]);b.cx(N[-1],nl)
 else:
  b.cx(N[-1],nl);b.x(bt[0]);dapply(b,E,N,BB,kappa,[]);correction();affine(b,E,u,t,NN,BB);mask=b.alloc_bits(24);b.hmr(km.array(list(u)+list(t)),mask);b.free(u);b.free(t)
  with evaluate(b,NN,BB,kappa) as (nv,bv):divide_phase(b,E,nv,bv,mask,kappa)
  b.push_condition(mask[12]);b.z(completion);b.pop_condition()
  for bit in range(12):b.push_condition(mask[12+bit]);b.cz(boundary,factor[bit]);b.pop_condition()
  b.bit_store(mask,False);b.free(mask)
 ab.release_recycle('affine_chain_emitter',sc);c=b.finish_circuit();return km.Circuit('\n'.join(headers)+'\n'+str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)))


def reference(u,t,N,B,kappa,comp):
 m=setup()[0];d=len(B);nn=[0]*(d+2);bb=[0]*(d+2);bp=B+[1]
 for i in range(d+1):
  nn[i]^=m.gmul(kappa,bp[i]);nn[i+1]^=N[i];bb[i]^=N[i]^bp[i];bb[i+1]^=bp[i]
 if comp:nn[d]^=kappa;bb[d]^=1
 nn=[m.gmul(t,a)^m.gmul(u,b) for a,b in zip(nn,bb)];assert bb[-1]==1;return nn,bb[:-1]


def test(d=3,shots=512):
 m,k,fw,bk,tr=setup();rng=np.random.default_rng(142063+d);uu=rng.integers(0,4096,shots);tt=rng.integers(1,4096,shots);nn=rng.integers(0,4096,(shots,d+1));bb=rng.integers(0,4096,(shots,d));choice=np.arange(shots)%3;ff=rng.integers(1,4096,shots);cc=choice==1;bound=choice==2;nn[cc|bound]=0;bb[cc|bound]=0;uu[cc]=0;tt[cc]=1;tt[bound]=ff[bound];want=[reference(int(u),int(t),n.tolist(),b.tolist(),k,bool(c or q)) for u,t,n,b,c,q in zip(uu,tt,nn,bb,cc,bound)];wn=np.array([v[0] for v in want]);wb=np.array([v[1] for v in want]);cf=build(d,k);ci=build(d,k,True);bits=lambda v,w=12:((v[None,:]>>np.arange(w)[:,None])&1).astype(bool);src={'u':uu,'t':tt,**{f'N{i}':nn[:,i] for i in range(d+1)},**{f'B{i}':bb[:,i] for i in range(d)},'Nlead':np.zeros(shots,dtype=int),'Btop':np.zeros(shots,dtype=int)};dst={'u':np.zeros(shots,dtype=int),'t':np.zeros(shots,dtype=int),**{f'N{i}':wn[:,i] for i in range(d+1)},**{f'B{i}':wb[:,i] for i in range(d)},'Nlead':wn[:,-1],'Btop':wb[:,-1]};rows=[]
 for inv,c,source,target in [(False,cf,src,dst),(True,ci,dst,src)]:
  sim=km.Simulator(batch_size=shots,seed=142063+inv);sim.use_same_registers_as(c)
  for key,v in source.items():sim.write_across_shots(c.register_data[key],bits(v))
  for key,v in [('completed',cc),('boundary',bound),('factor',ff)]:sim.write_across_shots(c.register_data[key],bits(v,len(c.register_data[key])))
  sim.do(c);errors={}
  for key,v in target.items():
   bad=np.flatnonzero((sim.read_across_shots(c.register_data[key])!=bits(v)).any(0))
   if len(bad):errors[key]=bad[:8].tolist()
  for key,v in [('completed',cc),('boundary',bound),('factor',ff)]:assert np.array_equal(sim.read_across_shots(c.register_data[key]),bits(v,len(c.register_data[key])))
  public={q for w in c.register_data.values() for q in w};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);junk=np.flatnonzero(sim.read_across_shots(scratch).any(0));phase=np.flatnonzero(sim.read_phase_flipped_across_shots());r=dict(d=d,inverse=inv,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,errors=errors,junk=junk[:8].tolist(),phase=phase[:8].tolist());rows.append(r);print(json.dumps(r),flush=True);assert not errors and not len(junk) and not len(phase)
  sim.do(cf if inv else ci)
  for key,v in source.items():assert np.array_equal(sim.read_across_shots(c.register_data[key]),bits(v)),key
  assert not sim.read_across_shots(scratch).any();assert not sim.read_phase_flipped_across_shots().any()
 result=dict(status='PASS',shots=shots*2,rows=rows,scope='affine-chain T D step, equal-degree generic polynomials, completed/boundary zero tails, exact coefficient model, both directions/compositions, restored flags/factor and clean phase/scratch');Path(f'results/schur_affine_step_d{d}_native.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--d',type=int,default=3);p.add_argument('--shots',type=int,default=512);a=p.parse_args();test(a.d,a.shots)
