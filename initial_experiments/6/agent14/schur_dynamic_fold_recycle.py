"""Schur synthesis with one or two mutually exclusive node-swap controls."""
from contextlib import contextmanager
from functools import lru_cache
import json
from pathlib import Path
import numpy as np
import kickmix as km
from extension_inverse import Adapter,setup
from schur_fold import axpy,xword,evaluate as fixed_evaluate,reference
from schur_dynamic_division import create as divide_create,phase as divide_phase
from schur_control import caxpy,cimul
from low_workspace_inverse import bind
from tower_inverse import map_ops,apply_map
from fast_dqi import FastBuilder


def dapply(b,E,nn,bb,kappa,controls,inverse=False):
 m=setup()[0];d=len(nn)
 def shear_u():
  for i in range(d):
   b.cx(bb[i],nn[i])
   if i:b.cx(bb[i-1],nn[i])
 def shear_v():
  for i in range(1,d):b.cx(nn[i-1],bb[i])
 def swap():
  for a,c in zip(nn,bb):b.swap(a,c)
 def scale(word,back=False):
  if not back:apply_map(b,word,map_ops([m.gmul(kappa,1<<i) for i in range(12)]))
  for flag,delta in (reversed(controls) if back else controls):cimul(m,E,flag,[q.id for q in word],m.gmul(kappa^delta,m.ginv(kappa)),back)
  if back:apply_map(b,word,map_ops([m.gmul(kappa,1<<i) for i in range(12)]),True)
 if inverse:swap();shear_v()
 else:shear_u()
 for i in (range(d) if inverse else reversed(range(d))):
  if not inverse:scale(bb[i])
  if i:b.cx(bb[i-1],bb[i])
  if i>1:b.cx(bb[i-2],bb[i])
  if inverse:scale(bb[i],True)
 if inverse:shear_u()
 else:shear_v();swap()


def powers(degree,kappa,m):
 out=[(1,0),(0,1)]
 for _ in range(2,degree+1):
  a,c=out[-1];out.append((m.gmul(kappa,c),a^c))
 return out


@contextmanager
def evaluate(b,E,nn,bb,kappa,controls):
 m=setup()[0];d=len(nn);base=powers(d,kappa,m)
 with fixed_evaluate(b,nn,bb,kappa) as (nv,dv):
  def correction(flag,delta):
   pp=powers(d,kappa^delta,m);dif=[(a^c,b^e) for (a,b),(c,e) in zip(base,pp)]
   for bank,target,monic in [(nn,nv,False),(bb,dv,True)]:
    if d>2:
     if d==3:
      # One coefficient supplies all changed output bits: rank <=12.
      rt=[]
      for bit,q in enumerate(bank[2]):
       lo=m.gmul(dif[2][0],1<<bit);hi=m.gmul(dif[2][1],1<<bit);dest=[target[j] for j in range(12) if lo>>j&1]+[target[12+j] for j in range(12) if hi>>j&1]
       if dest:rt.append(('P',(flag,),0,(q.id,),0,tuple(dest)))
      E.mac(rt)
     else:
      tmp=b.alloc_qubits(24)
      def linear():
       for j in range(2,d):axpy(b,bank[j],tmp[:12],dif[j][0],m);axpy(b,bank[j],tmp[12:],dif[j][1],m)
      linear()
      for q,t in zip(tmp,target):b.ccx(km.q(flag),q,km.q(t))
      linear();b.free(tmp)
    if monic:
     for bit in range(12):
      if dif[d][0]>>bit&1:b.cx(km.q(flag),km.q(target[bit]))
      if dif[d][1]>>bit&1:b.cx(km.q(flag),km.q(target[12+bit]))
  for flag,delta in controls:correction(flag,delta)
  yield nv,dv
  for flag,delta in reversed(controls):correction(flag,delta)


def arithmetic(b,E,aa,bs,nn,bb,nt,bt,kappa,controls,inverse=False):
 m=setup()[0];m.inv_word=bind(m);A=[q.id for q in aa];B=[q.id for q in bs];d=len(nn)
 with m.eq_test(E,B,0) as z:
  def original_linear():
   b.cx(aa,nt);axpy(b,bs,nt,kappa,m)
   for flag,delta in controls:caxpy(m,E,flag,B,[q.id for q in nt],delta)
  if not inverse:original_linear()
  E.cx(z,B[0])
  with m.inv_word(E,B) as bi,m.opened(E,'dynamic_fold_product',12) as temp,m.opened(E,'dynamic_fold_ratio',12) as r:
   rr=m.prod_rt(r,A,bi);E.mac(rr)
   def top():
    if d:b.cx(nn[-1],nt);b.cx(bb[-1],bt)
    E.mac(m.prod_rt([q.id for q in nt],A,r));b.cx(km.array([km.q(q) for q in r]),bt);b.x(bt[0]);E.cx(z,bt[0].id)
   if not inverse:top()
   for n,v in zip(nn,bb):
    N=[q.id for q in n];V=[q.id for q in v];E.mac(m.prod_rt(N,A,V));m.mon_slot(E,N,bi,B,temp)
   dapply(b,E,nn,bb,kappa,controls,inverse)
   for n,v in zip(nn,bb):
    N=[q.id for q in n];V=[q.id for q in v];m.mon_slot(E,N,B,bi,temp);E.mac(m.prod_rt(N,A,V))
   if inverse:top()
   E.free(rr)
  E.cx(z,B[0])
  if inverse:original_linear()


def build(d,kappa,deltas,inverse=False):
 m=setup()[0];b=FastBuilder(km.CircuitBuilder());headers=[];rid=0
 def reg(n,name):
  nonlocal rid
  qs=b.alloc_qubits(n)
  for q in qs:headers.append(f'APPEND_TO_REGISTER q{q.id} r{rid}')
  headers.append(f'REGISTER r{rid} "{name}"');rid+=1;return qs
 a=reg(12,'a');bv=reg(12,'b');nn=[reg(12,f'N{i}') for i in range(d)];bb=[reg(12,f'B{i}') for i in range(d)];nt=reg(12,'Nt');bt=reg(12,'Bt');flags=reg(len(deltas),'flags');controls=[(q.id,c) for q,c in zip(flags,deltas)];ab=Adapter(b);sc=ab.qubits(3);E=m.Em(ab,*sc,m.Alloc(ab))
 if inverse:
  with evaluate(b,E,nn+[nt],bb+[bt],kappa,controls) as (nv,dv):divide_create(b,E,nv,dv,[q.id for word in [a,bv] for q in word],kappa,controls)
  arithmetic(b,E,a,bv,nn,bb,nt,bt,kappa,controls,True)
 else:
  arithmetic(b,E,a,bv,nn,bb,nt,bt,kappa,controls);mask=b.alloc_bits(24);b.hmr(km.array(list(a)+list(bv)),mask);b.free(a);b.free(bv)
  with evaluate(b,E,nn+[nt],bb+[bt],kappa,controls) as (nv,dv):divide_phase(b,E,nv,dv,mask,kappa,controls)
  b.bit_store(mask,False);b.free(mask)
 ab.release_recycle('dynamic_fold_emitter',sc);c=b.finish_circuit();return km.Circuit('\n'.join(headers)+'\n'+str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)))


def test(d=3,shots=512):
 m,k,fw,bk,tr=setup();deltas=[tr(6^F.mul(6,6),fw),tr(14^F.mul(14,14),fw)];rng=np.random.default_rng(142050+d);aa=rng.integers(0,4096,shots);bs=rng.integers(1,4096,shots);nn=rng.integers(0,4096,(shots,d));bb=rng.integers(0,4096,(shots,d));choice=rng.integers(0,3,shots);ff=np.array([choice==1,choice==2]);aa[::11]=0;bs[::11]=0;nn[::11]=0;bb[::11]=0;want=[reference(int(a),int(b),n.tolist(),v.tolist(),k^(deltas[0] if f==1 else deltas[1] if f==2 else 0)) for a,b,n,v,f in zip(aa,bs,nn,bb,choice)];wn=np.array([v[0] for v in want]);wb=np.array([v[1] for v in want]);cf=build(d,k,deltas);ci=build(d,k,deltas,True);bits=lambda v,w=12:((v[None,:]>>np.arange(w)[:,None])&1).astype(bool);rows=[]
 src={'a':aa,'b':bs,'Nt':np.zeros(shots,dtype=int),'Bt':np.zeros(shots,dtype=int),**{f'N{i}':nn[:,i] for i in range(d)},**{f'B{i}':bb[:,i] for i in range(d)}};dst={'a':np.zeros(shots,dtype=int),'b':np.zeros(shots,dtype=int),'Nt':wn[:,-1],'Bt':wb[:,-1],**{f'N{i}':wn[:,i] for i in range(d)},**{f'B{i}':wb[:,i] for i in range(d)}}
 for inv,c,s,t in [(False,cf,src,dst),(True,ci,dst,src)]:
  sim=km.Simulator(batch_size=shots,seed=142050+inv);sim.use_same_registers_as(c)
  for name,v in s.items():sim.write_across_shots(c.register_data[name],bits(v))
  sim.write_across_shots(c.register_data['flags'],ff);sim.do(c);errors={}
  for name,v in t.items():
   bad=np.flatnonzero((sim.read_across_shots(c.register_data[name])!=bits(v)).any(0))
   if len(bad):errors[name]=bad[:8].tolist()
  assert np.array_equal(sim.read_across_shots(c.register_data['flags']),ff)
  public={q for w in c.register_data.values() for q in w};sc=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);junk=np.flatnonzero(sim.read_across_shots(sc).any(0));phase=np.flatnonzero(sim.read_phase_flipped_across_shots());r=dict(d=d,inverse=inv,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,errors=errors,junk=junk[:8].tolist(),phase=phase[:8].tolist());print(json.dumps(r),flush=True);rows.append(r);assert not errors and not len(junk) and not len(phase)
  sim.do(cf if inv else ci)
  for name,v in s.items():assert np.array_equal(sim.read_across_shots(c.register_data[name]),bits(v)),name
  assert np.array_equal(sim.read_across_shots(c.register_data['flags']),ff);assert not sim.read_phase_flipped_across_shots().any();assert not sim.read_across_shots(sc).any()
 Path(f'results/schur_dynamic_fold_recycle_d{d}_native.json').write_text(json.dumps(dict(status='PASS',shots=shots,rows=rows,scope='mutually exclusive coherent node controls; arbitrary pivots/tails and zero completions, both directions/compositions, clean scratch and phase'),indent=2)+'\n')
from quadratic_schur_model import F
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--d',type=int,default=3);p.add_argument('--shots',type=int,default=512);a=p.parse_args();test(a.d,a.shots)
