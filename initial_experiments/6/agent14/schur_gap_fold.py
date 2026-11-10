"""Degree-two scalar-chart synthesis with four consumed field words.

Input heads (a,0,u,h) synthesize the exact checked gap matrix. Numerator
and denominator coefficients grow by two; both leading coefficients stay
implicit. This is a component, not a coherent mixed-degree full synthesis.
"""
import argparse,json
from pathlib import Path
import numpy as np
import kickmix as km
from extension_inverse import Adapter,setup
from extension_division import create as divide_create,phase as divide_phase
from schur_fold import evaluate,axpy,xword
from schur_extract import imul_text
from tower_inverse import map_ops,apply_map
from fast_dqi import FastBuilder


def quadratic(b,bank,kappa,inverse=False):
 m=setup()[0];ops=map_ops([m.gmul(kappa,1<<i) for i in range(12)])
 for i in (range(len(bank)) if inverse else reversed(range(len(bank)))):
  if not inverse:apply_map(b,bank[i],ops)
  if i:b.cx(bank[i-1],bank[i])
  if i>1:b.cx(bank[i-2],bank[i])
  if inverse:apply_map(b,bank[i],ops,True)


def arithmetic(b,E,a,u,h,nn,bb,k1,k2,inverse=False):
 m=setup()[0];A=[q.id for q in a];U=[q.id for q in u];H=[q.id for q in h]
 def swap():
  for n,v in zip(nn,bb):b.swap(n,v)
 def bya():
  for n,v in zip(nn,bb):E.mac(m.prod_rt([q.id for q in v],A,[q.id for q in n]))
 def byv():
  for i,(n,v) in enumerate(zip(nn,bb)):
   E.mac(m.prod_rt([q.id for q in n],U,[q.id for q in v]))
   if i:E.mac(m.prod_rt([q.id for q in n],H,[q.id for q in bb[i-1]]))
 if inverse:swap();bya();quadratic(b,bb,k1,True);byv();quadratic(b,nn,k2,True);bya();swap()
 else:swap();bya();quadratic(b,nn,k2);byv();quadratic(b,bb,k1);bya();swap()


def head_tail_quotient(b,E,nn,bb,a,k2,delta,out=None,mask=None):
 m=setup()[0]
 with evaluate(b,nn,bb,k2) as (nv,bv):
  with m.opened(E,'gap_head_divisor',24) as den:
   rt=[('L',(q,),0,(t,)) for q,t in zip(nv,den)]
   rt+=m.prod_rt(den[:12],[q.id for q in a],bv[:12]);rt+=m.prod_rt(den[12:],[q.id for q in a],bv[12:]);E.mac(rt)
   nb=km.array([km.q(q) for q in bv]);ops=map_ops([m.gmul(delta,1<<i) for i in range(12)]);apply_map(b,nb[:12],ops);apply_map(b,nb[12:],ops)
   if out is None:divide_phase(b,E,bv,den,mask,k2)
   else:divide_create(b,E,bv,den,out,k2)
   apply_map(b,nb[12:],ops,True);apply_map(b,nb[:12],ops,True);E.free(rt)


def build(d,k1,k2,inverse=False):
 m=setup()[0];b=FastBuilder(km.CircuitBuilder());headers=[];rid=0
 def reg(n,name):
  nonlocal rid
  qs=b.alloc_qubits(n)
  for q in qs:headers.append(f'APPEND_TO_REGISTER q{q.id} r{rid}')
  headers.append(f'REGISTER r{rid} "{name}"');rid+=1;return qs
 a=reg(12,'a');zero=reg(12,'zero');u=reg(12,'u');h=reg(12,'h');nn=[reg(12,f'N{i}') for i in range(d)];bb=[reg(12,f'B{i}') for i in range(d)];nt=[reg(12,f'Nt{i}') for i in range(2)];bt=[reg(12,f'Bt{i}') for i in range(2)];N=nn+nt;B=bb+bt;ab=Adapter(b);sc=ab.qubits(3);E=m.Em(ab,*sc,m.Alloc(ab))
 if inverse:
  with evaluate(b,N,B,k1) as (nv,bv):divide_create(b,E,nv,bv,[q.id for w in [a,zero] for q in w],k1)
  head_tail_quotient(b,E,N,B,a,k2,k1^k2,out=[q.id for w in [u,h] for q in w]);arithmetic(b,E,a,u,h,N,B,k1,k2,True);b.x(bt[0][0])
 else:
  b.x(bt[0][0]);arithmetic(b,E,a,u,h,N,B,k1,k2);mask=b.alloc_bits(24);b.hmr(km.array(list(u)+list(h)),mask);b.free(u);b.free(h);head_tail_quotient(b,E,N,B,a,k2,k1^k2,mask=mask);b.bit_store(mask,False);b.free(mask)
  mask=b.alloc_bits(24);b.bit_store(mask,False);b.hmr(a,mask[:12]);b.free(a)
  with evaluate(b,N,B,k1) as (nv,bv):divide_phase(b,E,nv,bv,mask,k1)
  b.bit_store(mask,False);b.free(mask)
 ab.release_recycle('gap_fold_emitter',sc);c=b.finish_circuit();return km.Circuit('\n'.join(headers)+'\n'+str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)))


def reference(a,u,h,nn,bb,k1,k2,m):
 def add(a,b):
  out=a.copy()+[0]*max(0,len(b)-len(a))
  for i,v in enumerate(b):out[i]^=v
  return out
 def product(a,b):
  out=[0]*(len(a)+len(b)-1)
  for i,v in enumerate(a):
   for j,w in enumerate(b):out[i+j]^=m.gmul(v,w)
  return out
 B=bb+[1];N=nn;v=[u,h];av=[m.gmul(a,x) for x in v];M00=add([k1,1,1],av);M01=[m.gmul(a,k1^k2)^m.gmul(a,m.gmul(a,u)),m.gmul(a,m.gmul(a,h))];M11=add([k2,1,1],av);newn=add(product(M00,N),product(M01,B));newb=add(product(v,N),product(M11,B));return (newn+[0]*(len(nn)+2))[:len(nn)+2],(newb+[0]*(len(nn)+2))[:len(nn)+2]


def test(d=2,shots=512):
 m,k,fw,bk,tr=setup();from quadratic_schur_model import F
 k1=k;k2=k^tr(6^F.mul(6,6),fw);rng=np.random.default_rng(142062+d);aa,uu,hh=rng.integers(0,4096,(3,shots));n=rng.integers(0,4096,(shots,d));v=rng.integers(0,4096,(shots,d));aa[::11]=uu[::11]=hh[::11]=0;n[::11]=v[::11]=0;wn=[];wb=[]
 for a,u,h,x,y in zip(aa,uu,hh,n,v):
  xx,yy=reference(int(a),int(u),int(h),x.tolist(),y.tolist(),k1,k2,m);wn.append(xx);wb.append(yy)
 wn=np.array(wn);wb=np.array(wb);cf=build(d,k1,k2);ci=build(d,k1,k2,True);bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool);src={'a':aa,'u':uu,'h':hh,'zero':np.zeros(shots,dtype=int),**{f'N{i}':n[:,i] for i in range(d)},**{f'B{i}':v[:,i] for i in range(d)},**{f'Nt{i}':np.zeros(shots,dtype=int) for i in range(2)},**{f'Bt{i}':np.zeros(shots,dtype=int) for i in range(2)}};dst={name:np.zeros(shots,dtype=int) for name in ['a','u','h','zero']};dst.update({f'N{i}':wn[:,i] for i in range(d)});dst.update({f'B{i}':wb[:,i] for i in range(d)});dst.update({f'Nt{i}':wn[:,d+i] for i in range(2)});dst.update({f'Bt{i}':wb[:,d+i] for i in range(2)});rows=[]
 for inv,c,ss,tt in [(False,cf,src,dst),(True,ci,dst,src)]:
  sim=km.Simulator(batch_size=shots,seed=142062+inv);sim.use_same_registers_as(c)
  for name,values in ss.items():sim.write_across_shots(c.register_data[name],bits(values))
  sim.do(c);errors={}
  for name,values in tt.items():
   bad=np.flatnonzero((sim.read_across_shots(c.register_data[name])!=bits(values)).any(0))
   if len(bad):errors[name]=bad[:8].tolist()
  public={q for w in c.register_data.values() for q in w};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);junk=np.flatnonzero(sim.read_across_shots(scratch).any(0));phase=np.flatnonzero(sim.read_phase_flipped_across_shots());r=dict(d=d,inverse=inv,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,errors=errors,junk=junk[:8].tolist(),phase=phase[:8].tolist());rows.append(r);print(json.dumps(r),flush=True);assert not errors and not len(junk) and not len(phase)
  sim.do(cf if inv else ci)
  for name,values in ss.items():assert np.array_equal(sim.read_across_shots(c.register_data[name]),bits(values)),name
  assert not sim.read_across_shots(scratch).any();assert not sim.read_phase_flipped_across_shots().any()
 result=dict(status='PASS',shots=2*shots,rows=rows,scope='scalar-gap pair, arbitrary three head fields and monic proper tail, zero completion cases, exact coefficient model, both directions/compositions, scratch and phase');Path(f'results/schur_gap_fold_d{d}_native.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--d',type=int,default=2);p.add_argument('--shots',type=int,default=512);a=p.parse_args();test(a.d,a.shots)
