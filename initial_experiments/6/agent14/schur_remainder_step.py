"""Affine T D step in quadratic remainder/quotient coordinates.

No new polynomial slots are allocated. Old head slots hold the new numerator
remainder, and the old numerator slots hold the new denominator. Head HMR
phases remain explicit classical outputs until the matching inverse step.
"""
from contextlib import contextmanager
import argparse,json
from pathlib import Path
import numpy as np
import kickmix as km
from fast_dqi import FastBuilder
from extension_inverse import Adapter,setup,init_extension_inverse,inverse_phase
from schur_fold import axpy,xword
from schur_control import caxpy
from low_workspace_inverse import bind


def quotient(poly,kappa):
 m=setup()[0];p=list(poly);out=[0]*max(0,len(p)-2)
 for j in reversed(range(2,len(p))):
  v=p[j];out[j-2]=v;p[j-1]^=v;p[j-2]^=m.gmul(kappa,v)
 return out,p[:2]


def product_rt(m,a,b,out,kappa):
 a0,a1=a[:12],a[12:];b0,b1=b[:12],b[12:]
 rt=m.prod_rt(out[:12],a0,b0)+m.prod_rt(out[12:],a0,b1)+m.prod_rt(out[12:],a1,b0)+m.prod_rt(out[12:],a1,b1)
 # A three-product Karatsuba form; the mixed product carries both low/high.
 rt=[];aa=m.spec_plain(a0);ab=m.spec_plain(a1);ba=m.spec_plain(b0);bb=m.spec_plain(b1)
 m.bil_terms(rt,m.KT,aa,ba,tws_of=lambda js:[out[j] for j in js]+[out[12+j] for j in js])
 def targets(js):
  v=0
  for j in js:v^=m.gmul(kappa,1<<j)
  return [out[j] for j in range(12) if v>>j&1]
 m.bil_terms(rt,m.KT,ab,bb,tws_of=targets)
 m.bil_terms(rt,m.KT,[x^y for x,y in zip(aa,ab)],[x^y for x,y in zip(ba,bb)],tws_of=lambda js:[out[12+j] for j in js])
 return rt


def multiply(b,E,word,coefficient,kappa,ghost,reverse=False):
 m=setup()[0];src=[q.id for q in word];co=[q.id for q in coefficient]
 if reverse:
  inv=E.A.open('remainder_inverse',24);init_extension_inverse(m,E,co,inv,kappa)
  out=E.A.open('remainder_head_output',24);E.mac(product_rt(m,src,inv,out,kappa))
  # Erase the old remainder from the recovered head and live denominator.
  E.free(product_rt(m,out,co,src,kappa))
  for q,bit in zip(out,ghost):E.b.z(q,bit)
  b.bit_store(ghost,False)
  for q,r in zip(src,out):E.b.swap(q,r)
  E.A.close('remainder_head_output',out)
  mask=b.alloc_bits(24);b.hmr(km.array([km.q(q) for q in inv]),mask);inverse_phase(b,E,co,mask,kappa);b.bit_store(mask,False);b.free(mask);E.A.close('remainder_inverse',inv)
 else:
  out=E.A.open('remainder_head_output',24);E.mac(product_rt(m,src,co,out,kappa));b.bit_store(ghost,False);b.hmr(word,ghost)
  for q,r in zip(src,out):E.b.swap(q,r)
  E.A.close('remainder_head_output',out)


def build(d,kappa,inverse=False):
 m=setup()[0];m.inv_word=bind(m);b=FastBuilder(km.CircuitBuilder());headers=[];rid=1
 def reg(n,name):
  nonlocal rid
  v=b.alloc_qubits(n)
  for q in v:headers.append(f'APPEND_TO_REGISTER q{q.id} r{rid}')
  headers.append(f'REGISTER r{rid} "{name}"');rid+=1;return v
 u=reg(12,'u');t=reg(12,'t');N=[reg(12,f'N{i}') for i in range(d+1)];B=[reg(12,f'B{i}') for i in range(d)];C=reg(1,'completed')[0];Q=reg(1,'boundary')[0];factor=reg(12,'factor');ghost=b.create_classical_register(24,'head_ghost');ab=Adapter(b);sc=ab.qubits(3);E=m.Em(ab,*sc,m.Alloc(ab))
 def wd(back=False):
  if not d:return
  def top():b.x(N[d][0]);xword(b,N[d-1],kappa)
  if not back:top()
  for j in (range(2,d+1) if back else reversed(range(2,d+1))):
   axpy(b,N[j],N[j-2],kappa,m);b.cx(N[j],N[j-1])
  if back:top()
 def nd(back=False):
  words=[u,t]+B
  for j in (range(2,d+2) if back else reversed(range(2,d+2))):
   axpy(b,words[j],words[j-2],kappa,m);b.cx(words[j],words[j-1])
 def change_w():
  for i in range(d+1):
   if i<d:b.cx(B[i],N[i])
   else:b.x(N[i][0])
   if i:b.cx(B[i-1],N[i])
  b.cx(C,N[d][0]);b.cx(Q,N[d][0])
 qc=quotient([0]*d+[kappa,1],kappa)[0]
 def change_q():
  for i in range(d):
   b.cx(N[1] if not i else N[i+1],B[i])
   for bit in range(12):
    if qc[i]>>bit&1:b.cx(C,B[i][bit]);b.cx(Q,B[i][bit])
 def affine(back=False):
  if not d:return
  T=[q.id for q in t];U=[q.id for q in u]
  with m.inv_word(E,T) as ti,m.opened(E,'remainder_scale',12) as tmp:
   for i,w in enumerate(B):
    def shear():
     if i==d-1:b.cx(u,w)
     else:E.mac(m.prod_rt([q.id for q in w],U,[q.id for q in N[i+2]]))
    if back:shear()
    m.mon_slot(E,[q.id for q in w],ti if back else T,T if back else ti,tmp)
    if not back:shear()
 def shift_t():b.cx(C,t[0]);caxpy(m,E,Q.id,[q.id for q in factor],[q.id for q in t],1)
 if inverse:nd()
 if not inverse:change_w()
 wd()
 extra=None
 if not d:extra=b.alloc_qubits(12);b.x(extra[0]);den=list(N[0])+list(extra)
 else:den=list(N[0])+list(N[1])
 if inverse:
  multiply(b,E,km.array(list(u)+list(t)),km.array(den),kappa,ghost,True);shift_t();affine(True);change_q()
 else:
  change_q();affine();shift_t();multiply(b,E,km.array(list(u)+list(t)),km.array(den),kappa,ghost)
 if extra is not None:b.x(extra[0]);b.free(extra)
 wd(True)
 if inverse:change_w()
 else:nd(True)
 ab.release_recycle('remainder_emitter',sc);c=b.finish_circuit();held={v.id for v in ghost};return km.Circuit('\n'.join(headers)+'\n'+str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits) if i not in held))


def test(d=3,shots=512):
 from schur_affine_step import reference
 m,k,fw,bk,tr=setup();rng=np.random.default_rng(142071+d);uu=rng.integers(0,4096,shots);tt=rng.integers(1,4096,shots);nn=rng.integers(0,4096,(shots,d+1));bb=rng.integers(0,4096,(shots,d));cc=np.arange(shots)%3==1;qq=np.arange(shots)%3==2;ff=rng.integers(1,4096,shots);nn[cc|qq]=0;bb[cc|qq]=0;uu[cc]=0;tt[cc]=1;tt[qq]=ff[qq];wn=[];wb=[]
 for u,t,n,v,c,q in zip(uu,tt,nn,bb,cc,qq):
  a,bv=reference(int(u),int(t),n.tolist(),v.tolist(),k,bool(c or q));wn.append(a);wb.append(bv)
 wn=np.array(wn);wb=np.array(wb);cf=build(d,k);ci=build(d,k,True);bits=lambda v,w=12:((v[None,:]>>np.arange(w)[:,None])&1).astype(bool);src={'u':uu,'t':tt,**{f'N{i}':nn[:,i] for i in range(d+1)},**{f'B{i}':bb[:,i] for i in range(d)}};dst={'u':wn[:,0],'t':wn[:,1],**{f'B{i}':wn[:,i+2] for i in range(d)},**{f'N{i}':wb[:,i] for i in range(d+1)}};rows=[]
 for inv,c,ss,dd in [(False,cf,src,dst),(True,ci,dst,src)]:
  sim=km.Simulator(batch_size=shots,seed=142071+inv);sim.use_same_registers_as(c)
  for name,v in ss.items():sim.write_across_shots(c.register_data[name],bits(v))
  for name,v in [('completed',cc),('boundary',qq),('factor',ff)]:sim.write_across_shots(c.register_data[name],bits(v,len(c.register_data[name])))
  sim.do(c);errors={}
  for name,v in dd.items():
   bad=np.flatnonzero((sim.read_across_shots(c.register_data[name])!=bits(v)).any(0))
   if len(bad):errors[name]=bad[:8].tolist()
  public={q for name,r in c.register_data.items() if name!='head_ghost' for q in r};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);junk=np.flatnonzero(sim.read_across_shots(scratch).any(0));mask=sim.read_across_shots(c.register_data['head_ghost']);head=np.r_[bits(uu),bits(tt^cc.astype(int)^(ff*qq))];expect=np.zeros(shots,dtype=bool) if inv else ((mask&head).sum(0)%2).astype(bool);phase=np.flatnonzero(sim.read_phase_flipped_across_shots()!=expect);r=dict(inverse=inv,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,errors=errors,junk=junk[:8].tolist(),phase=phase[:8].tolist());rows.append(r);print(json.dumps(r),flush=True);assert not errors and not len(junk) and not len(phase)
  for name,v in [('completed',cc),('boundary',qq),('factor',ff)]:assert np.array_equal(sim.read_across_shots(c.register_data[name]),bits(v,len(c.register_data[name])))
  if not inv:
   sim.do(ci)
   for name,v in src.items():assert np.array_equal(sim.read_across_shots(c.register_data[name]),bits(v)),name
   assert not sim.read_phase_flipped_across_shots().any();assert not sim.read_across_shots(c.register_data['head_ghost']).any();assert not sim.read_across_shots(scratch).any()
 result=dict(status='PASS',d=d,shots=2*shots,rows=rows,scope='new quadratic remainder layout, exact coefficients, all completion/boundary modes, forward pending head phase, fresh inverse and paired cleanup');Path(f'results/schur_remainder_step_d{d}_native.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--d',type=int,default=3);p.add_argument('--shots',type=int,default=512);a=p.parse_args();test(a.d,a.shots)
