"""Demolish a polynomial remainder while recovering its extension quotient.

Inputs satisfy num=head*den on the supported polynomial-fold chart. Their
measurement phase is repaid by that bilinear identity after quotient creation.
"""
import kickmix as km
from extension_inverse import setup,inverse_output_rt
from low_workspace_inverse import norm_rt as norm12_rt
from schur_dynamic_division_compact import numerator,norm_rt
from cubic_mask_phase_lean import native_phase_inverse
from tower_inverse import apply_map,map_ops


def power5_rt(m,src,target):
 images={m._elt3(i):i for i in range(8)};values=[images[m.gpow(m._elt3(i),5)] for i in range(8)]
 for j in range(3):
  for i in range(8):
   if i>>j&1:values[i]^=values[i^(1<<j)]
 rt=[]
 for mask,c in enumerate(values):
  if not c:continue
  dst=tuple(target[j] for j in range(3) if c>>j&1);ii=[i for i in range(3) if mask>>i&1];assert len(ii) in [1,2]
  rt.append(('L',(src[ii[0]],),0,dst) if len(ii)==1 else ('P',(src[ii[0]],),0,(src[ii[1]],),0,dst))
 return rt


def inverse_into(m,E,src,out):
 n0=E.A.open('destroy_N0',3);n1=E.A.open('destroy_N1',3);n=n0+n1;mm=E.A.open('destroy_M3',3);s=E.A.open('destroy_S3',3)
 rN,rM,rS=m.norm_chain_terms(src,m.spec_rows(m.A0ROWS,src),m.spec_rows(m.A1ROWS,src),n,mm,s);E.macs([rN,rM,rS]);E.free(power5_rt(m,s,mm));E.A.close('destroy_M3',mm)
 def portion(rt,targets):
  allowed=set(targets);return [t[:-1]+(tuple(q for q in t[-1] if q in allowed),) for t in rt if any(q in allowed for q in t[-1])]
 si=m.spec_rows(m.SQ3R,s);yn0=E.A.open('destroy_Y0',3);rt=[];m.bil_terms(rt,m.KT3,[a^b for a,b in zip(m.spec_plain(n0),m.spec_plain(n1))],si,yn0);E.mac(rt);E.free(portion(rN,n0));E.A.close('destroy_N0',n0)
 yn1=E.A.open('destroy_Y1',3);rt=[];m.bil_terms(rt,m.KT3,m.spec_plain(n1),si,yn1);E.mac(rt);E.free(portion(rN,n1));E.A.close('destroy_N1',n1);y=yn0+yn1
 apply_map(E.b.raw,km.array([km.q(q) for q in s]),map_ops(m.SQ3R));rr=[];m.bil_terms(rr,m.KT3,m.spec_plain(y[:3]),m.spec_plain(y[3:]),s);m.lin_rows_terms(rr,m.LMROWS,y,s);E.free(rr);E.A.close('destroy_S3',s)
 E.mac(inverse_output_rt(m,src,y,out));E.free(norm12_rt(m,out,y));E.A.close('destroy_Y1',yn1);E.A.close('destroy_Y0',yn0)


def phase_terms(E,rt,mask_by_target):
 for t in rt:
  for tw in t[-1]:
   bit=mask_by_target[tw]
   if t[0]=='L':
    for q in t[1]:E.b.z(q,bit)
    if t[2]:E.b.neg(bit)
   else:E.pivot_pair(t[1],t[2],t[3],t[4],lambda a,b:E.b.cz(a,b,bit))


def product_rt(m,left,right,out,kappa):
 a,c=left[:12],left[12:];d,e=right[:12],right[12:];rt=m.prod_rt(out[:12],a,d)+m.prod_rt(out[12:],a,e)+m.prod_rt(out[12:],c,d)
 def targets(js):
  v=m.gmul(kappa,sum(1<<j for j in js));return [out[j] for j in range(12) if v>>j&1]+[out[12+j] for j in js]
 m.bil_terms(rt,m.KT,m.spec_plain(c),m.spec_plain(e),tws_of=targets);return rt


def create(b,E,num,den,out,kappa):
 m=setup()[0]
 with numerator(b,E,num,den,out,kappa,[],keep=True):pass
 mask=b.alloc_bits(24);b.hmr(km.array([km.q(q) for q in num]),mask);nn,ni=num[:12],num[12:];nr=norm_rt(m,den,nn,kappa,[]);E.mac(nr);inverse_into(m,E,nn,ni)
 with m.opened(E,'destroy_scale',12) as tmp:
  for word in [out[:12],out[12:]]:m.mon_slot(E,word,ni,nn,tmp)
 im=b.alloc_bits(12);b.hmr(km.array([km.q(q) for q in ni]),im);native_phase_inverse(b,nn,im);b.bit_store(im,False);b.free(im);E.free(nr)
 phase_terms(E,product_rt(m,out,den,num,kappa),{q:bit for q,bit in zip(num,mask)});b.bit_store(mask,False);b.free(mask)


def test(shots=16384):
 import json,numpy as np
 from pathlib import Path
 from extension_inverse import Adapter
 from lean_emitter import create as emit
 from fast_dqi import FastBuilder
 m,k,*_=setup();b=FastBuilder(km.CircuitBuilder());nn=b.create_quantum_register(24,'num');dd=b.create_quantum_register(24,'den');hh=b.create_quantum_register(24,'head');ab=Adapter(b);sc,E=emit(m,ab);create(b,E,[q.id for q in nn],[q.id for q in dd],[q.id for q in hh],k);ab.release_recycle('destroy_emitter',sc);c=b.finish_circuit();c=km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)));rng=np.random.default_rng(142147);num=rng.integers(0,1<<24,shots);den=rng.integers(1,1<<24,shots);den[:4096]=np.arange(4096);den[4096:8192]=np.arange(4096)<<12;num[den==0]=0
 def quotient(n,d):
  a,bb=int(n)&4095,int(n)>>12;u,v=int(d)&4095,int(d)>>12;iv=m.ginv(m.gmul(u,u)^m.gmul(u,v)^m.gmul(k,m.gmul(v,v)));lo=m.gmul(a,u^v)^m.gmul(k,m.gmul(bb,v));hi=m.gmul(a,v)^m.gmul(bb,u);return m.gmul(lo,iv)|(m.gmul(hi,iv)<<12)
 want=np.array([quotient(n,d) for n,d in zip(num,den)]);bits=lambda x:((np.asarray(x)[None,:]>>np.arange(24)[:,None])&1).astype(bool);sim=km.Simulator(batch_size=shots,seed=142147);sim.use_same_registers_as(c);sim.write_across_shots(nn,bits(num));sim.write_across_shots(dd,bits(den));sim.do(c);assert not sim.read_across_shots(nn).any();assert np.array_equal(sim.read_across_shots(dd),bits(den));assert np.array_equal(sim.read_across_shots(hh),bits(want));assert not sim.read_phase_flipped_across_shots().any();public=set(nn)|set(dd)|set(hh);junk=km.array([km.q(q) for q in range(c.num_qubits) if km.q(q) not in public]);assert not sim.read_across_shots(junk).any();r=dict(status='PASS',shots=shots,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,scope='nonzero denominators, both subfield axes, and zero/zero; numerator demolished, denominator preserved, quotient and clean phase/scratch');Path('results/destructive_quotient_lean_native.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
if __name__=='__main__':test()
