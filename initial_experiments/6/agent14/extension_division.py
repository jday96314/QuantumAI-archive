"""Native quadratic-extension division creation and measured-output phase.

The extension polynomial is x^2+x+kappa. Norm/inverse contraction retains
only base-field scratch. Existing attributed multiplier/tower helpers emit
all arithmetic. These are component kernels, not complete DQI decoders.
"""
import json
from pathlib import Path
import numpy as np
import kickmix as km
from extension_inverse import Adapter,setup,extension_norm_rt,classical_transpose
from low_workspace_inverse import bind
from tower_inverse import map_ops,apply_map
from recursive_forney import inverse_dot_data
from fast_dqi import FastBuilder


def numerator_rt(m,num,den,out,kappa):
 # num * conjugate(den), using three base-field multiplications.
 n0,n1=num[:12],num[12:];d0,d1=den[:12],den[12:]
 p=m.spec_plain(n0);q=m.spec_plain(n1);r=[a^b for a,b in zip(p,q)]
 s=m.spec_plain(d0);t=m.spec_plain(d1);u=[a^b for a,b in zip(s,t)]
 rt=[]
 m.bil_terms(rt,m.KT,p,u,tws_of=lambda js:[out[j] for j in js]+[out[12+j] for j in js])
 def ktargets(js):
  v=0
  for j in js:v^=m.gmul(kappa,1<<j)
  return [out[j] for j in range(12) if v>>j&1]
 m.bil_terms(rt,m.KT,q,t,tws_of=ktargets)
 m.bil_terms(rt,m.KT,r,s,tws_of=lambda js:[out[12+j] for j in js])
 return rt


def create(b,E,num,den,out,kappa):
 m=setup()[0];m.inv_word=bind(m)
 with m.opened(E,'division_norm',12) as nn,m.opened(E,'division_numerator',24) as pp:
  nr=extension_norm_rt(m,den,nn,kappa);pr=numerator_rt(m,num,den,pp,kappa)
  E.mac(nr);E.mac(pr)
  with m.inv_word(E,nn) as ni:
   E.mac(m.prod_rt(out[:12],pp[:12],ni));E.mac(m.prod_rt(out[12:],pp[12:],ni))
  E.free(pr);E.free(nr)


def phase(b,E,num,den,mask,kappa):
 m,_,fw,bk,tr=setup()
 with m.opened(E,'division_phase_norm',12) as nn,m.opened(E,'division_phase_numerator',24) as pp:
  nr=extension_norm_rt(m,den,nn,kappa);pr=numerator_rt(m,num,den,pp,kappa);E.mac(nr);E.mac(pr)
  ns=km.array([km.q(i) for i in nn]);ps=km.array([km.q(i) for i in pp]);bops=map_ops(bk)
  trans=tuple(sum(((fw[i]>>j)&1)<<i for i in range(12)) for j in range(12));mops=map_ops(trans)
  for word in [ns,ps[:12],ps[12:]]:apply_map(b,word,bops)
  for word in [mask[:12],mask[12:]]:apply_map(b,word,mops)
  u=b.alloc_qubits(12)
  def umap():
   classical_transpose(b,ps[:12],mask[:12],u);classical_transpose(b,ps[12:],mask[12:],u)
  umap();leaf,inv,cleanup=inverse_dot_data(b,ns,u);b.cz(leaf,inv);cleanup();umap();b.free(u)
  for word in [mask[:12],mask[12:]]:apply_map(b,word,mops,True)
  for word in [ns,ps[:12],ps[12:]]:apply_map(b,word,bops,True)
  E.free(pr);E.free(nr)


def build(phase_only=False):
 m,k,fw,bk,tr=setup();b=FastBuilder(km.CircuitBuilder());num=b.create_quantum_register(24,'num');den=b.create_quantum_register(24,'den')
 out=b.create_classical_register(24,'mask') if phase_only else b.create_quantum_register(24,'out')
 ab=Adapter(b);sc=ab.qubits(3);E=m.Em(ab,*sc,m.Alloc(ab))
 if phase_only:phase(b,E,[q.id for q in num],[q.id for q in den],out,k)
 else:create(b,E,[q.id for q in num],[q.id for q in den],[q.id for q in out],k)
 ab.release_recycle('division_emitter',sc);c=b.finish_circuit();public={v.id for v in out} if phase_only else set()
 return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits) if i not in public))


def test():
 m,k,fw,bk,tr=setup();rng=np.random.default_rng(142043);shots=4096;aa=rng.integers(0,1<<24,shots);dd=rng.integers(0,1<<24,shots);dd[0]=0;mm=rng.integers(0,1<<24,shots)
 def mul(x,y):
  a,b=x&4095,x>>12;c,d=y&4095,y>>12
  return (m.gmul(a,c)^m.gmul(k,m.gmul(b,d)))|((m.gmul(a,d)^m.gmul(b,c)^m.gmul(b,d))<<12)
 def inverse(x):
  a,b=x&4095,x>>12;n=m.gmul(a,a)^m.gmul(a,b)^m.gmul(k,m.gmul(b,b));v=m.ginv(n)
  return m.gmul(a^b,v)|(m.gmul(b,v)<<12)
 want=np.array([mul(int(x),inverse(int(y))) for x,y in zip(aa,dd)]);bits=lambda v:((v[None,:]>>np.arange(24)[:,None])&1).astype(bool);rows=[]
 for po in [False,True]:
  c=build(po);sim=km.Simulator(batch_size=shots,seed=142043+po);sim.use_same_registers_as(c);regs=c.register_data
  for name,v in [('num',aa),('den',dd)]+([('mask',mm)] if po else []):sim.write_across_shots(regs[name],bits(v))
  sim.do(c);assert np.array_equal(sim.read_across_shots(regs['num']),bits(aa));assert np.array_equal(sim.read_across_shots(regs['den']),bits(dd))
  if po:
   expected=np.array([(int(x)&int(y)).bit_count()%2 for x,y in zip(want,mm)],dtype=bool);assert np.array_equal(sim.read_across_shots(regs['mask']),bits(mm))
  else:expected=np.zeros(shots,dtype=bool);assert np.array_equal(sim.read_across_shots(regs['out']),bits(want))
  bad=np.flatnonzero(sim.read_phase_flipped_across_shots()!=expected);assert not len(bad),bad[:10]
  public={q for name,word in regs.items() if name!='mask' for q in word};sc=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);assert not sim.read_across_shots(sc).any()
  r=dict(phase_only=po,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,shots=shots);rows.append(r);print(json.dumps(r),flush=True)
 Path('results/extension_division_native.json').write_text(json.dumps(dict(status='PASS',rows=rows,scope='arbitrary extension numerator/divisor, zero divisor included; exact values, phase, scratch'),indent=2)+'\n')
if __name__=='__main__':test()
