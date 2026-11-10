"""Twelve-AND controlled field-linear maps with measured input erasure."""
import json
from pathlib import Path
import numpy as np
import kickmix as km
from extension_inverse import Adapter,setup
from fast_dqi import FastBuilder


def caxpy_rt(m,control,source,target,c):
 rt=[]
 for i,q in enumerate(source):
  v=m.gmul(c,1<<i);dest=[target[j] for j in range(12) if v>>j&1]
  if dest:rt.append(('P',(control,),0,(q,),0,tuple(dest)))
 return rt


def caxpy(m,E,control,source,target,c):E.mac(caxpy_rt(m,control,source,target,c))


def cimul(m,E,control,word,c,inverse=False):
 if c==1:return
 if inverse:c=m.ginv(c)
 with m.opened(E,'controlled_constant_output',12) as out:
  rt=[('L',(q,),0,(t,)) for q,t in zip(word,out)]+caxpy_rt(m,control,word,out,c^1);E.mac(rt)
  undo=[('L',(q,),0,(t,)) for q,t in zip(out,word)]+caxpy_rt(m,control,out,word,m.ginv(c)^1);E.free(undo)
  for a,b in zip(word,out):E.b.swap(a,b)


def test():
 m,k,fw,bk,tr=setup();b=FastBuilder(km.CircuitBuilder());flag=b.create_quantum_register(1,'flag');x=b.create_quantum_register(12,'x');y=b.create_quantum_register(12,'y');ab=Adapter(b);sc=ab.qubits(3);E=m.Em(ab,*sc,m.Alloc(ab));X=[q.id for q in x];Y=[q.id for q in y];caxpy(m,E,flag[0].id,X,Y,137);cimul(m,E,flag[0].id,X,2533);ab.release_recycle('control_emitter',sc);c=b.finish_circuit();c=km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)));shots=8192;xx=np.tile(np.arange(4096),2);ff=np.repeat([0,1],4096);rng=np.random.default_rng(142049);yy=rng.integers(0,4096,shots);wx=np.array([m.gmul(int(v),2533) if f else v for v,f in zip(xx,ff)]);wy=np.array([int(v)^(m.gmul(int(u),137) if f else 0) for u,v,f in zip(xx,yy,ff)]);bits=lambda v,w=12:((v[None,:]>>np.arange(w)[:,None])&1).astype(bool);sim=km.Simulator(batch_size=shots,seed=142049);sim.use_same_registers_as(c)
 for reg,v in [(flag,ff),(x,xx),(y,yy)]:sim.write_across_shots(reg,bits(v,len(reg)))
 sim.do(c)
 for reg,v in [(flag,ff),(x,wx),(y,wy)]:assert np.array_equal(sim.read_across_shots(reg),bits(v,len(reg)))
 public=set(flag)|set(x)|set(y);sc=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);assert not sim.read_across_shots(sc).any();assert not sim.read_phase_flipped_across_shots().any();r=dict(status='PASS',shots=shots,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,scope='controlled linear axpy and in-place scaling; both control values and all source elements, arbitrary target, scratch and phase');Path('results/schur_control_native.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
if __name__=='__main__':test()
