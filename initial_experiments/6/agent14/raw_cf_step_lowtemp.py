"""Unnormalized continued-fraction fold N <- (a+b*x)N+q*B, B <- N.

Two consumed heads provide the two added polynomial coefficients. This first
implementation uses two transient output fields and deferred head phases.
"""
from contextlib import contextmanager
import argparse,json
from pathlib import Path
import numpy as np
import kickmix as km
from extension_inverse import Adapter,setup
from fast_dqi import FastBuilder
from lean_emitter import create as lean_create
from schur_gap_fold import quadratic
from triangular_dp import DPArithmetic
from schur_fold import axpy
from schur_division_lowtemp import create as quotient


@contextmanager
def evaluate(b,N,B,kappa):
 m=setup()[0];power=[(1,0),(0,1)]
 for i in range(2,max(len(N),len(B))):
  a,c=power[-1];power.append((m.gmul(kappa,c),a^c))
 def apply():
  for bank in [N,B]:
   for i in range(2,len(bank)):
    axpy(b,bank[i],bank[0],power[i][0],m);axpy(b,bank[i],bank[1],power[i][1],m)
 apply();yield [q.id for w in N[:2] for q in w],[q.id for w in B[:2] for q in w];apply()


def build(size,kappa,inverse=False):
 assert size>=2;m=setup()[0];b=FastBuilder(km.CircuitBuilder());a=b.create_quantum_register(12,'a');h=b.create_quantum_register(12,'b');N=[b.create_quantum_register(12,f'N{i}') for i in range(size)];B=[b.create_quantum_register(12,f'B{i}') for i in range(size-1)];high=[b.create_quantum_register(12,f'Bt{i}') for i in range(2)];ghost=b.create_classical_register(24,'head_ghost');ab=Adapter(b);sc,E=lean_create(m,ab);dest=B+high;math=DPArithmetic(m,E);pts,I=m.vinv(size+1)
 def update():
  for col,x in enumerate(pts):
   left=math.spec([(1,[q.id for q in a],1),(x,[q.id for q in h],1)]);right=math.spec([(m.gpow(x,i),[q.id for q in w],1) for i,w in enumerate(N)])
   def targets(tb):return [dest[i][k].id for i in range(size+1) for k in m.bitlist(m.gmul(I[i][col],m.tmask(tb)),12)]
   math.product(left,right,targets)
 if inverse:
  with evaluate(b,dest,N,kappa) as (nv,dv):quotient(b,E,nv,dv,[q.id for q in list(a)+list(h)],kappa,[])
  for q,bit in zip(list(a)+list(h),ghost):ab.z(q.id,bit)
  b.bit_store(ghost,False);update();quadratic(b,dest,kappa,True)
 else:
  b.bit_store(ghost,False);quadratic(b,dest,kappa);update();b.hmr(km.array(list(a)+list(h)),ghost)
 ab.release_recycle('raw_cf_step_emitter',sc);c=b.finish_circuit();return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(24,c.num_bits)))


def test(size=3,shots=512):
 m,k,fw,bk,tr=setup();rng=np.random.default_rng(142114+size);aa,bb=rng.integers(0,4096,(2,shots));N=rng.integers(0,4096,(shots,size));B=rng.integers(0,4096,(shots,size-1));aa[::11]=bb[::11]=0;N[::11]=0;want=[]
 for a,h,nn,dd in zip(aa,bb,N,B):
  out=[0]*(size+1)
  for i,v in enumerate(nn):out[i]^=m.gmul(int(a),int(v));out[i+1]^=m.gmul(int(h),int(v))
  for i,v in enumerate(dd):out[i]^=m.gmul(k,int(v));out[i+1]^=int(v);out[i+2]^=int(v)
  want.append(out)
 want=np.array(want);cf=build(size,k);ci=build(size,k,True);bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool);src={'a':aa,'b':bb,**{f'N{i}':N[:,i] for i in range(size)},**{f'B{i}':B[:,i] for i in range(size-1)},'Bt0':np.zeros(shots,dtype=int),'Bt1':np.zeros(shots,dtype=int)};dst={'a':np.zeros(shots,dtype=int),'b':np.zeros(shots,dtype=int),**{f'N{i}':N[:,i] for i in range(size)},**{f'B{i}':want[:,i] for i in range(size-1)},**{f'Bt{i}':want[:,size-1+i] for i in range(2)}};rows=[]
 for inv,c,s,d in [(False,cf,src,dst),(True,ci,dst,src)]:
  sim=km.Simulator(batch_size=shots,seed=142114+inv);sim.use_same_registers_as(c)
  for key,val in s.items():sim.write_across_shots(c.register_data[key],bits(val))
  sim.do(c);errors=[key for key,val in d.items() if not np.array_equal(sim.read_across_shots(c.register_data[key]),bits(val))];pub={q for name,reg in c.register_data.items() if name!='head_ghost' for q in reg};junk=km.array([km.q(q) for q in range(c.num_qubits) if km.q(q) not in pub]);mask=sim.read_across_shots(c.register_data['head_ghost']);phase=np.zeros(shots,dtype=bool) if inv else ((mask&np.concatenate([bits(aa),bits(bb)])).sum(0)%2).astype(bool);r=dict(inverse=inv,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,errors=errors,scratch=bool(sim.read_across_shots(junk).any()),phase=bool((sim.read_phase_flipped_across_shots()!=phase).any()));rows.append(r);print(json.dumps(r),flush=True);assert not errors and not r['scratch'] and not r['phase']
  if not inv:
   sim.do(ci)
   for key,val in s.items():assert np.array_equal(sim.read_across_shots(c.register_data[key]),bits(val))
   assert not sim.read_across_shots(junk).any();assert not sim.read_across_shots(c.register_data['head_ghost']).any();assert not sim.read_phase_flipped_across_shots().any()
 Path(f'results/raw_cf_step_lowtemp_m{size}_native.json').write_text(json.dumps(dict(status='PASS',shots=shots,rows=rows),indent=2)+'\n')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--size',type=int,default=3);a=p.parse_args();test(a.size)
