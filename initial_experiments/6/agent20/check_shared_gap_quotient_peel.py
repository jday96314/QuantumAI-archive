"""Native independent long division, total inverse, phase and scratch checks."""
import argparse,hashlib,json,subprocess,time
from pathlib import Path
import numpy as np
import kickmix as km
from import_published import load
from projective_polynomial_inverse import bind
import shared_gap_quotient_peel as shared
import check_blocked_cfraction as ref
HERE=Path(__file__).resolve().parent

def build(n,both):
 m=load();m.inv_word=bind(m);b=m.Builder();F=[b.qubits(12) for _ in range(n+1)];G=[b.qubits(12) for _ in range(n)];flag=b.qubits(1)[0]
 E=m.Em(b,*b.qubits(3),m.Alloc(b));shared.peel(m,E,F,G,flag)
 if both:shared.peel(m,E,F,G,flag,True)
 tag=f'shared_gap_quotient_n{n}'+('_both' if both else '');tape=HERE/'parameters'/f'{tag}.npy';rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9);np.save(tape,rows)
 kb=km.CircuitBuilder()
 for name,bank in (('F',F),('G',G)):
  for j in range(len(bank)):kb.create_quantum_register(12,f'{name}{j}')
 kb.create_quantum_register(1,'gap');pre=HERE/'parameters'/f'{tag}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx');empty=HERE/'parameters/empty';out=HERE/'candidates'/f'{tag}.kmx'
 subprocess.run([str(ref.probe_translator()),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True,capture_output=True)
 c=km.Circuit.from_file(str(out),format='kmx');c.to_file(str(out.with_suffix('.kmb')),format='kmb');return c,out.with_suffix('.kmb')

def division(m,f,g):
 r=list(f);degree=len(g)-1;Q=[0]*(len(f)-degree);iv=m.ginv(g[-1])
 for i in reversed(range(len(Q))):
  Q[i]=m.gmul(r[degree+i],iv)
  for j,v in enumerate(g):r[i+j]^=m.gmul(Q[i],int(v))
 return Q,r[:degree]

def test(n,shots=256):
 m=load();rng=np.random.default_rng(2026102601+n);rows=[];started=time.time()
 for both in (False,True):
  F=rng.integers(0,4096,(shots,n+1));G=rng.integers(0,4096,(shots,n));flags=(np.arange(shots)%2).astype(bool)
  if not both:
   G[~flags,-1]=rng.integers(1,4096,(~flags).sum());G[flags,-1]=0;G[flags,-2]=rng.integers(1,4096,flags.sum());F[:16]=0
  else:
   G[:shots//4,-1]=0;G[:shots//4,-2]=0
  expected=F.copy()
  if not both:
   for s,(ff,gg,flag) in enumerate(zip(F,G,flags)):
    q,r=division(m,ff,gg[:-1] if flag else gg);expected[s]=np.r_[r,q]
  c,path=build(n,both);bits=lambda x:((np.asarray(x)[None,:]>>np.arange(12)[:,None])&1).astype(bool)
  for seed in (2026102602,2026102603):
   sim=km.Simulator(batch_size=shots,seed=seed);sim.use_same_registers_as(c)
   for name,val in (('F',F),('G',G)):
    for j in range(val.shape[1]):sim.write_across_shots(c.register_data[f'{name}{j}'],bits(val[:,j]))
   sim.write_across_shots(c.register_data['gap'],flags[None]);sim.do(c)
   for name,val in (('F',expected),('G',G)):
    for j in range(val.shape[1]):assert np.array_equal(sim.read_across_shots(c.register_data[f'{name}{j}']),bits(val[:,j])),(n,both,name,j)
   assert np.array_equal(sim.read_across_shots(c.register_data['gap']),flags[None]),(n,both,'flag')
   public=12*(2*n+1)+1;scratch=km.array([km.q(i) for i in range(public,c.num_qubits)])
   assert not sim.read_across_shots(scratch).any(),(n,both,'scratch');assert not sim.read_phase_flipped_across_shots().any(),(n,both,'phase')
  rows.append(dict(both=both,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,shots=2*shots,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
 out=dict(status='PASS_COMPONENT',n=n,circuits=rows,seconds=time.time()-started,scope='shared degree1/2 tail product; independent field long division, arbitrary paired zero-pivot inverse, phase and scratch; not a whole Euclid schedule')
 print(json.dumps(out),flush=True);return out
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--degrees',type=int,nargs='+',default=[3,7,17,35]);a=p.parse_args();rows=[]
 for n in a.degrees:
  rows.append(test(n));(HERE/'results/shared_gap_quotient_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
