"""Independent recurrence, padding, native phase/scratch/literal inverse checks."""
from pathlib import Path
import argparse,json,subprocess,time,hashlib
import numpy as np
import kickmix as km
from import_published import load
import projective_polynomial_inverse as projective
import completed_quotient_jacobi as complete
import check_blocked_cfraction as ref
HERE=Path(__file__).resolve().parent


def build(n,both=False):
 m=load();m.inv_word=projective.bind(m);b=m.Builder();W=[b.qubits(12) for _ in range(2*n)]
 E=m.Em(b,*b.qubits(3),m.Alloc(b));c=E.A.open('completed_quotient_terminal',12);flags=E.A.open('completed_quotient_flags',n)
 F,G=complete.synthesize(m,E,W,c,flags)
 if both:
  complete.synthesize(m,E,W,c,flags,True);E.A.close('completed_quotient_flags',flags);E.A.close('completed_quotient_terminal',c)
 tag=f'completed_quotient_jacobi_n{n}'+('_both' if both else '');tape=HERE/'parameters'/f'{tag}.npy'
 rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9);np.save(tape,rows);kb=km.CircuitBuilder()
 for j in range(2*n):kb.create_quantum_register(12,f'W{j}')
 pre=HERE/'parameters'/f'{tag}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx');empty=HERE/'parameters/empty';out=HERE/'candidates'/f'{tag}.kmx'
 subprocess.run([str(ref.probe_translator()),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True,capture_output=True)
 native=km.Circuit.from_file(str(out),format='kmx');native.to_file(str(out.with_suffix('.kmb')),format='kmb');return native,F,G,flags


def test(n,shots=256):
 m=load();rng=np.random.default_rng(2026102521+n);started=time.time();reports=[]
 def product(a,b):
  if not a or not b:return []
  out=[0]*(len(a)+len(b)-1)
  for i,x in enumerate(a):
   for j,y in enumerate(b):out[i+j]^=m.gmul(int(x),int(y))
  return out
 def trim(x):
  while len(x)>1 and x[-1]==0:x.pop()
  return x
 for both in (False,True):
  values=rng.integers(1,4096,(shots,2*n));degrees=np.arange(shots)%(n+1)
  for j,d in enumerate(degrees):values[j,2*d:]=0
  if both:values[shots//2:]=rng.integers(0,4096,(shots-shots//2,2*n))
  c,F,G,flags=build(n,both);expected=[complete.model(m,row) for row in values] if not both else None
  if not both:
   for raw,d,(ff,gg,fl) in zip(values,degrees,expected):
    P,Q=[],[1]
    for a,b in reversed(list(zip(raw[:2*d:2],raw[1:2*d:2]))):
     nq=product([1,int(b)],Q)
     if P:
      if len(nq)<len(P)+2:nq += [0]*(len(P)+2-len(nq))
      for i,v in enumerate(P):nq[i+2]^=v
     P,Q=[m.gmul(int(a),v) for v in Q],nq
    # Cross multiplication is invariant under the quantum common scale.
    lhs=trim(product(list(reversed(gg)),Q));rhs=trim(product(P,list(reversed(ff))))
    assert (lhs or [0])==(rhs or [0]),('rational identity',n,int(d),lhs,rhs)
    assert not any(ff[:n-d]) and not any(gg[:n-d]),('common padding',n,int(d))
  bits=lambda x:((np.asarray(x)[None,:]>>np.arange(12)[:,None])&1).astype(bool)
  for seed in (2026102522,2026102523):
   sim=km.Simulator(batch_size=shots,seed=seed);sim.use_same_registers_as(c)
   for j in range(2*n):sim.write_across_shots(c.register_data[f'W{j}'],bits(values[:,j]))
   sim.do(c)
   if both:
    for j in range(2*n):assert np.array_equal(sim.read_across_shots(c.register_data[f'W{j}']),bits(values[:,j])),(n,'inverse',j)
    live=set(range(24*n))
   else:
    for k,bank in enumerate((F,G)):
     for j,word in enumerate(bank):assert np.array_equal(sim.read_across_shots(km.array([km.q(i) for i in word])),bits([row[k][j] for row in expected])),(n,k,j)
    assert np.array_equal(sim.read_across_shots(km.array([km.q(i) for i in flags])),np.asarray([row[2] for row in expected],dtype=bool).T),(n,'flags')
    live={i for bank in (F,G) for word in bank for i in word}|set(flags)
   scratch=km.array([km.q(i) for i in range(c.num_qubits) if i not in live]);assert not sim.read_across_shots(scratch).any(),(n,both,'scratch')
   assert not sim.read_phase_flipped_across_shots().any(),(n,both,'phase')
  path=HERE/'candidates'/f'completed_quotient_jacobi_n{n}{"_both" if both else ""}.kmb'
  reports.append(dict(both=both,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,shots=2*shots,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
 r=dict(status='PASS_COMPONENT',n=n,circuits=reports,seconds=time.time()-started,scope='all completion degrees, independent Jacobi recurrence cross multiplication, monomial padding, arbitrary paired restoration, native phase and scratch; interior gaps unsupported')
 print(json.dumps(r),flush=True);return r
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--degrees',type=int,nargs='+',default=[3,7,17,35]);a=p.parse_args();rows=[]
 for n in a.degrees:
  rows.append(test(n));(HERE/'results/completed_quotient_jacobi_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
