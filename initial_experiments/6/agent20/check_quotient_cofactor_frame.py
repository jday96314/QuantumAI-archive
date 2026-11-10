"""Native parameter XOR with independent variable-degree continuants."""
from pathlib import Path
import argparse,json,hashlib,subprocess,time
import numpy as np
import kickmix as km
from import_published import load
from projective_polynomial_inverse import bind
import quotient_cofactor_frame as frame
import check_blocked_cfraction as ref
HERE=Path(__file__).resolve().parent

def build(n):
 m=load();m.inv_word=bind(m,limit=1024);b=m.Builder();W=[b.qubits(12) for _ in range(2*n)];out=[b.qubits(12) for _ in range(2*n)];E=m.Em(b,*b.qubits(3),m.Alloc(b))
 with frame.parameters(m,E,W) as (A,D,p):
  for src,dst in zip([*A[1:],*D[1:],p],out):
   for a,z in zip(src,dst):E.cx(a,z)
 rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9);tag=f'quotient_cofactor_n{n}';tape=HERE/'parameters'/f'{tag}.npy';np.save(tape,rows);kb=km.CircuitBuilder()
 for j in range(2*n):kb.create_quantum_register(12,f'W{j}')
 for j in range(2*n):kb.create_quantum_register(12,f'O{j}')
 pre=HERE/'parameters'/f'{tag}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx');empty=HERE/'parameters/empty';outpath=HERE/'candidates'/f'{tag}.kmx';subprocess.run([str(ref.probe_translator()),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(outpath),str(b.next_out),str(empty),str(len(rows))],check=True,capture_output=True);c=km.Circuit.from_file(str(outpath),format='kmx');c.to_file(str(outpath.with_suffix('.kmb')),format='kmb');return c,outpath.with_suffix('.kmb')

def model(m,H):
 packets=[]
 for a,b,h,e in zip(H[::4],H[1::4],H[2::4],H[3::4]):
  if not a:packets.append((int(b),[int(h),int(e),1]))
  else:
   d=int(e)^int(b);c=int(h)^m.gmul(int(b),d);packets.extend([(int(a),[int(b),1]),(c,[d,1])])
 P,Q=[0],[1];p=1
 for a,q in packets:
  out=[0]*(len(Q)+len(q)-1)
  for i,v in enumerate(Q):
   for j,w in enumerate(q):out[i+j]^=m.gmul(v,w)
  for i,v in enumerate(P):out[i]^=m.gmul(a,v)
  P,Q=Q,out;p=m.gmul(p,a)
 n=len(H)//2;return [*reversed(Q[:-1]),*reversed(P[:-1]),*([0]*(n-len(P))),p]

def test(n,shots=256):
 started=time.time();m=load();rng=np.random.default_rng(2026102635+n);values=rng.integers(0,4096,(shots,2*n))
 for s,H in enumerate(values):
  for i in range(0,2*n,4):
   a,b,c,d=map(int,rng.integers(1,4096,4))
   H[i:i+4]=[0,b,c,d] if (s+i)%3==0 else [a,b,c^m.gmul(b,d),b^d]
 expected=np.asarray([model(m,H) for H in values]);targets=rng.integers(0,4096,(shots,2*n));c,path=build(n);bits=lambda x:((np.asarray(x)[None,:]>>np.arange(12)[:,None])&1).astype(bool);scratch=km.array([km.q(i) for i in range(48*n,c.num_qubits)])
 for seed in (2026102636,2026102637):
  sim=km.Simulator(batch_size=shots,seed=seed);sim.use_same_registers_as(c)
  for j in range(2*n):sim.write_across_shots(c.register_data[f'W{j}'],bits(values[:,j]));sim.write_across_shots(c.register_data[f'O{j}'],bits(targets[:,j]))
  sim.do(c)
  for j in range(2*n):assert np.array_equal(sim.read_across_shots(c.register_data[f'W{j}']),bits(values[:,j])),(n,'source',j);assert np.array_equal(sim.read_across_shots(c.register_data[f'O{j}']),bits(targets[:,j]^expected[:,j])),(n,'output',j)
  assert not sim.read_across_shots(scratch).any(),(n,'scratch');assert not sim.read_phase_flipped_across_shots().any(),(n,'phase');sim.do(c)
  for j in range(2*n):assert np.array_equal(sim.read_across_shots(c.register_data[f'O{j}']),bits(targets[:,j]))
  assert not sim.read_across_shots(scratch).any();assert not sim.read_phase_flipped_across_shots().any()
 # Singular/completed heads promise restoration, not these generic formulas.
 arbitrary=rng.integers(0,4096,(shots,2*n));arbitrary[:shots//2,::2]=0;arbitrary[-1]=0;sim=km.Simulator(batch_size=shots,seed=2026102638);sim.use_same_registers_as(c)
 for j in range(2*n):sim.write_across_shots(c.register_data[f'W{j}'],bits(arbitrary[:,j]))
 sim.do(c);assert not sim.read_across_shots(scratch).any(),(n,'arbitrary scratch');assert not sim.read_phase_flipped_across_shots().any(),(n,'arbitrary phase');sim.do(c)
 for j in range(2*n):assert np.array_equal(sim.read_across_shots(c.register_data[f'W{j}']),bits(arbitrary[:,j]));assert not sim.read_across_shots(c.register_data[f'O{j}']).any()
 assert not sim.read_across_shots(scratch).any();assert not sim.read_phase_flipped_across_shots().any()
 r=dict(status='PASS_COFACTOR_COMPONENT',n=n,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,regular_gap_shots=2*shots,arbitrary_paired_shots=shots,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),seconds=time.time()-started,scope='reversed quotient cofactor A/D/p with retained gauge word; independent continuants, parameter XOR, arbitrary paired singular heads, phase/scratch; prefix completion decoding not claimed');print(json.dumps(r),flush=True);return r

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,nargs='+',default=[2,4,6]);a=p.parse_args();rows=[]
 for n in a.n:rows.append(test(n));(HERE/'results/quotient_cofactor_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
