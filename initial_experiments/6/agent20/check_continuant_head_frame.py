"""Native arbitrary-target parameter XOR, phase/scratch and complete zero suffixes."""
import argparse,hashlib,json,subprocess,time
from pathlib import Path
import numpy as np
import kickmix as km
from import_published import load
from projective_polynomial_inverse import bind
from cube_arithmetic import CubeArithmetic
import continuant_head_frame as frame
import check_general_jacobi as ref
import check_blocked_cfraction as buildref
HERE=Path(__file__).resolve().parent

def build(b):
 m=load();m.inv_word=bind(m);bb=m.Builder();W=[bb.qubits(12) for _ in range(2*b)];out=[bb.qubits(12) for _ in range(3*b)];E=m.Em(bb,*bb.qubits(3),m.Alloc(bb));math=CubeArithmetic(m,E)
 with frame.parameters(m,E,math,W) as (A,D,p,P):
  for src,dst in zip([*A[1:],*D[1:],p,*P],out):
   for a,z in zip(src,dst):E.cx(a,z)
 tag=f'continuant_head_frame_b{b}';tape=HERE/'parameters'/f'{tag}.npy';rows=np.asarray(bb._items._rows,dtype=np.uint32).reshape(-1,9);np.save(tape,rows);kb=km.CircuitBuilder()
 for i in range(2*b):kb.create_quantum_register(12,f'W{i}')
 for i in range(3*b):kb.create_quantum_register(12,f'O{i}')
 pre=HERE/'parameters'/f'{tag}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx');empty=HERE/'parameters/empty';path=HERE/'candidates'/f'{tag}.kmx'
 subprocess.run([str(buildref.probe_translator()),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(path),str(bb.next_out),str(empty),str(len(rows))],check=True,capture_output=True)
 c=km.Circuit.from_file(str(path),format='kmx');c.to_file(str(path.with_suffix('.kmb')),format='kmb');return c,path.with_suffix('.kmb')

def model(m,H):
 b=len(H)//2;p=int(H[0]);old=[1];cur=[int(H[1]),1]
 for ai,bi in zip(H[2::2],H[3::2]):
  if ai:
   new=[0]*(len(cur)+1)
   for i,v in enumerate(cur):new[i]^=m.gmul(int(bi),v);new[i+1]^=v
   for i,v in enumerate(old):new[i]^=m.gmul(int(ai),v)
   old,cur=cur,new;p=m.gmul(p,int(ai))
  else:old,cur=[int(bi),*old],[0,*cur]
 return [*reversed(cur[:-1]),*reversed(old[:-1]),p]

def test(b,shots=128):
 started=time.time();m=load();rng=np.random.default_rng(2026102605+b);heads=rng.integers(1,4096,(shots,2*b));heads[:,1::2]=rng.integers(0,4096,(shots,b));heads[:8,1::2]=0
 for i in range(shots//2,shots):heads[i,2*(i%(b+1)):]=0
 expected=[]
 for H in heads:
  A,D,p,P=ref.parameters(m,H.tolist());expectP=P
  expected.append([*model(m,H),*expectP])
 expected=np.asarray(expected);targets=rng.integers(0,4096,(shots,3*b));c,path=build(b);bits=lambda x:((np.asarray(x)[None,:]>>np.arange(12)[:,None])&1).astype(bool)
 for seed in (2026102606,2026102607):
  sim=km.Simulator(batch_size=shots,seed=seed);sim.use_same_registers_as(c)
  for i in range(2*b):sim.write_across_shots(c.register_data[f'W{i}'],bits(heads[:,i]))
  for i in range(3*b):sim.write_across_shots(c.register_data[f'O{i}'],bits(targets[:,i]))
  sim.do(c)
  for i in range(2*b):assert np.array_equal(sim.read_across_shots(c.register_data[f'W{i}']),bits(heads[:,i])),(b,'source',i)
  for i in range(3*b):assert np.array_equal(sim.read_across_shots(c.register_data[f'O{i}']),bits(targets[:,i]^expected[:,i])),(b,'output',i)
  scratch=km.array([km.q(i) for i in range(60*b,c.num_qubits)]);assert not sim.read_across_shots(scratch).any(),(b,'scratch');assert not sim.read_phase_flipped_across_shots().any(),(b,'phase')
  sim.do(c)
  for i in range(3*b):assert np.array_equal(sim.read_across_shots(c.register_data[f'O{i}']),bits(targets[:,i])),(b,'twice',i)
  assert not sim.read_phase_flipped_across_shots().any(),(b,'twice phase');assert not sim.read_across_shots(scratch).any(),(b,'twice scratch')
 # Arbitrary zero/head states have no affine correctness claim, but the
 # actual parameter XOR still restores sources/scratch/phase and is an involution.
 arbitrary=rng.integers(0,4096,(shots,2*b));arbitrary[:shots//2,2::2]=0;sim=km.Simulator(batch_size=shots,seed=2026102608);sim.use_same_registers_as(c)
 for i in range(2*b):sim.write_across_shots(c.register_data[f'W{i}'],bits(arbitrary[:,i]))
 sim.do(c);assert not sim.read_phase_flipped_across_shots().any();assert not sim.read_across_shots(scratch).any();sim.do(c)
 for i in range(2*b):assert np.array_equal(sim.read_across_shots(c.register_data[f'W{i}']),bits(arbitrary[:,i]))
 for i in range(3*b):assert not sim.read_across_shots(c.register_data[f'O{i}']).any()
 assert not sim.read_phase_flipped_across_shots().any();assert not sim.read_across_shots(scratch).any()
 out=dict(status='PASS_COMPONENT',pairs=b,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,regular_and_completion_shots=2*shots,arbitrary_paired_shots=shots,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),seconds=time.time()-started,scope='continuant frame and affine parameters, arbitrary-target XOR, all completion degrees, phase and scratch; interior gaps unsupported')
 print(json.dumps(out),flush=True);return out
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--pairs',type=int,nargs='+',default=[2,3,4]);a=p.parse_args();rows=[]
 for b in a.pairs:
  rows.append(test(b));(HERE/'results/continuant_head_frame_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
