"""Native parameter XOR, completion charts and arbitrary paired controls."""
from pathlib import Path
import argparse,json,hashlib,subprocess,time
import numpy as np
import kickmix as km
from import_published import load
from projective_polynomial_inverse import bind
import triple_gap_frame as triple
from check_continuant_head_frame import model as normal_model
import check_blocked_cfraction as ref
HERE=Path(__file__).resolve().parent

def build():
 m=load();m.inv_word=bind(m);b=m.Builder();W=[b.qubits(12) for _ in range(6)];flags=b.qubits(2);out=[b.qubits(12) for _ in range(6)];E=m.Em(b,*b.qubits(3),m.Alloc(b))
 with triple.parameters(m,E,W,flags) as (A,D,p):
  for src,dst in zip([*A[1:],*D[1:],p],out):
   for x,y in zip(src,dst):E.cx(x,y)
 rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9);tape=HERE/'parameters/triple_gap_frame.npy';np.save(tape,rows);kb=km.CircuitBuilder()
 for j in range(6):kb.create_quantum_register(12,f'W{j}')
 kb.create_quantum_register(2,'flags')
 for j in range(6):kb.create_quantum_register(12,f'O{j}')
 pre=HERE/'parameters/triple_gap_frame.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx');empty=HERE/'parameters/empty';outpath=HERE/'candidates/triple_gap_frame.kmx'
 subprocess.run([str(ref.probe_translator()),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(outpath),str(b.next_out),str(empty),str(len(rows))],check=True,capture_output=True);c=km.Circuit.from_file(str(outpath),format='kmx');c.to_file(str(outpath.with_suffix('.kmb')),format='kmb');return c,outpath.with_suffix('.kmb')

def test(shots=512):
 started=time.time();m=load();rng=np.random.default_rng(2026102623);heads=rng.integers(0,4096,(shots,6));controls=np.zeros((shots,2),dtype=bool);expected=[]
 for i,H in enumerate(heads):
  chart=i%3
  if chart==0:
   H[::2]=rng.integers(1,4096,3)
   if i%6==0:H[2*((i//6)%4):]=0
   expected.append(normal_model(m,H))
  elif chart==1:
   controls[i,0]=True;H[0]=0;z,u,h,e,a,b=map(int,H)
   expected.append([e^b,h^m.gmul(b,e),a^m.gmul(b,h),e,h,m.gmul(u,a or 1)])
  else:
   controls[i,1]=True;H[2]=0;a,b,z,u,h,e=map(int,H)
   expected.append([b^e,h^m.gmul(b,e),u^m.gmul(b,h),b,0,m.gmul(a,u or 1)])
 expected=np.asarray(expected);targets=rng.integers(0,4096,(shots,6));c,path=build();bits=lambda x:((np.asarray(x)[None,:]>>np.arange(12)[:,None])&1).astype(bool)
 for seed in (2026102624,2026102625):
  sim=km.Simulator(batch_size=shots,seed=seed);sim.use_same_registers_as(c)
  for j in range(6):sim.write_across_shots(c.register_data[f'W{j}'],bits(heads[:,j]));sim.write_across_shots(c.register_data[f'O{j}'],bits(targets[:,j]))
  sim.write_across_shots(c.register_data['flags'],controls.T);sim.do(c)
  for j in range(6):assert np.array_equal(sim.read_across_shots(c.register_data[f'W{j}']),bits(heads[:,j])),('input',j);assert np.array_equal(sim.read_across_shots(c.register_data[f'O{j}']),bits(targets[:,j]^expected[:,j])),('output',j)
  assert np.array_equal(sim.read_across_shots(c.register_data['flags']),controls.T)
  scratch=km.array([km.q(i) for i in range(146,c.num_qubits)]);assert not sim.read_across_shots(scratch).any(),'scratch';assert not sim.read_phase_flipped_across_shots().any(),'phase';sim.do(c)
  for j in range(6):assert np.array_equal(sim.read_across_shots(c.register_data[f'O{j}']),bits(targets[:,j])),('twice',j)
  assert not sim.read_across_shots(scratch).any();assert not sim.read_phase_flipped_across_shots().any()
 # Arbitrary marker/flag combinations assert only exact XOR involution.
 arbitrary=rng.integers(0,4096,(shots,6));ctrl=rng.integers(0,2,(shots,2),dtype=np.uint8).astype(bool);sim=km.Simulator(batch_size=shots,seed=2026102626);sim.use_same_registers_as(c)
 for j in range(6):sim.write_across_shots(c.register_data[f'W{j}'],bits(arbitrary[:,j]))
 sim.write_across_shots(c.register_data['flags'],ctrl.T);sim.do(c);assert not sim.read_across_shots(scratch).any();assert not sim.read_phase_flipped_across_shots().any();sim.do(c)
 for j in range(6):assert np.array_equal(sim.read_across_shots(c.register_data[f'W{j}']),bits(arbitrary[:,j]));assert not sim.read_across_shots(c.register_data[f'O{j}']).any()
 assert not sim.read_across_shots(scratch).any();assert not sim.read_phase_flipped_across_shots().any()
 r=dict(status='PASS_PARAMETER_COMPONENT',Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,regular_gap_completion_shots=2*shots,arbitrary_paired_shots=shots,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),seconds=time.time()-started,scope='three-slot cofactor parameter XOR; independent formulas, normal completion charts, both degree-two gap positions, arbitrary paired controls, phase/scratch; head extraction and full oracle unbuilt');(HERE/'results/triple_gap_frame_checks.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)

if __name__=='__main__':test()
