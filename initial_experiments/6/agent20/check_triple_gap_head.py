"""Native six-moment variable-degree head, with independent series model."""
from pathlib import Path
import inspect,json,hashlib,subprocess,time
import numpy as np
import kickmix as km
from import_published import load
from projective_polynomial_inverse import bind
import triple_gap_head as head
import conditional_ccz_builder
import check_lookahead_jacobi as look
import check_jacobi
from check_high_copy_jacobi import head as regular
HERE=Path(__file__).resolve().parent
source=inspect.getsource(look.model).replace('p=(m.gmul(a,c)^bg) or 1','p=(b if flag else m.gmul(a,c or 1)) or 1');namespace=dict(look.__dict__);exec(source,namespace);look_model=namespace['model']

def model(m,H,inverse=False):
 H=list(H)
 if not H[0]:
  if inverse:H[4:]=regular(m,H[4:],True);H=look_model(m,H,True)
  else:H=look_model(m,H);H[4:]=check_jacobi.pair(m,H[4:])
 else:
  if inverse:
   if H[2]:H[4]^=m.gmul(H[3],H[5]);H[5]^=H[3]
   H[2:]=look_model(m,H[2:],True);H=regular_one(m,H,True)
  else:
   H=check_jacobi.pair(m,H);H[2:]=look_model(m,H[2:])
   if H[2]:H[5]^=H[3];H[4]^=m.gmul(H[3],H[5])
 return H

def regular_one(m,H,inverse=False):
 if not inverse:return check_jacobi.pair(m,H)
 a,*tail=H;reciprocal=look.check.divide(m,[1],[1,*tail],len(H));return [a,*[m.gmul(a or 1,v) for v in reciprocal[1:]]]

def build(inverse=False):
 m=load();conditional_ccz_builder.install(m);m.inv_word=bind(m);b=m.Builder();W=[b.qubits(12) for _ in range(6)];E=m.Em(b,*b.qubits(3),m.Alloc(b));head.apply(m,E,W,inverse)
 tag='triple_gap_head_'+('inverse' if inverse else 'forward');rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9);tape=HERE/'parameters'/f'{tag}.npy';np.save(tape,rows);kb=km.CircuitBuilder()
 for j in range(6):kb.create_quantum_register(12,f'W{j}')
 pre=HERE/'parameters'/f'{tag}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx');empty=HERE/'parameters/empty';out=HERE/'candidates'/f'{tag}.kmx'
 subprocess.run([str(look.check.probe_translator()),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True,capture_output=True);c=km.Circuit.from_file(str(out),format='kmx');c.to_file(str(out.with_suffix('.kmb')),format='kmb');return c,out.with_suffix('.kmb')

def test(shots=512):
 started=time.time();m=load();rng=np.random.default_rng(2026102627);values=rng.integers(0,4096,(shots,6));values[:3]=[[0]*6,[1]*6,[1,0,0,0,0,0]]
 values[shots//3:2*shots//3,0]=0;expected=np.asarray([model(m,row) for row in values]);assert all(model(m,row,True)==old.tolist() for row,old in zip(expected,values));bits=lambda x:((np.asarray(x)[None,:]>>np.arange(12)[:,None])&1).astype(bool);reports=[]
 for inverse in (False,True):
  c,path=build(inverse);source,want=(expected,values) if inverse else (values,expected)
  for seed in (2026102628,2026102629):
   sim=km.Simulator(batch_size=shots,seed=seed);sim.use_same_registers_as(c)
   for j in range(6):sim.write_across_shots(c.register_data[f'W{j}'],bits(source[:,j]))
   sim.do(c)
   for j in range(6):assert np.array_equal(sim.read_across_shots(c.register_data[f'W{j}']),bits(want[:,j])),(inverse,j)
   scratch=km.array([km.q(i) for i in range(72,c.num_qubits)]);assert not sim.read_across_shots(scratch).any(),(inverse,'scratch');assert not sim.read_phase_flipped_across_shots().any(),(inverse,'phase')
  reports.append(dict(inverse=inverse,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,shots=2*shots,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
 r=dict(status='PASS_HEAD_COMPONENT',circuits=reports,seconds=time.time()-started,scope='six-moment head: single+degree2 or degree2+single, independent rational-series model, arbitrary and zero heads, native phase/scratch and both directions; full tail and decoder unbuilt');(HERE/'results/triple_gap_head_checks.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)

if __name__=='__main__':test()
