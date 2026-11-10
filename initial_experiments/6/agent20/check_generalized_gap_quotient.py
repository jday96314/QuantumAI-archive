"""Independent variable-degree continuants and arbitrary native inverse."""
from pathlib import Path
import argparse,json,hashlib,subprocess,time
import numpy as np
import kickmix as km
from import_published import load
from projective_polynomial_inverse import bind
import generalized_gap_quotient as general
import check_blocked_cfraction as ref
HERE=Path(__file__).resolve().parent
PREFIX='general_gap_quotient2_n'
SPARSE=False

def build(n,both=False):
 m=load();m.inv_word=bind(m,limit=1004);b=m.Builder();W=[b.qubits(12) for _ in range(2*n)];c=b.qubits(12);flags=b.qubits(n);F,G=general.base.layout(W,c);E=m.Em(b,*b.qubits(3),m.Alloc(b));gapwords=[None if SPARSE and i%3==2 else q for i,q in enumerate(flags)];general.synthesize(m,E,F,G,gapwords)
 if both:general.synthesize(m,E,F,G,gapwords,True)
 tag=f'{PREFIX}{n}'+('_both' if both else '');rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9);tape=HERE/'parameters'/f'{tag}.npy';np.save(tape,rows);kb=km.CircuitBuilder()
 for j in range(2*n):kb.create_quantum_register(12,f'W{j}')
 kb.create_quantum_register(12,'c');kb.create_quantum_register(n,'gap');pre=HERE/'parameters'/f'{tag}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx');empty=HERE/'parameters/empty';out=HERE/'candidates'/f'{tag}.kmx'
 subprocess.run([str(ref.probe_translator()),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True,capture_output=True);circ=km.Circuit.from_file(str(out),format='kmx');circ.to_file(str(out.with_suffix('.kmb')),format='kmb');return circ,F,G,out.with_suffix('.kmb')

def test(n,shots=128):
 started=time.time();m=load();rng=np.random.default_rng(2026102632+n);W=rng.integers(0,4096,(shots,2*n));constants=rng.integers(1,4096,shots);flags=np.zeros((shots,n),dtype=bool);expected=[]
 def prod(a,b):
  out=[0]*(len(a)+len(b)-1)
  for i,x in enumerate(a):
   for j,y in enumerate(b):out[i+j]^=m.gmul(int(x),int(y))
  return out
 for s in range(shots):
  packets=[];i=0
  while i<n:
   if i<n-1 and (s+i)%3==0 and not(SPARSE and i%3==2):
    q2=int(rng.integers(1,4096));q1=int(rng.integers(0,4096));q0=int(rng.integers(0,4096));W[s,2*i:2*i+4]=[q2,q1,0,q0];flags[s,i]=True;packets.append([q0,q1,q2]);i+=2
   else:
    q1=int(rng.integers(1,4096));W[s,2*i]=q1;packets.append([int(W[s,2*i+1]),q1]);i+=1
  F,G=[int(constants[s])],[]
  for q in reversed(packets):
   out=prod(q,F)
   for j,v in enumerate(G):out[j]^=v
   F,G=out,F
  expected.append((F+[0]*(n+1-len(F)),G+[0]*(n-len(G))))
 bits=lambda x:((np.asarray(x)[None,:]>>np.arange(12)[:,None])&1).astype(bool);reports=[]
 for both in (False,True):
  values=W.copy();cs=constants.copy();ctrl=flags.copy()
  if both:
   values[shots//2:]=rng.integers(0,4096,(shots-shots//2,2*n));values[shots//2::2,::2]=0;cs[shots//2::2]=0;ctrl[shots//2:]=rng.integers(0,2,(shots-shots//2,n)).astype(bool)
  circ,F,G,path=build(n,both)
  for seed in (2026102633,2026102634):
   sim=km.Simulator(batch_size=shots,seed=seed);sim.use_same_registers_as(circ)
   for j in range(2*n):sim.write_across_shots(circ.register_data[f'W{j}'],bits(values[:,j]))
   sim.write_across_shots(circ.register_data['c'],bits(cs));sim.write_across_shots(circ.register_data['gap'],ctrl.T);sim.do(circ)
   if both:
    for j in range(2*n):assert np.array_equal(sim.read_across_shots(circ.register_data[f'W{j}']),bits(values[:,j])),(n,'paired',j)
    assert np.array_equal(sim.read_across_shots(circ.register_data['c']),bits(cs)),(n,'paired c')
   else:
    for k,bank in enumerate((F,G)):
     for j,word in enumerate(bank):assert np.array_equal(sim.read_across_shots(km.array([km.q(i) for i in word])),bits([p[k][j] for p in expected])),(n,'continuant',k,j)
   assert np.array_equal(sim.read_across_shots(circ.register_data['gap']),ctrl.T)
   scratch=km.array([km.q(i) for i in range(24*n+12+n,circ.num_qubits)]);assert not sim.read_across_shots(scratch).any(),(n,both,'scratch');assert not sim.read_phase_flipped_across_shots().any(),(n,both,'phase')
  reports.append(dict(both=both,Q=circ.num_qubits,T=circ.operation_statistics().num_toffoli_like,shots=2*shots,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
 r=dict(status='PASS_GENERALIZED_CONTINUANT_COMPONENT',n=n,circuits=reports,seconds=time.time()-started,scope='arbitrary-position degree1/2 quotient continuants with shared products; independent ordinary recurrence; arbitrary paired adjacent controls/zero words; phase/scratch; not Jacobi conversion/extraction or full decoder');print(json.dumps(r),flush=True);return r

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--degrees',type=int,nargs='+',default=[3,7,17,35]);p.add_argument('--prefix',default=PREFIX);p.add_argument('--sparse-triples',action='store_true');a=p.parse_args();PREFIX=a.prefix;SPARSE=a.sparse_triples;rows=[]
 for n in a.degrees:
  rows.append(test(n));(HERE/'results'/(PREFIX+('_'.join(map(str,a.degrees)))+'_checks.json')).write_text(json.dumps(rows,indent=2)+'\n')
