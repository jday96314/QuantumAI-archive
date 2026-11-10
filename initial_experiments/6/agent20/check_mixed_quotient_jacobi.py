"""Native mixed-degree reconstruction, independent rational recurrence."""
from pathlib import Path
import argparse,hashlib,json,subprocess,time
import numpy as np
import kickmix as km
import mixed_quotient_jacobi as mixed
from import_published import load
from projective_polynomial_inverse import bind
import check_blocked_cfraction as ref
HERE=Path(__file__).resolve().parent

def build(n,both=False,through=0):
 m=load();m.inv_word=bind(m);b=m.Builder();W=[b.qubits(12) for _ in range(2*n)];E=m.Em(b,*b.qubits(3),m.Alloc(b))
 c=E.A.open('mixed_terminal',12);metadata=E.A.open('mixed_metadata',n.bit_length()+1);ng=sum(k>=through for k in range(0,n-1,2));flags=E.A.open('mixed_gaps',ng)
 F,G=mixed.synthesize(m,E,W,c,metadata,flags,through=through)
 if both:
  mixed.synthesize(m,E,W,c,metadata,flags,True,through);E.A.close('mixed_gaps',flags);E.A.close('mixed_metadata',metadata);E.A.close('mixed_terminal',c)
 tag=f'mixed_quotient_jacobi_n{n}_p{through}'+('_both' if both else '');tape=HERE/'parameters'/f'{tag}.npy';rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9);np.save(tape,rows);kb=km.CircuitBuilder()
 for j in range(2*n):kb.create_quantum_register(12,f'W{j}')
 pre=HERE/'parameters'/f'{tag}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx');empty=HERE/'parameters/empty';out=HERE/'candidates'/f'{tag}.kmx'
 subprocess.run([str(ref.probe_translator()),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True,capture_output=True);circuit=km.Circuit.from_file(str(out),format='kmx');circuit.to_file(str(out.with_suffix('.kmb')),format='kmb');return circuit,F,G,metadata,flags,out.with_suffix('.kmb')

def test(n,through=0,shots=128):
 m=load();rng=np.random.default_rng(2026102617+n);started=time.time();values=np.zeros((shots,2*n),dtype=np.int64);degrees=np.arange(shots)%(n+1);polys=[]
 def product(a,b):
  if len(a)==0 or len(b)==0:return []
  out=[0]*(len(a)+len(b)-1)
  for i,v in enumerate(a):
   for j,w in enumerate(b):out[i+j]^=m.gmul(int(v),int(w))
  return out
 def trim(a):
  while len(a)>1 and not a[-1]:a.pop()
  return a or [0]
 for s,d in enumerate(degrees):
  packets=[]
  for k in range(0,n-1,2):
   if k>=d:continue
   a,b,c,e=map(int,rng.integers(1,4096,4));d1=int(rng.integers(0,4096))
   if k+1==d:
    values[s,2*k:2*k+4]=[a,b,0,b] if k>=through else [a,b,0,0];packets.append((a,[b,1]))
   elif k>=through and (s+k)%3==0:
    values[s,2*k:2*k+4]=[0,b,c,e];packets.append((b,[c,e,1]))
   else:
    values[s,2*k:2*k+4]=[a,b,c^m.gmul(b,d1),b^d1] if k>=through else [a,b,c,d1];packets.extend([(a,[b,1]),(c,[d1,1])])
  if n%2 and n-1<d:
   a,b=map(int,rng.integers(1,4096,2));values[s,-2:]=[a,b];packets.append((a,[b,1]))
  F,G=[1],[]
  for a,q in reversed(packets):
   out=product(q,F)
   for j,v in enumerate(G):out[j]^=v
   F,G=out,[m.gmul(a,v) for v in F]
  polys.append((F,G))
 bits=lambda a:((np.asarray(a)[None,:]>>np.arange(12)[:,None])&1).astype(bool);reports=[]
 for both in (False,True):
  inputs=values.copy()
  if both:
   inputs[shots//2:]=rng.integers(0,4096,(shots-shots//2,2*n));inputs[shots//2::2,::2]=0;inputs[-1,:]=0
  circuit,F,G,metadata,flags,path=build(n,both,through)
  for seed in (2026102618,2026102619):
   sim=km.Simulator(batch_size=shots,seed=seed);sim.use_same_registers_as(circuit)
   for j in range(2*n):sim.write_across_shots(circuit.register_data[f'W{j}'],bits(inputs[:,j]))
   sim.do(circuit)
   if both:
    for j in range(2*n):assert np.array_equal(sim.read_across_shots(circuit.register_data[f'W{j}']),bits(inputs[:,j])),(n,'inverse',j)
    live=set(range(24*n))
   else:
    banks=[]
    for bank in (F,G):
     out=[]
     for word in bank:
      read=sim.read_across_shots(km.array([km.q(i) for i in word]));out.append(sum(read[j].astype(np.int64)<<j for j in range(12)))
     banks.append(np.asarray(out).T)
    for s,((f,g),(ff,gg),d) in enumerate(zip(zip(*banks),polys,degrees)):
     assert trim(product(f,gg))==trim(product(g,ff)),(n,'rational',s,int(d),f.tolist(),g.tolist(),ff,gg)
     assert not any(f[:n-d]) and not any(g[:n-d]),(n,'padding',s,int(d))
    count=sim.read_across_shots(km.array([km.q(i) for i in metadata[:-1]]));actual=sum(count[j].astype(np.int64)<<j for j in range(len(count)));assert np.array_equal(actual,n-degrees),(n,'completion count',actual.tolist(),(n-degrees).tolist())
    live={i for bank in (F,G) for word in bank for i in word}|set(metadata)|set(flags)
   scratch=km.array([km.q(i) for i in range(circuit.num_qubits) if i not in live]);assert not sim.read_across_shots(scratch).any(),(n,both,'scratch');assert not sim.read_phase_flipped_across_shots().any(),(n,both,'phase')
  reports.append(dict(both=both,Q=circuit.num_qubits,T=circuit.operation_statistics().num_toffoli_like,shots=2*shots,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
 r=dict(status='PASS_MIXED_COMPONENT',n=n,through=through,circuits=reports,seconds=time.time()-started,scope='shared degree1/2 reconstruction, independent variable-quotient rational identity, all completion degrees and common padding; arbitrary paired zero states, phase and scratch; extraction not included');print(json.dumps(r),flush=True);return r

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--degrees',type=int,nargs='+',default=[3,7,17]);p.add_argument('--through',type=int,default=0);p.add_argument('--out');a=p.parse_args();rows=[]
 for n in a.degrees:
  rows.append(test(n,a.through));(HERE/'results'/(a.out or f'mixed_quotient_jacobi_p{a.through}_n'+('_'.join(map(str,a.degrees)))+'_checks.json')).write_text(json.dumps(rows,indent=2)+'\n')
