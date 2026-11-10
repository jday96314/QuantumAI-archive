"""Coprime repad, direct amplitude conversion, scan, and literal inverse."""
import json
from pathlib import Path
import kickmix as km
from schur_repad import build as padding
from derivative_inverse_shared import build as convert
from schur_amplitude_scan import write_scan
from schur_extract import body,mapped,BRE


def write_middle(f,N,M,private_base,nn,bb,index,scratch):
 n=N//2;rows=[]
 def emit(c,ports,stage):
  qm={}
  for key,word in ports.items():
   assert len(c.register_data[key])==len(word),(key,len(c.register_data[key]),len(word));qm.update({q.id:v for q,v in zip(c.register_data[key],word)})
  used=set(qm.values());available=[q for q in scratch if q not in used];missing=[q for q in range(c.num_qubits) if q not in qm];assert len(missing)<=len(available);qm.update(dict(zip(missing,available)));f.write(BRE.sub(lambda x:f'b{private_base+int(x[1])}',mapped(body(c),qm)));rows.append(dict(stage=stage,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like))
 cp={f'B{i}':w for i,w in enumerate(bb)};cp['index']=index;pp={**cp,**{f'A{i}':w for i,w in enumerate(nn)}};conv={key:v for key,v in pp.items() if key!='index'}
 emit(padding(n,'count'),cp,'count');emit(padding(n),pp,'repad');emit(convert(n),conv,'amplitude');scan=write_scan(f,N,M,private_base,(nn,bb),scratch);emit(convert(n,True),conv,'amplitude_inverse');emit(padding(n,inverse=True),pp,'repad_inverse');emit(padding(n,'count',True),cp,'count_inverse');return dict(components=rows,**scan,middle_T=sum(r['T'] for r in rows)+M*scan['kernel_T'])


def test():
 from io import StringIO
 import numpy as np
 from quadratic_schur_model import F,node_maps
 from extension_inverse import setup
 N=6;n=N//2;M=63;shots=1024;m,_,fw,bk,tr=setup();b=km.CircuitBuilder();NN=[b.create_quantum_register(12,f'A{i}') for i in range(n)];BB=[b.create_quantum_register(12,f'B{i}') for i in range(n)]
 for j in range(M):b.create_classical_register(12,f'E{j}')
 nn=[[q.id for q in w] for w in NN];bb=[[q.id for q in w] for w in BB];base=24*n;index=list(range(base,base+(n+1).bit_length()));scratch=list(range(index[-1]+1,index[-1]+301));f=StringIO();f.write(str(b.finish_circuit())+'\n');r=write_middle(f,N,M,12*M,nn,bb,index,scratch);c=km.Circuit(f.getvalue());rng=np.random.default_rng(142089);_,_,P,_=node_maps(n);aa=[];bs=[];expected=[];mask=rng.integers(0,4096,(shots,M))
 for case in range(shots):
  d=case%(n+1);positions=rng.choice(M,d,replace=False).tolist();roots=[F.ex[j] for j in positions];values=rng.integers(1,4096,d).tolist();B=[1]
  for a in roots:B=F.product(B,[a,1])
  A=[]
  for a,e in zip(roots,values):
   q=[0]*(len(B)-1);q[-1]=B[-1]
   for i in reversed(range(len(q)-1)):q[i]=B[i+1]^F.mul(a,q[i+1])
   A=F.add(A,F.scale(q,F.mul(e,F.mul(a,F.evaluate(P,a)))))
  aa.append(([0]*(n-d)+A+[0]*n)[:n]);bs.append(([0]*(n-d)+B+[0]*n)[:n]);expected.append(sum((int(mask[case,j])&v).bit_count() for j,v in zip(positions,values))%2)
 bits=lambda v:((np.asarray(v)[None,:]>>np.arange(12)[:,None])&1).astype(bool);sim=km.Simulator(batch_size=shots,seed=142089);sim.use_same_registers_as(c)
 for key,values in [('A',aa),('B',bs)]:
  for i in range(n):sim.write_across_shots(c.register_data[f'{key}{i}'],bits([tr(v[i],fw) for v in values]))
 for j in range(M):sim.write_across_shots(c.register_data[f'E{j}'],bits(mask[:,j]))
 sim.do(c);bad=np.flatnonzero(sim.read_phase_flipped_across_shots()!=np.array(expected,dtype=bool));errors={}
 for key,values in [('A',aa),('B',bs)]:
  for i in range(n):
   if not np.array_equal(sim.read_across_shots(c.register_data[f'{key}{i}']),bits([tr(v[i],fw) for v in values])):errors[f'{key}{i}']=True
 for j in range(M):assert np.array_equal(sim.read_across_shots(c.register_data[f'E{j}']),bits(mask[:,j]))
 junk=km.array([km.q(q) for q in range(24*n,c.num_qubits)]);r.update(status='PASS' if not len(bad) and not errors and not sim.read_across_shots(junk).any() else 'FAIL',shots=shots,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,phase_errors=bad[:8].tolist(),errors=errors,scratch=bool(sim.read_across_shots(junk).any()));Path('results/schur_amplitude_shared_middle_n6_native.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True);assert r['status']=='PASS'
if __name__=='__main__':test()
