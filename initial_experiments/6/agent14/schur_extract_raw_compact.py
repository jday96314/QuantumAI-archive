"""Raw fixed-node Schur extraction without scalar normalization or swaps."""
import argparse,json,time
from pathlib import Path
import numpy as np
import kickmix as km
from extension_inverse import setup
from schur_extract import body,mapped,BRE,linear_map,imul_text,axpy_text
from spooky_extension_compact_phase import build as inverse
from quadratic_schur_model import node_maps
from raw_cf_model import extract


def build(N,reverse=False):
 start=time.time();n=N//2;m,k,fw,bk,tr=setup();nodes,kap,_,_=node_maps(n);cb=km.CircuitBuilder()
 for i in range(N):cb.create_quantum_register(12,f'S{i}')
 total=24*n*(n-1)//2;cb.create_classical_register(total,'ghost');heads=str(cb.finish_circuit())+'\n';data=[list(range(24*i,24*i+24)) for i in range(n)];scratch=list(range(12*N,12*N+80));c=inverse(reverse);txt=body(c);ops=linear_map(N);path=Path(f'candidates/schur_extract_raw_compact_n{N}_{"inverse" if reverse else "forward"}.kmx');pair={(i,j):sum(n-1-z for z in range(i))+j-i-1 for i in range(n) for j in range(i+1,n)}
 with path.open('w') as f:
  f.write(heads)
  if not reverse:
   for op,a,b in ops:f.write(f'{op.upper()} q{a} q{b}\n')
  for i in (reversed(range(n)) if reverse else range(n)):
   cur=data[i]
   if reverse:f.write(axpy_text(cur[12:],cur[:12],tr(nodes[i]&4095,fw),m))
   for j in (reversed(range(i+1,n)) if reverse else range(i+1,n)):
    target=data[j];constant=tr(kap[i]^kap[j],fw)
    def translate():
     for q,r in zip(cur,target):f.write(f'CX q{q} q{r}\n')
     f.write(axpy_text(cur[12:],target[:12],tr((nodes[i]^nodes[j])&4095,fw),m))
    def scale(back=False):
     for w in [target[:12],target[12:]]:f.write(imul_text(w,constant,m,back))
    if reverse:scale(True)
    else:translate()
    qm=dict(enumerate(target));qm.update({24+j:q for j,q in enumerate(scratch[:c.num_qubits-24])});base=24*pair[i,j];f.write(BRE.sub(lambda x:f'b{base+int(x[1]) if int(x[1])<24 else total+int(x[1])-24}',mapped(txt,qm)))
    if reverse:translate()
    else:scale()
   if not reverse:f.write(axpy_text(cur[12:],cur[:12],tr(nodes[i]&4095,fw),m))
  if reverse:
   for op,a,b in reversed(ops):f.write(f'{op.upper()} q{a} q{b}\n')
 c=km.Circuit.from_file(str(path),format='kmx');c.to_file(str(path.with_suffix('.kmb')),format='kmb');r=dict(N=N,reverse=reverse,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,seconds=time.time()-start);Path(f'results/schur_extract_raw_compact_n{N}_{"inverse" if reverse else "forward"}_build.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True);return c


def test(N=6,shots=512):
 from quadratic_schur_model import F
 cf=build(N);ci=build(N,True);m,k,fw,bk,tr=setup();rng=np.random.default_rng(142112+N);n=N//2;ss=rng.integers(0,4096,(shots,N));ss[0]=0
 for case in range(0,shots,2):
  d=(case//2)%(n+1);pos=rng.choice(np.arange(1,4096),d,replace=False).tolist();val=rng.integers(1,4096,d).tolist();ss[case]=F.syndrome(pos,val,N)
 models=[extract(s.tolist()) for s in ss];want=np.array([[tr(x,fw) for h in row[0] for x in h] for row in models]);norms=[]
 for row in models:
  rr=[]
  for den in row[1]:
   a=tr(den&4095,fw);b=tr(den>>12,fw);rr.append(m.gmul(a,a)^m.gmul(a,b)^m.gmul(k,m.gmul(b,b)))
  norms.append(rr)
 norms=np.array(norms).T;bits=lambda v:((np.asarray(v)[None,:]>>np.arange(12)[:,None])&1).astype(bool);rows=[]
 for reverse,c,src,dst in [(False,cf,ss,want),(True,ci,want,ss)]:
  sim=km.Simulator(batch_size=shots,seed=142112+reverse);sim.use_same_registers_as(c)
  for i in range(N):sim.write_across_shots(c.register_data[f'S{i}'],bits(src[:,i]))
  sim.do(c);bad=[]
  for i in range(N):
   if not np.array_equal(sim.read_across_shots(c.register_data[f'S{i}']),bits(dst[:,i])):bad.append(i)
  scratch=km.array([km.q(q) for q in range(12*N,c.num_qubits)]);mask=sim.read_across_shots(c.register_data['ghost']).reshape(-1,24,shots);expect=np.zeros(shots,dtype=bool) if reverse else ((mask[:,:12,:]&((norms[:,None,:]>>np.arange(12)[None,:,None])&1)).sum((0,1))%2).astype(bool);r=dict(reverse=reverse,errors=bad,scratch=bool(sim.read_across_shots(scratch).any()),phase=bool((sim.read_phase_flipped_across_shots()!=expect).any()));rows.append(r);print(json.dumps(r),flush=True);assert not bad and not r['scratch'] and not r['phase']
  if not reverse:
   sim.do(ci)
   for i in range(N):assert np.array_equal(sim.read_across_shots(c.register_data[f'S{i}']),bits(src[:,i]))
   assert not sim.read_across_shots(scratch).any();assert not sim.read_across_shots(c.register_data['ghost']).any();assert not sim.read_phase_flipped_across_shots().any()
 r=dict(status='PASS',N=N,shots=shots,rows=rows);Path(f'results/schur_extract_raw_compact_n{N}_native.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=6);p.add_argument('--shots',type=int,default=512);a=p.parse_args();test(a.n,a.shots)
