"""Independent native-field forward/inverse tests for a compact gap block."""
import argparse,json,subprocess
from pathlib import Path
import numpy as np,kickmix as km
from import_published import load
from low_workspace_inverse import bind
from mobius_series import BoundedArithmetic
from gap_block import apply
from classical_repair import mul,inv,product,plus,repair
from published_dqi import isomorphism

def sp(a,b,n):
 r=[0]*n
 for i,x in enumerate(a[:n]):
  for j,y in enumerate(b[:n-i]):r[i+j]^=mul(x,y)
 return r
def sd(a,b,n):
 a=list(a)+[0]*max(0,n-len(a));r=[0]*n
 for i in range(n):
  z=a[i]
  for j in range(1,min(i+1,len(b))):z^=mul(b[j],r[i-j])
  r[i]=mul(z,inv(b[0]))
 return r
def reference(bp,b,c,pi):
 n=6;t=list(reversed(c))[:6]
 F=sd(sp(t,b[::-1],6),[inv(pi),0]+sp(t,bp[::-1],4),6)
 if F[0]:raise ValueError('first pivot nonzero')
 U,V,pinv=repair(F);e,d,a,bb,_,_=F
 A2=[d^mul(e,U[1]),e];A3=[a^mul(d,V[2])^mul(e,V[1]),d^mul(e,V[2]),e]
 newp=plus(product(A2,bp),product(U,b));newb=plus(product(A3,bp),product(V,b))
 newc=sd(sp(newb[::-1],c[::-1],len(c))[6:],b[::-1],len(c)-6)[::-1]
 return newp,newb,newc,mul(pi,pinv)

def build(N,k,which):
 m=load();m.inv_word=bind(m);bb=m.Builder();W=[bb.qubits(12) for _ in range(N)];E=m.Em(bb,*bb.qubits(3),m.Alloc(bb));math=BoundedArithmetic(m,E);math.virtual_zero=True
 bp=W[:k-1];b=W[k-1:2*k-1];c=W[2*k-1:-1];pi=W[-1]
 for back in ([False,True] if which=='both' else [True,False] if which=='inverse_both' else [which=='inverse']):apply(m,E,math,bp,b,c,pi,back)
 rows=np.asarray(bb._items._rows,dtype=np.uint32).reshape(-1,9);stem=f'gap_block_n{N}_k{k}_{which}';tape=Path('parameters')/(stem+'.npy');np.save(tape,rows)
 kb=km.CircuitBuilder()
 for i in range(N):kb.create_quantum_register(12,f'W{i}')
 pre=Path('parameters')/(stem+'.pre.kmx');kb.finish_circuit().to_file(str(pre),format='kmx');empty=Path('parameters/empty');empty.write_text('');out=Path('candidates')/(stem+'.kmx')
 subprocess.run(['./translate_prefix',str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(bb.next_out),str(empty),str(len(rows))],check=True)
 c=km.Circuit.from_file(str(out),format='kmx');c.to_file(str(out.with_suffix('.kmb')),format='kmb');return c

def test(N,k,shots):
 rng=np.random.default_rng(816919+k);m=load();fw,bw,tr=isomorphism(m);starts=[];wants=[]
 ci=list(range(2*k-1,N-1));raw=ci[::-1][:6];Q=[raw[i] for i in [1,3,2,4,5,0]];newpi=list(range(k-1))+Q[:3];newbi=list(range(k-1,2*k-1))+[Q[3],Q[5],Q[4]]
 while len(starts)<shots:
  row=rng.integers(0,4096,N).tolist();row[-1]=int(rng.integers(1,4096));row[-2]=0;row[-3]=int(rng.integers(1,4096));bp=row[:k-1]+[1];b=row[k-1:2*k-1]+[1];c=row[2*k-1:-1];pi=row[-1]
  try:np_,nb,nc,npi=reference(bp,b,c,pi)
  except ValueError:continue
  out=row.copy()
  for i,v in zip(newpi,np_[:-1]):out[i]=v
  for i,v in zip(newbi,nb[:-1]):out[i]=v
  for i,v in zip(ci[:-6],nc):out[i]=v
  out[-1]=npi;starts.append(row);wants.append(out)
 convert=lambda a:np.array([[tr(int(v),fw) for v in row] for row in np.array(a).T]);starts=convert(starts);wants=convert(wants);bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool);results=[]
 for which in ['forward','inverse','both','inverse_both']:
  c=build(N,k,which);sim=km.Simulator(batch_size=shots,seed=912917);sim.use_same_registers_as(c);src=wants if which in ['inverse','inverse_both'] else starts;dst=wants if which in ['forward','inverse_both'] else starts
  for i,v in enumerate(src):sim.write_across_shots(c.register_data[f'W{i}'],bits(v))
  sim.do(c);errors={}
  for i,v in enumerate(dst):
   bad=np.flatnonzero((sim.read_across_shots(c.register_data[f'W{i}'])!=bits(v)).any(0))
   if len(bad):errors[i]=bad[:8].tolist()
  public={q for reg in c.register_data.values() for q in reg};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);junk=int(sim.read_across_shots(scratch).any(0).sum());phase=int(sim.read_phase_flipped_across_shots().sum());results.append(dict(which=which,status='PASS' if not errors and not junk and not phase else 'FAIL',errors=errors,junk=junk,phase=phase,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like))
 result=dict(N=N,k=k,shots=shots,status='PASS' if all(r['status']=='PASS' for r in results) else 'FAIL',scope='complete local compact degree-gap block: high remainder reconstruction and global cofactors, inverse and both compositions; no global dispatch/decoder score',circuits=results)
 Path(f'results/gap_block_n{N}_k{k}.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=22);p.add_argument('--k',type=int,default=3);p.add_argument('--shots',type=int,default=512);a=p.parse_args();test(a.n,a.k,a.shots)
