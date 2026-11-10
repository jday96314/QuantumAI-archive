"""Quantum independent-value and arbitrary-input composition terminal-block checks."""
import argparse,json,subprocess,time
from pathlib import Path
import numpy as np
import kickmix as km
from import_published import load
from low_workspace_inverse import bind
from mobius_series import BoundedArithmetic
from blocked_monic_euclid import schedule
from series_euclid_block import apply
from check_deferred_padding import monomial_reference
ROOT=Path(__file__).resolve().parent


def build(N,width,which,k):
 m=load();m.inv_word=bind(m);bb=m.Builder();words=[bb.qubits(12) for _ in range(N)];flags=bb.qubits(width);Bp=words[:k-1];B=words[k-1:2*k-1];C=words[2*k-1:-1];p=words[-1];E=m.Em(bb,*bb.qubits(3),m.Alloc(bb));math=BoundedArithmetic(m,E);math.virtual_zero=True
 for inv in ([False,True] if which=='both' else [which=='inverse']):apply(m,E,math,Bp,B,C,p,width,inv,flags)
 rows=np.asarray(bb._items._rows,dtype=np.uint32).reshape(-1,9);stem=f'series_block_n{N}_k{k}_b{width}_{which}';tape=ROOT/f'parameters/{stem}.npy';np.save(tape,rows);kb=km.CircuitBuilder()
 for i in range(N):kb.create_quantum_register(12,f'S{i}')
 kb.create_quantum_register(width,'flags');pre=ROOT/f'parameters/{stem}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx');empty=ROOT/f'parameters/{stem}.empty';empty.write_text('');out=ROOT/f'candidates/{stem}.kmx'
 subprocess.run([str(ROOT/'translate_prefix'),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(bb.next_out),str(empty),str(len(rows))],check=True)
 c=km.Circuit.from_file(str(out),format='kmx');c.to_file(str(out.with_suffix('.kmb')),format='kmb');return c


def main():
 ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,default=22);ap.add_argument('--block',type=int,default=4);ap.add_argument('--shots',type=int,default=512);ap.add_argument('--k',type=int,default=7);a=ap.parse_args();rng=np.random.default_rng(65123);N=a.n;k=a.k;inputs=[];outputs=[];ff=[];skipped=0;begin=time.time();bpidx=list(range(k-1));bidx=list(range(k-1,2*k-1));cidx=list(range(2*k-1,N-1));_,(npidx,nbidx)=schedule(bpidx,bidx,list(reversed(cidx))[:2*a.block])
 while len(inputs)<a.shots:
  row=rng.integers(0,4096,N).tolist();row[-1]=int(rng.integers(1,4096))
  if len(inputs)%4==0:
   for j in cidx:row[j]=0
  try:(bp,bb,cc,p),flags=monomial_reference(k,a.block,row[:k-1]+[1],row[k-1:2*k-1]+[1],row[2*k-1:-1],row[-1])
  except ValueError:skipped+=1;continue
  out=row.copy()
  for j,v in zip(npidx,bp[:-1]):out[j]=v
  for j,v in zip(nbidx,bb[:-1]):out[j]=v
  
  for j,v in zip(cidx[:-2*a.block],cc):out[j]=v
  out[-1]=p;inputs.append(row);outputs.append(out);ff.append(flags)
 inputs=np.array(inputs).T;outputs=np.array(outputs).T;flags=np.array(ff,dtype=bool).T;bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool);reports=[]
 for which in ['forward','inverse','both']:
  c=build(N,a.block,which,k);start=outputs if which=='inverse' else inputs;want=outputs if which=='forward' else inputs
  if which=='both':
   start=rng.integers(0,4096,(N,a.shots));start[:,0]=0;start[-2,1::4]=0;start[-1]=rng.integers(1,4096,a.shots);want=start
  sim=km.Simulator(batch_size=a.shots,seed=814);sim.use_same_registers_as(c)
  for i,values in enumerate(start):sim.write_across_shots(c.register_data[f'S{i}'],bits(values))
  if which=='inverse':sim.write_across_shots(c.register_data['flags'],flags)
  sim.do(c);failures={}
  for i,values in enumerate(want):
   bad=np.flatnonzero((sim.read_across_shots(c.register_data[f'S{i}'])!=bits(values)).any(0))
   if len(bad):failures[i]=bad[:8].tolist()
  assert np.array_equal(sim.read_across_shots(c.register_data['flags']),flags if which=='forward' else np.zeros_like(flags)),which
  public={q for reg in c.register_data.values() for q in reg};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);junk=np.flatnonzero(sim.read_across_shots(scratch).any(0));phase=np.flatnonzero(sim.read_phase_flipped_across_shots());r=dict(which=which,status='PASS' if not failures and not len(junk) and not len(phase) else 'FAIL',Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,failures=failures,scratch_errors=junk[:8].tolist(),phase_errors=phase[:8].tolist());reports.append(r);print(json.dumps(r),flush=True)
 result=dict(status='PASS' if all(r['status']=='PASS' for r in reports) else 'FAIL',scope='full-series block independent monomial reference, forward/inverse, fresh arbitrary-input composition, flags/scratch/phase; not full decoding',N=N,k=k,block=a.block,shots=a.shots,skipped_nongeneric_value_cases=skipped,circuits=reports,seconds=time.time()-begin);Path(f'results/series_block_n{N}_k{k}_b{a.block}_check.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
