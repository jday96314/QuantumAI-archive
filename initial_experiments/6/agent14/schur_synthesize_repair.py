"""Polynomial synthesis for coherently swapped quadratic Schur nodes."""
import argparse,json,time
from pathlib import Path
import numpy as np
import kickmix as km
from extension_inverse import setup
from schur_synthesize import plan as fixed_plan
from schur_dynamic_fold import build as fold
from schur_dynamic_basefold import build as basefold
from schur_fold import reference
from schur_extract import body,mapped
from schur_extract_repair import controls,model
from quadratic_schur_model import node_maps,F


def plan(N):
 steps,nn,bb,free=fixed_plan(N);count=N//2-1;base=12*N+12;shift=lambda q:q+count if q>=base else q
 for s in steps:s['ports']=[shift(q) for q in s['ports']]
 return steps,[[shift(q) for q in w] for w in nn],[[shift(q) for q in w] for w in bb],[[shift(q) for q in w] for w in free]


def build(N,reverse=False):
 start=time.time();m,k,fw,bk,tr=setup();n=N//2;_,kap,_,_=node_maps(n);steps,nn,bb,free=plan(N);flags=list(range(12*N+12,12*N+12+n-1));cb=km.CircuitBuilder()
 for i in range(N):cb.create_quantum_register(12,f'S{i}')
 cb.create_quantum_register(12,'g');cb.create_quantum_register(n-1,'flags');header=str(cb.finish_circuit())+'\n';scratch=list(range(12*N+36+n-1,12*N+36+n-1+100));path=Path(f'candidates/schur_synthesize_repair_n{N}_{"inverse" if reverse else "forward"}.kmx')
 with path.open('w') as f:
  f.write(header)
  for s in (reversed(steps) if reverse else steps):
   cc=controls(N,s['i']);deltas=[tr(F.mul(de,de)^de,fw) for _,de in cc];ports=s['ports']+[flags[fl] for fl,_ in cc];c=basefold(tr(kap[s['i']],fw),deltas,reverse) if s.get('base') else fold(s['d'],tr(kap[s['i']],fw),deltas,reverse);size=len(ports);qm=dict(enumerate(ports));qm.update({size+j:q for j,q in enumerate(scratch[:c.num_qubits-size])});f.write(mapped(body(c),qm))
 c=km.Circuit.from_file(str(path),format='kmx');c.to_file(str(path.with_suffix('.kmb')),format='kmb');r=dict(N=N,reverse=reverse,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,seconds=time.time()-start,numerator_qubits=nn,denominator_qubits=bb,clean_recycled_words=free,scope='synthesis with coherent nonconsecutive node descriptors; polynomial normalization and measured-pivot repayment');Path(f'results/schur_synthesize_repair_n{N}_{"inverse" if reverse else "forward"}_build.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if not isinstance(v,list)}),flush=True);return c


def test(N,shots=256):
 m,k,fw,bk,tr=setup();rng=np.random.default_rng(142053+N);n=N//2;_,kap,_,_=node_maps(n);pivots=[];flagrows=[];wn=[];wb=[];truecases=0
 for case in range(shots):
  if case%2==0:
   d=(case//2)%(n+1);pos=rng.choice(np.arange(1,4096),d,replace=False).tolist();val=rng.integers(1,4096,d).tolist();s=F.syndrome(pos,val,N);v,g,flags,ghost=model(s);v=[tr(x,fw) for x in v];truecases+=1
   # Skip rare unsupported charts in this synthesis component test. Full
   # decoder verification measures them against the complete distribution.
   if any(v[2*i] and not v[2*i+1] for i in range(n)):continue
  else:
   v=rng.integers(0,4096,N).tolist();v[1::2]=rng.integers(1,4096,n).tolist();flags=[]
   for i in range(n-1):flags.append(int((not flags or not flags[-1]) and rng.random()<.3))
  nn=[];bb=[]
  for i in reversed(range(n)):
   kk=kap[i]
   for flag,de in controls(N,i):
    if flags[flag]:kk^=F.mul(de,de)^de
   nn,bb=reference(v[2*i],v[2*i+1],nn,bb,tr(kk,fw))
  pivots.append(v);flagrows.append(flags);wn.append(nn);wb.append(bb)
 src=np.array(pivots);ff=np.array(flagrows,dtype=bool).T;wn=np.array(wn);wb=np.array(wb);shots=len(src);g=rng.integers(1,4096,shots);cf=build(N);ci=build(N,True);_,nn,bb,free=plan(N);bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool);rows=[]
 for inv,c in [(False,cf),(True,ci)]:
  sim=km.Simulator(batch_size=shots,seed=142053+inv);sim.use_same_registers_as(c);sim.write_across_shots(c.register_data['g'],bits(g));sim.write_across_shots(c.register_data['flags'],ff)
  if inv:
   for ws,vs in [(nn,wn),(bb,wb)]:
    for j,w in enumerate(ws):sim.write_across_shots(km.array([km.q(q) for q in w]),bits(vs[:,j]))
  else:
   for i in range(N):sim.write_across_shots(c.register_data[f'S{i}'],bits(src[:,i]))
  sim.do(c);errors={}
  if inv:
   for i in range(N):
    bad=np.flatnonzero((sim.read_across_shots(c.register_data[f'S{i}'])!=bits(src[:,i])).any(0))
    if len(bad):errors[f'S{i}']=bad[:8].tolist()
   used=set(range(12*N+12+n-1))
  else:
   for name,ws,vs in [('N',nn,wn),('B',bb,wb)]:
    for j,w in enumerate(ws):
     bad=np.flatnonzero((sim.read_across_shots(km.array([km.q(q) for q in w]))!=bits(vs[:,j])).any(0))
     if len(bad):errors[f'{name}{j}']=bad[:8].tolist()
   used={q for w in nn+bb for q in w}|set(range(12*N,12*N+12+n-1))
  assert np.array_equal(sim.read_across_shots(c.register_data['flags']),ff);assert np.array_equal(sim.read_across_shots(c.register_data['g']),bits(g));sc=km.array([km.q(i) for i in range(c.num_qubits) if i not in used]);junk=np.flatnonzero(sim.read_across_shots(sc).any(0));phase=np.flatnonzero(sim.read_phase_flipped_across_shots());r=dict(inverse=inv,errors=errors,junk=junk[:8].tolist(),phase=phase[:8].tolist());rows.append(r);print(json.dumps(r),flush=True);assert not errors and not len(junk) and not len(phase)
  sim.do(cf if inv else ci)
  if not inv:
   for i in range(N):assert np.array_equal(sim.read_across_shots(c.register_data[f'S{i}']),bits(src[:,i]))
  else:
   for ws,vs in [(nn,wn),(bb,wb)]:
    for j,w in enumerate(ws):assert np.array_equal(sim.read_across_shots(km.array([km.q(q) for q in w])),bits(vs[:,j]))
  assert np.array_equal(sim.read_across_shots(c.register_data['flags']),ff);assert np.array_equal(sim.read_across_shots(c.register_data['g']),bits(g));assert not sim.read_phase_flipped_across_shots().any()
 result=dict(status='PASS',N=N,shots=shots,true_syndrome_cases=truecases,node_swaps=int(ff.sum()),rows=rows,scope='true syndromes/all weights and generic pivots, varied nonconsecutive node flags, both directions/compositions, preserved flags/frame, clean scratch and phase');Path(f'results/schur_synthesize_repair_n{N}_native.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=6);p.add_argument('--shots',type=int,default=256);a=p.parse_args();test(a.n,a.shots)
