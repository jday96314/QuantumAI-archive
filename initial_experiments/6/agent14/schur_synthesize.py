"""Compose 35 tested Schur folds using two recycled coefficient words."""
import argparse,json,time
from pathlib import Path
import numpy as np
import kickmix as km
from extension_inverse import setup
from schur_fold import build as fold,reference
from schur_basefold import build as basefold
from schur_extract import body,mapped
from quadratic_schur_model import node_maps,F,transform,extract


def plan(N):
 n=N//2;words=[list(range(12*i,12*(i+1))) for i in range(N)];free=[list(range(12*N+12,12*N+24)),list(range(12*N+24,12*N+36))];nn=[];bb=[];steps=[]
 for i in reversed(range(n)):
  a,b=words[2*i:2*i+2]
  if i==n-1:
   steps.append(dict(i=i,d=0,ports=a+b,base=True));nn=[b];bb=[a];continue
  nt,bt=free
  ports=a+b+[q for w in nn for q in w]+[q for w in bb for q in w]+nt+bt
  steps.append(dict(i=i,d=len(nn),ports=ports));nn=nn+[nt];bb=bb+[bt];free=[a,b]
 return steps,nn,bb,free


def build(N,reverse=False):
 start=time.time();m,k,fw,bk,tr=setup();_,kap,_,_=node_maps(N//2);steps,nn,bb,free=plan(N);cb=km.CircuitBuilder()
 for i in range(N):cb.create_quantum_register(12,f'S{i}')
 cb.create_quantum_register(12,'g');headers='\n'.join(line for line in str(cb.finish_circuit()).splitlines() if line.startswith(('REGISTER','APPEND_TO_REGISTER')))+'\n';scratch=list(range(12*N+36,12*N+36+97))
 path=Path(f'candidates/schur_synthesize_n{N}_{"inverse" if reverse else "forward"}.kmx')
 with path.open('w') as f:
  f.write(headers)
  for step in (reversed(steps) if reverse else steps):
   d=step['d'];c=(basefold(tr(kap[step['i']],fw),reverse) if step.get('base') else fold(d,tr(kap[step['i']],fw),reverse));ports=step['ports'];count=len(ports);qmap=dict(enumerate(ports));qmap.update({count+j:scratch[j] for j in range(c.num_qubits-count)});f.write(mapped(body(c),qmap))
 c=km.Circuit.from_file(str(path),format='kmx');c.to_file(str(path.with_suffix('.kmb')),format='kmb')
 r=dict(N=N,reverse=reverse,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,seconds=time.time()-start,numerator_qubits=nn,denominator_qubits=bb,clean_recycled_words=free,scope='synthesis only, assumes generic nonzero pivots or completed zero tails')
 Path(f'results/schur_synthesize_n{N}_{"inverse" if reverse else "forward"}_build.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if not isinstance(v,list)}),flush=True);return c


def test(N,shots=256):
 m,k,fw,bk,tr=setup();rng=np.random.default_rng(142045+N);n=N//2;nodes,kap,P,maps=node_maps(n);rows=[];pivots=[];wantn=[];wantb=[];truecases=0
 for case in range(shots):
  if case%2==0:
   d=(case//2)%(n+1);pos=rng.choice(np.arange(1,4096),d,replace=False).tolist();val=rng.integers(1,4096,d).tolist();s=F.syndrome(pos,val,N);stages,exc,cc=extract(nodes,kap,transform(s,maps),False)
   if exc:continue
   v=[tr(x,fw) for ss in stages for x in [ss['a'],ss['b']]]+[0]*(N-2*len(stages));truecases+=1
  else:
   v=rng.integers(0,4096,N).tolist();v[1::2]=rng.integers(1,4096,n).tolist()
  nn=[];bb=[]
  for i in reversed(range(n)):nn,bb=reference(v[2*i],v[2*i+1],nn,bb,tr(kap[i],fw))
  pivots.append(v);wantn.append(nn);wantb.append(bb)
 shots=len(pivots);src=np.array(pivots);wn=np.array(wantn);wb=np.array(wantb);cf=build(N);ci=build(N,True);steps,nn,bb,free=plan(N);bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool);g=rng.integers(1,4096,shots)
 for inv,c in [(False,cf),(True,ci)]:
  sim=km.Simulator(batch_size=shots,seed=142045+inv);sim.use_same_registers_as(c);sim.write_across_shots(c.register_data['g'],bits(g))
  if not inv:
   for i in range(N):sim.write_across_shots(c.register_data[f'S{i}'],bits(src[:,i]))
  else:
   for ws,vals in [(nn,wn),(bb,wb)]:
    for j,w in enumerate(ws):sim.write_across_shots(km.array([km.q(q) for q in w]),bits(vals[:,j]))
  sim.do(c);errors={}
  if inv:
   for i in range(N):
    bad=np.flatnonzero((sim.read_across_shots(c.register_data[f'S{i}'])!=bits(src[:,i])).any(0))
    if len(bad):errors[f'S{i}']=bad[:8].tolist()
   used=set(range(12*N+12))
  else:
   for name,ws,vals in [('N',nn,wn),('B',bb,wb)]:
    for j,w in enumerate(ws):
     bad=np.flatnonzero((sim.read_across_shots(km.array([km.q(q) for q in w]))!=bits(vals[:,j])).any(0))
     if len(bad):errors[f'{name}{j}']=bad[:8].tolist()
   used={q for w in nn+bb for q in w}|set(range(12*N,12*N+12))
  bad=np.flatnonzero((sim.read_across_shots(c.register_data['g'])!=bits(g)).any(0))
  if len(bad):errors['g']=bad[:8].tolist()
  sc=km.array([km.q(q) for q in range(c.num_qubits) if q not in used]);junk=np.flatnonzero(sim.read_across_shots(sc).any(0));phase=np.flatnonzero(sim.read_phase_flipped_across_shots());r=dict(inverse=inv,errors=errors,junk=junk[:8].tolist(),phase=phase[:8].tolist());rows.append(r);print(json.dumps(r),flush=True);assert not errors and not len(junk) and not len(phase)
  sim.do(cf if inv else ci)
  if not inv:
   for i in range(N):assert np.array_equal(sim.read_across_shots(c.register_data[f'S{i}']),bits(src[:,i])),i
  else:
   for ws,vals in [(nn,wn),(bb,wb)]:
    for j,w in enumerate(ws):assert np.array_equal(sim.read_across_shots(km.array([km.q(q) for q in w])),bits(vals[:,j]))
  assert not sim.read_phase_flipped_across_shots().any()
 Path(f'results/schur_synthesize_n{N}_native.json').write_text(json.dumps(dict(status='PASS',shots=shots,true_syndrome_cases=truecases,rows=rows,scope='all weights, generic pivots, full coefficient model, both directions/compositions, frame preserved, clean scratch and phase'),indent=2)+'\n')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=6);p.add_argument('--shots',type=int,default=256);a=p.parse_args();test(a.n,a.shots)
