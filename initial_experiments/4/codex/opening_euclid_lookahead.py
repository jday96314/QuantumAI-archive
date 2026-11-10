"""Initial three-step Euclid transform using the pivot-robust Padé map.

Copies only six head words while updating the remainder tail. The head map
has a total singular extension, but its Padé meaning requires D2,p nonzero.
This component does not arm the returned reciprocal normalizer.
"""
import argparse,json,subprocess,time
from pathlib import Path
import numpy as np
import kickmix as km
from import_published import load
from low_workspace_inverse import bind
from lookahead_six_quantum import apply as six,reference
from mobius_series import BoundedArithmetic
from monic_euclid_local import first_step
from blocked_monic_euclid import apply as block
from monomial_serial_euclid import apply as serial
from block_euclid_ratio import series_product

ROOT=Path(__file__).resolve().parent


def layout(S):
 e,d,a,b,c,f=list(reversed(S))[:6]
 return [a,b],[d,e,f],S[:-6],c


def apply(m,E,math,S,inverse=False):
 head=list(reversed(S))[:6]
 with m.opened(E,'opening_six_head',72) as scratch:
  Q=[scratch[12*i:12*(i+1)] for i in range(6)]
  def copy(reverse=False):
   for source,target in (zip(Q,head) if reverse else zip(head,Q)):
    for q,t in zip(source,target):E.cx(q,t)
  rev=[Q[1],Q[5],Q[0],Q[1]]
  def swap():
   for a,b in zip(Q,head):E.swaps(None,list(zip(a,b)))
  if not inverse:
   copy();six(m,E,*Q)
   math.triangular_suffix(rev,1,list(reversed(S)),6)
   six(m,E,*head);copy(True);swap()
  else:
   swap();copy(True);six(m,E,*head,True)
   math.triangular_suffix(rev,1,list(reversed(S)),6,True)
   six(m,E,*Q,True);copy()
 return layout(S)


def build(N,which,baseline=False,prefix4=False):
 m=load();m.inv_word=bind(m);bb=m.Builder();words=[bb.qubits(12) for _ in range(N)];flags=bb.qubits(4) if prefix4 or baseline else []
 E=m.Em(bb,*bb.qubits(3),m.Alloc(bb));math=BoundedArithmetic(m,E);math.virtual_zero=True
 for inv in ([False,True] if which=='both' else [True,False] if which=='inverse_both' else [which=='inverse']):
  if baseline:
   if not inv:
    with m.eq_test(E,words[-1],0) as zero:E.cx(zero,flags[0])
    E.cx(flags[0],words[-1][0]);first_step(m,E,math,words)
    block(m,E,math,[],[words[-2]],words[:-2],words[-1],3,flags=flags[1:],monomial=True,suffix=True,linear_padding=True)
   else:
    block(m,E,math,[],[words[-2]],words[:-2],words[-1],3,True,flags=flags[1:],monomial=True,suffix=True,linear_padding=True)
    first_step(m,E,math,words,True);E.cx(flags[0],words[-1][0])
    with m.eq_test(E,words[-1],0) as zero:E.cx(zero,flags[0])
  else:
   if not inv:state=apply(m,E,math,words)
   else:state=layout(words)
   if prefix4:
    bp,b,c,p=state
    # A persistent guard prevents division by a zero normalizer on the
    # deliberately extended singular-head branch. It is separate from
    # completion flags and is not interpreted as a monomial padding flag.
    if not inv:
     with m.eq_test(E,p,0) as zero:E.cx(zero,flags[0])
     E.cx(flags[0],p[0])
    serial(m,E,math,bp,b,c,p,flags[3],inv,linear_padding=True)
    if inv:
     E.cx(flags[0],p[0])
     with m.eq_test(E,p,0) as zero:E.cx(zero,flags[0])
   if inv:apply(m,E,math,words,True)
 rows=np.asarray(bb._items._rows,dtype=np.uint32).reshape(-1,9);stem=f'opening_lookahead_n{N}_{which}'+('_baseline' if baseline else '')+('_prefix4' if prefix4 else '');tape=ROOT/f'parameters/{stem}.npy';np.save(tape,rows)
 kb=km.CircuitBuilder()
 for i in range(N):kb.create_quantum_register(12,f'S{i}')
 if flags:kb.create_quantum_register(len(flags),'flags')
 pre=ROOT/f'parameters/{stem}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx');empty=ROOT/f'parameters/{stem}.empty';empty.write_text('');out=ROOT/f'candidates/{stem}.kmx'
 subprocess.run([str(ROOT/'translate_prefix'),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(bb.next_out),str(empty),str(len(rows))],check=True)
 c=km.Circuit.from_file(str(out),format='kmx');c.to_file(str(out.with_suffix('.kmb')),format='kmb');return c


def main():
 ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,default=22);ap.add_argument('--shots',type=int,default=512);ap.add_argument('--profile',action='store_true');a=ap.parse_args();start=time.time()
 if a.profile:
  rows=[]
  for baseline in [True,False]:
   c=build(a.n,'forward',baseline=baseline,prefix4=True)
   rows.append(dict(baseline=baseline,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like))
  report=dict(scope='resource-only initial four-degree prefix; no full decoder score',N=a.n,rows=rows,seconds=time.time()-start)
  Path('results/opening_lookahead_prefix_profile.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));return
 rng=np.random.default_rng(831915);m=load();inputs=rng.integers(0,4096,(a.n,a.shots));outputs=inputs.copy()
 for i in range(a.shots):
  if i%5==0:inputs[-1,i]=0
  elif i%5==1:inputs[-2:,i]=0
  elif i%5==2:inputs[:,i]=0
  F=list(reversed(inputs[:,i].tolist()));coords=reference(m,F[:6]);v1,v0,u,v,pi,v2=coords
  prod=series_product(F,[1,v2,v1,v0],a.n)
  outputs[:,i]=inputs[:,i];outputs[:-6,i]=list(reversed(prod[6:]));outputs[-6:,i]=list(reversed(coords))
 bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool);reports=[]
 for which in ['forward','inverse','both','inverse_both']:
  c=build(a.n,which);sim=km.Simulator(batch_size=a.shots,seed=728);sim.use_same_registers_as(c);initial=outputs if which=='inverse' else inputs;want=outputs if which=='forward' else inputs
  for i,values in enumerate(initial):sim.write_across_shots(c.register_data[f'S{i}'],bits(values))
  sim.do(c);failures={}
  for i,values in enumerate(want):
   bad=np.flatnonzero((sim.read_across_shots(c.register_data[f'S{i}'])!=bits(values)).any(0))
   if len(bad):failures[f'S{i}']=bad[:8].tolist()
  public={q for reg in c.register_data.values() for q in reg};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);junk=np.flatnonzero(sim.read_across_shots(scratch).any(0));phase=np.flatnonzero(sim.read_phase_flipped_across_shots());report=dict(which=which,status='PASS' if not failures and not len(junk) and not len(phase) else 'FAIL',Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,failures=failures,scratch_errors=junk[:8].tolist(),phase_errors=phase[:8].tolist());reports.append(report);print(json.dumps(report),flush=True)
 result=dict(status='PASS' if all(r['status']=='PASS' for r in reports) else 'FAIL',scope='initial six-moment transform plus remainder-tail update, independent series-product model, forward/inverse/fresh arbitrary-input compositions, scratch and phase; no complete decoder',N=a.n,shots=a.shots,circuits=reports,seconds=time.time()-start);Path('results/opening_lookahead_check.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
