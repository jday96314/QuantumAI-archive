"""In-place two-moment Padé solve with a pivot selected before division.

Kept moments e,d and live moments a,b map to e,d,u,v with
u=(d*b+a*a)/(d*d+e*a), v=(e*b+d*a)/(d*d+e*a).
The first pivot e may be zero; then d supplies it. A reversible extension
handles a singular second minor, too, while only the nonzero-minor branch
has the stated Padé meaning. This is not a complete decoder.
"""
import argparse,json,subprocess,time
from pathlib import Path
import numpy as np
import kickmix as km
from import_published import load
from low_workspace_inverse import bind
from tower_inverse import map_ops
ROOT=Path(__file__).resolve().parent


def apply(m,E,e,d,a,b,inverse=False):
 def mac(x,y,target):E.xor(m.prod_rt(target,x,y))
 def square_xor(source,target):
  for bit,q in enumerate(source):
   value=m.gmul(1<<bit,1<<bit)
   for j in range(12):
    if value>>j&1:E.cx(q,target[j])
 def square_ops(word,inv=False):
  ops=map_ops([m.gmul(1<<j,1<<j) for j in range(12)])
  for name,i,j in (reversed(ops) if inv else ops):getattr(E.b,name)(word[i],word[j])
 def masked_square(flag,source,target):
  square_ops(source)
  for q,t in zip(source,target):E.b.ccx(flag,q,t)
  square_ops(source,True)
 with m.opened(E,'lookahead_D2',12) as D2,m.opened(E,'lookahead_pivot_flag',3) as fw:
  flag,minor_bad,pivot_bad=fw
  if inverse:
   with m.opened(E,'lookahead_original_a',12) as olda:
    rt=m.prod_rt(olda,e,a)+m.prod_rt(olda,d,b);E.mac(rt);E.mac(m.prod_rt(D2,e,olda));E.free(rt)
  else:E.mac(m.prod_rt(D2,e,a))
  square_xor(d,D2)
  with m.eq_test(E,D2,0) as zero:E.cx(zero,minor_bad)
  E.cx(minor_bad,D2[0])
  with m.eq_test(E,e,0) as zero:E.cx(zero,flag)
  E.swaps(flag,list(zip(e,d)))
  with m.eq_test(E,e,0) as zero:E.cx(zero,pivot_bad)
  E.cx(pivot_bad,e[0])
  with m.inv_word(E,e) as pi:
   if inverse:E.swaps(flag,list(zip(a,b)))
   with m.opened(E,'lookahead_product',12) as prod:
    if inverse:
     m.mon_slot(E,a,e,pi,prod);mac(d,b,a)
    else:
     m.mon_slot(E,b,e,pi,prod);mac(d,a,b);masked_square(flag,a,b)
   with m.inv_word(E,D2) as di,m.opened(E,'lookahead_product',12) as prod:
    m.mon_slot(E,b,D2 if inverse else di,di if inverse else D2,prod)
   with m.opened(E,'lookahead_product',12) as prod:
    if inverse:
     masked_square(flag,a,b);mac(d,a,b);m.mon_slot(E,b,pi,e,prod)
    else:
     mac(d,b,a);m.mon_slot(E,a,pi,e,prod)
   if not inverse:E.swaps(flag,list(zip(a,b)))
  E.cx(pivot_bad,e[0])
  with m.eq_test(E,e,0) as zero:E.cx(zero,pivot_bad)
  E.swaps(flag,list(zip(e,d)))
  with m.eq_test(E,e,0) as zero:E.cx(zero,flag)
  E.cx(minor_bad,D2[0])
  with m.eq_test(E,D2,0) as zero:E.cx(zero,minor_bad)
  square_xor(d,D2)
  if inverse:E.free(m.prod_rt(D2,e,a))
  else:
   with m.opened(E,'lookahead_original_a',12) as olda:
    rt=m.prod_rt(olda,e,a)+m.prod_rt(olda,d,b);E.mac(rt);E.free(m.prod_rt(D2,e,olda));E.free(rt)


def build(which):
 m=load();m.inv_word=bind(m);bb=m.Builder();words=[bb.qubits(12) for _ in range(4)];E=m.Em(bb,*bb.qubits(3),m.Alloc(bb))
 for inv in ([False,True] if which=='both' else [which=='inverse']):apply(m,E,*words,inv)
 rows=np.asarray(bb._items._rows,dtype=np.uint32).reshape(-1,9);stem='lookahead_pair_'+which;tape=ROOT/f'parameters/{stem}.npy';np.save(tape,rows)
 kb=km.CircuitBuilder()
 for name in ['e','d','a','b']:kb.create_quantum_register(12,name)
 pre=ROOT/f'parameters/{stem}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx');empty=ROOT/f'parameters/{stem}.empty';empty.write_text('');out=ROOT/f'candidates/{stem}.kmx'
 subprocess.run([str(ROOT/'translate_prefix'),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(bb.next_out),str(empty),str(len(rows))],check=True)
 c=km.Circuit.from_file(str(out),format='kmx');c.to_file(str(out.with_suffix('.kmb')),format='kmb');return c


def main():
 ap=argparse.ArgumentParser();ap.add_argument('--shots',type=int,default=4096);args=ap.parse_args();rng=np.random.default_rng(613715);m=load();inputs=[];outputs=[];skipped=0;begin=time.time()
 while len(inputs)<args.shots:
  e,d,a,b=map(int,rng.integers(0,4096,4));kind=len(inputs)%6
  if kind==0:e=0
  elif kind==1:d=0
  elif kind==2:
   e=e or 1;a=m.gmul(m.gmul(d,d),m.ginv(e))
  elif kind==3:e=d=0
  elif kind==4:d=a=0
  D2=m.gmul(d,d)^m.gmul(e,a)
  if D2:
   di=m.ginv(D2);u=m.gmul(m.gmul(d,b)^m.gmul(a,a),di);v=m.gmul(m.gmul(e,b)^m.gmul(d,a),di)
  elif e:
   v=m.gmul(e,b)^m.gmul(d,a);u=m.gmul(a^m.gmul(d,v),m.ginv(e))
  else:u=b^m.gmul(a,a);v=a
  inputs.append([e,d,a,b]);outputs.append([e,d,u,v])
 inputs=np.array(inputs).T;outputs=np.array(outputs).T;bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool);reports=[]
 for which in ['forward','inverse','both']:
  c=build(which);sim=km.Simulator(batch_size=args.shots,seed=714);sim.use_same_registers_as(c);initial=outputs if which=='inverse' else inputs;want=outputs if which=='forward' else inputs
  for name,values in zip(['e','d','a','b'],initial):sim.write_across_shots(c.register_data[name],bits(values))
  sim.do(c);failures={}
  for name,values in zip(['e','d','a','b'],want):
   bad=np.flatnonzero((sim.read_across_shots(c.register_data[name])!=bits(values)).any(0))
   if len(bad):failures[name]=bad[:8].tolist()
  public={q for reg in c.register_data.values() for q in reg};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);junk=np.flatnonzero(sim.read_across_shots(scratch).any(0));phase=np.flatnonzero(sim.read_phase_flipped_across_shots());report=dict(which=which,status='PASS' if not failures and not len(junk) and not len(phase) else 'FAIL',Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,failures=failures,scratch_errors=junk[:8].tolist(),phase_errors=phase[:8].tolist());reports.append(report);print(json.dumps(report),flush=True)
 result=dict(status='PASS' if all(r['status']=='PASS' for r in reports) else 'FAIL',scope='in-place first-Pade solve on nonzero second minor, total reversible extension on singular minors and both zero kept moments, forward/inverse/composition/scratch/phase; no complete decoder claim',shots=args.shots,skipped_singular_second_minor=skipped,circuits=reports,seconds=time.time()-begin);Path('results/lookahead_pair_quantum_check.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
