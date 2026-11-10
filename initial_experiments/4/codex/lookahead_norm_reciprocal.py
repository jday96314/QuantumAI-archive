"""Reciprocal of a binary quadratic norm, retaining its coefficients.

For Q=d^2+e*d*v+e^2*u, map (e,d) to p*(e,d)/Q.
This involution supplies a reversible extension at Q=0 and (e,d)=(0,0).
The multiplier p must be nonzero. No complete decoder claim is made.
"""
import argparse,json,subprocess,time
from pathlib import Path
import numpy as np
import kickmix as km
from import_published import load
from low_workspace_inverse import bind
from monic_euclid_step import inverse_keep,replace_inverse

ROOT=Path(__file__).resolve().parent


def square_xor(m,E,source,target):
 for bit,q in enumerate(source):
  value=m.gmul(1<<bit,1<<bit)
  for j in range(12):
   if value>>j&1:E.cx(q,target[j])


def invert(m,E,word):
 with m.opened(E,'norm_reciprocal_inverse',12) as wi:
  inverse_keep(m,E,word,wi);replace_inverse(m,E,word,wi)


def apply(m,E,e,d,u,v,p,p_inverse=None):
 def mac(x,y,target):E.xor(m.prod_rt(target,x,y))
 def scale(word,factor,reciprocal):
  with m.opened(E,'norm_reciprocal_product',12) as product:
   m.mon_slot(E,word,factor,reciprocal,product)
 with m.opened(E,'norm_reciprocal_flags',2) as flags:
  swap,bad=flags
  with m.eq_test(E,e,0) as zero:E.cx(zero,swap)
  E.swaps(swap,list(zip(e,d)))
  with m.inv_word(E,e) as ei:scale(d,ei,e)
  with m.opened(E,'norm_reciprocal_H',12) as H:
   def emit_H():
    square_xor(m,E,d,H);mac(d,v,H)
    for q,t in zip(u,H):
     E.cx(q,t);E.b.ccx(swap,q,t)
    E.cx(swap,H[0])
   emit_H()
   with m.eq_test(E,H,0) as zero:E.cx(zero,bad)
   E.cx(bad,H[0])
   with m.inv_word(E,H) as hi:scale(e,H,hi)
   E.cx(bad,H[0])
   with m.eq_test(E,H,0) as zero:E.cx(zero,bad)
   emit_H()
  invert(m,E,e)
  if p_inverse is None:
   with m.inv_word(E,p) as pi:scale(e,p,pi)
  else:scale(e,p,p_inverse)
  with m.inv_word(E,e) as ei:scale(d,e,ei)
  E.swaps(swap,list(zip(e,d)))
  with m.eq_test(E,e,0) as zero:E.cx(zero,swap)


def reference(m,e,d,u,v,p):
 if not(e or d):return [0,0,u,v,p]
 Q=m.gmul(d,d)^m.gmul(m.gmul(e,d),v)^m.gmul(m.gmul(e,e),u)
 if Q:s=m.gmul(p,m.ginv(Q))
 else:
  pivot=e or d;s=m.gmul(p,m.ginv(m.gmul(pivot,pivot)))
 return [m.gmul(e,s),m.gmul(d,s),u,v,p]


def build(which):
 m=load();m.inv_word=bind(m);bb=m.Builder();words=[bb.qubits(12) for _ in range(5)];E=m.Em(bb,*bb.qubits(3),m.Alloc(bb))
 for _ in range(2 if which=='both' else 1):apply(m,E,*words)
 rows=np.asarray(bb._items._rows,dtype=np.uint32).reshape(-1,9);stem='lookahead_norm_reciprocal_'+which;tape=ROOT/f'parameters/{stem}.npy';np.save(tape,rows)
 kb=km.CircuitBuilder()
 for name in ['e','d','u','v','p']:kb.create_quantum_register(12,name)
 pre=ROOT/f'parameters/{stem}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx');empty=ROOT/f'parameters/{stem}.empty';empty.write_text('');out=ROOT/f'candidates/{stem}.kmx'
 subprocess.run([str(ROOT/'translate_prefix'),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(bb.next_out),str(empty),str(len(rows))],check=True)
 c=km.Circuit.from_file(str(out),format='kmx');c.to_file(str(out.with_suffix('.kmb')),format='kmb');return c


def main():
 ap=argparse.ArgumentParser();ap.add_argument('--shots',type=int,default=4096);args=ap.parse_args();rng=np.random.default_rng(781415);m=load();inputs=[];outputs=[];begin=time.time()
 for i in range(args.shots):
  e,d,u,v,p=map(int,rng.integers(0,4096,5));p=p or 1
  if i%6==0:e=0
  elif i%6==1:e=d=0
  elif i%6==2:
   e=e or 1;u=m.gmul(m.gmul(d,d)^m.gmul(m.gmul(e,d),v),m.ginv(m.gmul(e,e)))
  elif i%6==3:d=0
  row=[e,d,u,v,p];out=reference(m,*row);assert reference(m,*out)==row
  inputs.append(row);outputs.append(out)
 inputs=np.array(inputs).T;outputs=np.array(outputs).T;bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool);reports=[]
 for which in ['forward','inverse','both']:
  c=build(which);sim=km.Simulator(batch_size=args.shots,seed=715);sim.use_same_registers_as(c);initial=outputs if which=='inverse' else inputs;want=outputs if which=='forward' else inputs
  for name,values in zip(['e','d','u','v','p'],initial):sim.write_across_shots(c.register_data[name],bits(values))
  sim.do(c);failures={}
  for name,values in zip(['e','d','u','v','p'],want):
   bad=np.flatnonzero((sim.read_across_shots(c.register_data[name])!=bits(values)).any(0))
   if len(bad):failures[name]=bad[:8].tolist()
  public={q for reg in c.register_data.values() for q in reg};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);junk=np.flatnonzero(sim.read_across_shots(scratch).any(0));phase=np.flatnonzero(sim.read_phase_flipped_across_shots());report=dict(which=which,status='PASS' if not failures and not len(junk) and not len(phase) else 'FAIL',Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,failures=failures,scratch_errors=junk[:8].tolist(),phase_errors=phase[:8].tolist());reports.append(report);print(json.dumps(report),flush=True)
 result=dict(status='PASS' if all(r['status']=='PASS' for r in reports) else 'FAIL',scope='quadratic norm reciprocal with nonzero p, including zero norm and zero pair; forward/inverse/composition, scratch and phase; no decoder acceptance',shots=args.shots,circuits=reports,seconds=time.time()-begin);Path('results/lookahead_norm_reciprocal_check.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
