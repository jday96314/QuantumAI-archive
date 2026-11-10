"""In-place six-moment Padé coordinates with a total reversible extension.

The nonzero D2,p branch agrees with lookahead_moment_map.forward, including
e=0. Output physical order is v1,v0,u0,u1,pinv,v2. This is a research
component; no complete global-cofactor adapter or decoder is claimed.
"""
import argparse,json,subprocess,time
from pathlib import Path
import numpy as np
import kickmix as km
from import_published import load
from low_workspace_inverse import bind
from lookahead_pair_quantum import apply as pair
from lookahead_norm_reciprocal import apply as norm_reciprocal,invert
from tower_inverse import map_ops
from lookahead_moment_map import forward as pade_forward

ROOT=Path(__file__).resolve().parent


def apply(m,E,e,d,u,v,p,alpha,inverse=False):
 def mac(x,y,target):E.xor(m.prod_rt(target,x,y))
 def square(word,undo=False):
  ops=map_ops([m.gmul(1<<j,1<<j) for j in range(12)])
  for name,i,j in (reversed(ops) if undo else ops):getattr(E.b,name)(word[i],word[j])
 def shear_p():
  with m.opened(E,'six_original_a',12) as a:
   ra=m.prod_rt(a,e,u)+m.prod_rt(a,d,v);E.mac(ra)
   square(v);mac(a,v,p);square(v,True)
   square(u);mac(e,u,p);square(u,True)
   E.free(ra)
 def shear_alpha():
  with m.opened(E,'six_alpha_T',12) as T,m.opened(E,'six_alpha_K',12) as K,m.opened(E,'six_alpha_H',12) as H:
   rt=m.prod_rt(T,e,v);E.mac(rt)
   for q,t in zip(d,T):E.cx(q,t)
   square(v);rk=m.prod_rt(K,u,v);E.mac(rk)
   for q,t in zip(u,H):E.cx(q,t)
   for q,t in zip(v,H):E.cx(q,t)
   square(H);mac(d,H,alpha);square(H,True);mac(T,K,alpha)
   for q,t in zip(v,H):E.cx(q,t)
   for q,t in zip(u,H):E.cx(q,t)
   E.free(rk);square(v,True)
   for q,t in zip(d,T):E.cx(q,t)
   E.free(rt)
 def scale_alpha(pi,undo=False):
  with m.opened(E,'six_scale_product',12) as product:
   m.mon_slot(E,alpha,p if undo else pi,pi if undo else p,product)
 def output_shear(undo=False):
  if undo:
   for q,t in zip(v,alpha):E.cx(q,t)
  mac(alpha,u,d);mac(alpha,v,e)
  for q,t in zip(u,e):E.cx(q,t)
  if not undo:
   for q,t in zip(v,alpha):E.cx(q,t)
 if not inverse:
  pair(m,E,e,d,u,v);shear_p()
 with m.opened(E,'six_zero_p',1) as flagword:
  flag=flagword[0]
  with m.eq_test(E,p,0) as zero:E.cx(zero,flag)
  E.cx(flag,p[0])
  if inverse:
   invert(m,E,p);output_shear(True)
   with m.inv_word(E,p) as pi:
    norm_reciprocal(m,E,e,d,u,v,p,pi);scale_alpha(pi,True)
   shear_alpha()
  else:
   shear_alpha()
   with m.inv_word(E,p) as pi:
    scale_alpha(pi);norm_reciprocal(m,E,e,d,u,v,p,pi)
   output_shear();invert(m,E,p)
  E.cx(flag,p[0])
  with m.eq_test(E,p,0) as zero:E.cx(zero,flag)
 if inverse:
  shear_p();pair(m,E,e,d,u,v,True)


def reference(m,row):
 e,d,a,b,c,f=row;mul=m.gmul;inv=m.ginv
 D2=mul(d,d)^mul(e,a)
 if D2:
  u=mul(mul(d,b)^mul(a,a),inv(D2));v=mul(mul(e,b)^mul(d,a),inv(D2))
 elif e:
  v=mul(e,b)^mul(d,a);u=mul(a^mul(d,v),inv(e))
 else:u=b^mul(a,a);v=a
 # Reconstruct using the polynomial formulas even on singular extensions.
 aa=mul(e,u)^mul(d,v);bb=mul(d,u)^mul(aa,v)
 p=c^mul(aa,u)^mul(bb,v);armed=p or 1
 alpha=mul(f^mul(d,mul(u,u))^mul(bb,mul(v,v)),inv(armed))
 from lookahead_norm_reciprocal import reference as norm_ref
 re,rd,_,_,_=norm_ref(m,e,d,u,v,armed)
 return [re^u^mul(alpha,v),rd^mul(alpha,u),u,v,inv(p) if p else 0,alpha^v]


def build(which):
 m=load();m.inv_word=bind(m);bb=m.Builder();words=[bb.qubits(12) for _ in range(6)];E=m.Em(bb,*bb.qubits(3),m.Alloc(bb))
 for direction in ([False,True] if which=='both' else [which=='inverse']):apply(m,E,*words,direction)
 rows=np.asarray(bb._items._rows,dtype=np.uint32).reshape(-1,9);stem='lookahead_six_'+which;tape=ROOT/f'parameters/{stem}.npy';np.save(tape,rows)
 kb=km.CircuitBuilder()
 for name in ['e','d','a','b','c','f']:kb.create_quantum_register(12,name)
 pre=ROOT/f'parameters/{stem}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx');empty=ROOT/f'parameters/{stem}.empty';empty.write_text('');out=ROOT/f'candidates/{stem}.kmx'
 subprocess.run([str(ROOT/'translate_prefix'),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(bb.next_out),str(empty),str(len(rows))],check=True)
 c=km.Circuit.from_file(str(out),format='kmx');c.to_file(str(out.with_suffix('.kmb')),format='kmb');return c


def main():
 ap=argparse.ArgumentParser();ap.add_argument('--shots',type=int,default=8192);args=ap.parse_args();rng=np.random.default_rng(815913);m=load();inputs=[];outputs=[];begin=time.time();pade_checks=0
 for i in range(args.shots):
  row=list(map(int,rng.integers(0,4096,6)));kind=i%7
  if kind==0:row[0]=0
  elif kind==1:row[0]=row[1]=0
  elif kind==2:
   row[0]=row[0] or 1;row[2]=m.gmul(m.gmul(row[1],row[1]),m.ginv(row[0]))
  elif kind==3:row=[0]*6
  elif kind==4:
   preliminary=reference(m,row);u,v=preliminary[2:4];e,d=row[:2]
   aa=m.gmul(e,u)^m.gmul(d,v);bb=m.gmul(d,u)^m.gmul(aa,v)
   row[4]=m.gmul(aa,u)^m.gmul(bb,v)
  out=reference(m,row)
  try:
   u,v,v0,v1,v2,pi=pade_forward(row)
   assert out==[v1,v0,u,v,pi,v2],(row,out)
   pade_checks+=1
  except ValueError:pass
  inputs.append(row);outputs.append(out)
 inputs=np.array(inputs).T;outputs=np.array(outputs).T;bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool);reports=[]
 for which in ['forward','inverse','both']:
  c=build(which);sim=km.Simulator(batch_size=args.shots,seed=719);sim.use_same_registers_as(c);initial=outputs if which=='inverse' else inputs;want=outputs if which=='forward' else inputs
  for name,values in zip(['e','d','a','b','c','f'],initial):sim.write_across_shots(c.register_data[name],bits(values))
  sim.do(c);failures={}
  for name,values in zip(['e','d','a','b','c','f'],want):
   bad=np.flatnonzero((sim.read_across_shots(c.register_data[name])!=bits(values)).any(0))
   if len(bad):failures[name]=bad[:8].tolist()
  public={q for reg in c.register_data.values() for q in reg};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);junk=np.flatnonzero(sim.read_across_shots(scratch).any(0));phase=np.flatnonzero(sim.read_phase_flipped_across_shots());report=dict(which=which,status='PASS' if not failures and not len(junk) and not len(phase) else 'FAIL',Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,failures=failures,scratch_errors=junk[:8].tolist(),phase_errors=phase[:8].tolist());reports.append(report);print(json.dumps(report),flush=True)
 result=dict(status='PASS' if all(r['status']=='PASS' for r in reports) else 'FAIL',scope='six-moment in-place transform, total singular extension, independent Pade reference on nonsingular D2,p, forward/inverse/composition, scratch and phase; no full decoder',shots=args.shots,pade_checks=pade_checks,circuits=reports,seconds=time.time()-begin);Path('results/lookahead_six_quantum_check.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
