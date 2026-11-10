"""Reversible degree-(2,1) cofactor realignment for one missing Euclid pivot.

New round-5 component: consumes six words, including one retained exceptional
coordinate, and restores the ordinary consecutive degree layout after three
degree slots. Every product is emitted and charged on every circuit execution.
The full decoder still needs a total guarded local-ratio adapter and dispatch.
"""
import json,subprocess
from pathlib import Path
import numpy as np,kickmix as km
from import_published import load
from low_workspace_inverse import bind
from classical_repair import product,plus,mul,inv
ROOT=Path(__file__).resolve().parent

def apply(m,E,Bp,B,Q,pinv,inverse=False,norm=True):
 beta,q0,q1,gamma,alpha,extra=Q
 mid=Bp+[beta,q0,q1];final=B+[gamma,extra,alpha]
 def xor(a,b):
  for q,t in zip(a,b):E.cx(q,t)
 def mac(a,b,t):E.xor(m.prod_rt(t,a,b))
 def quadratic(back=False):
  k=len(B)
  def top():
   def change_beta():
    if len(B)>2:xor(B[-3],beta)
    if len(B)>1:mac(q1,B[-2],beta)
    mac(q0,B[-1],beta)
   if not back:change_beta()
   if back:xor(B[-1],q1)
   if len(B)>1:xor(B[-2],q0)
   mac(q1,B[-1],q0)
   if not back:xor(B[-1],q1)
   if back:change_beta()
  if back:top()
  with m.inv_word(E,beta) as bi,m.opened(E,'gap_quadratic_product',12) as t:
   for i,w in enumerate(Bp):
    if not back:m.mon_slot(E,w,beta,bi,t)
    mac(q0,B[i],w)
    if i:mac(q1,B[i-1],w)
    if i>1:xor(B[i-2],w)
    if back:m.mon_slot(E,w,bi,beta,t)
   if norm:m.mon_slot(E,pinv,beta if back else bi,bi if back else beta,t)
  if not back:top()
 def linear(back=False):
  def top():
   if back:xor(mid[-1],alpha)
   mac(alpha,mid[len(B)],gamma)
   xor(mid[len(B)-1],gamma)
   mac(alpha,mid[-1],extra);xor(mid[-2],extra)
   if not back:xor(mid[-1],alpha)
  if back:top()
  with m.inv_word(E,gamma) as gi,m.opened(E,'gap_linear_product',12) as t:
   for i,w in enumerate(B):
    if not back:m.mon_slot(E,w,gamma,gi,t)
    mac(alpha,mid[i],w)
    if i:xor(mid[i-1],w)
    if back:m.mon_slot(E,w,gi,gamma,t)
   if norm:m.mon_slot(E,pinv,gamma if back else gi,gi if back else gamma,t)
  if not back:top()
 if inverse:linear(True);quadratic(True)
 else:quadratic();linear()
 return mid,final

def build(k,which):
 m=load();m.inv_word=bind(m);b=m.Builder();N=2*k+6;words=[b.qubits(12) for _ in range(N)];E=m.Em(b,*b.qubits(3),m.Alloc(b))
 bp=words[:k-1];bb=words[k-1:2*k-1];Q=words[2*k-1:2*k+5];pi=words[-1]
 for back in ([False,True] if which=='both' else [True,False] if which=='inverse_both' else [which=='inverse']):apply(m,E,bp,bb,Q,pi,back)
 rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9);stem=f'gap_cofactor_k{k}_{which}';tape=ROOT/'parameters'/f'{stem}.npy';np.save(tape,rows)
 kb=km.CircuitBuilder()
 for i in range(N):kb.create_quantum_register(12,f'W{i}')
 pre=ROOT/'parameters'/f'{stem}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx');empty=ROOT/'parameters/empty';empty.write_text('');out=ROOT/'candidates'/f'{stem}.kmx'
 subprocess.run([str(ROOT/'translate_prefix'),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True)
 c=km.Circuit.from_file(str(out),format='kmx');c.to_file(str(out.with_suffix('.kmb')),format='kmb');return c

def test(k,shots=4096):
 # Independent model uses the host native field, circuit uses public field.
 from published_dqi import isomorphism
 m=load();fw,bw,tr=isomorphism(m);rng=np.random.default_rng(16512+k);N=2*k+6;start=rng.integers(0,4096,(N,shots));want=start.copy()
 start[2*k-1]=rng.integers(1,4096,shots);start[2*k+2]=rng.integers(1,4096,shots);start[-1]=rng.integers(1,4096,shots)
 # Half the cases carry arbitrary incoming exceptional coordinate; half zero.
 start[2*k+4,::2]=0
 for s in range(shots):
  bp=start[:k-1,s].tolist()+[1];bb=start[k-1:2*k-1,s].tolist()+[1];beta,q0,q1,gam,al,extra=start[2*k-1:2*k+5,s].tolist();pi=int(start[-1,s])
  mid=plus([mul(beta,v) for v in bp],product([q0,q1,1],bb))
  final=plus([mul(gam,v) for v in bb],product([al,1],mid));final += [0]*max(0,k+4-len(final));final[k+1]^=extra
  out=start[:,s].copy();mi=list(range(k-1))+list(range(2*k-1,2*k+2));fi=list(range(k-1,2*k-1))+[2*k+2,2*k+4,2*k+3]
  for i,v in zip(mi,mid[:-1]):out[i]=v
  for i,v in zip(fi,final[:-1]):out[i]=v
  out[-1]=mul(pi,inv(mul(beta,gam)));want[:,s]=out
 convert=lambda a:np.array([[tr(int(v),fw) for v in row] for row in a])
 start=convert(start);want=convert(want);bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool);reports=[]
 for which in ['forward','inverse','both','inverse_both']:
  c=build(k,which);sim=km.Simulator(batch_size=shots,seed=819126);sim.use_same_registers_as(c);source=want if which=='inverse' else start;target=want if which=='forward' else start
  if which=='inverse_both':
   # A fresh arbitrary output checks the other composition order.
   source=want[:,rng.permutation(shots)].copy();target=source
  for i,v in enumerate(source):sim.write_across_shots(c.register_data[f'W{i}'],bits(v))
  sim.do(c);errors={}
  for i,v in enumerate(target):
   bad=np.flatnonzero((sim.read_across_shots(c.register_data[f'W{i}'])!=bits(v)).any(0))
   if len(bad):errors[i]=bad[:8].tolist()
  public={q for reg in c.register_data.values() for q in reg};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);junk=int(sim.read_across_shots(scratch).any(0).sum());phase=int(sim.read_phase_flipped_across_shots().sum())
  reports.append(dict(which=which,status='PASS' if not errors and not junk and not phase else 'FAIL',errors=errors,junk=junk,phase=phase,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like))
 result=dict(k=k,shots=shots,status='PASS' if all(r['status']=='PASS' for r in reports) else 'FAIL',scope='degree-(2,1) cofactor repair only; exact native-field reference; arbitrary incoming exceptional coordinate, forward, inverse, both composition orders, scratch, phase; not complete decoder',circuits=reports)
 Path(f'results/gap_cofactor_k{k}.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--k',type=int,default=10);p.add_argument('--shots',type=int,default=4096);a=p.parse_args();test(a.k,a.shots)
