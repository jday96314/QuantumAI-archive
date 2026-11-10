"""Reversible conversion of a consumed extension sample to a Schur pivot."""
import json
from pathlib import Path
import numpy as np
import kickmix as km
from extension_inverse import Adapter,setup
from low_workspace_inverse import bind
from base_inverse import inplace as invert12
from tower_inverse import map_ops,apply_map
from fast_dqi import FastBuilder


def build(node,inverse=False,first=False):
 m,k,fw,bk,tr=setup();m.inv_word=bind(m);b=FastBuilder(km.CircuitBuilder())
 lo=b.create_quantum_register(12,'a');hi=b.create_quantum_register(12,'b');g=b.create_quantum_register(12,'g')
 prev_a=b.create_quantum_register(12,'prev_a');prev_b=b.create_quantum_register(12,'prev_b')
 ab=Adapter(b);scratch=ab.qubits(3);E=m.Em(ab,*scratch,m.Alloc(ab));A=[q.id for q in lo];B=[q.id for q in hi];G=[q.id for q in g]
 def axpy(src,dst,c):
  for i,q in enumerate(src):
   value=m.gmul(c,1<<i)
   for j in range(12):
    if value>>j&1:b.cx(q,dst[j])
 def gupdate():
  invert12(b,g)
  with m.eq_test(E,B,0) as z:
   E.cx(z,B[0]);ops=map_ops([m.gpow(1<<i,2) for i in range(12)]);apply_map(b,hi,ops)
   with m.inv_word(E,B) as bi,m.opened(E,'pivot_scale_g',12) as temp:m.mon_slot(E,G,B,bi,temp)
   apply_map(b,hi,ops,True);E.cx(z,B[0])
 def offset():
  if not first:b.cx(prev_a,lo);b.cx(prev_b,lo);b.cx(prev_b,hi)
 def scale(back=False):
  with m.inv_word(E,G) as gi,m.opened(E,'pivot_scale_values',12) as temp:
   for word in [A,B]:m.mon_slot(E,word,gi if back else G,G if back else gi,temp)
 if inverse:
  gupdate();offset();scale(True);axpy(hi,lo,node)
 else:
  axpy(hi,lo,node);scale();offset();gupdate()
 ab.release_recycle('pivot_emitter',scratch);c=b.finish_circuit()
 return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)))

def test():
 m,k,fw,bk,tr=setup();node=tr(26,fw);rng=np.random.default_rng(142038);n=2048
 aa,bb,pa,pb=rng.integers(0,4096,(4,n));gg=rng.integers(1,4096,n)
 # Force the armed zero branch, without selecting a favorable verifier sample.
 bb[::7]=np.array([m.gmul(int(p),m.ginv(int(g))) for p,g in zip(pb[::7],gg[::7])])
 a=np.array([m.gmul(int(g),int(x)^m.gmul(node,int(y)))^int(p)^int(q) for g,x,y,p,q in zip(gg,aa,bb,pa,pb)])
 d=np.array([m.gmul(int(g),int(y))^int(q) for g,y,q in zip(gg,bb,pb)])
 gout=np.array([m.gmul(m.gpow(int(v) or 1,2),m.ginv(int(g))) for g,v in zip(gg,d)])
 bits=lambda x:((x[None,:]>>np.arange(12)[:,None])&1).astype(bool)
 rows=[]
 for inv in [False,True]:
  c=build(node,inv);sim=km.Simulator(batch_size=n,seed=142038+inv);sim.use_same_registers_as(c);regs=c.register_data;src=[a,d,gout,pa,pb] if inv else [aa,bb,gg,pa,pb];want=[aa,bb,gg,pa,pb] if inv else [a,d,gout,pa,pb]
  for name,x in zip(['a','b','g','prev_a','prev_b'],src):sim.write_across_shots(regs[name],bits(x))
  sim.do(c);errors={}
  for name,x in zip(['a','b','g','prev_a','prev_b'],want):
   bad=np.flatnonzero((sim.read_across_shots(regs[name])!=bits(x)).any(0))
   if len(bad):errors[name]=bad[:8].tolist()
  public={q for word in regs.values() for q in word};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);junk=np.flatnonzero(sim.read_across_shots(scratch).any(0));phase=np.flatnonzero(sim.read_phase_flipped_across_shots())
  r=dict(inverse=inv,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,errors=errors,junk=junk[:8].tolist(),phase=phase[:8].tolist());rows.append(r);print(json.dumps(r),flush=True);assert not errors and not len(junk) and not len(phase)
  c.to_file(f'candidates/schur_pivot_{"inverse" if inv else "forward"}.kmb',format='kmb')
 Path('results/schur_pivot_native.json').write_text(json.dumps(dict(status='PASS',shots=4096,rows=rows,scope='arbitrary pivot/frame values, nonzero frame g, forced zero b pivot, both directions'),indent=2)+'\n')
if __name__=='__main__':test()
