"""Tower inverse context with direct classical-mask parity cleanup."""
from contextlib import contextmanager
import kickmix as km
from low_workspace_inverse import norm_rt,inverse_output_rt
from inverse_mask_kernel import phase_inverse
from tower_inverse import apply_map,map_ops


def bind(m):
 @contextmanager
 def inv_word(E,src):
  assert not E.J.recording
  names=[('low_inv_N',6),('low_inv_M',3),('low_inv_S',3),('low_inv_Y',6)];wins=[E.A.open(name,n) for name,n in names];n,mm,s,y=wins
  rN,rM,rS=m.norm_chain_terms(src,m.spec_rows(m.A0ROWS,src),m.spec_rows(m.A1ROWS,src),n,mm,s);rY=m.chain_rts(src,wins,src)[-2]
  E.macs([rN,rM,rS,rY]);E.frees([rN,rM,rS])
  for (name,_),word in reversed(list(zip(names[:3],wins[:3]))):E.A.close(name,word)
  out=E.A.open('low_inv_out',12);E.mac(inverse_output_rt(m,src,y,out));E.free(norm_rt(m,out,y));E.A.close(names[3][0],y)
  yield out
  b=E.b.raw;mask=b.alloc_bits(12);b.hmr(km.array([km.q(q) for q in out]),mask);E.A.close('low_inv_out',out)
  from extension_inverse import setup
  _,_,fw,bk,_=setup();source=km.array([km.q(q) for q in src]);bops=map_ops(bk);trans=tuple(sum(((fw[i]>>j)&1)<<i for i in range(12)) for j in range(12));mops=map_ops(trans)
  apply_map(b,source,bops);apply_map(b,mask,mops);phase_inverse(b,source,mask);apply_map(b,mask,mops,True);apply_map(b,source,bops,True);b.bit_store(mask,False);b.free(mask)
 return inv_word


def test():
 import json,numpy as np
 from pathlib import Path
 from extension_inverse import Adapter,setup
 from fast_dqi import FastBuilder
 from low_workspace_inverse import bind as old_bind
 m=setup()[0];rows=[]
 for tag,context in [('old',old_bind),('fast',bind)]:
  b=FastBuilder(km.CircuitBuilder());src=b.create_quantum_register(12,'src');target=b.create_quantum_register(12,'target');ab=Adapter(b);sc=ab.qubits(3);E=m.Em(ab,*sc,m.Alloc(ab))
  with context(m)(E,[q.id for q in src]) as inv:
   for q,r in zip(inv,target):E.cx(q,r.id)
  ab.release_recycle('inverse_context_emitter',sc);c=b.finish_circuit();c=km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)));shots=8192;x=np.tile(np.arange(4096),2);rng=np.random.default_rng(142092);y=rng.integers(0,4096,shots);want=y^np.array([m.ginv(int(v)) for v in x]);bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool);sim=km.Simulator(batch_size=shots,seed=142092);sim.use_same_registers_as(c);sim.write_across_shots(src,bits(x));sim.write_across_shots(target,bits(y));sim.do(c);public=set(src)|set(target);junk=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);r=dict(tag=tag,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,source=bool((sim.read_across_shots(src)!=bits(x)).any()),value=bool((sim.read_across_shots(target)!=bits(want)).any()),scratch=bool(sim.read_across_shots(junk).any()),phase=bool(sim.read_phase_flipped_across_shots().any()));print(json.dumps(r),flush=True);rows.append(r);assert not any(r[k] for k in ['source','value','scratch','phase'])
 Path('results/low_workspace_inverse_fast_native.json').write_text(json.dumps(dict(status='PASS',shots=shots,rows=rows,scope='every native field input twice, arbitrary XOR target, source/scratch/phase restoration'),indent=2)+'\n')
if __name__=='__main__':test()
