"""Classical-mask inverse parity with a 31-AND cubic trace kernel.

A classical Gaussian elimination factors multiplication by the trace-dual
mask. Its recorded row operations scale the quantum source in place, avoiding
a quantum scaling register. All classical matrices depend only on the mask.
"""
import inspect,json
from functools import lru_cache
from pathlib import Path
import numpy as np
import kickmix as km
import cubic_inverse_product_kernel as base
from cubic_trace_field import data
from tower_inverse import map_ops,apply_map
from fast_dqi import FastBuilder

# Retain the attributed cubic Kummer norm/inversion construction; replace its
# root-gated parity with the ungated trace dot product.
s=inspect.getsource(base.trace_inverse)
a=s.index('  ops=map_ops(');z=s.index('\n nm=',a)
s=s[:a]+"  ops=map_ops(z['trace4_columns']);apply_map(b,Q,ops);b.cz(Q,iv);apply_map(b,Q,ops,True)"+s[z:]
ns=dict(vars(base));exec(s,ns);trace_inverse=ns['trace_inverse']


@lru_cache(None)
def matrix_masks():
 f=km.GF2Field(12)
 return tuple(tuple(sum(((f.mul(1<<j,1<<k)>>i)&1)<<k for k in range(12)) for j in range(12)) for i in range(12))


def classical_scale_inverse(b,source,coefficient):
 masks=matrix_masks();matrix=[b.alloc_bits(12) for _ in range(12)];ops=[];cb=base.ClassicalBuilder(b)
 for row,rr in zip(matrix,masks):
  for cell,mask in zip(row,rr):
   b.bit_store(cell,False);b.cx(km.array([coefficient[k] for k in range(12) if mask>>k&1]),cell)
 for j in range(12):
  for i in range(j+1,12):
   flag=b.alloc_bits(1)[0];b.x(matrix[j][j]);cb.init_and(matrix[j][j],matrix[i][j],flag);b.x(matrix[j][j]);b.push_condition(flag);b.swap(matrix[j],matrix[i]);b.swap(source[j],source[i]);b.pop_condition();ops.append(('swap',j,i,flag))
  for i in range(12):
   if i==j:continue
   flag=b.alloc_bits(1)[0];b.bit_store(flag,False);b.cx(matrix[i][j],flag);b.push_condition(flag);b.cx(matrix[j],matrix[i]);b.cx(source[j],source[i]);b.pop_condition();ops.append(('cx',j,i,flag))
 for row in matrix:b.bit_store(row,False);b.free(row)
 return ops


def phase_inverse(b,source,mask):
 active=b.alloc_bits(1)[0];b.bit_store(active,False);b.bit_store(active,True,control=mask);b.push_condition(active);ops=map_ops(data()['dual_columns']);apply_map(b,mask,ops);record=classical_scale_inverse(b,source,mask);trace_inverse(b,None,source)
 for kind,j,i,flag in reversed(record):
  b.push_condition(flag);getattr(b,kind)(source[j],source[i]);b.pop_condition();b.bit_store(flag,False);b.free(km.array([flag]))
 apply_map(b,mask,ops,True);b.pop_condition();b.bit_store(active,False);b.free(km.array([active]))


def test():
 b=FastBuilder(km.CircuitBuilder());source=b.create_quantum_register(12,'source');mask=b.create_classical_register(12,'mask');phase_inverse(b,source,mask);c=b.finish_circuit();f=km.GF2Field(12);values=np.arange(4096);iv=np.array([f.pow(int(v),4094) for v in values]);rng=np.random.default_rng(142107);bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool);public=set(source);scratch=km.array([km.q(q) for q in range(c.num_qubits) if km.q(q) not in public]);rows=[]
 for j,mm in enumerate([np.full(4096,1<<i,dtype=int) for i in range(12)]+[np.zeros(4096,dtype=int),rng.integers(0,4096,4096)]):
  sim=km.Simulator(batch_size=4096,seed=142107+j);sim.use_same_registers_as(c);sim.write_across_shots(source,bits(values));sim.write_across_shots(mask,bits(mm));sim.do(c);want=np.array([(int(v)&int(w)).bit_count()%2 for v,w in zip(iv,mm)],dtype=bool);r=dict(case=j,value=bool((sim.read_across_shots(source)!=bits(values)).any()),mask=bool((sim.read_across_shots(mask)!=bits(mm)).any()),scratch=bool(sim.read_across_shots(scratch).any()),phase=bool((sim.read_phase_flipped_across_shots()!=want).any()));rows.append(r);print(json.dumps(r),flush=True);assert not any(r[k] for k in ['value','mask','scratch','phase'])
 r=dict(status='PASS',shots=4096*len(rows),Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,rows=rows);Path('results/cubic_mask_phase_native.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
if __name__=='__main__':test()


def native_phase_inverse(b,source,mask,fw_unused=None,bk_unused=None):
 from extension_inverse import setup
 m,k,fw,bk,tr=setup();word=km.array([km.q(q) for q in source]);bops=map_ops(bk);trans=tuple(sum(((fw[i]>>j)&1)<<i for i in range(12)) for j in range(12));mops=map_ops(trans);apply_map(b,word,bops);apply_map(b,mask,mops);phase_inverse(b,word,mask);apply_map(b,mask,mops,True);apply_map(b,word,bops,True)
