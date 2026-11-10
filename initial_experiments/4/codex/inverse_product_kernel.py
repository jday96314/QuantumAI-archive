"""Root-gated classical parity of inverse(a*b), preserving both operands."""
import json
from pathlib import Path
import numpy as np
import kickmix as km
from inverse_mask_kernel import apply as phase_inverse
from tower_inverse import product
from fast_dqi import FastBuilder


def apply(b,root,a,c,mask):
    value=b.alloc_qubits(12);product(b,a,c,value)
    phase_inverse(b,root,value,mask)
    b.del_gf2_mul(a,c,target=value,field=km.GF2Field(12));b.free(value)


def test():
    b=FastBuilder(km.CircuitBuilder());root=b.create_quantum_register(12,'root');a=b.create_quantum_register(12,'a');c=b.create_quantum_register(12,'c');mask=b.create_classical_register(12,'mask');apply(b,root,a,c,mask);circ=b.finish_circuit()
    n=4096;rng=np.random.default_rng(89127);av=np.arange(n);cv=rng.integers(0,n,n);rv=rng.integers(0,n,n);rv[::2]=0;mv=rng.integers(0,n,n)
    bits=lambda a:((a[None,:]>>np.arange(12)[:,None])&1).astype(bool)
    sim=km.Simulator(batch_size=n,seed=51);sim.use_same_registers_as(circ)
    for reg,v in [(root,rv),(a,av),(c,cv),(mask,mv)]:sim.write_across_shots(reg,bits(v))
    sim.do(circ);f=km.GF2Field(12);want=np.array([r==0 and (f.pow(f.mul(int(x),int(y)),4094)&int(m)).bit_count()%2==1 for r,x,y,m in zip(rv,av,cv,mv)])
    assert np.array_equal(sim.read_phase_flipped_across_shots(),want)
    for reg,v in [(root,rv),(a,av),(c,cv),(mask,mv)]:assert np.array_equal(sim.read_across_shots(reg),bits(v))
    public=set(root)|set(a)|set(c);scratch=km.array([km.q(i) for i in range(circ.num_qubits) if km.q(i) not in public]);assert not sim.read_across_shots(scratch).any()
    out=dict(status='PASS',scope='all a field elements, random b and masks, root/non-root paths, source restoration, scratch and phase',shots=n,Q=circ.num_qubits,T=circ.operation_statistics().num_toffoli_like)
    Path('results/inverse_product_kernel_check.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
if __name__=='__main__':test()
