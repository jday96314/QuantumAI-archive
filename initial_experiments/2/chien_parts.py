import kickmix as km, numpy as np
from fast_dqi import FastBuilder
gf=km.GF2Field(12)
def cost(fn,shots=512):
    b=FastBuilder(km.CircuitBuilder())
    s=b.create_quantum_register(12,'s');sp=b.create_quantum_register(12,'sp');om=b.create_quantum_register(12,'om')
    M=b.create_classical_register(12,'M')
    fn(b,s,sp,om,M);c=b.finish_circuit()
    sim=km.Simulator(batch_size=shots,seed=3,counting_operations=True);sim.use_same_registers_as(c)
    rng=np.random.default_rng(0)
    for r in (s,sp,om,M):
        v=rng.integers(1,4096,shots); 
        if r is s: v[::2]=0
        sim.write_across_shots(r,((v[None,:]>>np.arange(12)[:,None])&1)==1)
    sim.do(c);return sim.operation_statistics().num_toffoli_like/shots
def inv(b,s,sp,om,M):
    t=b.alloc_qubits(12);_,ch=b.init_gf2_inverse_with_scaffold(sp,target=t,field=gf)
    b.del_gf2_inverse_with_scaffold(sp,target=t,scaffold=ch,field=gf);b.free(ch);b.free(t)
def initinv(b,s,sp,om,M):
    t=b.alloc_qubits(12);_,ch=b.init_gf2_inverse_with_scaffold(sp,target=t,field=gf)
def zero(b,s,sp,om,M):
    z=b.alloc_qubits(1)[0];b.flip_if_equal(s,0,target=z);b.flip_if_equal(s,0,target=z);b.free(z)
def phase(b,s,sp,om,M):
    z=b.alloc_qubits(1)[0];loc=b.alloc_qubits(12);b.x(z)
    b.init_and(z,om,loc);b.gf2_phase_by_product(loc,sp,mask=M,field=gf);b.del_and(z,om,loc);b.x(z);b.free(loc);b.free(z)
for f in (inv,initinv,zero,phase):print(f.__name__,cost(f))
