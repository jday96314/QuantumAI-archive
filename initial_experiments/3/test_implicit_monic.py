"""Check the sweep with its implicit unit coefficient and retained scale."""
import numpy as np
import kickmix as km
from fast_dqi import FastBuilder
from stream_forney import chien,GF

rng=np.random.default_rng(20261023);N=256;M=71
b=FastBuilder(km.CircuitBuilder())
sigma=[b.create_quantum_register(12,f's{i}') for i in range(36)]
omega=[b.create_quantum_register(12,f'o{i}') for i in range(35)]
masks=[b.create_classical_register(12,f'm{i}') for i in range(M)]
chien(b,sigma,omega,masks);c=b.finish_circuit()
ss=np.zeros((36,N),dtype=np.int64);ss[0]=rng.integers(1,4096,N)
oo=rng.integers(0,4096,(35,N));mm=rng.integers(0,4096,(M,N))
for shot in range(N):
    poly=[1]
    for j in rng.choice(M,size=1+shot%35,replace=False):
        root=GF.pow(GF.primitive_element,int(j));iv=GF.invert(root)
        nxt=poly+[0]
        for i,v in enumerate(poly):nxt[i+1]^=GF.mul(int(v),iv)
        poly=nxt
    ss[1:len(poly),shot]=poly[1:]
expected=np.zeros(N,dtype=bool)
for j in range(M):
    point=GF.pow(GF.primitive_element,j)
    powers=[GF.pow(point,i) for i in range(36)]
    for shot in range(N):
        sv=1;dv=0;ov=0
        for i in range(1,36):
            v=GF.mul(int(ss[i,shot]),powers[i]);sv^=v
            if i%2:dv^=v
        if sv:continue
        for i in range(35):ov^=GF.mul(int(oo[i,shot]),powers[i])
        value=GF.mul(GF.mul(ov,point),GF.pow(dv,4094))
        expected[shot]^=bool((value&int(mm[j,shot])).bit_count()%2)
bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool)
sim=km.Simulator(batch_size=N,seed=103);sim.use_same_registers_as(c)
for regs,vals in [(sigma,ss),(omega,oo),(masks,mm)]:
    for reg,val in zip(regs,vals):sim.write_across_shots(reg,bits(val))
sim.do(c)
assert np.array_equal(sim.read_phase_flipped_across_shots(),expected)
for regs,vals in [(sigma,ss),(omega,oo),(masks,mm)]:
    for reg,val in zip(regs,vals):assert np.array_equal(sim.read_across_shots(reg),bits(val))
public=set(q for r in sigma+omega for q in r)
assert not sim.read_across_shots(km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])).any()
print('PASS implicit monic sweep:',N,'shots,',M,'positions,',c.num_qubits,'qubits',flush=True)
