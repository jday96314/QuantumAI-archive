"""Exhaustive reduced-width test of signed row phasing and exponent routing."""
import hashlib,json
import numpy as np
import kickmix as km
from invsqrt.verify_submission_invsqrt import half_turns_to_units
from direct_rows import apply
from model import HERE,product

def circuit(ns,L,e):
    b=km.CircuitBuilder();a=b.create_quantum_register(ns,'a');z=b.create_quantum_register(L,'z')
    flags=b.create_quantum_register(5,'flags');apply(b,a,z,flags,e)
    return b.finish_circuit(),a,z,flags

def main():
    c,a,z,flags=circuit(4,6,-3)
    values=np.array([(aa,zz,ss) for aa in range(16) for zz in range(64) for ss in range(11)],dtype=np.int64)
    expected=(-product(values[:,0],values[:,1],6,-3)<<values[:,2]).astype(np.uint64)<<np.uint64(27)
    public=set(a)|set(z)|set(flags);scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
    bits=lambda v,n:((v[None,:]>>np.arange(n)[:,None])&1).astype(bool)
    for salt in ('row-exhaust-v1','row-exhaust-v2','row-exhaust-v3','row-exhaust-v4'):
        sim=km.Simulator(batch_size=256,seed=int(hashlib.sha256(salt.encode()).hexdigest(),16));sim.use_same_registers_as(c)
        for first in range(0,len(values),256):
            aa,zz,ss=values[first:first+256].T;sim.clear_for_shot()
            for reg,v in ((a,aa),(z,zz),(flags,2*ss)):sim.write_across_shots(reg,bits(v,len(reg)))
            sim.do(c)
            for reg,v in ((a,aa),(z,zz),(flags,2*ss)):assert np.array_equal(sim.read_across_shots(reg),bits(v,len(reg)))
            assert not sim.read_across_shots(scratch).any()
            actual=np.array([half_turns_to_units(sim.read_shot_phase(j)) for j in range(256)],dtype=np.uint64)
            assert np.array_equal(actual,expected[first:first+256]),(salt,first)
    full,*_=circuit(17,14,-12)
    stats=full.operation_statistics()
    report=dict(status='PASS',scope='every reduced-width source, signed coordinate and exponent 0..10; four independent measurement seeds',
        cases=len(values),shots=4*len(values),small_Q=c.num_qubits,
        full_width_component=dict(Q=full.num_qubits,Toffolis=stats.num_toffoli_like,
            rotations=stats.num_z_pow(min_clifford_hierarchy_level=4)))
    full.to_file(str(HERE/'candidates/direct_rows_component.kmb'),format='kmb')
    (HERE/'results/components.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))

if __name__=='__main__':main()
