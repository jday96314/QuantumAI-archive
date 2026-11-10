"""Independent dense-table arbitrary-target XOR controls, not submissions."""
import numpy as np
import kickmix as km
from core import *
from resources import build,bound

def bits(v,n):return ((v[None,:]>>np.arange(n)[:,None])&1).astype(bool)

def main():
    rows=[];rng=np.random.default_rng(924152)
    for A,brot,J in [(796,7,8),(796,8,7),(972,8,7),(796,9,6)]:
        delta=rng.integers(0,1<<(brot-1),(A,J));delta|=(delta>>(brot-2)&1)<<(brot-1)
        bnd=bound(A,brot,J);c=build(delta,brot,bnd['lowbits'])
        assert c.num_qubits<=128
        for seed in [134932,716582]:
            addresses=np.repeat(np.arange(A),64);N=len(addresses)
            sim=km.Simulator(batch_size=N,seed=seed);sim.use_same_registers_as(c)
            public=set();values={}
            for name,reg in c.register_data.items():
                public.update(reg);v=addresses if name=='selection' else rng.integers(0,1<<brot,N)
                values[name]=v;sim.write_across_shots(reg,bits(v,len(reg)))
            sim.do(c)
            for name,reg in c.register_data.items():
                want=values[name] if name=='selection' else values[name]^delta[addresses,int(name.split('_')[1])]
                assert np.array_equal(sim.read_across_shots(reg),bits(want,len(reg)))
            scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
            assert not sim.read_across_shots(scratch).any() and not sim.read_phase_flipped_across_shots().any()
            sim.do(c)
            for name,reg in c.register_data.items():assert np.array_equal(sim.read_across_shots(reg),bits(values[name],len(reg)))
            assert not sim.read_across_shots(scratch).any() and not sim.read_phase_flipped_across_shots().any()
            row=dict(status='PASS',A=A,brot=brot,J=J,seed=seed,shots_per_address=64,shots=N,Q=c.num_qubits,static_T_per_pass=c.operation_statistics().num_toffoli_like)
            rows.append(row);emit(row)
    report('oracle_controls.json',dict(status='PASS',scope='Synthetic dense-table single-oracle XOR, phase and workspace controls on arbitrary incoming targets, forward and repeated application; not a chemistry submission or physical fit.',rows=rows))

if __name__=='__main__':main()
