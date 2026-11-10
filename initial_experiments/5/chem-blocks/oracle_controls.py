"""Native arbitrary-target/phase/clean-scratch controls, independent RNG seeds."""
import numpy as np
import kickmix as km
from signed_oracle import bound,build
from common import HERE,report

def check(c,delta,bits,seed,shots=2048):
    rng=np.random.default_rng(seed);A=len(delta);addr=np.arange(shots)%A
    sim=km.Simulator(batch_size=shots,seed=seed^647912);sim.use_same_registers_as(c)
    public=set();values={}
    def binary(v,n):return ((v[None,:]>>np.arange(n)[:,None])&1).astype(bool)
    for name,reg in c.register_data.items():
        public.update(reg);v=addr if name=='selection' else rng.integers(0,1<<bits,shots)
        values[name]=v;sim.write_across_shots(reg,binary(v,len(reg)))
    sim.do(c)
    for name,reg in c.register_data.items():
        target=values[name] if name=='selection' else values[name]^delta[addr,int(name.split('_')[-1])]
        assert np.array_equal(sim.read_across_shots(reg),binary(target,len(reg))),(seed,name)
    scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
    assert not sim.read_across_shots(scratch).any()
    assert not sim.read_phase_flipped_across_shots().any()
    return dict(status='PASS',shots=shots,Q=c.num_qubits,T=int(c.operation_statistics().num_toffoli_like),seed=seed)

def main():
    records=[];rng=np.random.default_rng(748192)
    for A,bits,J in [(795,8,7),(796,7,8),(778,8,7),(1100,9,6),(844,7,8)]:
        ob=bound(A,bits,J);d=rng.integers(0,1<<(bits-1),(A,J),dtype=np.int64)
        d|=((d>>(bits-2))&1)<<(bits-1)
        c=build(d,bits,ob['lowbits'])
        r=check(c,d,bits,917123+bits+A)
        assert r['T']<=ob['T_per_oracle'] and r['Q']<=ob['Q_oracle']
        records.append(dict(A=A,brot=bits,J=J,**r,bound=ob))
    report(HERE/'results/oracle_controls.json',dict(status='PASS',scope='random synthetic XOR tables, not a physical candidate',controls=records))

if __name__=='__main__':main()
