"""Independent native phase-grid decode and arbitrary-target ANGLE tests."""
import argparse
import kickmix as km
from core import *
from phase_grid import phase_encode
from oracles import build
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas

def check(c,delta,seed):
    rng=np.random.default_rng(seed);shots=2*len(delta);addresses=np.arange(shots)%len(delta)
    sim=km.Simulator(batch_size=shots,seed=seed);sim.use_same_registers_as(c);public=set();values={}
    for name,reg in c.register_data.items():
        public.update(reg);v=addresses if name=='selection' else rng.integers(0,1<<len(reg),shots)
        values[name]=v;sim.write_across_shots(reg,((v[None,:]>>np.arange(len(reg))[:,None])&1).astype(bool))
    scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
    for reverse in [False,True]:
        sim.do(c)
        for name,reg in c.register_data.items():
            want=values[name] if name=='selection' or reverse else values[name]^delta[addresses,int(name.split('_')[-1])]
            assert np.array_equal(sim.read_across_shots(reg),((want[None,:]>>np.arange(len(reg))[:,None])&1).astype(bool))
        assert not sim.read_across_shots(scratch).any() and not sim.read_phase_flipped_across_shots().any()
    return dict(Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,shots_per_direction=shots,status='PASS',seed=seed)

def main():
    p=argparse.ArgumentParser();p.add_argument('--submission',nargs='*');a=p.parse_args();records=[]
    if a.submission:
        for path in a.submission:
            z=np.load(Path(path)/'submission.npz');words=z['angle_words'];bits=int(z.get('brot',8))
            from dfthcsos.scoreinfra.circuits import load_circuits
            for k,(c,delta) in enumerate(zip(load_circuits(path),batch_deltas(words,7 if bits==8 else 6))):
                assert c.num_qubits<=128
                records.append(dict(submission=path,k=k,**check(c,delta,743391+k)))
        report(HERE/'results/actual_oracle_components.json',records)
    else:
        rng=np.random.default_rng(3331);u=rng.normal(size=(100,N));u/=np.linalg.norm(u,axis=-1,keepdims=True)
        tr=ordered(np.r_[0,rng.permutation(np.arange(1,N))]);phase=rng.integers(0,2,75)
        words,v=phase_encode(u,tr,phase);host=reconstruct_u_vectors(words/256,tr)
        assert np.max(abs(v-host))<1e-14
        assert np.all((words&1)==phase[None,:])
        w0,v0=phase_encode(u,tr,np.zeros(75,dtype=np.int64));w8,v8=quantize(u,tr)
        assert np.array_equal(w0,2*w8) and np.max(abs(v0-v8))<1e-14
        for A,bits in [(778,9),(1034,8),(1034,9)]:
            J=6 if bits>8 else 7
            coarse=rng.integers(-64,64,(A,J));delta=(coarse*(1<<(bits-8)))%(1<<bits)
            if bits>8:delta+=rng.integers(0,2,J)[None,:]
            c=build(delta,bits,bits-8);assert c.num_qubits<=128
            for seed in [743391,999841]:records.append(dict(A=A,bits=bits,J=J,**check(c,delta,seed)))
        report(HERE/'results/grid_oracle_controls.json',dict(native_decoder_max=float(np.max(abs(v-host))),
            zero_phase_exact_word_match=True,tests=records))

if __name__=='__main__':main()
