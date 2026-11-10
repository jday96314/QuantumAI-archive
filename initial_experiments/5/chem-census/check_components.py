"""All-address, arbitrary-target reversible XOR/clean-scratch/phase checks."""
from census_common import *
import kickmix as km
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas
from dfthcsos.scoreinfra.circuits import load_circuits
path=HERE/'candidates/structured';words=np.load(path/'submission.npz')['angle_words'];cs=load_circuits(path)
rows=[];A=len(words)
for seed in [2026102406,2026102407]:
    rng=np.random.default_rng(seed);addresses=np.tile(np.arange(A),2);N=len(addresses)
    bits=lambda v,n:((v[None,:]>>np.arange(n)[:,None])&1).astype(bool)
    for k,(c,d) in enumerate(zip(cs,batch_deltas(words,4))):
        sim=km.Simulator(batch_size=N,seed=seed+k);sim.use_same_registers_as(c);values={};public=set()
        for name,reg in c.register_data.items():
            public.update(reg)
            values[name]=addresses if name=='selection' else rng.integers(0,256,N)
            sim.write_across_shots(reg,bits(values[name],len(reg)))
        scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
        for reverse in [False,True]:
            sim.do(c)
            for name,reg in c.register_data.items():
                want=values[name] if name=='selection' or reverse else values[name]^d[addresses,int(name[5:])]
                assert np.array_equal(sim.read_across_shots(reg),bits(want,len(reg))),(seed,k,name,reverse)
            assert not sim.read_across_shots(scratch).any(),(seed,k,'scratch',reverse)
            assert not sim.read_phase_flipped_across_shots().any(),(seed,k,'phase',reverse)
        rows.append(dict(seed=seed,k=k,shots_per_direction=N,status='PASS'))
        print(json.dumps(rows[-1]),flush=True)
write('structured_components.json',dict(status='PASS',rows=rows,arbitrary_incoming_words=True,all_valid_addresses=True,
    forward_and_reverse=True,clean_scratch=True,no_phase=True))
