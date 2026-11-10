"""All-address arbitrary-target forward/reverse tests, including 11-bit padding."""
import json
import kickmix as km
from core import *
from oracles import build
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas
rng=np.random.default_rng(359211);records=[]
for source in [INC/'submission.npz',HERE/'candidates/qd_b479_C175.npz']:
    z,u,w,tree,words=load(source)
    words=encode(reconstruct_u_vectors(words/128,tree),tree,8,True)
    if len(words)==u.shape[0]*u.shape[1]:words=np.concatenate([np.load(INC/'submission.npz')['angle_words'][:76],words])
    for k,delta in enumerate(batch_deltas(words,7)):
        c=build(delta);assert c.num_qubits<=128
        shots=2*len(words);addresses=np.arange(shots)%len(words);sim=km.Simulator(batch_size=shots,seed=359211+k);sim.use_same_registers_as(c)
        public=set();values={}
        for name,reg in c.register_data.items():
            public.update(reg);v=addresses if name=='selection' else rng.integers(0,256,shots);values[name]=v
            sim.write_across_shots(reg,((v[None,:]>>np.arange(len(reg))[:,None])&1).astype(bool))
        scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
        for reverse in [False,True]:
            sim.do(c)
            for name,reg in c.register_data.items():
                want=values[name] if name=='selection' or reverse else values[name]^delta[addresses,int(name.split('_')[-1])]
                assert np.array_equal(sim.read_across_shots(reg),((want[None,:]>>np.arange(len(reg))[:,None])&1).astype(bool)),(k,name,reverse)
            assert not sim.read_across_shots(scratch).any()
            assert not sim.read_phase_flipped_across_shots().any()
        records.append(dict(source=str(source),k=k,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,shots_per_direction=shots,status='PASS'))
    print(json.dumps(records[-11:]),flush=True)
(HERE/'results/oracle_forward_reverse.json').write_text(json.dumps(records,indent=2)+'\n')
