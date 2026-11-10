"""All-address arbitrary incoming target XOR and reverse tests."""
import argparse,json
import kickmix as km
from core import *
from oracles import build
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas
p=argparse.ArgumentParser();p.add_argument('fits',nargs='+');a=p.parse_args();rng=np.random.default_rng(915342);records=[]
for source in a.fits:
    z,u,w,tree,words=load(source);R,B,C=w.shape
    allwords=z.get('angle_words')
    if allwords is None or len(allwords)!=76+R*B:
        one=z.get('onebody_angle_words',np.load(INC/'submission.npz')['angle_words'][:76]);allwords=np.concatenate([one,words])
    allwords=encode(reconstruct_u_vectors(allwords/128,tree),tree,8,True)
    for k,delta in enumerate(batch_deltas(allwords,7)):
        c=build(delta);assert c.num_qubits<=128;shots=2*len(allwords);addresses=np.arange(shots)%len(allwords)
        sim=km.Simulator(batch_size=shots,seed=915342+k);sim.use_same_registers_as(c);public=set();values={}
        for name,reg in c.register_data.items():
            public.update(reg);v=addresses if name=='selection' else rng.integers(0,256,shots);values[name]=v
            sim.write_across_shots(reg,((v[None,:]>>np.arange(len(reg))[:,None])&1).astype(bool))
        scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
        for reverse in [False,True]:
            sim.do(c)
            for name,reg in c.register_data.items():
                want=values[name] if name=='selection' or reverse else values[name]^delta[addresses,int(name.split('_')[-1])]
                assert np.array_equal(sim.read_across_shots(reg),((want[None,:]>>np.arange(len(reg))[:,None])&1).astype(bool)),(k,name,reverse)
            assert not sim.read_across_shots(scratch).any();assert not sim.read_phase_flipped_across_shots().any()
        row=dict(source=source,k=k,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,shots_per_direction=shots,status='PASS')
        records.append(row)
    print(json.dumps(records[-11:]),flush=True)
(HERE/'results/oracle_forward_reverse.json').write_text(json.dumps(records,indent=2))
