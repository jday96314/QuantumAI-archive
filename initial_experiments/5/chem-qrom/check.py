"""Exhaustive-address arbitrary-target oracle and complete-family tests."""
import argparse,time
import kickmix as km
from common import *
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas
from dfthcsos.scoreinfra.circuits import load_circuits,recover_words_from_circuits

def component(c,d,seed=2026102405):
    A,J=d.shape;addresses=np.arange(A);rng=np.random.default_rng(seed)
    sim=km.Simulator(batch_size=A,seed=seed);sim.use_same_registers_as(c);public=set();values={}
    bits=lambda v,n:((v[None,:]>>np.arange(n)[:,None])&1).astype(bool)
    for name,reg in c.register_data.items():
        public.update(reg);v=addresses if name=='selection' else rng.integers(0,256,A)
        values[name]=v;sim.write_across_shots(reg,bits(v,len(reg)))
    scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
    for reverse in [False,True]:
        sim.do(c)
        for name,reg in c.register_data.items():
            want=values[name] if name=='selection' or reverse else values[name]^d[addresses,int(name[5:])]
            assert np.array_equal(sim.read_across_shots(reg),bits(want,len(reg))),(name,reverse)
        assert not sim.read_across_shots(scratch).any()
        assert not sim.read_phase_flipped_across_shots().any()
    return dict(status='PASS',shots_each_direction=A,seed=seed)

def check(path,repeats=4,seed=2026102401):
    path=Path(path);words=np.load(path/'submission.npz')['angle_words'];circuits=load_circuits(path)
    J=len([name for name in circuits[0].register_data if name.startswith('word_')])
    rng=np.random.default_rng(seed);N=repeats*len(words);addresses=np.tile(np.arange(len(words)),repeats);rows=[]
    bits=lambda v,n:((v[None,:]>>np.arange(n)[:,None])&1).astype(bool)
    for k,(c,d) in enumerate(zip(circuits,batch_deltas(words,J))):
        sim=km.Simulator(batch_size=N,seed=seed+k);sim.use_same_registers_as(c);public=set();values={}
        for name,reg in c.register_data.items():
            public.update(reg);v=addresses if name=='selection' else rng.integers(0,256,N)
            values[name]=v;sim.write_across_shots(reg,bits(v,len(reg)))
        scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
        for reverse in [False,True]:
            sim.do(c)
            for name,reg in c.register_data.items():
                want=values[name] if name=='selection' or reverse else values[name]^d[addresses,int(name[5:])]
                assert np.array_equal(sim.read_across_shots(reg),bits(want,len(reg))),(k,name,reverse)
            assert not sim.read_across_shots(scratch).any(),(k,reverse,'scratch')
            assert not sim.read_phase_flipped_across_shots().any(),(k,reverse,'phase')
        rows.append(dict(k=k,status='PASS',Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,shots_each_direction=N,sha256=sha(path/f'angle_k_{k}.kmb')))
    recovered,layout,t=recover_words_from_circuits(circuits,num_vectors=len(words),num_angles=words.shape[1],shots_per_address=repeats,seed=seed+1000)
    assert np.array_equal(recovered,words)
    result=dict(status='PASS',seed=seed,rows=rows,family_mean_toffolis=t,batch=J,
                exact_npz_sha256=sha(path/'submission.npz'),predicted=predicted(J,rows))
    save(path.name+'_check.json',result);print(json.dumps(result),flush=True)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('path',type=Path);p.add_argument('--repeats',type=int,default=4);p.add_argument('--seed',type=int,default=2026102401);a=p.parse_args();check(a.path,a.repeats,a.seed)
