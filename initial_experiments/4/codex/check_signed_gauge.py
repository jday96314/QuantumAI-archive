"""Check every sign-gauge XOR oracle on arbitrary incoming public words."""
import argparse,json
from pathlib import Path
import numpy as np
import kickmix as km
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas

p=argparse.ArgumentParser();p.add_argument('--candidate',default='candidates/signed_gauge7');a=p.parse_args()
rng=np.random.default_rng(983741);root=Path(a.candidate);words=np.load(root/'submission.npz')['angle_words'];N=2048;addresses=np.arange(N)%len(words);records=[]
for k,delta in enumerate(batch_deltas(words,7)):
    circuit=km.Circuit.from_file(str(root/f'angle_k_{k}.kmb'),format='kmb');sim=km.Simulator(batch_size=N,seed=917+k);sim.use_same_registers_as(circuit)
    public=set();values={}
    for name,reg in circuit.register_data.items():
        public.update(reg);v=addresses if name=='selection' else rng.integers(0,256,N);values[name]=v
        sim.write_across_shots(reg,((v[None,:]>>np.arange(len(reg))[:,None])&1).astype(bool))
    sim.do(circuit)
    for name,reg in circuit.register_data.items():
        want=values[name] if name=='selection' else values[name]^delta[addresses,int(name.split('_')[-1])]
        assert np.array_equal(sim.read_across_shots(reg),((want[None,:]>>np.arange(len(reg))[:,None])&1).astype(bool)),(k,name)
    scratch=km.array([km.q(i) for i in range(circuit.num_qubits) if km.q(i) not in public])
    assert not sim.read_across_shots(scratch).any(),(k,'scratch');assert not sim.read_phase_flipped_across_shots().any(),(k,'phase')
    records.append(dict(k=k,status='PASS',shots=N,Q=circuit.num_qubits,T=circuit.operation_statistics().num_toffoli_like))
Path(f'results/{root.name}_random_targets.json').write_text(json.dumps(records,indent=2)+'\n');print(json.dumps(records),flush=True)
