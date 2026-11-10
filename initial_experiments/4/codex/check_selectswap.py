"""Check XOR-oracle behavior on random incoming words, beyond the host stream."""
import io,json,zipfile
import numpy as np
import kickmix as km
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas
from build_selectswap_chem import build
rng=np.random.default_rng(90172)
with zipfile.ZipFile('baseline/submission.zip') as archive:
    with np.load(io.BytesIO(archive.read('dfthcsos/submission.npz'))) as z:words=z['angle_words']
N=1024;addr=np.arange(N)%len(words);records=[]
for k,delta in enumerate(batch_deltas(words,7)):
    c=build(delta,8,monomial=True);sim=km.Simulator(batch_size=N,seed=983+k);sim.use_same_registers_as(c)
    public=set();allq=[]
    for name in ['selection']+[f'word_{j}' for j in range(7)]:
        q=c.register_data[name];public.update(q)
        if name=='selection':values=addr
        else:values=rng.integers(0,256,N)
        bits=((values[None,:]>>np.arange(len(q))[:,None])&1).astype(bool)
        sim.write_across_shots(q,bits);allq.append((q,values))
    sim.do(c)
    for j,(q,values) in enumerate(allq):
        expected=values if j==0 else values^delta[addr,j-1]
        bits=((expected[None,:]>>np.arange(len(q))[:,None])&1).astype(bool)
        assert np.array_equal(sim.read_across_shots(q),bits),(k,j)
    scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
    assert not sim.read_across_shots(scratch).any(),k
    assert not sim.read_phase_flipped_across_shots().any(),k
    records.append(dict(k=k,status='PASS',shots=N))
print(json.dumps(records,indent=2))
