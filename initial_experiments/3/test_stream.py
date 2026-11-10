"""Check restoration and the complete locator's phase parity."""
import json,sys
import numpy as np
import kickmix as km
from dqi.verify_submission_dqi import ProblemSpecDQI,FieldTables,sample_shots,to_bits
from stream_model import decode,pub
from published_dqi import isomorphism
path=sys.argv[1];sigma_test='sigma' in path
c=km.Circuit.from_file(path);rng=np.random.default_rng(231);N=1024
spec=ProblemSpecDQI.benchmark_instance(m=71,n=70,num_shots=N)
syn,mask,ph=sample_shots(spec,FieldTables(km.GF2Field(12)),rng,N)
sim=km.Simulator(batch_size=N,seed=78);sim.use_same_registers_as(c)
qs=km.array([km.q(i) for i in range(840)]);inputbits=to_bits(syn,12)
sim.write_across_shots(qs,inputbits);sim.do(c)
values=int(np.count_nonzero(np.any(sim.read_across_shots(qs)^inputbits,axis=0)))
junk=int(np.count_nonzero(np.any(sim.read_across_shots(km.array([km.q(i) for i in range(840,c.num_qubits)])),axis=0)))
phases=sim.read_phase_flipped_across_shots()
if sigma_test:
    f,back,tr=isomorphism(pub)
    expected=[]
    for row in syn:
        sig,why=decode([tr(int(x),f) for x in row])
        expected.append(sum(x.bit_count() for x in sig)%2==1)
    errors=int(np.count_nonzero(phases^expected))
else:errors=int(np.count_nonzero(phases))
print(json.dumps(dict(value_errors=values,junk_errors=junk,phase_errors=errors,qubits=c.num_qubits,toffoli=c.operation_statistics().num_toffoli_like)),flush=True)
assert values==0 and junk==0
assert errors<=3
