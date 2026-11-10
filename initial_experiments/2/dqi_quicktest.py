import sys,time,numpy as np,kickmix as km
from dqi.verify_submission_dqi import ProblemSpecDQI,FieldTables,sample_shots,verify
import importlib;generic_dqi=importlib.import_module(sys.argv[3] if len(sys.argv)>3 else "generic_dqi")
t0=time.time();c=generic_dqi.build_circuit();print('built',c.num_qubits,c.operation_statistics().num_toffoli_like,'%.0fs'%(time.time()-t0),flush=True)
if len(sys.argv)>1:c.to_file(sys.argv[1],format='kmb')
spec=ProblemSpecDQI.benchmark_instance(num_shots=int(sys.argv[2]) if len(sys.argv)>2 else 512)
rng=np.random.default_rng(7);tables=FieldTables(km.GF2Field(12))
shots=sample_shots(spec,tables,rng,spec.num_shots)
r=verify(c,spec,shots=shots,max_workers=14)
print({k:r[k] for k in ('value_errors','phase_errors','toffoli_count','num_qubits','score','repetition_factor')},'%.0fs'%(time.time()-t0))
