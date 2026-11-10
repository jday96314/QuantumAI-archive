"""Bind all native PASS JSONs to their delivered exact circuit/NPZ bytes."""
import json
from core import *
from dfthcsos.scoreinfra.circuits import load_circuits,submission_sha256
records=[]
for p in sorted((HERE/'results').glob('*_native.json')):
    result=json.loads(p.read_text());tag=p.name.removesuffix('_native.json')
    directory=INC if tag=='incumbent' else HERE/'candidates'/tag
    z=dict(np.load(directory/'submission.npz'));family=load_circuits(directory)
    digest=submission_sha256(family,directory/'submission.npz',z['tree_edges'])
    assert result['status']=='PASS' and digest==result['submission_sha256']
    assert max(c.num_qubits for c in family)<=128
    score=result['lambda_eff']*result['toffoli_count']*result['num_qubits'];assert abs(score-result['score'])<1e-6
    metrics=result['metrics'];assert metrics['error_1']<=.020 and metrics['error_2']<=.020 and metrics['Lambda']>49 and metrics['E_gap']>3.914898
    records.append(dict(tag=tag,directory=str(directory),status='PASS',submission_sha256=digest,score=score,verified_shots=result['verified_shots']))
(HERE/'results/native_digest_audit.json').write_text(json.dumps(records,indent=2));print(json.dumps(records),flush=True)
