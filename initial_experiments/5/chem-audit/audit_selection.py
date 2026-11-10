"""Select only a full-default native-passing chemistry task-cost improvement."""
import argparse,hashlib,json,math
from audit_core import *
from dfthcsos.scoreinfra.circuits import submission_sha256,load_circuits
p=argparse.ArgumentParser();p.add_argument('--candidate',required=True);p.add_argument('--result',required=True);a=p.parse_args()
path=Path(a.candidate).resolve();result=json.loads(Path(a.result).read_text());z=dict(np.load(path/'submission.npz'))
family=load_circuits(path);actual=submission_sha256(family,z,z['tree_edges'])
assert result['status']=='PASS'
assert result['submission_sha256']==actual
assert result['verified_shots']==(76+np.prod(z['W'].shape[:2]))*64
assert passes_accuracy(result['metrics'])
assert math.isclose(result['score'],result['lambda_eff']*result['toffoli_count']*result['num_qubits'],rel_tol=1e-13)
assert result['score']<708547388.9553173
if (HERE/'selection.json').exists():
    current=json.loads((HERE/'selection.json').read_text())
    assert result['score']<=current['score']
    if result['score']==current['score']:assert actual==current['canonical_digest']
rec=dict(status='SELECTED_NATIVE_COST_IMPROVEMENT',candidate=str(path),native_result=str(Path(a.result).resolve()),
    score=result['score'],previous_score=708547388.9553173,improvement=708547388.9553173/result['score'],
    combined_reference_rescaled=9.3677221531*(708547388.9553173/result['score'])**(1/3),
    inverse=182956.5272,dqi=1598204229385,canonical_digest=actual,verified_shots=result['verified_shots'],
    metrics=result['metrics'],resources=result['resources'],
    artifact_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in path.iterdir() if p.is_file()})
write_json(HERE/'selection.json',rec);print(json.dumps(rec),flush=True)
