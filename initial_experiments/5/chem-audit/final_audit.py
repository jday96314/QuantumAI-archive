"""Bind selected measured results, preserved incumbent and research artifacts."""
import hashlib,json,math
from audit_core import *
from dfthcsos.scoreinfra.circuits import load_circuits,submission_sha256

selected=json.loads((HERE/'selection.json').read_text());path=Path(selected['candidate'])
z=dict(np.load(path/'submission.npz'));circuits=load_circuits(path)
assert submission_sha256(circuits,z,z['tree_edges'])==selected['canonical_digest']
for name,digest in selected['artifact_sha256'].items():
    assert hashlib.sha256((path/name).read_bytes()).hexdigest()==digest
base=json.loads((HERE/'results/calibration.json').read_text())
for name,digest in base['source_hashes'].items():
    assert hashlib.sha256((INC/name).read_bytes()).hexdigest()==digest
native=json.loads(Path(selected['native_result']).read_text())
assert native['status']=='PASS' and native['score']==selected['score']
assert native['verified_shots']==(76+np.prod(z['W'].shape[:2]))*64
assert passes_accuracy(native['metrics'])
assert math.isclose(native['score'],native['lambda_eff']*native['toffoli_count']*native['num_qubits'],rel_tol=1e-13)
independent=json.loads((HERE/'results/banked_native_independent.json').read_text())
assert independent['status']=='PASS' and independent['verified_shots']==50880
assert independent['submission_sha256']==json.loads((HERE/'results/banked_native.json').read_text())['submission_sha256']
package=json.loads((HERE/'results/package_independent_audit.json').read_text())
metric=json.loads((HERE/'results/package_metric.json').read_text())
summary=json.loads((HERE/'submission.summary.json').read_text())
assert package['status']==metric['status']=='PASS'
archive_digest=hashlib.sha256((HERE/'submission.zip').read_bytes()).hexdigest()
assert archive_digest==metric['archive_sha256']==package['archive_sha256']
assert summary['overall_score']==metric['local_combined_score']==package['local_combined_score']
assert metric['verified_shots']==dict(invsqrt=200540928,dqi=102400,dfthcsos=50880)
reused=[ROOT/'from_host/verifiers/LICENSE.txt',ROOT/'from_host/verifiers/dfthcsos/verify.py',
    ROOT/'from_host/verifiers/dfthcsos/scoreinfra/dfthc_math.py',ROOT/'from_host/verifiers/dfthcsos/scoreinfra/costs_prepare.py',
    ROOT/'from_host/verifiers/dfthcsos/scoreinfra/reference.py',ROOT/'from_host/verifiers/dfthcsos/scoreinfra/circuits.py',
    ROOT/'from_host/verifiers/dfthcsos/FeMoCo76o113e.fixed_shifts.npz']
reused += [ROOT/'initial_experiments/5/chem-qfit'/name for name in ['core.py','physical.py','discrete.py','oracles.py']]
reused += [ROOT/'initial_experiments/4/codex'/name for name in ['angle_gauge_probe.py','build_selectswap_chem.py','monomial_unlookup.py','jacobi_l1','NORMALIZATION_BOUND.md']]
reused.append(ROOT/'initial_experiments/4/agy/scripts/fast_optimize_offsets.py')
reused += [ROOT/'from_host/metric.py',ROOT/'initial_experiments/4/codex/package_submission.py',
    ROOT/'initial_experiments/4/codex/audit_submission.py']
write_json(HERE/'results/reused_sources.json',{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in reused})
manifest={}
for folder,pattern in [(HERE,'*.py'),(HERE,'*.sh'),(HERE,'RESULTS.md'),(HERE,'selection.json'),
    (HERE,'submission.zip'),(HERE,'submission.summary.json'),(HERE/'parameters','*.npz'),
    (HERE/'results','*.json'),(HERE/'results','*.jsonl'),(HERE/'results','*.txt'),
    (HERE/'results','*.log'),(HERE/'results','*.csv'),(HERE/'results','*.png'),
    (HERE/'results','*.svg'),(HERE/'results/code','*.py'),
    (HERE/'package_validation','*.csv'),(HERE/'package_validation/eval_inputs','*.csv'),
    (HERE/'package_validation/eval_inputs','*.json')]:
    if not folder.exists():continue
    for p in sorted(folder.glob(pattern)):
        if p.name in ['final_audit.json','artifact_manifest.json']:continue
        manifest[str(p.relative_to(HERE))]=hashlib.sha256(p.read_bytes()).hexdigest()
write_json(HERE/'artifact_manifest.json',manifest)
rec=dict(status='PASS',selected=selected['candidate'],score=selected['score'],canonical_digest=selected['canonical_digest'],
    preserved_round4_artifacts=len(base['source_hashes']),selected_artifacts=len(selected['artifact_sha256']),
    default_native_shots=native['verified_shots'],independent_banked_native_shots=independent['verified_shots'],
    archive_sha256=archive_digest,package_status=package['status'],
    local_combined_score=metric['local_combined_score'],full_package_shots=metric['verified_shots'],
    research_artifacts_hashed=len(manifest),manifest_sha256=hashlib.sha256((HERE/'artifact_manifest.json').read_bytes()).hexdigest())
write_json(HERE/'results/final_audit.json',rec);print(json.dumps(rec),flush=True)
