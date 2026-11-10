"""Independently reopen complete ZIP and bind full host-metric measurements."""
from pathlib import Path
import hashlib,json,io,math,zipfile
import numpy as np
import kickmix as km
from dfthcsos.scoreinfra.circuits import submission_sha256

HERE=Path(__file__).resolve().parent
def sha_file(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8<<20),b''):h.update(b)
    return h.hexdigest()
assembly=json.loads((HERE/'results/package_assembly.json').read_text())
metric=json.loads((HERE/'results/package_metric.json').read_text())
summary=json.loads((HERE/'submission.summary.json').read_text())
archive=HERE/'submission.zip';digest=sha_file(archive)
assert metric['status']=='PASS' and summary['status']=='PACKAGE_FULL_METRIC_PASS'
assert digest==assembly['sha256']==metric['archive_sha256']==summary['sha256']
assert metric['configuration']['verification_shot_divisor']==1
assert metric['configuration']['process_cpu_count']==2
native={k:json.loads(Path(p).read_text()) for k,p in metric['native_reports'].items()}
expected=dict(invsqrt=200540928,dqi=102400,dfthcsos=50880)
for task,r in native.items():
    assert r['status']=='PASS' and r['verified_shots']==expected[task]
    assert math.isclose(r['score'],metric['task_details'][task]['cost'],rel_tol=1e-13)
with zipfile.ZipFile(archive) as z:
    assert len(z.namelist())==14 and set(z.namelist())==set(assembly['members'])
    member_hashes={}
    for name in z.namelist():
        h=hashlib.sha256()
        with z.open(name) as f:
            for chunk in iter(lambda:f.read(8<<20),b''):h.update(chunk)
        member_hashes[name]=h.hexdigest();assert h.hexdigest()==assembly['members'][name]
    assert member_hashes['invsqrt.kmb']==native['invsqrt']['circuit_sha256']
    assert member_hashes['dqi.kmb']==native['dqi']['circuit_sha256']
    circuits=[km.Circuit.from_file(io.BytesIO(z.read(f'dfthcsos/angle_k_{k}.kmb')),format='kmb') for k in range(11)]
    with np.load(io.BytesIO(z.read('dfthcsos/submission.npz'))) as q:
        chem_digest=submission_sha256(circuits,{k:q[k] for k in ['U','W','H2_I','H1_I','Hsym_lower']},q['tree_edges'])
    assert chem_digest==native['dfthcsos']['submission_sha256']==assembly['canonical_chemistry_digest']
actual=math.prod(v['improvement_over_baseline'] for v in metric['task_details'].values())**(1/3)
assert math.isclose(actual,metric['local_combined_score'],rel_tol=1e-14)
assert summary['overall_score']==metric['local_combined_score']
summary['scope']='Complete archive passed unmodified host metric readers, full-default native verification and metric.score; independently reopened and bound to those reports.'
(HERE/'submission.summary.json').write_text(json.dumps(summary,indent=2)+'\n')
record=dict(status='PASS',archive=str(archive),archive_sha256=digest,archive_bytes=archive.stat().st_size,members=member_hashes,
    chemistry_digest=chem_digest,verified_shots=expected,local_combined_score=actual,accepted_reference_score=assembly['accepted_reference_score'],
    selected_dqi=assembly['dqi_decision']['selected'],dqi_default_and_independent_evidence=assembly['dqi_decision']['proofs'],
    scope='ZIP read to EOF/CRC checked and SHA256 matched for all 14 members; native reports and metric.score bound without reusing estimated task costs.')
(HERE/'results/package_independent_audit.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({k:v for k,v in record.items() if k!='members'},indent=2))
