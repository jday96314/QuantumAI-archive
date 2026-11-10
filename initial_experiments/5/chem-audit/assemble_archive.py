"""Assemble accepted artifacts and independently reopen/hash every ZIP member.

Layout and accepted-reference calculation follow round-4 Codex packaging;
native reader/verification/score validation is performed by package_metric.py.
"""
from pathlib import Path
import hashlib,json,math,zipfile,shutil
import numpy as np
from dfthcsos.scoreinfra.circuits import load_circuits,submission_sha256

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];R4=ROOT/'initial_experiments/4/codex'
BASE=dict(invsqrt=2076774.233715,dqi=24439463583672.929688,dfthcsos=3355599588.315537)
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(8<<20),b''):h.update(chunk)
    return h.hexdigest()
def dump(path,data):Path(path).write_text(json.dumps(data,indent=2)+'\n')

probe=json.loads((HERE/'results/subtree_probe.json').read_text())
assert probe['optimized_projected_total_saving']<probe['build_threshold_T'] and not probe['pruned_family_built']
decision_path=ROOT/'initial_experiments/5/dqi-core/RESULTS.md'
decision=decision_path.read_text()
trimmed='ACCEPTED: trimmed decoder' in decision
core=ROOT/'initial_experiments/5/dqi-core'
selection=json.loads((HERE/'selection.json').read_text());chem=Path(selection['candidate'])
assert chem==HERE/'candidates/banked'
sources=dict(invsqrt=R4/'candidates/zero_head_affine.kmb',dqi=R4/'candidates/cubic_trace_product_n70_m4095.kmb')
reports=dict(invsqrt=R4/'results/zero_head_affine_verify.json',dqi=R4/'results/cubic_trace_product_n70_m4095_verify.json',dfthcsos=Path(selection['native_result']))
proofs={}
if trimmed:
    acceptance=json.loads((core/'results/accepted_decoder.json').read_text())
    assert acceptance['status']=='ACCEPTED' and acceptance['name']=='trimmed decoder'
    sources['dqi']=core/'candidates/trimmed_decoder.kmb'
    reports['dqi']=core/'results/trimmed_default_verify.json'
    for name in ['incumbent_seed_2026102405','trimmed_seed_2026102405','trimmed_default_verify','trimmed_independent_verify']:
        p=core/f'results/{name}.json';r=json.loads(p.read_text())
        assert r['status']=='PASS' and r['verified_shots']==102400 and r['num_qubits']==961 and r['log2_delta']==-128.0
        if name.startswith('trimmed'):assert r['circuit_sha256']==acceptance['circuit_sha256']
        proofs[name]=dict(path=str(p),sha256=sha(p),report=r)
    old=proofs['incumbent_seed_2026102405']['report'];new=proofs['trimmed_seed_2026102405']['report']
    for field in ['value_errors','phase_errors','verified_shots','num_qubits','effective_error_rate','repetition_factor']:assert old[field]==new[field]
    assert old['toffoli_count']-new['toffoli_count']==25596.0
    assert acceptance['native_verifier_sha256']==sha(ROOT/'from_host/verifiers/dqi/verify_submission_dqi.py')
    assert acceptance['default']['score']==proofs['trimmed_default_verify']['report']['score']
    dump(HERE/'results/package_trimmed_dqi_proofs.json',dict(acceptance=acceptance,runs=proofs))
native={k:json.loads(p.read_text()) for k,p in reports.items()}
expected=dict(invsqrt=200540928,dqi=102400,dfthcsos=50880)
for k,r in native.items():assert r['status']=='PASS' and r['verified_shots']==expected[k]
for k,p in sources.items():assert sha(p)==native[k]['circuit_sha256']
z=dict(np.load(chem/'submission.npz'))
assert submission_sha256(load_circuits(chem),z,z['tree_edges'])==native['dfthcsos']['submission_sha256']==selection['canonical_digest']
members={'invsqrt.kmb':sources['invsqrt'],'dqi.kmb':sources['dqi'],'dfthcsos/submission.npz':chem/'submission.npz'}
members.update({f'dfthcsos/angle_k_{k}.kmb':chem/f'angle_k_{k}.kmb' for k in range(11)})
assert len(members)==14
hashes={name:sha(path) for name,path in members.items()}
temporary=HERE/'submission.zip.tmp';archive=HERE/'submission.zip'
with zipfile.ZipFile(temporary,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as out:
    for name,path in members.items():out.write(path,arcname=name)
with zipfile.ZipFile(temporary) as check:
    assert len(check.namelist())==14 and set(check.namelist())==set(members)
    for name in members:
        h=hashlib.sha256()
        with check.open(name) as f:
            for chunk in iter(lambda:f.read(8<<20),b''):h.update(chunk)
        assert h.hexdigest()==hashes[name] # Reads to EOF also validate ZIP CRC.
if archive.exists() and sha(archive)!=sha(temporary):
    backup=HERE/'revisions';backup.mkdir(exist_ok=True)
    shutil.copy2(archive,backup/f'package_{sha(archive)}.zip')
temporary.replace(archive)
accepted_score=math.prod(BASE[k]/native[k]['score'] for k in BASE)**(1/3)
record=dict(status='ASSEMBLED_BOUND_TO_ACCEPTED_RESULTS',archive=str(archive),sha256=sha(archive),bytes=archive.stat().st_size,
    accepted_reference_score=accepted_score,baseline_costs=BASE,accepted_reference_task_scores={k:r['score'] for k,r in native.items()},
    source_reports={k:str(p) for k,p in reports.items()},source_artifacts={k:str(p) for k,p in sources.items()},chemistry=str(chem),
    canonical_chemistry_digest=selection['canonical_digest'],members=hashes,verified_reference_shots=expected,
    dqi_decision=dict(path=str(decision_path),sha256=sha(decision_path),accepted_trimmed_decoder=trimmed,full_trimmed_native_runs=bool(proofs),
        selected='accepted trimmed decoder' if trimmed else 'round-4 Codex',proofs=str(HERE/'results/package_trimmed_dqi_proofs.json') if trimmed else None),
    scope='ZIP CRC/content hashes and full accepted native reports bound; fresh host-metric validation pending')
dump(HERE/'results/package_assembly.json',record)
dump(HERE/'submission.summary.json',record)
(HERE/'results/dqi_core_selection_evidence.txt').write_text(decision)
print(json.dumps({k:v for k,v in record.items() if k not in ['members']},indent=2),flush=True)
