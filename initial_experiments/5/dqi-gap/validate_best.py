"""Unchanged 102400-shot host verification, fixed candidate and independent salts."""
import concurrent.futures,json,subprocess,sys
from pathlib import Path
from build_hash import digest

manifest=json.loads(Path('parameters/best_candidate_manifest.json').read_text())
assert digest(manifest['circuit'])==manifest['sha256']
for source,expected in manifest['source_sha256'].items():assert digest(source)==expected,source
adversarial=[sys.executable,'verify_adversarial.py',manifest['circuit'],'--shots','1024','--output','results/best_adversarial.json']
existing=Path('results/best_adversarial.json')
cached=json.loads(existing.read_text()) if existing.exists() else {}
if not(cached.get('status')=='PASS' and cached.get('sha256')==manifest['sha256'] and cached.get('shots')==1024):
    with Path('results/best_adversarial.log').open('w') as log:subprocess.run(adversarial,stdout=log,stderr=subprocess.STDOUT,check=True)
print('Adversarial phase/restoration/composition checks passed.',flush=True)
runs=[('default',None),('salt_a','dqi-gap-validation-20261024-a'),('salt_b','dqi-gap-validation-20261024-b')]
def verify(item):
    name,salt=item;cmd=[sys.executable,'../../../from_host/verifiers/dqi/verify_submission_dqi.py',manifest['circuit'],'--max-workers','4']
    if salt:cmd += ['--salt',salt]
    with Path(f'results/best_verify_{name}.json').open('w') as out,Path(f'results/best_verify_{name}.err').open('w') as err:
        done=subprocess.run(cmd,stdout=out,stderr=err)
    report=json.loads(Path(f'results/best_verify_{name}.json').read_text());report['validation_name']=name;report['exit_code']=done.returncode
    print(json.dumps({k:report[k] for k in ['validation_name','status','value_errors','phase_errors','score'] if k in report}),flush=True)
    return report
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:reports=list(pool.map(verify,runs))
for r in reports:
    assert r['status']=='PASS' and r['verified_shots']==102400,r
    assert r['circuit_sha256']==manifest['sha256'],r
incumbent=1598204229384.9265
summary=dict(status='PASS',scope='fixed resource-selected candidate, unchanged host confidence/scoring, default hash seed and two independent salts; no salt selection',circuit=manifest['circuit'],sha256=manifest['sha256'],runs=[dict(name=r['validation_name'],Q=r['num_qubits'],T=r['toffoli_count'],value_errors=r['value_errors'],phase_errors=r['phase_errors'],R=r['repetition_factor'],effective_error_rate=r['effective_error_rate'],score=r['score'],improvement=incumbent/r['score']) for r in reports])
Path('results/best_validation_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)
