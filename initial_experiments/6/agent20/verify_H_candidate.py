"""Run both frozen full-default host checks; never vary salts or circuit bytes."""
import hashlib,json,subprocess,time,sys
from pathlib import Path

HERE=Path(__file__).resolve().parent
manifest=json.loads((HERE/'results/repaired_H_build.json').read_text())
frozen=json.loads((HERE/'parameters/repaired_H_pre_native.json').read_text())
path=Path(manifest['circuit']);digest=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert manifest['Q']==frozen['expected_Q'] and manifest['static_T']==frozen['expected_static_T']
assert digest(path)==manifest['sha256']
assert all(digest(source)==sha for source,sha in manifest['source_sha256'].items())
assert all(digest(source)==sha for source,sha in frozen['sources'].items())
verifier=HERE/'../../../from_host/verifiers/dqi/verify_submission_dqi.py'
runs=[]
for name,salt in [('default',''),('independent',frozen['independent_salt'])]:
    command=[sys.executable,str(verifier),str(path),'--max-workers','2']
    if salt:command+=['--salt',salt]
    print(json.dumps(dict(stage='FULL_NATIVE_START',run=name,salt=salt,sha256=manifest['sha256'])),flush=True)
    started=time.time()
    with (HERE/f'results/dqi_H_{name}_verify.json').open('w') as out,(HERE/f'results/dqi_H_{name}_verify.err').open('w') as err:
        result=subprocess.run(command,stdout=out,stderr=err,cwd=HERE)
    report=json.loads((HERE/f'results/dqi_H_{name}_verify.json').read_text())
    row=dict(run=name,salt=salt,seconds=time.time()-started,returncode=result.returncode,status=report['status'],score=report.get('score'),value_errors=report.get('value_errors'),phase_errors=report.get('phase_errors'))
    runs.append(row);(HERE/'results/dqi_H_verification_runs.json').write_text(json.dumps(runs,indent=2)+'\n');print(json.dumps(row),flush=True)
    assert result.returncode==0 and report['status']=='PASS'
    assert report['verified_shots']==102400 and report['value_errors']==0 and report['circuit_sha256']==manifest['sha256']
assert digest(path)==manifest['sha256']
print(json.dumps(dict(status='TWO_FULL_NATIVE_RUNS_COMPLETE',sha256=manifest['sha256'],runs=runs)),flush=True)
