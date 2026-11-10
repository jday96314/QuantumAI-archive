"""Verify frozen Apache-2.0 public components with unchanged host CLIs.

Fixed salts are declared before verification. No sampling or cost changes.
See public/ahmed_v3/download_manifest.json and NOTICE.md for attribution.
"""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
FROZEN = {
    'invsqrt': ('58c7ccb60f61088c9c3d14a562f6c1d795bcb8ac12b364cb16941b54b74f5cf1', 200540928),
    'dqi': ('a91c66ac9d0a720e112976ae13a0c4b855223c6605c42b0809de2ea7e6ad5a25', 102400),
}


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    reports = []
    for task, (sha, shots) in FROZEN.items():
        circuit = HERE / 'public/ahmed_v3/candidate' / f'{task}.kmb'
        assert digest(circuit) == sha
        verifier = (HERE / '../../../from_host/verifiers' / task / f'verify_submission_{task}.py').resolve()
        verifier_sha = digest(verifier)
        for run, salt in [('default', ''), ('independent', f'agent20-public-ahmed-v3-{task}-independent-20261025')]:
            destination = HERE / 'results' / f'public_{task}_{run}_verify.json'
            err = destination.with_suffix('.err')
            command = [sys.executable, str(verifier), str(circuit), '--max-workers', '2']
            if salt:
                command += ['--salt', salt]
            print(json.dumps(dict(stage='FULL_NATIVE_START', task=task, run=run, salt=salt, circuit_sha256=sha, verifier_sha256=verifier_sha)), flush=True)
            started = time.time()
            with destination.open('w') as out, err.open('w') as error:
                result = subprocess.run(command, stdout=out, stderr=error, cwd=HERE)
            report = json.loads(destination.read_text())
            row = dict(task=task, run=run, salt=salt, seconds=time.time()-started, returncode=result.returncode,
                       status=report['status'], score=report.get('score'), verified_shots=report.get('verified_shots'),
                       phase_errors=report.get('phase_errors'), value_errors=report.get('value_errors'),
                       circuit_sha256=sha, verifier_sha256=verifier_sha)
            reports.append(row)
            (HERE / 'results/public_component_verification_runs.json').write_text(json.dumps(reports, indent=2)+'\n')
            print(json.dumps(row), flush=True)
            assert result.returncode == 0 and report['status'] == 'PASS'
            assert report['verified_shots'] == shots and report['circuit_sha256'] == sha
            if task == 'dqi':
                assert report['value_errors'] == 0
            assert digest(circuit) == sha and digest(verifier) == verifier_sha
    print(json.dumps(dict(status='FOUR_FULL_NATIVE_CHECKS_COMPLETE', runs=reports)), flush=True)


if __name__ == '__main__':
    main()
