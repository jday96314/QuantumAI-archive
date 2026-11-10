"""Full, unchanged host CLI checks for a predeclared local release manifest."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('manifest', type=Path)
    args = parser.parse_args()
    manifest_sha = digest(args.manifest)
    frozen = json.loads(args.manifest.read_text())
    runs = []
    for task, spec in frozen['tasks'].items():
        circuit = HERE/spec['path']
        verifier = (HERE/'../../../from_host/verifiers'/task/f'verify_submission_{task}.py').resolve()
        for run, salt in [('default', ''), ('independent', spec['independent_salt'])]:
            assert digest(circuit) == spec['sha256'] and digest(verifier) == spec['verifier_sha256']
            assert digest(args.manifest) == manifest_sha
            destination = HERE/'results'/f"{frozen['report_prefix']}_{task}_{run}_verify.json"
            command = [sys.executable, str(verifier), str(circuit), '--max-workers', '2']
            if salt:
                command += ['--salt', salt]
            print(json.dumps(dict(stage='FULL_NATIVE_START', task=task, run=run, salt=salt,
                                  circuit_sha256=spec['sha256'], manifest_sha256=manifest_sha)), flush=True)
            start = time.time()
            with destination.open('w') as out, destination.with_suffix('.err').open('w') as err:
                result = subprocess.run(command, stdout=out, stderr=err, cwd=HERE)
            report = json.loads(destination.read_text())
            row = dict(task=task, run=run, salt=salt, seconds=time.time()-start, returncode=result.returncode,
                       status=report['status'], score=report.get('score'), verified_shots=report.get('verified_shots'),
                       phase_errors=report.get('phase_errors'), value_errors=report.get('value_errors'),
                       circuit_sha256=spec['sha256'], manifest_sha256=manifest_sha)
            runs.append(row)
            (HERE/'results'/f"{frozen['report_prefix']}_verification_runs.json").write_text(json.dumps(runs, indent=2)+'\n')
            print(json.dumps(row), flush=True)
            assert result.returncode == 0 and report['status'] == 'PASS'
            assert report['verified_shots'] == spec['shots'] and report['circuit_sha256'] == spec['sha256']
            assert report['num_qubits'] == spec['Q']
            if task == 'dqi' and spec.get('require_zero_value_errors', True):
                assert report['value_errors'] == 0
            assert digest(circuit) == spec['sha256'] and digest(verifier) == spec['verifier_sha256']
    print(json.dumps(dict(status='FULL_NATIVE_CHECKS_COMPLETE', runs=runs)), flush=True)


if __name__ == '__main__':
    main()
