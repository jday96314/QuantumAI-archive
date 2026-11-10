"""Independently reopen the final ZIP and bind it to full native results."""
import hashlib
import io
import json
import math
import zipfile
from pathlib import Path

import kickmix as km
import numpy as np
from dfthcsos.scoreinfra.circuits import submission_sha256


def main():
    sha = lambda data: hashlib.sha256(data).hexdigest()
    summary = json.loads(Path('submission.summary.json').read_text())
    archive = Path('submission.zip')
    assert sha(archive.read_bytes()) == summary['sha256']
    native = {task: json.loads(Path(path).read_text()) for task, path in summary['sources'].items()}
    expected_shots = dict(invsqrt=200540928, dqi=102400, dfthcsos=50880)
    for task, result in native.items():
        assert result['status'] == 'PASS'
        assert result['verified_shots'] == expected_shots[task]
        assert result['score'] == summary['task_scores'][task]
    with zipfile.ZipFile(archive) as z:
        expected = {'invsqrt.kmb', 'dqi.kmb', 'dfthcsos/submission.npz'}
        expected.update(f'dfthcsos/angle_k_{k}.kmb' for k in range(native['dfthcsos']['num_circuits']))
        assert len(z.namelist()) == len(expected) and set(z.namelist()) == expected
        assert z.testzip() is None
        hashes = {name: sha(z.read(name)) for name in sorted(expected)}
        assert hashes['invsqrt.kmb'] == native['invsqrt']['circuit_sha256']
        assert hashes['dqi.kmb'] == native['dqi']['circuit_sha256']
        circuits = [km.Circuit.from_file(io.BytesIO(z.read(f'dfthcsos/angle_k_{k}.kmb')), format='kmb')
                    for k in range(native['dfthcsos']['num_circuits'])]
        with np.load(io.BytesIO(z.read('dfthcsos/submission.npz'))) as arrays:
            params = {key: arrays[key] for key in ['U', 'W', 'H2_I', 'H1_I', 'Hsym_lower']}
            digest = submission_sha256(circuits, params, arrays['tree_edges'])
        assert digest == native['dfthcsos']['submission_sha256']
    baseline = dict(invsqrt=2076774.233715, dqi=24439463583672.929688, dfthcsos=3355599588.315537)
    score = math.prod(baseline[t] / native[t]['score'] for t in baseline)**(1/3)
    assert score == summary['overall_score']
    record = dict(status='PASS', archive_sha256=summary['sha256'], members=hashes,
                  chemistry_digest=digest, verified_shots=expected_shots, overall_score=score,
                  scope='ZIP reopened independently; circuit hashes, canonical array digest and full native results bound')
    Path('results/archive_independent_audit.json').write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps({k: v for k, v in record.items() if k != 'members'}, indent=2))


if __name__ == '__main__':
    main()
