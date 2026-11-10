"""Package one independently verified candidate and record its local score.

Run from the repository root. A milestone archive is written only when its
integer threshold is met. The ordinary submission.zip always holds the best
fully verified candidate built in this experiment.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import zipfile


BASELINE = {
    'invsqrt': 2076774.233715,
    'dqi': 24439463583672.929688,
    'dfthcsos': 3355599588.315537,
}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--invsqrt', type=Path, required=True)
    p.add_argument('--invsqrt-result', type=Path, required=True)
    p.add_argument('--dqi', type=Path, required=True)
    p.add_argument('--dqi-result', type=Path, required=True)
    p.add_argument('--chemistry', type=Path, required=True)
    p.add_argument('--chemistry-result', type=Path, required=True)
    p.add_argument('--output', type=Path, default=Path('initial_experiments/3/submission.zip'))
    p.add_argument('--milestone', type=int)
    a = p.parse_args()
    sources = {'invsqrt': a.invsqrt_result, 'dqi': a.dqi_result,
               'dfthcsos': a.chemistry_result}
    results = {k: json.loads(v.read_text()) for k, v in sources.items()}
    assert all(v['status'] == 'PASS' for v in results.values())
    improvements = {k: BASELINE[k] / results[k]['score'] for k in BASELINE}
    overall = math.prod(improvements.values()) ** (1 / 3)
    if a.milestone is not None and overall < a.milestone:
        raise ValueError(f'candidate score {overall:.6f} does not reach milestone {a.milestone}')
    files = [('invsqrt.kmb', a.invsqrt), ('dqi.kmb', a.dqi)]
    files += [(f'dfthcsos/{path.name}', path)
              for path in sorted(a.chemistry.glob('angle_k_*.kmb'))]
    files += [('dfthcsos/submission.npz', a.chemistry / 'submission.npz')]
    expected_circuits = results['dfthcsos']['num_circuits']
    assert len(files) == expected_circuits + 3 and all(path.is_file() for _, path in files)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    tmp = a.output.with_suffix('.zip.tmp')
    with zipfile.ZipFile(tmp, 'w', compression=zipfile.ZIP_DEFLATED,
                         compresslevel=6) as archive:
        for name, path in files:
            archive.write(path, name)
    with zipfile.ZipFile(tmp) as archive:
        assert archive.testzip() is None
    tmp.replace(a.output)
    summary = dict(overall_score=overall,
                   task_scores={k: results[k]['score'] for k in BASELINE},
                   improvements=improvements,
                   sources={k: str(v) for k, v in sources.items()},
                   archive=str(a.output),
                   sha256=hashlib.sha256(a.output.read_bytes()).hexdigest())
    a.output.with_suffix('.summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
