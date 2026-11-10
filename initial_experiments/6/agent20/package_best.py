"""Package hash-bound full-verifier evidence and save crossed milestones."""
import argparse
import hashlib
import importlib.util
import inspect
import json
import shutil
from pathlib import Path

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parents[2]
ROUND5=PROJECT/'initial_experiments/5'


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--inverse',type=Path,default=HERE/'candidates/rows_odd_best.kmb')
    p.add_argument('--inverse-result',type=Path,default=HERE/'results/rows_odd_default_verify.json')
    p.add_argument('--dqi',type=Path,default=HERE/'candidates/safe_retained_dqi.kmb')
    p.add_argument('--dqi-result',type=Path,default=HERE/'results/dqi_safe_default_verify.json')
    p.add_argument('--output',type=Path,default=HERE/'best')
    p.add_argument('--notice',type=Path)
    p.add_argument('--license',type=Path)
    p.add_argument('--chemistry',type=Path)
    p.add_argument('--chemistry-result',type=Path)
    args=p.parse_args()
    spec=importlib.util.spec_from_file_location('round5_packager',ROUND5/'integration/package_submission.py')
    package=importlib.util.module_from_spec(spec);spec.loader.exec_module(package)
    package.HERE=args.output.resolve()
    package.RESULTS=package.HERE/'results';package.EXTRACTED=package.HERE/'submission'
    package.ARCHIVE=package.HERE/'submission.zip'
    package.INVERSE=args.inverse.resolve();package.DQI=args.dqi.resolve()
    package.CHEMISTRY=ROUND5/'chem-census/candidates/structured'
    if args.chemistry:
        package.CHEMISTRY=args.chemistry.resolve()
    package.PREVIOUS_SUMMARY=ROUND5/'integration/wave3/submission.summary.json'
    package.EVIDENCE=dict(invsqrt=args.inverse_result.resolve(),dqi=args.dqi_result.resolve(),
                          dfthcsos=ROUND5/'integration/wave3/results/chemistry_verify_default.json')
    if args.chemistry_result:
        package.EVIDENCE['dfthcsos']=args.chemistry_result.resolve()
    if args.notice or args.license:
        if not (args.notice and args.license):
            p.error('--notice and --license must be supplied together')
        package.ADDITIONAL_FILES={'NOTICE.md':args.notice.resolve(),'LICENSE':args.license.resolve()}
        assert all(path.is_file() for path in package.ADDITIONAL_FILES.values())
        source=inspect.getsource(package.build)
        marker='    for task, source in EVIDENCE.items():'
        assert source.count(marker)==1
        source=source.replace(marker,'    files.update(ADDITIONAL_FILES)\n'+marker)
        exec(source,package.__dict__)
    package.build()
    package.verify(3,reuse_tasks=('invsqrt','dqi','dfthcsos'))
    package.score_and_audit()
    summary=json.loads((package.HERE/'submission.summary.json').read_text())
    for threshold in range(10,21):
        if summary['overall_score']<threshold:continue
        milestone=HERE/'milestones'/str(threshold)
        if (milestone/'submission.zip').exists():continue
        milestone.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(package.ARCHIVE,milestone/'submission.zip')
        shutil.copyfile(package.HERE/'submission.summary.json',milestone/'submission.summary.json')
        shutil.copytree(package.RESULTS,milestone/'results',dirs_exist_ok=True)
        with (milestone/'submission.zip').open('rb') as f:
            assert hashlib.file_digest(f,'sha256').hexdigest()==summary['sha256']
        (milestone/'milestone.json').write_text(json.dumps(dict(threshold=threshold,
            verified_local_score=summary['overall_score'],archive_sha256=summary['sha256'],
            scope='full default-shot native evidence bound to delivered archive members'),indent=2)+'\n')
        print(f"Saved milestone {threshold}: {summary['overall_score']:.12f}",flush=True)


if __name__=='__main__':main()
