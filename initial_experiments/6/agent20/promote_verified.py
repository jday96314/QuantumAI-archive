"""Promote an audited archive while preserving the previous best directory."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import shutil

HERE = Path(__file__).resolve().parent


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('candidate',type=Path)
    parser.add_argument('--label',required=True)
    parser.add_argument('--independent-invsqrt',type=Path,required=True)
    parser.add_argument('--independent-dqi',type=Path,required=True)
    parser.add_argument('--independent-chemistry',type=Path,required=True)
    args = parser.parse_args()
    assert args.label.replace('_','').isalnum()
    candidate = args.candidate.resolve()
    summary = json.loads((candidate/'submission.summary.json').read_text())
    audit = json.loads((candidate/'results/metric_preflight_and_score.json').read_text())
    previous = json.loads((HERE/'best/submission.summary.json').read_text())
    assert audit['status']==summary['status']=='PASS'
    assert abs(summary['overall_score']-audit['overall_score'])<1e-12
    assert summary['overall_score']>previous['overall_score']
    assert digest(candidate/'submission.zip')==summary['sha256']==audit['archive_sha256']
    independent = {}
    for task,path in [('invsqrt',args.independent_invsqrt),('dqi',args.independent_dqi),('dfthcsos',args.independent_chemistry)]:
        report = json.loads(path.read_text())
        default = json.loads((candidate/'results'/f'{task}_verify.json').read_text())
        assert report['status']=='PASS' and report['verified_shots']==default['verified_shots']
        key = 'submission_sha256' if task=='dfthcsos' else 'circuit_sha256'
        assert report[key]==default[key]
        independent[task] = report['score']
    independent_score = math.prod(summary['baseline_costs'][k]/v for k,v in independent.items())**(1/3)
    destination = HERE/f'verified_before_{args.label}'
    pending = HERE/'best.pending'
    assert not destination.exists() and not pending.exists()
    shutil.copytree(candidate,pending)
    (HERE/'best').rename(destination)
    pending.rename(HERE/'best')
    assert digest(destination/'submission.zip')==previous['sha256']
    assert digest(HERE/'best/submission.zip')==summary['sha256']
    out = dict(status='PROMOTED_VERIFIED_IMPROVEMENT',label=args.label,
               default_score=summary['overall_score'],independent_combination_score=independent_score,
               archive_sha256=summary['sha256'],previous_score=previous['overall_score'],
               previous_archive_sha256=previous['sha256'],goal20_complete=False,
               milestones={str(k):digest(HERE/f'milestones/{k}/submission.zip') for k in range(10,21)
                           if (HERE/f'milestones/{k}/submission.zip').exists()})
    (HERE/'results'/f'accepted_{args.label}_improvement.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2),flush=True)


if __name__=='__main__':
    main()
