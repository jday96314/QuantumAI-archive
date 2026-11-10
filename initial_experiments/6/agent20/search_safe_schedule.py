"""Reoptimize the decoder's one/two/three-slot schedule after arithmetic changes.

Only resource counts select candidates. These profiles do not establish
decoding accuracy; selected schedules still need whole-circuit verification.
"""
import json
import multiprocessing as mp
import argparse
from functools import partial
from pathlib import Path

import dqi_budget_probe as probe
import profile_schedule as original
from safe_retained_inverse import bind

HERE=Path(__file__).resolve().parent
RETENTION_CAP=1004


def stage(task):
    original.bind=partial(bind,limit=RETENTION_CAP)
    original.BoundedArithmetic=probe.trimmed.TrimmedArithmetic
    return original.stage(task)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--retention-cap',type=int,default=1004)
    args=parser.parse_args();RETENTION_CAP=args.retention_cap
    suffix='' if RETENTION_CAP==1004 else f'_{RETENTION_CAP}'
    path=HERE/f'results/safe_schedule_stages{suffix}.json'
    rows=json.loads(path.read_text()) if path.exists() else []
    done={(r['k'],r['count']) for r in rows}
    tasks=[(k,c) for k in range(1,35) for c in (1,2,3) if k+c<=35 and (k,c) not in done]
    with mp.Pool(2) as pool:
        for row in pool.imap_unordered(stage,tasks):
            rows.append(row);path.write_text(json.dumps(rows,indent=2)+'\n')
            print(json.dumps(row),flush=True)
    costs={(r['k'],r['count']):r for r in rows}
    reference=json.loads((HERE/'results/safe_inverse_profile.json').read_text())
    fixed=reference['T']-sum(r['T'] for r in reference['rows'] if r['label'].startswith(('forward_','inverse_')))
    frontier=[]
    for cap in sorted({975,*[r['Q'] for r in rows if r['Q']>=975]}):
        dp={35:(0,[])}
        for k in reversed(range(1,35)):
            candidates=[]
            for count in (1,2,3):
                row=costs.get((k,count))
                if row and row['Q']<=cap and k+count in dp:
                    candidates.append((row['T']+dp[k+count][0],[count]+dp[k+count][1]))
            if candidates:dp[k]=min(candidates)
        if 1 not in dp:continue
        total=dp[1][0]+fixed+371150
        frontier.append(dict(Q=cap,static_Toffolis=total,core_Toffolis=total-371150,
                             objective_without_repetition=cap*cap*total,counts=dp[1][1]))
    frontier.sort(key=lambda r:r['objective_without_repetition'])
    (HERE/f'results/safe_schedule_frontier{suffix}.json').write_text(json.dumps(frontier,indent=2)+'\n')
    print(json.dumps(dict(status='RESOURCE_SEARCH_COMPLETE',best=frontier[:4])),flush=True)
