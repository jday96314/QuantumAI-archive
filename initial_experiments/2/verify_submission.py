"""Run the unchanged host verifiers at full competition sizes and report the local score."""
import argparse
import json
import math
from pathlib import Path
import sys
import time
ROOT=Path(__file__).resolve().parents[1]
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(ROOT/'from_host/kickmix'),str(ROOT/'from_host/verifiers'),str(HERE)]
BASELINES={'invsqrt':2076774.233715282,'dqi':24439463583672.93,'dfthcsos':3355599588.3155375}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--submission',type=Path,default=HERE/'submission')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--salt',default='')
    p.add_argument('--workers',type=int,default=4)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    import kickmix as km
    from invsqrt.verify_submission_invsqrt import verify as inv_verify,ProblemSpecInvSqrt
    from dqi.verify_submission_dqi import verify as dqi_verify,ProblemSpecDQI
    from dfthcsos.verify import verify as chem_verify
    results={}
    for task in ['invsqrt','dqi','dfthcsos']:
        start=time.monotonic()
        if task=='dfthcsos':result=chem_verify(a.submission/task,salt=a.salt)
        else:
            c=km.Circuit.from_file(a.submission/f'{task}.kmb',format='kmb')
            if task=='invsqrt':result=inv_verify(c,ProblemSpecInvSqrt.for_n_p(10),salt=a.salt,max_workers=a.workers)
            else:result=dqi_verify(c,ProblemSpecDQI.benchmark_instance(),salt=a.salt,max_workers=a.workers)
        assert result['status']=='PASS'
        result['wall_seconds']=time.monotonic()-start
        results[task]=result
        (a.output/f'{task}.json').write_text(json.dumps(result,indent=2)+'\n')
        print(task,result['status'],'cost',result['score'],'seconds',result['wall_seconds'],flush=True)
    ratios={k:BASELINES[k]/v['score'] for k,v in results.items()}
    summary={'status':'PASS','salt':a.salt,'baseline_costs':BASELINES,'costs':{k:v['score'] for k,v in results.items()},'improvement_ratios':ratios,'overall_score':math.prod(ratios.values())**(1/3),'note':'Local normalization against the recorded empty-salt host baselines; official split salts and reference costs are private.'}
    (a.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':main()
