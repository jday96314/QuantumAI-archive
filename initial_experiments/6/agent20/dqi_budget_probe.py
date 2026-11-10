"""Resource-only port of retained inversion to the accepted gap-aware core."""
import importlib.util
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROUND5 = HERE.parent.parent/'5'
GAP = ROUND5/'dqi-gap'
sys.path.insert(0,str(GAP))
os.chdir(GAP)
import hybrid_decoder as decoder
from blocked_monic_euclid import schedule
from profile_hybrid import module as count_module
import low_workspace_inverse
os.chdir(HERE)


def module_from_file(name, path):
    spec = importlib.util.spec_from_file_location(name,path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


trimmed = module_from_file('agent20_trimmed_middle',ROUND5/'dqi-core/components/trimmed_middle.py')
budget = module_from_file('agent20_budget_inverse',ROUND5/'dqi-core/components/budget_inverse.py')
counts = json.loads((GAP/'parameters/best_schedule.json').read_text())['counts']


def frames(words,width):
    previous,current,tail,pivot = [],[words[-2]],words[:-2],words[-1]
    k,stages=1,[]
    for count in counts:
        stages.append((k,count,previous,current,tail,pivot))
        _,(previous,current) = schedule(previous,current,list(reversed(tail))[:2*count])
        tail=tail[:-2*count];k+=count
    assert k==35 and not tail
    return stages,(previous,current,tail,pivot)


if __name__ == '__main__':
    decoder.frames=frames
    decoder.BoundedArithmetic=trimmed.TrimmedArithmetic
    decoder.load=count_module
    decoder.build.profile_only=True
    decoder.build.slot_counts=counts
    original_bind=low_workspace_inverse.bind
    reports=[]
    for limit in [0,940,960,980,1004]:
        low_workspace_inverse.bind=lambda m,limit=limit:budget.bind(m,limit)
        report=decoder.build(70,3,100,4095,Path(f'candidates/budget_profile_{limit}.kmb'),
                             linear_padding=True,virtual_zero=True,virtual_lead=True,
                             low_inverse=True,lazy_flags=True)
        report.update(retention_limit=limit,saved_core_Toffolis=801388-report['T'])
        reports.append(report)
        print(json.dumps({k:v for k,v in report.items() if k!='rows'}),flush=True)
    low_workspace_inverse.bind=original_bind
    (HERE/'results/dqi_budget_profiles.json').write_text(json.dumps(reports,indent=2)+'\n')
