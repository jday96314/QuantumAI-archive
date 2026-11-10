"""Compare complete C-fraction arithmetic with the new cube primitive.

These are generic-chart resource profiles, not accepted DQI scores. Zero-head
chart handling remains required before any full-size submission could use it.
"""
import importlib
import json
import sys
from pathlib import Path

import dqi_budget_probe as probe
from cube_arithmetic import CubeArithmetic

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(probe.ROUND5/'dqi-core-phase'))
cf=importlib.import_module('cfraction_decoder')


if __name__=='__main__':
    reports=[]
    for name,arithmetic in [('trimmed',probe.trimmed.TrimmedArithmetic),('cube',CubeArithmetic)]:
        cf.BoundedArithmetic=arithmetic
        _,_,_,_,report=cf.emit(70,count_only=True)
        report.update(architecture=name,phase_scan_Toffolis=80*4095,
                      whole_oracle_static_Toffolis=report['T']+80*4095)
        reports.append(report);print(json.dumps(report),flush=True)
        (HERE/'results/cfraction_cube_profiles.json').write_text(json.dumps(reports,indent=2)+'\n')
