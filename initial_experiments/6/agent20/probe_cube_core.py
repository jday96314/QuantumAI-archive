"""Exact resource profile for the cube-specialized gap-aware decoder core."""
import json
from pathlib import Path
import dqi_budget_probe as probe
from cube_arithmetic import CubeArithmetic
import low_workspace_inverse


if __name__=='__main__':
    decoder=probe.decoder
    decoder.frames=probe.frames
    decoder.BoundedArithmetic=CubeArithmetic
    decoder.load=probe.count_module
    decoder.build.profile_only=True;decoder.build.slot_counts=probe.counts
    low_workspace_inverse.bind=lambda m:probe.budget.bind(m,980)
    report=decoder.build(70,3,100,4095,Path('candidates/cube_core_profile.kmb'),
                         linear_padding=True,virtual_zero=True,virtual_lead=True,
                         low_inverse=True,lazy_flags=True)
    report.update(saved_core_Toffolis=801388-report['T'])
    Path('results/cube_core_profile.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='rows'}),flush=True)
