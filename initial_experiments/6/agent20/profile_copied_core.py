"""Exact core profiles with budgeted reciprocal-prefix storage, no tape retained."""
from pathlib import Path
import json
import dqi_budget_probe as probe
import copied_prefix_inverse as copied
import low_workspace_inverse
import fused_gap_block
import gap_serial
from safe_retained_inverse import bind
HERE=Path(__file__).resolve().parent

if __name__=='__main__':
    rows=[];instances=[]
    probe.decoder.frames=probe.frames;probe.decoder.load=probe.count_module
    probe.decoder.build.profile_only=True;probe.decoder.build.slot_counts=probe.counts
    low_workspace_inverse.bind=bind
    for limit in (940,980,1004,1080,1200):
        class Arithmetic(copied.CopiedArithmetic):
            def __init__(self,*args,**kwargs):
                super().__init__(*args,limit=limit,**kwargs);instances.append(self)
        probe.decoder.BoundedArithmetic=Arithmetic
        fused_gap_block.BoundedArithmetic=Arithmetic
        gap_serial.BoundedArithmetic=Arithmetic
        instances.clear()
        r=probe.decoder.build(70,3,100,4095,HERE/'candidates'/f'copied_core_limit{limit}.kmb',
            linear_padding=True,virtual_zero=True,virtual_lead=True,low_inverse=True,lazy_flags=True)
        r.update(copy_limit=limit,copied_calls=[call for ar in instances for call in ar.copied_calls],
                 baseline_core_T=791596,scope='exact emitted accepted-core arithmetic profile; no native full circuit')
        rows.append(r);(HERE/'results/copied_core_profiles_nested.json').write_text(json.dumps(rows,indent=2)+'\n')
        print(json.dumps({k:v for k,v in r.items() if k!='rows'}),flush=True)
