"""Count source-preserving retained inversions, then build a small oracle."""
import argparse
import json
from pathlib import Path

import dqi_budget_probe as probe
import low_workspace_inverse
from safe_retained_inverse import bind


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--small',action='store_true');args=p.parse_args()
    decoder=probe.decoder
    decoder.frames=probe.frames;decoder.BoundedArithmetic=probe.trimmed.TrimmedArithmetic
    decoder.load=probe.count_module
    decoder.build.profile_only=True;decoder.build.slot_counts=probe.counts
    low_workspace_inverse.bind=bind
    if args.small:
        from import_published import load
        from blocked_euclid_sequence import frames
        decoder.load=load;decoder.frames=frames;decoder.build.profile_only=False
        decoder.build.slot_counts=None
        decoder.build(14,2,100,127,Path('candidates/safe_small_n14.kmb'),
                      linear_padding=True,virtual_zero=True,virtual_lead=True,
                      low_inverse=True,lazy_flags=True)
    else:
        report=decoder.build(70,3,100,4095,Path('candidates/safe_inverse_profile.kmb'),
                             linear_padding=True,virtual_zero=True,virtual_lead=True,
                             low_inverse=True,lazy_flags=True)
        report.update(saved_core_Toffolis=801388-report['T'])
        Path('results/safe_inverse_profile.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps({k:v for k,v in report.items() if k!='rows'}),flush=True)
