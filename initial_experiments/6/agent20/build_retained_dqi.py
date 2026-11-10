"""Compile the resource-selected, fixed-schedule retained-inverse decoder."""
import hashlib
import json
from pathlib import Path

import dqi_budget_probe as probe
import low_workspace_inverse
from import_published import load

HERE=Path(__file__).resolve().parent


def main():
    reports=json.loads((HERE/'results/dqi_budget_profiles.json').read_text())
    feasible=[r for r in reports if r['Q']<=1004]
    best=min(feasible,key=lambda r:r['T'])
    assert best['retention_limit']==980
    manifest=dict(status='FROZEN_BEFORE_NATIVE_VERIFICATION',
                  design='retained field-inversion intermediates under a measured lifetime budget; fixed accepted gap schedule and trimmed middle products',
                  selected_profile={k:v for k,v in best.items() if k!='rows'},
                  expected_static_Toffolis=1172538-best['saved_core_Toffolis'],
                  independent_salt='agent20-retained-dqi-independent-20261025',
                  sources={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [
                      Path(__file__),HERE/'dqi_budget_probe.py',
                      probe.ROUND5/'dqi-core/components/budget_inverse.py',
                      probe.ROUND5/'dqi-core/components/trimmed_middle.py',
                      probe.GAP/'hybrid_decoder.py',probe.GAP/'parameters/best_schedule.json']})
    (HERE/'parameters/retained_dqi_pre_native.json').write_text(json.dumps(manifest,indent=2)+'\n')
    decoder=probe.decoder
    decoder.load=load;decoder.frames=probe.frames
    decoder.BoundedArithmetic=probe.trimmed.TrimmedArithmetic
    decoder.build.profile_only=False;decoder.build.slot_counts=probe.counts
    low_workspace_inverse.bind=lambda m:probe.budget.bind(m,best['retention_limit'])
    path=HERE/'candidates/retained_inverse_dqi.kmb'
    circuit=decoder.build(70,3,100,4095,path,linear_padding=True,virtual_zero=True,
                          virtual_lead=True,low_inverse=True,lazy_flags=True)
    assert circuit.num_qubits<=1004
    assert circuit.operation_statistics().num_toffoli_like==manifest['expected_static_Toffolis']
    manifest.update(status='BUILT_REQUIRES_NATIVE_VERIFICATION',Q=circuit.num_qubits,
                    static_Toffolis=circuit.operation_statistics().num_toffoli_like,
                    circuit=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    (HERE/'results/retained_dqi_build.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest,indent=2),flush=True)


if __name__=='__main__':main()
