"""Build the source-preserving retained-inversion decoder after small validation."""
import hashlib
import json
from pathlib import Path

import dqi_budget_probe as probe
import low_workspace_inverse
from import_published import load
from safe_retained_inverse import bind

HERE=Path(__file__).resolve().parent


if __name__=='__main__':
    small=json.loads((HERE/'results/safe_small_verify.json').read_text())
    assert small['status']=='PASS' and small['value_errors']==small['phase_errors']==0
    profile=json.loads((HERE/'results/safe_inverse_profile.json').read_text())
    assert profile['Q']<=1004 and profile['T']<801388
    manifest=dict(status='FROZEN_BEFORE_NATIVE_VERIFICATION',
                  design='retained inverse only when the source is preserved; unchanged accepted gap schedule, trimmed middle products, and phase scan',
                  profile={k:v for k,v in profile.items() if k!='rows'},
                  expected_static_Toffolis=profile['T']+371150,
                  independent_salt='agent20-safe-retained-dqi-independent-20261025',
                  sources={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [
                      Path(__file__),HERE/'safe_retained_inverse.py',HERE/'dqi_budget_probe.py',
                      probe.ROUND5/'dqi-core/components/trimmed_middle.py',
                      probe.GAP/'hybrid_decoder.py',probe.GAP/'parameters/best_schedule.json',
                      probe.GAP/'serial_opening.py',probe.GAP/'gap_serial.py']})
    (HERE/'parameters/safe_dqi_pre_native.json').write_text(json.dumps(manifest,indent=2)+'\n')
    decoder=probe.decoder
    decoder.load=load;decoder.frames=probe.frames
    decoder.BoundedArithmetic=probe.trimmed.TrimmedArithmetic
    decoder.build.profile_only=False;decoder.build.slot_counts=probe.counts
    low_workspace_inverse.bind=bind
    path=HERE/'candidates/safe_retained_dqi.kmb'
    circuit=decoder.build(70,3,100,4095,path,linear_padding=True,virtual_zero=True,
                          virtual_lead=True,low_inverse=True,lazy_flags=True)
    assert circuit.num_qubits<=1004
    assert circuit.operation_statistics().num_toffoli_like==manifest['expected_static_Toffolis']
    manifest.update(status='BUILT_REQUIRES_NATIVE_VERIFICATION',Q=circuit.num_qubits,
                    static_Toffolis=circuit.operation_statistics().num_toffoli_like,
                    circuit=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    (HERE/'results/safe_dqi_build.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest,indent=2),flush=True)
