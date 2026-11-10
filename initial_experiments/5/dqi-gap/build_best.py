"""Build the resource-selected candidate, frozen before any native samples."""
import hashlib,json
from pathlib import Path
import hybrid_decoder as decoder
from blocked_monic_euclid import schedule


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8<<20),b''):h.update(block)
    return h.hexdigest()


selection=json.loads(Path('parameters/best_schedule.json').read_text());counts=selection['counts']
assert sum(counts)==34
sources=['build_best.py','hybrid_decoder.py','fused_gap_block.py','gap_mobius.py','onepass_gap.py','laurent_cofactor.py','gap_serial.py','quantum_control.py','serial_opening.py','recycled_scan.py','allocation_profile.py','env.sh','translate_prefix']
manifest=dict(scope='candidate fixed by exact resource profiles before any hash-derived sample',resource_selection=selection,post_selection_optimization='reuse one history bit for previous-shift control and the completion bit for padding control; unchanged emitted T',source_sha256={p:digest(p) for p in sources},host_verifier_sha256=digest('../../../from_host/verifiers/dqi/verify_submission_dqi.py'))
Path('parameters/best_candidate_pre_native.json').write_text(json.dumps(manifest,indent=2)+'\n')


def frames(S,width):
    bp=[];b=[S[-2]];c=S[:-2];p=S[-1];k=1;stages=[]
    for count in counts:
        stages.append((k,count,bp,b,c,p));_,(bp,b)=schedule(bp,b,list(reversed(c))[:2*count]);c=c[:-2*count];k+=count
    assert k==35 and not c
    return stages,(bp,b,c,p)


decoder.frames=frames;decoder.build.slot_counts=counts
out=Path('candidates/gap_best_n70_m4095.kmb')
c=decoder.build(70,3,100,4095,out,linear_padding=True,virtual_zero=True,virtual_lead=True,low_inverse=True,lazy_flags=True)
manifest.update(circuit=str(out),sha256=digest(out),Q=c.num_qubits,static_Toffolis=c.operation_statistics().num_toffoli_like)
Path('parameters/best_candidate_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps({k:manifest[k] for k in ['circuit','sha256','Q','static_Toffolis']}),flush=True)
