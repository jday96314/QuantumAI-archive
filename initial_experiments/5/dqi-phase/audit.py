"""Bind component evidence to source/artifact hashes; no full acceptance claim."""
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent
OLD=(ROOT/'../../4/codex').resolve()


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        while chunk:=f.read(8<<20):h.update(chunk)
    return h.hexdigest()


def main():
    selected=json.loads((OLD/'parameters/cubic_trace_product_n70_m4095_source_manifest.json').read_text())
    core=json.loads((OLD/'parameters/compact_product_slim_source_manifest.json').read_text())
    checked={}
    for name,digest in [*selected['source_sha256'].items(),*core['sources'].items()]:
        path=(OLD/name).resolve()
        assert sha(path)==digest,('read-only selected source changed',name)
        checked[str(path)]=digest
    evidence={}
    for name in ['baseline_audit','oracle','control_identity','point_kernels','psi_schedule','modular_schedule','fast_modular','mask_scan']:
        p=ROOT/'results'/f'{name}.json';data=json.loads(p.read_text())
        reports=data if isinstance(data,list) else [data]
        assert all(row['status']=='PASS' for row in reports),name
        evidence[p.name]=sha(p)
    exceptional=json.loads((ROOT/'results/psi_exceptions.json').read_text())
    assert exceptional['all_sources_restored'] and exceptional['all_scratch_clean'] and exceptional['all_phases_clean']
    evidence['psi_exceptions.json']=sha(ROOT/'results/psi_exceptions.json')
    budget=json.loads((ROOT/'results/budget.json').read_text())
    assert budget['status']=='NO_FULL_CANDIDATE_PROMOTED'
    assert budget['control']['circuit_sha256']==selected['circuit_sha256']
    evidence['budget.json']=sha(ROOT/'results/budget.json')
    artifacts={str(p.relative_to(ROOT)):dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted((ROOT/'candidates').glob('*.kmb'))}
    local={str(p.relative_to(ROOT)):sha(p) for p in [*sorted(ROOT.glob('*.py')),ROOT/'env.sh',ROOT/'RESULTS.md',ROOT/'parameters/cubic_joint_quadratic_16.txt',*sorted((ROOT/'reference').glob('*'))]}
    # Record additional read-only direct dependencies beyond the selected
    # manifest (the joint cover is only a component control).
    extras=[OLD/'joint_cubic_inverse_product_kernel.py',OLD/'monic_compact_euclid.py',OLD/'compact_euclid.py',
            OLD/'inverse_product_kernel.py',ROOT/'../../../from_host/verifiers/dqi/verify_submission_dqi.py']
    for p in extras:checked[str(p.resolve())]=sha(p)
    manifest=dict(scope='round-5 phase components, no full-task acceptance',local_sources=local,read_only_dependencies=checked,
                  artifacts=artifacts,evidence_sha256=evidence,selected_circuit_sha256=selected['circuit_sha256'])
    (ROOT/'parameters/source_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    out=dict(status='PASS',scope='source/evidence/artifact integrity audit; not native full-decoder verification',
             artifact_count=len(artifacts),read_only_source_pins=len(checked),evidence_count=len(evidence),
             source_manifest_sha256=sha(ROOT/'parameters/source_manifest.json'),
             selected_circuit_sha256=selected['circuit_sha256'],new_full_candidate_accepted=False,
             combined_reference=9.367722153134604,selected_DQI_cost=1598204229384.9265)
    (ROOT/'results/audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)


if __name__=='__main__':main()
