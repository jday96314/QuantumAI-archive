"""Native artifact/count audit and read-only incumbent preservation check."""
import ast, hashlib, json
from collections import Counter
from pathlib import Path
import kickmix as km
ROOT=Path(__file__).resolve().parent
PROJECT=ROOT.parents[2]

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        while block:=f.read(8<<20):h.update(block)
    return h.hexdigest()

def main():
    cases=[]
    for path in sorted((ROOT/'results').glob('*_r*.json')):
        r=json.loads(path.read_text())
        if 'kmb' not in r:continue
        assert r['status']=='PASS' and r['mismatching_bits']==r['scratch_bits']==r['phase_errors']==0
        p=ROOT/r['kmb'];assert sha(p)==r['sha256']
        c=km.Circuit.from_file(str(p),format='kmb')
        assert c.num_qubits==r['Q']
        assert c.operation_statistics().num_toffoli_like==r['T']
        cases.append(dict(path=str(path.relative_to(ROOT)),shots=r['shots'],Q=r['Q'],T=r['T'],sha256=r['sha256']))
        del c
    assert len(cases)==10
    upstream=PROJECT/'initial_experiments/3/vendor/uncompute/tasks/dqi-reed-solomon/circuits/discovered/liveness.json'
    live=json.loads(upstream.read_text());peak=live['totals']['peak_qubits']
    t=next(t for t,q in live['skyline'] if q==peak)
    alive=[v for v in live['intervals'] if (v['born'] or 0)<=t<(v['died'] if v['died'] is not None else float('inf'))]
    assert len({v['q'] for v in alive})==peak
    partition=dict(T=t,live_intervals=len(alive),live_unique_qubits=len({v['q'] for v in alive}),owner_counts=Counter('/'.join(v['sets']) for v in alive),source=str(upstream.relative_to(PROJECT)),method='interval sample at first peak; same-Toffoli boundary ordering is unavailable; owner sum cross-checks exact published peak')
    (ROOT/'results/published_peak.json').write_text(json.dumps(partition,indent=2)+'\n')
    identities=json.loads((ROOT/'results/fprime_identity.json').read_text())
    assert identities['status']=='PASS' and all(r['failures']==0 for r in identities['runs'])
    for p in ROOT.glob('*.py'):ast.parse(p.read_text())
    inputs=[PROJECT/'initial_experiments/4/REPORT.md',PROJECT/'initial_experiments/4/codex/CONTINUATION.md',
            PROJECT/'from_host/verifiers/dqi/verify_submission_dqi.py',PROJECT/'initial_experiments/3/import_published.py',
            PROJECT/'initial_experiments/3/published_dqi.py',PROJECT/'initial_experiments/4/codex/mobius_series.py',
            PROJECT/'initial_experiments/4/codex/series_inverse.py',PROJECT/'initial_experiments/4/codex/low_workspace_inverse.py',
            PROJECT/'initial_experiments/4/codex/translate_prefix',PROJECT/'initial_experiments/5/dqi-phase/oracle.py',
            PROJECT/'initial_experiments/5/dqi-phase/fast_modular.py',PROJECT/'initial_experiments/5/dqi-phase/modular_schedule.py',
            PROJECT/'initial_experiments/3/vendor/uncompute/tasks/dqi-reed-solomon/circuits/discovered/circuit.py',
            PROJECT/'initial_experiments/3/vendor/uncompute/tasks/dqi-reed-solomon/circuits/discovered/constants.json',
            PROJECT/'initial_experiments/3/vendor/uncompute/tasks/dqi-reed-solomon/circuits/discovered/liveness.json',
            PROJECT/'initial_experiments/5/dqi-hybrid/results/round3_liveness.json',PROJECT/'initial_experiments/5/dqi-hybrid/results/round4_liveness.json']
    for family in ['dqi-repair','dqi-phase','dqi-hybrid','dqi-core','dqi-dispatch']:
        p=PROJECT/'initial_experiments/5'/family/'RESULTS.md';inputs.append(p)
        snap=ROOT/'reference'/f'{family}_RESULTS_read.md'
        if not snap.exists():snap.write_bytes(p.read_bytes())
    reference=PROJECT/'initial_experiments/4/codex/candidates/cubic_trace_product_n70_m4095.kmb'
    digest=sha(reference)
    assert digest=='4989da5cb633529e8851b9aa0c8c52abeb8990423a3a9463a70f673661a909c0'
    manifest=dict(vendor_revision='64b5a13aa249f1360846657392e6c66817b26d0f',incumbent=dict(path=str(reference.relative_to(PROJECT)),sha256=digest),read_only_sources={str(p.relative_to(PROJECT)):sha(p) for p in inputs},external_snapshots={p.name:sha(p) for p in (ROOT/'reference').iterdir() if p.is_file()})
    (ROOT/'results/source_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    excluded={'results/audit.json','results/artifact_manifest.json'}
    artifacts={str(p.relative_to(ROOT)):sha(p) for p in sorted(ROOT.rglob('*')) if p.is_file() and str(p.relative_to(ROOT)) not in excluded}
    (ROOT/'results/artifact_manifest.json').write_text(json.dumps(artifacts,indent=2)+'\n')
    report=dict(status='PASS',native_circuits=len(cases),final_suite_shots=sum(c['shots'] for c in cases),native_checks=cases,identity_runs=identities['runs'],incumbent_sha256=digest,incumbent_preserved=True,artifact_count=len(artifacts),full_DQI_verifications=0,replacement_accepted=False,scope='exact KMB counts/digests + retained native component results; no full decoder or native failure-rate claim')
    (ROOT/'results/audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['native_checks','identity_runs']}),flush=True)
if __name__=='__main__':main()
