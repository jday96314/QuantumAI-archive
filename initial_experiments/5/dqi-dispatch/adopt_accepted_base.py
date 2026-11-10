"""Record the accepted dqi-core base, leaving its builder and artifact read-only."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parent;CORE=ROOT/'../dqi-core'
assert (CORE/'RESULTS.md').read_text().startswith('ACCEPTED: trimmed decoder')
r=json.loads((CORE/'results/accepted_decoder.json').read_text());assert r['status']=='ACCEPTED'
assert r['Q']==961 and r['static_T']==994358
assert r['default']['verified_shots']==r['independent']['verified_shots']==102400
assert r['paired_comparison']['mean_T_reduction']==25596
h=hashlib.sha256()
with (CORE/r['circuit']).open('rb') as f:
    for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
assert h.hexdigest()==r['circuit_sha256']
report=dict(status='ADOPTED_READ_ONLY_ACCEPTED_BASE',builder=str((CORE/r['builder']).resolve()),
            circuit=str((CORE/r['circuit']).resolve()),sha256=h.hexdigest(),Q=r['Q'],static_T=r['static_T'],
            default_mean_T=r['default']['toffoli_count'],default_phase_errors=r['default']['phase_errors'],
            default_value_errors=r['default']['value_errors'],default_cost=r['default']['score'],
            acceptance_evidence=str((CORE/'results/accepted_decoder.json').resolve()),
            full_builder_executed_here=False,exception_correct_replacement_accepted_here=False)
(ROOT/'parameters'/'base_adoption.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
