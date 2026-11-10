"""Audit experiment evidence and fingerprint the retained research artifacts."""
import hashlib,json
from pathlib import Path

root=Path(__file__).resolve().parent
reports={p.name:json.loads(p.read_text()) for p in (root/'results').glob('*.json')}
expected={'host_102400_2026102405_g3.json':0,'law_1000000_2026102405_g3.json':1,
          'law_1000000_2026102406_g3.json':0,'controls_4608_2026102423_g3.json':0}
for name,bad in expected.items():assert reports[name]['wrong_decodes']==bad
components=['matrix_transfer_n3.json','matrix_transfer_n18.json','matrix_coefficients.json',
            'completion_counter_b2.json','completion_counter_b3.json']
for name in components:
    assert reports[name]['status']=='PASS'
    assert all(c['status']=='PASS' for c in reports[name]['circuits'])
screen=reports['transfer_resource_screen.json']
assert screen['status']=='KILL_THIS_CONSTRUCTION'
assert screen['partial_projected_delta']>90000 and screen['optimistic_embedded_bridge_Q']>965
assert reports['whole_block_selector_screen.json']['optimistic_added_T']>90000
assert reports['serial_failure_diagnosis.json']['gaps']==[{'degree_stage':24,'quotient_degree':3}]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
sources=list(root.glob('*.py'))+list(root.glob('*.cpp'))+[root/'env.sh',root/'RESULTS.md']
inputs=[root/'../dqi-repair/parameters'/n for n in ['host_cases.bin','antilog.bin','log.bin','probs.txt']]
prior=[root/'../../4/codex'/n for n in ['build_cubic_compact_decoder.py','compact_product_decoder.py']]
report=dict(status='PASS',disposition='KILLED_CONCRETE_NATIVE_CONSTRUCTION_NO_REPLACEMENT',
            baseline_sha256_rechecked_in_session='4989da5cb633529e8851b9aa0c8c52abeb8990423a3a9463a70f673661a909c0',
            full_native_decoder_emitted=False,full_native_verifier_run=False,
            classical_wrong_decodes=expected,components=components,
            partial_projected_added_T=screen['partial_projected_delta'],
            source_sha256={str(p.relative_to(root)):sha(p) for p in sources},
            readonly_input_sha256={str(p.relative_to(root)):sha(p) for p in inputs+prior},
            component_sha256={str(p.relative_to(root)):sha(p) for p in (root/'candidates').glob('*.kmb')},
            evidence_sha256={str(p.relative_to(root)):sha(p) for p in (root/'results').glob('*.json') if p.name!='final_audit.json'})
(root/'results'/'final_audit.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if not k.endswith('sha256')}))
