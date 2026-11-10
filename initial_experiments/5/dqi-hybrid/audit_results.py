"""Bind research evidence to exact control and component bytes.

No inherited artifact is modified. This is a research audit, not a replacement
for the unchanged host's full DQI verifier.
"""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from dqi.verify_submission_dqi import clopper_pearson_upper_bound,repetition_overhead
from profile_selected import HERE,R3,R4,TAG,digest


def main():
    summary=json.loads((R4/'submission.summary.json').read_text())
    assert digest(R4/'submission.zip')==summary['sha256']
    sources={}
    for manifest,field in [(R4/'parameters/compact_product_slim_source_manifest.json','sources'),
                           (R4/'parameters/cubic_trace_product_n70_m4095_source_manifest.json','source_sha256')]:
        data=json.loads(manifest.read_text())
        for path,want in data[field].items():
            source=R4/path;actual=digest(source)
            assert actual==want,(source,actual,want)
            sources[str(source)]=actual
    controls=[]
    for family,path,result in [('round3',R3/f'candidates/{TAG}_full.kmb',R3/f'results/{TAG}_verify.json'),
                               ('round4',R4/'candidates/cubic_trace_product_n70_m4095.kmb',R4/'results/cubic_trace_product_n70_m4095_verify.json')]:
        data=json.loads(result.read_text())
        assert digest(path)==data['circuit_sha256'],path
        assert data['verified_shots']==102400 and data['status']=='PASS'
        assert abs(data['score']-data['repetition_factor']*data['num_qubits']**2*data['toffoli_count'])<.01
        (HERE/f'results/control_{family}_native.json').write_bytes(result.read_bytes())
        controls.append(dict(family=family,circuit=str(path),circuit_sha256=data['circuit_sha256'],native_result=str(result),
                             Q=data['num_qubits'],mean_T=data['toffoli_count'],R=data['repetition_factor'],cost=data['score']))
    for family,tape in [('round3',R3/f'parameters/{TAG}_tape.npy'),('round4',R4/'candidates/cubic_trace_product_n70_m4095.npy')]:
        data=json.loads((HERE/f'results/{family}_liveness.json').read_text())
        assert data['tape_comparison']['equal']
        rows=np.load(tape,mmap_mode='r')
        with (HERE/f'results/{family}_timeline.csv').open('w') as f:
            writer=csv.DictWriter(f,fieldnames=['region','begin_op','end_op','Toffolis','peak_live','end_live','physical_peak_roles','logical_layout'])
            writer.writeheader()
            total=0
            for row in data['timeline']:
                a,b=row['begin'],row['end'];T=int(np.count_nonzero(np.isin(rows[a:b,0]&255,[6,7])))
                total+=T
                writer.writerow(dict(region=row['region'],begin_op=a,end_op=b,Toffolis=T,peak_live=row['peak']['live'],end_live=row['end_state']['live'],
                                     physical_peak_roles=json.dumps(row['peak']['physical_by_role'],sort_keys=True),logical_layout=json.dumps(row['logical_layout'],sort_keys=True)))
            assert total==data['core_forward_T']+data['core_reverse_T']
    components=[]
    for path in sorted((HERE/'results').glob('consume_l*_independent-*.json')):
        data=json.loads(path.read_text());assert data['status']=='PASS' and data['shots_per_circuit']==102400
        for row in data['circuits']:
            circuit=HERE/f"components/consume_l{row['L']}_{row['which']}.kmb"
            assert digest(circuit)==row['circuit_sha256']
            assert digest(HERE/'consume_component.py')==row['source_sha256']
            assert row['mean_native_T']==row['T']
            assert row['value_errors']==row['scratch_errors']==row['phase_errors']==0
        components.append(dict(result=str(path.relative_to(HERE)),shots_per_circuit=data['shots_per_circuit'],L=data['L'],salt=data['salt'],circuits=len(data['circuits'])))
    law=json.loads((HERE/'results/law_law-independent-a.json').read_text())
    assert law['status']=='PASS' and law['shots']==1000000 and law['value_roundtrip_errors']==0
    breakdown=json.loads((HERE/'results/native_breakdown.json').read_text())
    assert all(row['core_plus_scan_matches_selected'] for row in breakdown)
    minimum_R=repetition_overhead(clopper_pearson_upper_bound(0,102400))
    audit=dict(status='PASS',outcome='NO_REPLACEMENT; research stopped under the register-design scope limit',
        control_archive_sha256=summary['sha256'],preserved_combined_score=summary['overall_score'],controls=controls,
        selected_source_hashes=sources,component_checks=components,
        total_component_native_shots=sum(x['shots_per_circuit']*x['circuits'] for x in components),
        law_component_cases=law['shots'],
        price_examples=dict(scope='hypothetical arithmetic examples, not verified candidates',
            target_Q1050_T600000_R123=1.23*1050**2*600000,
            minimum_R_at_102400_zero_errors=minimum_R,
            published_Q1837_T442000_with_optimistic_host_pricing=minimum_R*1837**2*442000),
        limitations=['No new complete DQI circuit or task-cost win.',
                    'Component tests are independent seeded native tests, not host DQI verifier acceptance.',
                    'The million-case check samples host-law syndromes but independently varies locator operands; it is not a million full decoder cases.',
                    'The zero-pivot information bound is local and conditional on fixing all other operands; no general width lower bound is claimed.'])
    (HERE/'results/audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    manifest={str(path.relative_to(HERE)):digest(path) for path in sorted(HERE.rglob('*'))
              if path.is_file() and path.suffix!='.log' and path.name!='artifact_manifest.json'}
    (HERE/'results/artifact_manifest.json').write_text(json.dumps(dict(
        scope='all experiment files except transient logs and this manifest; no inherited artifact is copied or modified',
        sha256=manifest),indent=2)+'\n')
    print(json.dumps({key:audit[key] for key in ['status','outcome','preserved_combined_score','total_component_native_shots','law_component_cases']}),flush=True)


if __name__=='__main__':main()
