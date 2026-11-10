"""Bind frozen development evidence, rebuilds, model margins and retained reference."""
import hashlib,json,math
from pathlib import Path
import kickmix as km
from build_experiment import build
H=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    exact=json.loads((H/'results/prototypes.json').read_text())
    matches=json.loads((H/'results/complete_model_matches.json').read_text())
    reference=json.loads((H/'reference/zero_head_affine_verify.json').read_text())
    assert reference['status']=='PASS' and reference['verified_shots']==200540928
    assert sha(H/'reference/zero_head_affine.kmb')==reference['circuit_sha256']
    assert exact['integer']['inputs']==3133452 and exact['integer']['min_margin_units']==19995185
    rows=[]
    for mode in ('borrowed','circular','walsh'):
        p=H/f'candidates/{mode}_best.kmb';record=json.loads((H/f'results/{mode}_best.json').read_text())
        record['path']=str(p.relative_to(H));(H/f'results/{mode}_best.json').write_text(json.dumps(record,indent=2)+'\n')
        match=matches[mode];assert sha(p)==match['sha256']
        c,rebuild=build(mode,*[record[k] for k in ('limit','keep','delay','product')])
        tmp=H/f'results/{mode}_rebuilt.kmb';c.to_file(str(tmp),format='kmb');assert sha(tmp)==sha(p)
        cstats=c.operation_statistics()
        assert cstats.num_toffoli_like==record['static_Toffolis']
        assert cstats.num_z_pow(min_clifford_hierarchy_level=4)==record['rotations']
        rows.append(dict(mode=mode,sha256=sha(p),Q=c.num_qubits,static_Toffolis=record['static_Toffolis'],
            sampled_mean_Toffolis=match['mean_Toffolis'],sample_shots=match['shots'],
            rotations=record['rotations'],synthesis_work=record['synthesis_work'],
            exact_model_margin_units=19995185,
            projection_with_sampled_mean=c.num_qubits*(match['mean_Toffolis']+record['synthesis_work']),
            byte_identical_rebuild=True,accepted=False))
    grids={}
    for mode in ('borrowed','circular','walsh'):
        data=[json.loads(s) for s in (H/f'results/{mode}_grid.jsonl').read_text().splitlines()]
        ok=[x for x in data if x['status']=='BUILT']
        grids[mode]=dict(attempts=len(data),built=len(ok),infeasible=len(data)-len(ok),minimum_compiled_Q=min(x['Q'] for x in ok))
    report=dict(status='DEVELOPMENT_AUDIT_PASS',selection='preserved incumbent',new_full_verifier_runs=0,
        reference=reference,reference_sha256=sha(H/'reference/zero_head_affine.kmb'),
        table_sha256=sha(H/'candidates/margin_curve8_t5.9e-08.npy'),
        reference_cache_sha256=sha(H.parent/'inv-phase/reference/exact_target.npz'),
        prototypes=rows,schedules=grids,
        source_sha256={str(p.relative_to(H)):sha(p) for p in sorted(H.glob('*.py'))})
    (H/'results/audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(rows,indent=2))

if __name__=='__main__':main()
