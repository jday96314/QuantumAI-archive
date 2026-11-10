"""Collect measured frontiers and bind native evidence to frozen candidates."""
import argparse,csv,hashlib,json
from pathlib import Path
from build import HERE,OLD
from reproduce import SPECS

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def collect(require_native=True):
    approximations=[]
    for name,*spec in SPECS:
        m=json.loads((HERE/'models'/f'{name}.json').read_text())
        c=json.loads((HERE/'results'/f'{name}_certificate.json').read_text())
        approximations.append(dict(model=name,mode=m['mode'],widths='/'.join(map(str,m['widths'])),
          leaves=len(m['leaves']),depth_histogram=json.dumps(c['depth_histogram'],sort_keys=True),
          status=c['status'],worst_error_over_tolerance=c['worst_error_over_tolerance'],
          min_margin_units=c['min_margin_units'],synthesis_work=c['synthesis_toffoli_equivalent'],
          accumulator_bound_pass=c.get('signed_accumulator_magnitudes_fit'),
          product_sign_violations=json.dumps(c.get('product_sign_violations'))))
    with (HERE/'results/approximation_frontier.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=approximations[0].keys());writer.writeheader();writer.writerows(approximations)
    resources=[]
    for path in sorted((HERE/'results').glob('*_build.json')):
        d=json.loads(path.read_text())
        if 'Q' not in d or 'certificate' not in d:continue
        cert=d['certificate'];q=d['Q'];t=d['T_static']
        resources.append(dict(report=path.name,model=Path(d['model']).stem,Q=q,T_static=t,
          synthesis_work=cert['synthesis_toffoli_equivalent'],
          compiled_model_upper_cost=q*(t+cert['synthesis_toffoli_equivalent']),
          phase_entry_live=next((x['live'] for x in d['stages'] if x['stage']=='phase'),None)))
    with (HERE/'results/resource_frontier.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=resources[0].keys());writer.writeheader();writer.writerows(resources)
    schedules=[]
    for path in sorted((HERE/'results').glob('*_schedules.jsonl')):
        rows=[json.loads(s) for s in path.read_text().splitlines()]
        schedules.append(dict(model=path.name.removesuffix('_schedules.jsonl'),
          feasible=sum(r.get('status')!='INFEASIBLE' for r in rows),infeasible=sum(r.get('status')=='INFEASIBLE' for r in rows)))
    audit=dict(status='PENDING',numerical_models=len(approximations),
               feasible_schedules=sum(r['feasible'] for r in schedules),schedules=schedules,native={})
    for name in ['final_shared','final_local']:
        path=HERE/'candidates'/f'{name}.kmb';result=HERE/'results'/f'{name}_verify.json'
        if not result.exists() or not result.stat().st_size:
            if require_native:raise AssertionError(f'{name} native verifier incomplete')
            continue
        d=json.loads(result.read_text());assert d['status']=='PASS'
        assert d['verified_shots']==200540928
        assert d['circuit_sha256']==digest(path)
        build=json.loads((HERE/'results'/f'{name}_build.json').read_text())
        assert d['num_qubits']==build['Q'] and d['max_toffoli_count']==build['T_static']
        cert=json.loads((HERE/'results'/f'{Path(build["model"]).stem}_certificate.json').read_text())
        assert abs(d['worst_error_over_tolerance']-cert['worst_error_over_tolerance'])<1e-12
        assert abs(d['rotation_eps']/cert['rotation_eps']-1)<1e-12
        rebuilt=HERE/'candidates'/f'{name}_rebuilt.kmb'
        assert rebuilt.read_bytes()==path.read_bytes()
        audit['native'][name]=d
    old=HERE/'reference/incumbent.kmb';oldresult=json.loads((HERE/'reference/incumbent_verify.json').read_text())
    assert digest(old)==oldresult['circuit_sha256']==digest(OLD/'candidates/zero_head_affine.kmb')
    audit['incumbent_unchanged']=True
    audit['native_shots_total']=sum(d['verified_shots'] for d in audit['native'].values())
    if len(audit['native'])==2:audit['status']='PASS'
    (HERE/'results/audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    sources=[*HERE.glob('*.py'),HERE/'fit.cpp',HERE/'env.sh',
       OLD/'binary_centered_invsqrt.py',OLD/'booth_adder.py',OLD/'short_adder.py',OLD/'monomial_unlookup.py',
       OLD/'compact_wires.cpp',OLD/'compact_wires',HERE.parents[2]/'from_host/verifiers/invsqrt/verify_submission_invsqrt.py',
       HERE.parents[2]/'from_host/verifiers/invsqrt/nr_reference_circuit.py']
    manifest=dict(sources={str(p):digest(p) for p in sources},
                  models={p.name:digest(p) for p in (HERE/'models').glob('*.json')},
                  final_candidates={name:digest(HERE/'candidates'/f'{name}.kmb') for name in ['final_shared','final_local']})
    (HERE/'results/source_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(audit),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--pending',action='store_true');a=p.parse_args();collect(not a.pending)
