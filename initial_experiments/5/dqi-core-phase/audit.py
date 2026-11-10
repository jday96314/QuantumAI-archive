"""Bind final research evidence to exact bytes; prior experiments are read-only."""
from pathlib import Path
import hashlib,json
import numpy as np
import kickmix as km
from dqi.verify_submission_dqi import optimal_dqi_weights

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parents[2]


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(4*1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def main():
    previous=HERE/'../dqi-gap';manifest=json.loads((previous/'parameters/best_candidate_manifest.json').read_text())
    baseline=previous/manifest['circuit'];baseline_hash=digest(baseline)
    assert baseline_hash==manifest['sha256']
    verifier=PROJECT/'from_host/verifiers/dqi/verify_submission_dqi.py'
    verifier_hash=digest(verifier);assert verifier_hash==manifest['host_verifier_sha256']
    preserved_sources={}
    for name,expected in manifest['source_sha256'].items():
        actual=digest(previous/name);assert actual==expected,name;preserved_sources[name]=actual
    weights=np.loadtxt(HERE/'parameters/host_weights.txt')
    assert np.allclose(weights,optimal_dqi_weights(4095,35,12,2016)**2,rtol=1e-15,atol=0)
    checks=[]
    def bind(report_name,artifact_name):
        r=json.loads((HERE/'results'/report_name).read_text());p=HERE/'candidates'/artifact_name
        actual=digest(p);assert actual==r['sha256'],(report_name,actual,r['sha256'])
        checks.append(dict(report=report_name,circuit=artifact_name,sha256=actual,status=r['status']))
    for suffix,stem in [('a','cfraction_n6_m31_allocated_scan'),('b','cfraction_n6_m31'),('c','cfraction_n6_m63'),('d','cfraction_n14_m127')]:
        report_stem=stem.replace('_allocated_scan','')
        bind(f'{report_stem}_core-phase-cfraction-{suffix}.json',stem+'.kmb')
    for d,salt in [(6,'core-phase-H-even-a'),(7,'core-phase-H-a'),(35,'core-phase-H-b'),(35,'core-phase-H-c')]:
        bind(f'h_degree{d}_{salt}.json',f'h_degree{d}_phase.kmb')
    for d,salt in [(7,'core-phase-H-multipoint-a'),(35,'core-phase-H-multipoint-b')]:
        bind(f'h_degree{d}_m9_{salt}.json',f'h_degree{d}_m9_phase.kmb')
    bind('inverse_only_point.json','inverse_only_point.kmb')
    for kind,k,salt in [('moment',3,'core-phase-recurrence-a'),('moment',35,'core-phase-recurrence-b'),('transport',3,'core-phase-recurrence-a')]:
        bind(f'{kind}_k{k}_{salt}.json',f'{kind}_k{k}_phase.kmb')
    for name,t in [('sylvester_t17_core-phase-sylvester-b.json',17),('sylvester_t3_restoration_core-phase-restoration-a.json',3)]:
        for r in json.loads((HERE/'results'/name).read_text()):
            p=HERE/'candidates'/f'sylvester_t{t}_{r["which"]}.kmb';actual=digest(p);assert actual==r['sha256'],name
            checks.append(dict(report=name,circuit=p.name,sha256=actual,status=r['status']))
    for kind in ['pair','H']:
        name=f'padding_{kind}_n35.json'
        for r in json.loads((HERE/'results'/name).read_text()):
            p=HERE/'candidates'/f'padding_{kind}_n35_{r["which"]}.kmb'
            c=km.Circuit.from_file(str(p),format='kmb')
            assert c.num_qubits==r['Q'] and c.operation_statistics().num_toffoli_like==r['T']
            checks.append(dict(report=name,circuit=p.name,sha256=digest(p),status=r['status']))
    algebra=json.loads((HERE/'results/algebra.json').read_text());assert algebra['status']=='PASS'
    budget=json.loads((HERE/'results/resource_model.json').read_text());assert budget['status']=='NO_FULL_CANDIDATE_PROMOTED'
    for label in ['a','b']:
        census=json.loads((HERE/f'results/census_complete_{label}.json').read_text())
        assert census['status']=='PASS_CLASSICAL_MODEL' and census['cf_wrong_supported_pairs']==0
    # Hash sources, native bytes and result JSONs. Large construction tapes
    # remain preserved; list their lengths rather than rereading them twice.
    files={};construction=[]
    for p in sorted(HERE.rglob('*')):
        if not p.is_file() or p.name in ['artifact_manifest.json','audit.json']:continue
        rel=str(p.relative_to(HERE))
        if p.suffix in ['.kmx','.npy']:
            construction.append(dict(path=rel,bytes=p.stat().st_size));continue
        if p.suffix in ['.py','.cpp','.sh','.md','.json','.txt','.kmb'] or p.name=='census':files[rel]=digest(p)
    imported=[
        'initial_experiments/3/import_published.py','initial_experiments/3/published_dqi.py',
        'initial_experiments/3/recursive_forney.py','initial_experiments/2/fast_dqi.py',
        'initial_experiments/2/tower_inverse.py','initial_experiments/4/codex/mobius_series.py',
        'initial_experiments/4/codex/series_inverse.py','initial_experiments/4/codex/low_workspace_inverse.py',
        'initial_experiments/4/codex/cubic_inverse_product_kernel.py',
        'initial_experiments/4/codex/cubic_trace_field.py','initial_experiments/4/codex/controlled_linear.py',
        'initial_experiments/4/codex/coprime_padding.py','initial_experiments/4/codex/quantum_coprime_padding.py',
        'initial_experiments/4/codex/translate_prefix',
        'initial_experiments/3/vendor/uncompute/tasks/dqi-reed-solomon/circuits/discovered/circuit.py']
    dependencies={name:digest(PROJECT/name) for name in imported}
    out=dict(status='PASS',baseline_sha256=baseline_hash,host_verifier_sha256=verifier_hash,
             unchanged_gap_sources=preserved_sources,bound_component_checks=checks,
             native_full_benchmark_runs=0,verified_task_cost_improvement=False,
             scope='final component byte bindings, unchanged baseline/verifier/source checks, host weights and model consistency; no full-task acceptance')
    (HERE/'results/audit.json').write_text(json.dumps(out,indent=2)+'\n')
    (HERE/'results/artifact_manifest.json').write_text(json.dumps(dict(files=files,dependencies=dependencies,construction_tapes=construction),indent=2)+'\n')
    print(json.dumps(dict(status='PASS',bound_component_checks=len(checks),hashed_files=len(files),baseline_sha256=baseline_hash,
                          native_full_benchmark_runs=0,verified_task_cost_improvement=False)),flush=True)


if __name__=='__main__':main()
