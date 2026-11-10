"""Rebuild frozen bytes and audit immutable references and native evidence."""
import argparse
import kickmix as km
from common import *
from build import build as linear_build
from quadratic import build as nonlinear_build,fit
from hybrid import search
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas
from dfthcsos.scoreinfra.circuits import load_circuits,submission_sha256

def rebuild(refit=False):
    final=HERE/'candidates/final';meta=json.loads((HERE/'results/final_build.json').read_text());trans=json.loads((HERE/'results/final_transform.json').read_text())
    z=dict(np.load(INC/'submission.npz'));original=z['angle_words'];order=np.array(trans['order']);matrix=trans['matrix']
    if refit:
        oo,mm,history=search(original)
        assert np.array_equal(oo,order) and mm==matrix and history==trans['features']
    z['angle_words']=original[order];two=order[76:]-76;z['U']=z['U'][:,two,:];z['W']=z['W'][:,two,:]
    out=HERE/'tmp/rebuild';out.mkdir(parents=True,exist_ok=True);np.savez(out/'submission.npz',**z)
    assert sha(out/'submission.npz')==sha(final/'submission.npz')
    rows=[]
    for k,d in enumerate(batch_deltas(z['angle_words'],4)):
        if meta['rows'][k]['compiler']=='linear':c,info=linear_build(d,select=2,lowbits=7,mix=matrix)
        else:
            frozen_model=meta['rows'][k]['info']['model']
            c,info=nonlinear_build(d,k,model=None if refit else frozen_model)
            if refit:assert info['model']==frozen_model
        path=out/f'angle_k_{k}.kmb';c.to_file(str(path),format='kmb');assert sha(path)==sha(final/path.name)
        rows.append(dict(k=k,sha256=sha(path),identical=True))
    result=dict(status='PASS',refitted=refit,rows=rows,npz_identical=True);save('final_refit_rebuild.json' if refit else 'final_rebuild.json',result);print(json.dumps(result))

def audit():
    final=HERE/'candidates/final';freeze=json.loads((HERE/'results/final_freeze.json').read_text())
    assert {f.name:sha(f) for f in sorted(final.iterdir())}==freeze['files']
    baseline=json.loads((HERE/'results/baseline_reproduction.json').read_text())
    assert all(r['identical_bytes'] for r in baseline['rows'])
    assert all(sha(INC/f'angle_k_{r["k"]}.kmb')==r['sha256'] for r in baseline['rows'])
    assert sha(INC/'submission.npz')==baseline['npz_sha256']
    words=np.load(final/'submission.npz')['angle_words'];cs=load_circuits(final)
    digest=submission_sha256(cs,final/'submission.npz',np.load(final/'submission.npz')['tree_edges'])
    records=[]
    for letter in ['A','B']:
        r=json.loads((HERE/'results'/f'final_verify_{letter}.json').read_text())
        assert r['status']=='PASS' and r['verified_shots']==795*64
        assert r['submission_sha256']==digest
        assert r['score']<BASE_COST and r['num_qubits']==454 and r['toffoli_count']==31608
        assert r['resources']['angles_qubits']<=128
        assert r['resources']['temporary_qubits']==337
        records.append(r)
    checks=json.loads((HERE/'results/final_check.json').read_text());assert checks['status']=='PASS'
    assert checks['exact_npz_sha256']==sha(final/'submission.npz')
    for k,c in enumerate(cs):
        assert c.num_qubits<=128
        assert checks['rows'][k]['sha256']==sha(final/f'angle_k_{k}.kmb')
    paths=list(HERE.glob('*.py'))+[HERE/'env.sh']
    paths += [PROJECT/p for p in [
        'initial_experiments/4/codex/build_signed_gauge_chem.py',
        'initial_experiments/4/codex/build_selectswap_chem.py',
        'initial_experiments/4/codex/monomial_unlookup.py',
        'initial_experiments/5/inv-dyadic/stream_erase.py',
        'initial_experiments/5/chem-qfit/oracles.py',
        'from_host/verifiers/dfthcsos/verify.py',
        'from_host/verifiers/dfthcsos/scoreinfra/circuits.py',
        'from_host/verifiers/dfthcsos/scoreinfra/costs_prepare.py',
        'from_host/verifiers/dfthcsos/scoreinfra/dfthc_math.py']]
    save('source_manifest.json',{str(p.relative_to(PROJECT)):sha(p) for p in sorted(paths)})
    result=dict(status='PASS',submission_sha256=digest,independent_full_native_runs=2,
                shots_per_run=795*64,total_native_shots=2*795*64,
                incumbent_unchanged=True,candidate_files_match_freeze=True,
                score=records[0]['score'],relative_cost_reduction=1-records[0]['score']/BASE_COST,
                multiplicative_improvement=BASE_COST/records[0]['score'])
    save('audit.json',result);print(json.dumps(result))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--rebuild',action='store_true');p.add_argument('--refit',action='store_true');p.add_argument('--audit',action='store_true');a=p.parse_args()
    if a.rebuild or a.refit:rebuild(a.refit)
    if a.audit:audit()
