"""Audit exact artifacts, full native evidence, and scoped dictionary bounds."""
from census_common import *
from scipy.linalg import qr,eigh
from dfthcsos.scoreinfra.circuits import load_circuits,submission_sha256

def main():
    original=sha(BANK/'submission.npz')
    assert original==json.loads((HERE/'results/calibration.json').read_text())['source_npz_sha256']
    out=HERE/'candidates/structured';z=dict(np.load(out/'submission.npz'));bank=dict(np.load(BANK/'submission.npz'))
    transform=json.loads((HERE/'results/structured_transform.json').read_text());order=np.array(transform['order'])
    assert np.array_equal(z['angle_words'],bank['angle_words'][order])
    assert np.array_equal(z['U'],bank['U'][:,order[76:]-76,:])
    assert np.array_equal(z['W'],bank['W'][:,order[76:]-76,:])
    for key in ['H2_I','H1_I','Hsym_lower','tree_edges']:assert np.array_equal(z[key],bank[key])
    build=json.loads((HERE/'results/structured_build.json').read_text())
    files={p.name:sha(p) for p in sorted(out.iterdir())};assert files==build['files']
    cs=load_circuits(out);digest=submission_sha256(cs,out/'submission.npz',z['tree_edges'])
    assert digest==build['digest']
    checks=json.loads((HERE/'results/structured_components.json').read_text());assert checks['status']=='PASS'
    native=[]
    for file in ['structured_native.json','structured_native_independent.json']:
        r=json.loads((HERE/'results'/file).read_text());assert r['status']=='PASS'
        assert r['submission_sha256']==digest and r['verified_shots']==795*64
        assert r['score']==r['lambda_eff']*r['toffoli_count']*r['num_qubits']
        assert r['score']<BASE_COST and r['resources']['angles_qubits']<=128
        native.append(r)
    assert native[0]['score']==native[1]['score']
    research=[]
    for name in ['continuation703','growth511','growth671']:
        r=json.loads((HERE/'results'/f'{name}_native.json').read_text())
        path=HERE/'candidates'/name;zz=np.load(path/'submission.npz')
        dd=submission_sha256(load_circuits(path),path/'submission.npz',zz['tree_edges'])
        assert r['submission_sha256']==dd and r['recovered_words_match_saved']
        assert r['verified_shots']==len(zz['angle_words'])*64
        assert r['status']==('PASS' if passes_accuracy(r['metrics']) else 'REJECT')
        if r['status']=='PASS':assert r['score']>=native[0]['score'],'New winner needs independent native acceptance before selection'
        research.append(dict(candidate=name,status=r['status'],metrics=r['metrics'],
            resources=r['resources'],conditional_objective=r['conditional_objective'],
            digest=dd,verified_shots=r['verified_shots'],reason=r['rejection_reason']))
    # Independent CPU SciPy QR/eigensolve, using exact native-decoded vectors.
    m=Metric(z['tree_edges'],device='cpu');G=m.G.numpy();bounds=[]
    # A dictionary-independent obstruction: every native model is F.T@F,
    # PSD and rank <= R*C. Retain the full negative native spectrum.
    global_ev=eigh(G,eigvals_only=True)
    def global_bound(D):
        tail=global_ev.copy()
        if D:tail[-D:]=np.minimum(tail[-D:],0)
        return float(np.linalg.norm(tail))
    minimum_rank=next(D for D in range(len(global_ev)+1) if global_bound(D)<=.020)
    rank_rows=[dict(RC=D,error2_lower_bound=global_bound(D)) for D in sorted({224,256,280,minimum_rank-1,minimum_rank,300,350})]
    h=heuristic_costs(76,1,719,minimum_rank,8)
    global_rank=dict(minimum_RC_for_20mHa=minimum_rank,rows=rank_rows,
        proof='For native M=F.T@F with rank(M)<=R*C, PSD rank-D approximation error is sqrt(sum negative eigenvalues squared plus discarded positive eigenvalues squared).',
        conditional_R1_B719_minimum_C=dict(C=minimum_rank,T_J7=h['toffoli']+10604,Q_J7=h['qubits']+56,
            cost_at_banked_lambda=(h['toffoli']+10604)*(h['qubits']+56)*BASE_LAMBDA),
        caveat='Native float64 spectrum; analytic PSD/rank theorem, numerical eigenvalues, no interval arithmetic')
    write('global_rank.json',global_rank)
    for file in ['repair_B703_extended.npz','growth511_sequential_residual_B511.npz',
                 'growth511_sequential_reinsertion_B511.npz','growth671_residual_B671.npz',
                 'growth671_reinsertion_B671.npz']:
        zz,uu,ww,tree,words=load(HERE/'parameters'/file)
        decoded=reconstruct_u_vectors(words/128,tree).reshape(uu.shape)
        P=m.feats(m.tensor(decoded))[0].numpy();Q,R=qr(P.T,mode='economic')
        H=Q.T@G@Q;ev=eigh((H+H.T)/2,eigvals_only=True)
        outside=np.linalg.norm(G-Q@H@Q.T)
        tail=ev.copy();C=ww.shape[-1];tail[-C:]=np.minimum(tail[-C:],0)
        opt=float(np.sqrt(outside**2+np.sum(tail**2)))
        X=ww[0].T@P;actual=float(np.linalg.norm(X.T@X-G))
        assert actual>=opt-1e-10
        bounds.append(dict(file=file,B=uu.shape[1],C=C,span_bound=float(outside),
            rank_psd_bound=opt,saved_error2=actual,
            coefficient_rescue_excluded=opt>.020+1e-8))
    write('independent_bounds.json',dict(rows=bounds,method='CPU SciPy QR and symmetric eigensolve; full signed native target',
        caveat='Floating-point certificates with explicit 1e-8 Ha exclusion slack, not interval arithmetic; fixed words only'))
    selected=native[0]
    selection=dict(status='SELECTED_NATIVE_COST_IMPROVEMENT',candidate=str(out),score=selected['score'],
        previous_score=BASE_COST,previous_candidate=str(BANK),improvement=BASE_COST/selected['score'],target=TARGET,
        target_achieved=selected['score']<=TARGET,submission_sha256=digest,
        full_native_runs=2,shots_per_run=795*64,native_results=['results/structured_native.json','results/structured_native_independent.json'],
        metrics=selected['metrics'],resources=selected['resources'],files=files,
        previous_external_selection_unchanged=True)
    write('../selection.json',selection)
    source_paths=list(HERE.glob('*.py'))+list(HERE.glob('*.sh'))+[BANK/'submission.npz']
    for relative in ['from_host/verifiers/dfthcsos/verify.py','from_host/verifiers/dfthcsos/scoreinfra/costs_prepare.py',
        'from_host/verifiers/dfthcsos/scoreinfra/dfthc_math.py','from_host/verifiers/dfthcsos/scoreinfra/circuits.py',
        'initial_experiments/5/chem-audit/audit_core.py','initial_experiments/5/chem-audit/refit.py',
        'initial_experiments/5/chem-audit/integer_refine.py','initial_experiments/5/chem-qfit/core.py',
        'initial_experiments/5/chem-qfit/physical.py','initial_experiments/5/chem-qfit/oracles.py',
        'initial_experiments/4/codex/angle_gauge_probe.py','initial_experiments/3/parameters/reference_spectrum.npz']:
        source_paths.append(ROOT/relative)
    source_paths+=list((ROOT/'initial_experiments/5/chem-qrom').glob('*.py'))
    write('source_manifest.json',{str(p.relative_to(ROOT)):sha(p) for p in sorted(source_paths)})
    artifacts={str(p.relative_to(HERE)):sha(p) for base in ['candidates','parameters'] for p in sorted((HERE/base).rglob('*')) if p.is_file()}
    (HERE/'artifact_manifest.json').write_text(json.dumps(artifacts,indent=2)+'\n')
    record=dict(status='PASS',baseline_npz_unchanged=True,selected_digest=digest,selected_score=selected['score'],
        improvement=BASE_COST/selected['score'],target_achieved=selected['score']<=TARGET,
        complete_native_passes=2,research_evaluations=research,independent_bounds=bounds,
        global_rank=global_rank,exact_physical_permutation=True,arbitrary_target_component_checks=len(checks['rows']))
    write('audit.json',record);print(json.dumps(record),flush=True)
if __name__=='__main__':main()
