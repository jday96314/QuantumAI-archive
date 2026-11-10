"""Bind research outputs to the unchanged incumbent; never select invalid fits."""
import ast
import hashlib
import json
from pathlib import Path
import numpy as np
from common import HERE,ROUND4,SOURCE,SPECTRUM,report

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    def read(name):return json.loads((HERE/'results'/name).read_text())
    incumbent=read('incumbent_native.json')
    prior=json.loads((ROUND4/'results/norm_cycle2_signed_verify.json').read_text())
    assert incumbent['status']=='PASS' and incumbent['verified_shots']==50880
    for key in ['score','lambda_eff','toffoli_count','num_qubits','submission_sha256','metrics']:
        assert incumbent[key]==prior[key],key
    assert incumbent['score']==708547388.9553173
    frontier=read('frontier.json');oracles=read('oracle_controls.json')
    assert oracles['status']=='PASS'
    assert all(c['status']=='PASS' and c['Q']<=128 and c['T']<=c['bound']['T_per_oracle'] for c in oracles['controls'])
    fits=[]
    for name in ['r4b180c88_adam','r4b180c88_lbfgs','r4b180c88_pca_adam','r4b180c88_pca_lbfgs','r2b351c175_adam','r2b351c175_lbfgs']:
        z=read(name+'_relax.json')
        fits.append(dict(name=name,R=z['R'],B=z['B'],C=z['C'],optimizer=z['optimizer'],steps=z['steps'],best_error_2=z['best_error_2'],evaluations=z['records'][-1]['evaluations']))
    quantized=[]
    for name in ['r4b180c88_pca_adam','r4b180c88_pca_lbfgs','r2b351c175_adam','r2b351c175_lbfgs']:
        z=read(name+'_quantized.json')
        for r in z['records']:
            m=r['native_physical_metrics']
            valid=m['error_1']<=.02 and m['error_2']<=.02 and m['Lambda']>49 and m['E_gap']>3.914898 and m['gap_in_domain']
            assert valid==r['passes_physical_screen']
            assert not valid,'A newly promising candidate requires native oracle integration before finalizing.'
            quantized.append(dict(name=name,brot=r['bits'],metrics=m,margins=r['margins'],passes_physical_screen=valid,
                                  rounded_error_2=r['rounded_error_2'],accepted_integer_moves=r['integer_search']['accepted_moves']))
    # Distance to the rank-k PSD set is 1-Lipschitz in Frobenius norm. The
    # measured cache discrepancy cannot change the frontier rank-290 screen.
    delta=z['cached_target_vs_host_frobenius_discrepancy']
    eig=np.load(SPECTRUM)['eigenvalues']
    def rank_error(k):return float(np.linalg.norm(np.r_[np.minimum(eig[:k],0),eig[k:]]))
    assert rank_error(289)-delta>.02 and rank_error(290)+delta<.02
    reference_audit=dict(cached_target_vs_host_frobenius_discrepancy=delta,
                         cached_rank289_error=rank_error(289),cached_rank290_error=rank_error(290),
                         certified_native_minimum_unrestricted_PSD_rank=290,
                         continuous_error_uncertainty_bound=delta,
                         quantized_metrics_use_exact_native_target=True)
    scripts=sorted(HERE.glob('*.py'))
    for path in scripts:ast.parse(path.read_text(),filename=str(path))
    artifacts=sorted((SOURCE.parent).glob('angle_k_*.kmb'))+[SOURCE]
    hashes={str(p.relative_to(ROUND4)):digest(p) for p in artifacts}
    evidence={str(p.relative_to(HERE)):digest(p) for p in sorted((HERE/'results').glob('*.json')) if p.name!='final_audit.json'}
    code={p.name:digest(p) for p in scripts}
    parameter_hashes={str(p.relative_to(HERE)):digest(p) for p in sorted((HERE/'parameters').glob('*.npz'))}
    report(HERE/'results/final_audit.json',dict(status='COMPLETE_NO_REPLACEMENT',
        selected_candidate=str(SOURCE),incumbent_native=incumbent,
        local_combined_reference=9.3677221531,chemistry_improvement=1.0,target_improvement=1.5,target_achieved=False,
        frontier_configurations=frontier['rows'],frontier_pareto_points=len(frontier['pareto_T_Q']),
        continuous_fits=fits,quantized_screens=quantized,reference_audit=reference_audit,
        oracle_controls=oracles,incumbent_file_sha256=hashes,evidence_sha256=evidence,code_sha256=code,parameter_sha256=parameter_hashes,
        scope='Two shape families, bounded multi-start continuous fits, one eight-bit fallback, no new accepted submission.'))

if __name__=='__main__':main()
