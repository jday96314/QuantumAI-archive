"""Rank curves on the three previously optimized dictionaries."""
import time
from core import *
from solve import solve

def main():
    inputs=[('prior_r4b180',ROOT/'initial_experiments/5/chem-blocks/parameters/r4b180c88_pca_lbfgs.npz'),
            ('prior_r4b224',ROOT/'initial_experiments/5/chem-blocks/jad_search/candidates/r4_b224_c88_float.npz'),
            ('prior_r2b351',ROOT/'initial_experiments/5/chem-blocks/parameters/r2b351c175_lbfgs.npz')]
    summary=[]
    for name,path in inputs:
        z=dict(np.load(path));m=Covariance(z['U']);start=time.time()
        S0=m.from_weights(z['w']);initial=m.diagnostics(S0,True)
        saved=HERE/f'parameters/{name}_psd.npz'
        if saved.exists() and (HERE/f'results/{name}_psd.json').exists():
            S=tensor(np.load(saved)['S']);full=m.diagnostics(S,True)
        else:
            S,hist=solve(m,S0,2400,None,'fista',200)
            full=m.diagnostics(S,True)
            np.savez(saved,U=z['U'],S=S.cpu().numpy(),w=m.weights(S).cpu().numpy(),tree_edges=z['tree_edges'])
            report(HERE/f'results/{name}_psd.json',dict(input=str(path),initial=initial,diagnostics=full,history=hist))
        curves=[]
        ranks=[88,128,160,175,180] if m.R==4 else [88,128,160,175,192,224,256,351]
        for C in ranks:
            if C>m.B:continue
            old=HERE/f'results/{name}_c{C}.json'
            if old.exists():
                curves.append(dict(C=C,**json.loads(old.read_text())['diagnostics']));continue
            truncated=project(S,C)
            if float((truncated-S).norm())<1e-10:
                Sc,hc=truncated,[]
            else:Sc,hc=solve(m,S,600,C,'coordinate',200)
            dc=m.diagnostics(Sc,True,C);curves.append(dict(C=C,**dc))
            np.savez(HERE/f'parameters/{name}_c{C}.npz',U=z['U'],w=m.weights(Sc,C).cpu().numpy(),S=Sc.cpu().numpy(),tree_edges=z['tree_edges'])
            report(HERE/f'results/{name}_c{C}.json',dict(rank=C,diagnostics=dc,history=hc))
        rec=dict(name=name,full_PSD=full,rank_curve=curves,seconds=time.time()-start);summary.append(rec)
        report(HERE/'results/experiment_a.json',summary)
if __name__=='__main__':main()
