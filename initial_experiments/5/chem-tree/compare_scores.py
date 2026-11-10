"""Evaluate surrogate rankings against actual decoded rank-C tensor fits."""
from core import *
from scipy.stats import spearmanr

def expected_isotropic(u,w,tr):
    P=features(u);F=w.transpose(0,2,1)@P;Z=w@F
    weight=np.sum(Z*Z,axis=-1).reshape(-1);mass=(u*u).reshape(-1,N).copy();total=np.zeros(len(mass))
    for p,c in tr[::-1]:
        mass[:,p]+=mass[:,c];total+=mass[:,p];mass[:,c]=0
    return float(np.sqrt((np.pi/128)**2/3*np.dot(weight,total)))

def main():
    output=[]
    for name in ['incumbent','b639_float','b351_float','b479_exact']:
        _,u,w,_,_=load(SOURCES[name]);rows=json.loads((HERE/f'results/screen_{name}.json').read_text())
        unique={}
        for r in rows:
            tr=load(r['checkpoint'])[3];sig=binary_signature(tr)
            r['binary_signature']=sig;r['expected_isotropic_tensor_rms']=expected_isotropic(u,w,tr)
            unique[sig]=r
        # Tensor Huffman is an exact optimizer of this one isotropic term,
        # under independent uniform small angle errors. The anisotropic,
        # coherent rounding and nonlinear refit terms remain empirical.
        tensor_huffman=next(r for r in rows if r['family']=='huffman_tensor')
        assert tensor_huffman['expected_isotropic_tensor_rms']<=min(r['expected_isotropic_tensor_rms'] for r in rows)+1e-12
        for broad in [True,False]:
            used=[r for r in unique.values() if r['family']!='incumbent' and
                   (broad or r['family'] not in ['star_power','chain_power'])]
            targetval=[r['error_refit'] for r in used]
            correlations={key:float(spearmanr([r[key] for r in used],targetval).statistic)
                 for key in ['vector_rms','projector_rms','linear_tensor_delta','predicted_error','frozen_error','expected_isotropic_tensor_rms']}
            output.append(dict(control=name,include_chain_and_star=broad,n=len(used),spearman=correlations))
        report(HERE/f'results/screen_{name}.json',rows,False)
    report(HERE/'results/surrogate_comparison.json',output)
if __name__=='__main__':main()
