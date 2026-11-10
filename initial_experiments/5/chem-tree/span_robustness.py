"""Quantify fragile continuous dictionary directions, without global claims."""
from core import *

def main():
    source=SOURCES['b351_float'];_,u,w,tr,_=load(source);P=tensor(features(u));Q0,_=torch.linalg.qr(P.transpose(1,2))
    sv0=torch.linalg.svdvals(P);records=[]
    for name in ['b351_float__incumbent','b351_float__balanced_power','b351_float__chain_power',
                 'b351_float__star_power','b351_float__huffman_tensor','topology_b351_polished','grid9_b351']:
        path=HERE/f'candidates/{name}.npz';_,v,ww,tt,words=load(path);PQ=tensor(features(v));delta=PQ-P
        Q1,_=torch.linalg.qr(PQ.transpose(1,2));cosines=torch.linalg.svdvals(Q0.transpose(1,2)@Q1)
        dsv=torch.linalg.svdvals(delta);ratio=dsv[:,0]/sv0[:,-1]
        row=dict(candidate=name,source=str(source),perturbation_operator_norm=dsv[:,0].cpu().tolist(),
                 minimum_source_singular_value=sv0[:,-1].cpu().tolist(),
                 perturbation_over_minimum_singular=ratio.cpu().tolist(),
                 minimum_span_cosine=cosines[:,-1].cpu().tolist(),
                 maximum_principal_angle_degrees=(torch.arccos(cosines[:,-1].clamp(-1,1))*180/np.pi).cpu().tolist(),
                 median_span_cosine=cosines.median(-1).values.cpu().tolist(),
                 physical_error=physical_error(v,ww),feature_frobenius_perturbation=float(delta.norm()))
        records.append(row);print(json.dumps(row),flush=True)
    report(HERE/'results/span_robustness.json',records,False)
if __name__=='__main__':main()
