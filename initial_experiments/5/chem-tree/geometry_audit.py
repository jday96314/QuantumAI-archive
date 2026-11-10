"""Span obstructions, conditioning and physical sensitivity of saved dictionaries."""
import argparse
from core import *

def main():
    p=argparse.ArgumentParser();p.add_argument('fits',nargs='+');a=p.parse_args();H=target()
    output=HERE/'results/geometry_audit.json';records=json.loads(output.read_text()) if output.exists() else []
    for path in a.fits:
        z,u,w,tree,words=load(path);P=tensor(features(u));F=tensor(w).transpose(1,2)@P
        Z=tensor(w)@F;sv=torch.linalg.svdvals(P)
        flat=P.reshape(-1,P.shape[-1]);Q,Rm=torch.linalg.qr(flat.T)
        diag=Rm.diagonal().abs();assert float(diag.min())>1e-12
        HC=Q.T@tensor(H)@Q
        residual=tensor(H)-Q@HC@Q.T;lower=float(residual.norm())
        rec=dict(path=str(path),sha256=sha(path),shape=list(w.shape),exact_words=words is not None,
                 error=physical_error(u,w,H),joint_span_lower_bound=lower,
                 span_excludes_20mHa=lower>.020,
                 block_feature_min_singular_values=sv[:,-1].cpu().tolist(),
                 block_feature_condition=(sv[:,0]/sv[:,-1]).cpu().tolist(),
                 weight_frobenius=float(np.linalg.norm(w)),weight_l1=float(np.abs(w).sum()),
                 tensor_row_sensitivity_norm=float(Z.norm()),max_row_sensitivity=float(Z.norm(dim=-1).max()),
                 smallest_union_QR_diagonal=float(diag.min()),tree=summary(tree))
        records=[r for r in records if r['path']!=str(path)]+[rec];print(json.dumps(rec),flush=True)
    report(HERE/'results/geometry_audit.json',records,False)
if __name__=='__main__':main()
