"""Independently reconstruct primal tensors and explicit feasible dual tensors.

This deliberately constructs dense matrices only for final certificates, not
inside the optimization. Numeric safety slack is checked independently of
cached block contractions. No interval-arithmetic/global-shape claim is made.
"""
import argparse
from core import *

@torch.no_grad()
def main():
    p=argparse.ArgumentParser();p.add_argument('inputs',nargs='+');a=p.parse_args();records=[]
    for path in a.inputs:
        z=np.load(path);m=Covariance(z['U']);S=tensor(z['S']) if 'S' in z else m.from_weights(z['w'])
        approximation,primal=m.direct(S);Y=approximation-m.H
        U,s,_=torch.linalg.svd(m.Q.flatten(0,1).T,full_matrices=False);U=U[:,s>float(s.max())*1e-10];Pi=U@U.T
        restrictions=m.Q@Y@m.Q.transpose(1,2)
        minimum=float(torch.linalg.eigvalsh((restrictions+restrictions.transpose(1,2))/2).min())
        delta=max(0.,-minimum)+1e-10;dual=Y+delta*Pi
        dblocks=m.Q@dual@m.Q.transpose(1,2)
        dual_min=float(torch.linalg.eigvalsh((dblocks+dblocks.transpose(1,2))/2).min())
        value=float(-dual.square().sum()-2*(dual*m.H).sum())
        assert dual_min>0 and value<=primal*primal+1e-8
        identity=float((m.Q@Pi@m.Q.transpose(1,2)-torch.eye(m.B,**D)).abs().max())
        w=m.weights(S);X=(w.transpose(1,2)@m.P).flatten(0,1);weight_error=float((X.T@X-m.H).norm())
        assert abs(weight_error-primal)<1e-7
        rec=dict(input=path,primal_error_2=primal,explicit_dual_error_lower_bound=max(0.,value)**.5,
                 explicit_dual_squared_raw=value,dual_min_eigenvalue=dual_min,dual_shift=delta,
                 union_rank=U.shape[1],projector_restriction_identity_error=identity,
                 factor_weight_reconstruction_error_2=weight_error,excludes_20mHa=value>.0004)
        records.append(rec);print(json.dumps(rec),flush=True)
    report(HERE/'results/explicit_dual_certificates.json',dict(scope='Numerical certificates for exact saved dictionaries only',records=records))
if __name__=='__main__':main()
