"""Native signed target and metric-correct covariance algebra.

Host DFTHC imports are Apache-2.0, from from_host/verifiers (read-only).
The target is constructed from native float32 integrals and fixed shifts.
"""
from pathlib import Path
import json
import os
import numpy as np
import torch
from dfthcsos.verify import load_inputs
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
INC=ROOT/'initial_experiments/4/codex/candidates/norm_cycle2_signed/submission.npz'
torch.set_num_threads(2)
torch.backends.cuda.matmul.allow_tf32=False
D=dict(device=os.environ.get('PSD_DEVICE','cuda'),dtype=torch.float64)
ii,jj=np.tril_indices(76)
fac=np.where(ii==jj,1.,np.sqrt(2.))

def tensor(x):return torch.as_tensor(x,**D)
def features(u):return u[...,ii]*u[...,jj]*fac
def unpack(f):
    out=np.zeros((*f.shape[:-1],76,76));out[...,ii,jj]=f/fac;out[...,jj,ii]=f/fac
    return out
def report(path,data):
    Path(path).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    print(json.dumps(data,allow_nan=False),flush=True)

def target():
    path=HERE/'parameters/native_target.npz'
    if path.exists():return np.load(path)['H']
    native=load_inputs(INC);s=native.S[ii,jj]*fac;ident=(ii==jj).astype(float)
    H=native.ref['tei'][ii[:,None],jj[:,None],ii[None,:],jj[None,:]]*fac[:,None]*fac[None,:]
    H+=.5*(s[:,None]*ident[None,:]+ident[:,None]*s[None,:])
    np.savez(path,H=H)
    return H

def native_spectrum():
    path=HERE/'parameters/native_eigenvalues.npy'
    if path.exists():return np.load(path)
    ev=np.linalg.eigvalsh(target())[::-1].copy();np.save(path,ev);return ev

class Covariance:
    """H_r=Q_r.T S_r Q_r, Q_r Q_r.T=I; S_r >=0.

    Original atom covariance K=Rm^-1 S Rm^-T for P.T=Q.T Rm.
    No jitter changes the feature metric. Inputs require independent feature
    rows inside each block; saved solutions undergo explicit factor-weight
    reconstruction checks. Sharing atoms across blocks is allowed.
    """
    def __init__(self,u,H=None):
        self.u=np.asarray(u);self.R,self.B,_=self.u.shape
        self.P=tensor(features(self.u));self.H=tensor(target() if H is None else H)
        Qt,self.Rm=torch.linalg.qr(self.P.transpose(1,2),mode='reduced')
        self.Q=Qt.transpose(1,2).contiguous()
        self.orth_error=float((self.Q@Qt-torch.eye(self.B,**D)).abs().max())
        self.A=self.Q@self.H@Qt
        flat=self.Q.flatten(0,1)
        self.cross=(flat@flat.T).reshape(self.R,self.B,self.R,self.B).permute(0,2,1,3).contiguous()
        self.h2=self.H.square().sum();self.dim=min(self.R*self.B,self.H.shape[0])
    def apply(self,S):
        return ((self.cross@S[None])@self.cross.transpose(-1,-2)).sum(1)
    def objective(self,S):
        return self.h2-2*(S*self.A).sum()+(S*self.apply(S)).sum()
    def direct(self,S):
        approx=(self.Q.transpose(1,2)@S@self.Q).sum(0)
        return approx, float((approx-self.H).norm())
    def from_weights(self,w):
        F=self.Rm@tensor(w);return F@F.transpose(1,2)
    def weights(self,S,C=None):
        C=self.B if C is None else C
        ev,V=torch.linalg.eigh((S+S.transpose(1,2))/2)
        F=V[:,:,-C:]*ev[:,-C:].clamp(min=0).sqrt()[:,None,:]
        return torch.linalg.solve_triangular(self.Rm,F,upper=True)
    def diagnostics(self,S,direct=False,C=None):
        g=self.apply(S)-self.A
        mineig=float(torch.linalg.eigvalsh((g+g.transpose(1,2))/2).min())
        # Y=A(S)-H is a dual variable. Y+delta*Pi_union is dual feasible
        # since Q_r Pi_union Q_r.T=I. Use dim<=R*B for a conservative norm
        # penalty, so union rank need not be computed. Safety covers roundoff.
        delta=max(0.,-mineig)+1e-11
        f=float(self.objective(S));comp=float((g*S).sum())
        lower=f-2*comp-2*delta*float(S.diagonal(dim1=1,dim2=2).sum())-delta*delta*self.dim
        pg=S-project(S-g)
        out=dict(error_2=max(0.,f)**.5,dual_error_lower_bound=max(0.,lower)**.5,
                 dual_squared_bound_raw=lower,dual_shift=delta,gradient_min_eigenvalue=mineig,
                 complementarity=comp,projected_gradient_norm=float(pg.norm()),
                 primal_dual_gap_squared=f-lower,orthogonality_error=self.orth_error)
        if C is not None:out['rank_projected_gradient_norm']=float((S-project(S-g,C)).norm())
        if direct:
            _,err=self.direct(S);out['direct_error_2']=err;out['objective_squared_discrepancy']=abs(err*err-f)
            assert out['objective_squared_discrepancy']<2e-8
        return out

def project(S,C=None):
    ev,V=torch.linalg.eigh((S+S.transpose(-1,-2))/2)
    if C is not None:ev,V=ev[:,-C:],V[:,:,-C:]
    return (V*ev.clamp(min=0)[:,None,:])@V.transpose(-1,-2)
