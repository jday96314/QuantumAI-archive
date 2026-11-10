"""Exact native pair metric and reusable compact rank-constrained refits.

Metric/QR algebra follows completed chem-psd and chem-discrete (read-only).
No inherited spectrum is used as the physical target.
"""
from pathlib import Path
import json,os,time,hashlib
import numpy as np
import torch
from dfthcsos.verify import load_inputs
from trees import *
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
INC=ROOT/'initial_experiments/4/codex/candidates/norm_cycle2_signed/submission.npz'
torch.set_num_threads(2)
torch.backends.cuda.matmul.allow_tf32=False
DEVICE=os.environ.get('TREE_DEVICE','cpu')
D=dict(device=DEVICE,dtype=torch.float64)
ii,jj=np.tril_indices(N);fac=np.where(ii==jj,1.,np.sqrt(2.))
def tensor(x):return torch.as_tensor(x,**D)
def features(u):return u[...,ii]*u[...,jj]*fac
def unpack(p):
    out=np.zeros((*p.shape[:-1],N,N));out[...,ii,jj]=p/fac;out[...,jj,ii]=p/fac
    return out
def report(path,data,verbose=True):
    def convert(x):
        if isinstance(x,np.generic):return x.item()
        if isinstance(x,Path):return str(x)
        raise TypeError(type(x).__name__)
    Path(path).write_text(json.dumps(data,indent=2,allow_nan=False,default=convert)+'\n')
    if verbose:print(json.dumps(data,allow_nan=False,default=convert),flush=True)
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def native():return load_inputs(INC)
def target():
    path=HERE/'candidates/native_target.npz'
    if path.exists():return np.load(path)['H']
    inp=native();s=inp.S[ii,jj]*fac;ident=(ii==jj).astype(float)
    H=inp.ref['tei'][ii[:,None],jj[:,None],ii[None,:],jj[None,:]]*fac[:,None]*fac[None,:]
    H+=.5*(s[:,None]*ident[None,:]+ident[:,None]*s[None,:])
    np.savez(path,H=H)
    return H
SOURCES={
 'incumbent':INC,
 'b639_float':ROOT/'initial_experiments/5/chem-qfit/candidates/matched_B639_warmup_float.npz',
 'b351_float':ROOT/'initial_experiments/5/chem-psd/parameters/claude_r2b351_threshold_c189.npz',
 'b479_exact':ROOT/'initial_experiments/5/chem-qfit/candidates/legacy_exact_split_B479.npz',
 'b639_exact':ROOT/'initial_experiments/5/chem-discrete/candidates/native_norm_B639.npz',
}
def load(path):
    z=dict(np.load(path));u=z['U'];w=z['w'] if 'w' in z else z['W'];tree=z.get('tree_edges',np.load(INC)['tree_edges'])
    words=z.get('angle_words');words=None if words is None else words.reshape(-1,75)[-u.shape[0]*u.shape[1]:]
    if words is not None:u=decode(words*np.pi/(1<<(int(z.get('brot',8))-1)),tree).reshape(u.shape)
    u=u/np.linalg.norm(u,axis=-1,keepdims=True)
    return z,u,w,tree,words

def physical_error(u,w,H=None):
    H=target() if H is None else H
    P=features(u);F=np.matmul(w.transpose(0,2,1),P).reshape(-1,len(ii))
    return float(np.linalg.norm(F.T@F-H))

class Refit:
    """Rank-C PSD cyclic exact minimization in block-orthonormal coordinates.

    QR and target compressions are cached. Inner sweeps work only on BxB
    matrices; objective is independently checked against dense physical error.
    Columns are aligned to warm weights to keep normalization comparisons fair.
    """
    def __init__(self,u,H=None):
        self.u=np.asarray(u);self.R,self.B,_=u.shape;self.P=tensor(features(u))
        self.H=tensor(target() if H is None else H)
        Qt,self.Rm=torch.linalg.qr(self.P.transpose(1,2),mode='reduced')
        self.Q=Qt.transpose(1,2).contiguous()
        self.A=self.Q@self.H@Qt
        flat=self.Q.flatten(0,1)
        self.cross=(flat@flat.T).reshape(self.R,self.B,self.R,self.B).permute(0,2,1,3).contiguous()
        self.h2=self.H.square().sum()
    def from_weights(self,w):
        F=self.Rm@tensor(w);return F@F.transpose(1,2)
    def objective(self,S):
        applied=((self.cross@S[None])@self.cross.transpose(-1,-2)).sum(1)
        return self.h2-2*(S*self.A).sum()+(S*applied).sum()
    @torch.no_grad()
    def fit(self,w,cycles=8,rank=None,align=True):
        C=w.shape[-1] if rank is None else rank;S=self.from_weights(w)
        old=float(self.objective(S));history=[max(0,old)**.5]
        for k in range(cycles):
            for r in range(self.R):
                H=self.A[r].clone()
                for s in range(self.R):
                    if r!=s:H-=self.cross[r,s]@S[s]@self.cross[s,r]
                ev,V=torch.linalg.eigh((H+H.T)/2)
                F=V[:,-C:]*ev[-C:].clamp(min=0).sqrt()[None]
                S[r]=F@F.T
            f=float(self.objective(S));history.append(max(0,f)**.5)
            if abs(f-old)<1e-15:break
            old=f
        ev,V=torch.linalg.eigh((S+S.transpose(1,2))/2)
        F=V[:,:,-C:]*ev[:,-C:].clamp(min=0).sqrt()[:,None,:]
        out=torch.linalg.solve_triangular(self.Rm,F,upper=True)
        if align and C==w.shape[-1]:
            left,_,right=torch.linalg.svd(out.transpose(1,2)@tensor(w),full_matrices=False)
            out=out@(left@right)
        ww=out.cpu().numpy();err=physical_error(self.u,ww,self.H.cpu().numpy())
        assert abs(err**2-float(self.objective(S)))<2e-8
        return ww,err,history

def save(path,u,w,tree,words=None,**kw):
    data=dict(U=u,w=w,tree_edges=tree,**kw)
    if words is not None:data.update(angle_words=words.reshape(-1,75).astype(np.int64)%256,brot=8)
    np.savez(path,**data)

def surrogate(u,v,w,H=None):
    """Sensitivity to decoded projector increments, not mean vector L2.

    dF=w.T(P(v)-P(u)); first-order dM=F.T dF+dF.T F. Exact contractions
    include leverage, factor-column coherence, pair metric and interblock
    cancellation. Also report exact frozen-weight tensor error.
    """
    H=target() if H is None else H
    P=features(u);DP=features(v)-P
    F=(w.transpose(0,2,1)@P).reshape(-1,len(ii))
    dF=(w.transpose(0,2,1)@DP).reshape(-1,len(ii))
    FF=F@F.T;DD=dF@dF.T;FD=F@dF.T
    linear2=2*np.sum(FF*DD)+2*np.trace(FD@FD)
    E=F.T@F-H
    cross=4*np.sum((F@E)*dF)
    frozen=F+dF
    return dict(vector_rms=float(np.sqrt(np.mean(np.minimum(np.sum((u-v)**2,axis=-1),np.sum((u+v)**2,axis=-1))))),
                projector_rms=float(np.sqrt(np.mean(DP*DP)*len(ii))),
                linear_tensor_delta=float(np.sqrt(max(0,linear2))),
                predicted_error=float(np.sqrt(max(0,np.sum(E*E)+cross+linear2))),
                frozen_error=float(np.linalg.norm(frozen.T@frozen-H)))
