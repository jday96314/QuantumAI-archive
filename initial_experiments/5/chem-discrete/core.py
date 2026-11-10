"""Shared exact word, pair metric and unregularized block projection operations.

Givens conventions follow the native verifier. Signed gauge encoder and circuit
builder are imported read-only from round-4 Codex; pair metric follows the
Apache-2.0 host DFTHC code. No score rules are modified.
"""
from pathlib import Path
import os
import numpy as np
import torch
from angle_gauge_probe import encode
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors

ROOT=Path('/media/shared/ExpansionDrive2/QuantumAI')
HERE=Path(__file__).resolve().parent
INC=ROOT/'initial_experiments/4/codex/candidates/norm_cycle2_signed'
SP=ROOT/'initial_experiments/3/parameters/reference_spectrum.npz'
torch.set_num_threads(3)
torch.backends.cuda.matmul.allow_tf32=False
torch.backends.cudnn.allow_tf32=False

def load(path):
    z=dict(np.load(path)); u=z['U']
    if 'w' in z:w=z['w']
    elif 'W' in z and z['W'].ndim==3:w=z['W']
    else:
        ws=[]
        Ks=z['K'][None] if z['K'].ndim==2 else z['K']
        for K in Ks:
            ev,V=np.linalg.eigh(K);ix=np.flatnonzero(ev>max(ev.max()*1e-9,1e-12))
            ws.append(V[:,ix]*np.sqrt(ev[ix]))
        w=np.stack(ws)
    tree=z.get('tree_edges',np.load(INC/'submission.npz')['tree_edges'])
    words=z.get('angle_words')
    if words is not None: words=words[-u.shape[0]*u.shape[1]:]
    return z,u,w,tree,words

class Metric:
    def __init__(self,tree,dtype=torch.float64,align_columns=True):
        self.align_columns=align_columns
        self.d=dict(device='cuda',dtype=dtype);sp=np.load(SP)
        self.ii=torch.tensor(sp['ii'],device='cuda');self.jj=torch.tensor(sp['jj'],device='cuda')
        self.fac=torch.tensor(sp['fac'],**self.d)
        self.ev=sp['eigenvalues'];native=HERE/'candidates/native_pair_target.npz'
        if not native.exists():
            from native_target import build
            build()
        with np.load(native) as reference:self.G=torch.tensor(reference['G'],**self.d)
        self.tree=tree;self.stages=[]
        for t,(p,c) in enumerate(tree):
            if not self.stages or any(p in (pp,cc) or c in (pp,cc) for _,pp,cc in self.stages[-1]):self.stages.append([])
            self.stages[-1].append((t,int(p),int(c)))
        self.stages=[tuple(torch.tensor(x,device='cuda') for x in zip(*s)) for s in self.stages]
    def tensor(self,x):return torch.tensor(x,**self.d)
    def decode(self,t):
        u=torch.zeros((*t.shape[:-1],76),**self.d);u[...,0]=1
        for ts,ps,cs in self.stages:
            an=t[...,ts]*(np.pi/128);co,sn=an.cos(),an.sin();up,uc=u[...,ps],u[...,cs]
            new=u.clone();new[...,ps]=co*up-sn*uc;new[...,cs]=sn*up+co*uc;u=new
        return u
    def feats(self,u):
        u=u/u.norm(dim=-1,keepdim=True)
        return u[...,self.ii]*u[...,self.jj]*self.fac
    def loss(self,u,w):
        P=self.feats(u);X=torch.bmm(w.transpose(1,2),P).flatten(0,1)
        return (X.T@X-self.G).square().sum()
    @torch.no_grad()
    def project(self,u,w,cycles=2,ridge=0):
        """Exact PSD rank-C minimizer per block, QR, no ridge by default.

        Every accepted block update decreases the actual residual. Warm weights
        remain feasible, including when QR/eigensolver roundoff would worsen it.
        """
        P=self.feats(u);R,B,C=w.shape;w=w.clone()
        cont=[(w[r].T@P[r]).T@(w[r].T@P[r]) for r in range(R)];total=sum(cont)
        for _ in range(cycles):
            for r in range(R):
                target=self.G-(total-cont[r])
                if ridge:
                    S=P[r]@P[r].T;L=torch.linalg.cholesky(S+ridge*torch.eye(B,**self.d))
                    H=P[r]@target@P[r].T
                    H=torch.linalg.solve_triangular(L,H,upper=False)
                    H=torch.linalg.solve_triangular(L,H.T,upper=False).T
                    ev,V=torch.linalg.eigh((H+H.T)/2)
                    wr=torch.linalg.solve_triangular(L.T,V[:,-C:]*ev[-C:].clamp(min=0).sqrt(),upper=True)
                else:
                    Q,Rm=torch.linalg.qr(P[r].T)
                    H=Q.T@target@Q;ev,V=torch.linalg.eigh((H+H.T)/2)
                    wr=torch.linalg.solve(Rm,V[:,-C:]*ev[-C:].clamp(min=0).sqrt())
                X=wr.T@P[r];new=X.T@X;nt=total-cont[r]+new
                if (nt-self.G).square().sum() < (total-self.G).square().sum():
                    # Keep the closest column gauge to the warm weights. This
                    # preserves the tensor and prevents each eigensolve from
                    # scrambling the optimizer state and LCU normalization.
                    if self.align_columns:
                        left,_,right=torch.linalg.svd(wr.T@w[r],full_matrices=False)
                        wr=wr@(left@right)
                    total=nt;cont[r]=new;w[r]=wr
        return w,float((total-self.G).norm())

def save(path,u,w,tree,words=None,**kw):
    def arr(t):return t.detach().cpu().numpy() if isinstance(t,torch.Tensor) else t
    data=dict(U=arr(u),w=arr(w),tree_edges=tree,**kw)
    if words is not None:data.update(angle_words=np.rint(arr(words)).astype(np.int64).reshape(-1,75)%256,brot=8)
    path=Path(path);tmp=path.with_name(path.name+f'.{os.getpid()}.tmp')
    with tmp.open('wb') as f:np.savez(f,**data)
    tmp.replace(path)
