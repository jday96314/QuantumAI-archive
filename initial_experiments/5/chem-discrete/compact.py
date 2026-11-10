"""Exact full symmetric-pair metric in a fixed dictionary's joint span.

The QR basis includes the union of all block spans; PSD factors remain
separate, with exactly C columns per block. The outside residual is measured
explicitly rather than by subtracting two large norms.
"""
from core import *

class Compact:
    def __init__(self,m,u):
        self.m=m;self.P=m.feats(u);self.R,self.B,_=self.P.shape
        self.q=[];self.tr=[]
        for p in self.P:
            q,tr=torch.linalg.qr(p.T);self.q.append(q);self.tr.append(tr)
        self.q=torch.stack(self.q);self.tr=torch.stack(self.tr)
        self.Q=torch.linalg.qr(self.q.transpose(1,2).reshape(-1,2926).T).Q
        self.H=self.Q.T@m.G@self.Q
        self.A=torch.stack([self.Q.T@q for q in self.q])
        self.outside=(m.G-self.Q@self.H@self.Q.T).square().sum()
    def factors(self,w):return torch.bmm(self.tr,w)
    def weights(self,f):return torch.linalg.solve(self.tr,f)
    def loss(self,f):
        x=torch.bmm(self.A,f).transpose(1,2).flatten(0,1)
        return (x.T@x-self.H).square().sum()+self.outside
    @torch.no_grad()
    def project(self,f,cycles=1):
        C=f.shape[-1];f=f.clone();x=torch.bmm(self.A,f)
        for _ in range(cycles):
            for r in range(self.R):
                target=self.H-sum(x[s]@x[s].T for s in range(self.R) if s!=r)
                h=self.A[r].T@target@self.A[r];ev,v=torch.linalg.eigh((h+h.T)/2)
                fr=v[:,-C:]*ev[-C:].clamp(min=0).sqrt()
                left,_,right=torch.linalg.svd(fr.T@f[r],full_matrices=False)
                f[r]=fr@(left@right);x[r]=self.A[r]@f[r]
        return f,float(self.loss(f).sqrt())

def derivatives(m,t):
    """Exact decoded vector and word derivatives in the native tree order."""
    shape=t.shape;tf=t.reshape(-1,75);v=torch.zeros((len(tf),76),**m.d);v[:,0]=1
    J=torch.zeros((len(tf),75,76),**m.d);alpha=np.pi/128
    for j,(p,c) in enumerate(m.tree):
        co=(tf[:,j]*alpha).cos();sn=(tf[:,j]*alpha).sin()
        jp,jc=J[:,:,p].clone(),J[:,:,c].clone()
        J[:,:,p]=co[:,None]*jp-sn[:,None]*jc;J[:,:,c]=sn[:,None]*jp+co[:,None]*jc
        vp,vc=v[:,p].clone(),v[:,c].clone();v[:,p]=co*vp-sn*vc;v[:,c]=sn*vp+co*vc
        J[:,j,p]=-alpha*v[:,c];J[:,j,c]=alpha*v[:,p]
    return v.reshape(*shape[:-1],76),J.reshape(*shape,76)
