"""Exact per-vector Riemannian Hessian of the full signed tensor residual.

Off-vector Hessian couplings are omitted for a block-Jacobi proposal only;
every proposed move is still accepted against the complete physical tensor.
Negative eigenvalues are damped, not interpreted as an optimality certificate.
"""
from core import *

@torch.no_grad()
def hessians(m,S):
    u=tensor(m.u).flatten(0,1);approx,_=m.direct(S);E=approx-m.H
    Y=torch.linalg.solve_triangular(m.Rm,S@m.Q,upper=True).flatten(0,1)
    w=m.weights(S);alpha=w.square().sum(2).flatten();yy=Y.square().sum(1)
    def matrix(f):
        out=torch.zeros((len(f),76,76),**D);out[:,ii,jj]=f/tensor(fac);out[:,jj,ii]=f/tensor(fac);return out
    A=matrix(Y@E);a=(matrix(Y)@u[:,:,None]).squeeze(-1)
    # Symmetric-pair reshuffling avoids expanding an N^4 tensor. The
    # symmetrized (ps,qt) kernel contracts all u[p]u[s] in one pair-space GEMM.
    p,q=np.indices((76,76));idx=np.maximum(p,q)*(np.maximum(p,q)+1)//2+np.minimum(p,q)
    index=torch.tensor(idx,device=D['device']);fmat=tensor(np.where(p==q,1.,np.sqrt(2.)))
    ti=torch.tensor(ii,device=D['device']);tj=torch.tensor(jj,device=D['device']);ff=tensor(fac)
    E=(E+E.T)/2;reshuffled=torch.empty_like(E)
    for lo in range(0,len(ii),128):
        hi=min(lo+128,len(ii));pa,sa=ti[lo:hi,None],tj[lo:hi,None];qb,tb=ti[None],tj[None]
        pq,st,sq,pt=index[pa,qb],index[sa,tb],index[sa,qb],index[pa,tb]
        term1=E[pq,st]/(fmat[pa,qb]*fmat[sa,tb]);term2=E[sq,pt]/(fmat[sa,qb]*fmat[pa,tb])
        reshuffled[lo:hi]=.5*(term1+term2)*ff[lo:hi,None]*ff[None]
    contraction=matrix(tensor(features(m.u)).flatten(0,1)@reshuffled)
    eye=torch.eye(76,**D)[None];radial=8*(u*(A@u[:,:,None]).squeeze(-1)).sum(1)
    H=(8*yy-radial)[:,None,None]*eye+16*a[:,:,None]*a[:,None,:]+16*alpha[:,None,None]*contraction+8*A
    P=eye-u[:,:,None]*u[:,None,:];H=P@H@P+8*yy[:,None,None]*u[:,:,None]*u[:,None,:]
    rawg=8*(A@u[:,:,None]).squeeze(-1);g=rawg-u*(u*rawg).sum(1,keepdim=True)
    return (H+H.transpose(1,2))/2,g,yy

@torch.no_grad()
def direction(m,S):
    H,g,yy=hessians(m,S);ev,V=torch.linalg.eigh(H)
    floor=(8*yy*.05).clamp(min=1e-7)
    d=-(V@((V.transpose(1,2)@g[:,:,None]).squeeze(-1)/ev.clamp(min=floor[:,None]))[:,:,None]).squeeze(-1)
    u=tensor(m.u).flatten(0,1);d-=u*(d*u).sum(1,keepdim=True)
    return d.reshape(m.u.shape),dict(hessian_min_eigenvalue=float(ev.min()),damped_eigenvalues=int((ev<floor[:,None]).sum()))

def main():
    z=np.load(HERE/'parameters/claude_r2b351_q8_psd.npz');m=Covariance(z['U']);S=tensor(z['S']);H,g,_=hessians(m,S)
    w=m.weights(S);K=w@w.transpose(1,2);u=tensor(z['U']);r,b=0,17;row=r*m.B+b
    rng=np.random.default_rng(1801);v=tensor(rng.normal(size=76));v-=u[r,b]*(u[r,b]@v);v/=v.norm();eps=1e-5
    def grad_at(sign):
        un=u.clone();un[r,b]+=sign*eps*v;un[r,b]/=un[r,b].norm()
        P=un[:,:,ii]*un[:,:,jj]*tensor(fac);Y=K@P;E=(P.transpose(1,2)@Y).sum(0)-m.H;f=Y[r,b]@E
        A=torch.zeros((76,76),**D);A[ii,jj]=f/tensor(fac);A[jj,ii]=f/tensor(fac)
        raw=8*A@un[r,b];return raw-un[r,b]*(raw@un[r,b])
    fd=(grad_at(1)-grad_at(-1))/(2*eps);fd-=u[r,b]*(fd@u[r,b]);analytic=H[row]@v
    err=float((fd-analytic).norm());assert err<1e-5,(err,float(fd.norm()))
    report(HERE/'results/hessian_control.json',dict(absolute_error=err,finite_difference_norm=float(fd.norm()),analytic_norm=float(analytic.norm())))
if __name__=='__main__':main()
