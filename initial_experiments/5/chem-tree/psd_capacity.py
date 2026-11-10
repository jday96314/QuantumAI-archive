"""Convex full-PSD block capacity and independent feasible residual dual.

Whitened covariance/FISTA follows completed chem-psd; sequential dual cone
repair follows completed chem-discrete/psd_certificate.py. Numerical certificates
apply only to the exact decoded dictionaries, not to tree families or shapes.
"""
import argparse
from core import *

def psd(S):
    ev,V=torch.linalg.eigh((S+S.transpose(-1,-2))/2)
    return (V*ev.clamp(min=0)[:,None,:])@V.transpose(-1,-2)

def main():
    p=argparse.ArgumentParser();p.add_argument('--init',required=True);p.add_argument('--tag',required=True)
    p.add_argument('--steps',type=int,default=2500);a=p.parse_args();z,u,w,tr,words=load(a.init)
    m=Refit(u);R,B,C=w.shape;S=m.from_weights(w);Y=S.clone();tau=1.;start=time.time();history=[]
    def apply(S):return ((m.cross@S[None])@m.cross.transpose(-1,-2)).sum(1)
    best=S.clone();bestf=float(m.objective(S))
    for k in range(a.steps):
        Sn=psd(Y-(apply(Y)-m.A)/R);f=float(m.objective(Sn))
        if f<bestf:best,bestf=Sn.clone(),f
        taun=(1+np.sqrt(1+4*tau*tau))/2
        if float(((Sn-S)*(Y-Sn)).sum())>0:tau=1.;Y=Sn.clone()
        else:Y=Sn+(tau-1)/taun*(Sn-S);tau=taun
        S=Sn
        if k%100==0 or k==a.steps-1:
            g=apply(S)-m.A;gamma=max(0.,-float(torch.linalg.eigvalsh((g+g.transpose(-1,-2))/2).min()))+1e-11
            comp=float((S*g).sum());lower=f-2*comp-2*gamma*float(S.diagonal(dim1=1,dim2=2).sum())-gamma*gamma*m.H.shape[0]
            row=dict(step=k,primal=max(0,f)**.5,scalar_dual=max(0,lower)**.5,dual_shift=gamma,seconds=time.time()-start)
            history.append(row);report(HERE/f'results/{a.tag}_history.json',history,False);print(json.dumps(row),flush=True)
    S=best;ev,V=torch.linalg.eigh((S+S.transpose(1,2))/2)
    ww=torch.linalg.solve_triangular(m.Rm,V*ev.clamp(min=0).sqrt()[:,None,:],upper=True)
    ww=ww.cpu().numpy();M=(m.Q.transpose(1,2)@S@m.Q).sum(0);D0=m.H-M;Ddual=D0.clone()
    for Q in m.Q:
        compressed=Q@Ddual@Q.T;ev,V=torch.linalg.eigh((compressed+compressed.T)/2)
        F=Q.T@(V*ev.clamp(min=0).sqrt()[None]);Ddual-=F@F.T
    maxev=[float(torch.linalg.eigvalsh(Q@Ddual@Q.T)[-1]) for Q in m.Q]
    shift=max(0,max(maxev))+1e-11;Ddual-=shift*torch.eye(len(m.H),**D)
    finalev=[float(torch.linalg.eigvalsh(Q@Ddual@Q.T)[-1]) for Q in m.Q]
    lower2=float(2*(m.H*Ddual).sum()-Ddual.square().sum());direct=physical_error(u,ww)
    assert abs(direct-float(D0.norm()))<1e-9 and max(finalev)<0
    save(HERE/f'candidates/{a.tag}.npz',u,ww,tr,words,err2=direct)
    np.savez_compressed(HERE/f'candidates/{a.tag}_dual.npz',D=Ddual.cpu().numpy(),angle_words=words,tree_edges=tr)
    report(HERE/f'results/{a.tag}_final.json',dict(source=a.init,source_sha256=sha(a.init),shape=list(w.shape),
           full_PSD_error=direct,dual_lower_bound=max(0,lower2)**.5,dual_squared_bound=lower2,
           block_max_dual_eigenvalues=finalev,dual_safety_shift=shift,steps=a.steps,seconds=time.time()-start,
           fixed_dictionary_excluded=bool(lower2>.020**2),precision='float64 numerical safety slack, not interval arithmetic'))
if __name__=='__main__':main()
