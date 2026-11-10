"""Cheap exact native-semantics normalization on a decoded pair dictionary."""
from core import *

class Normalization:
    def __init__(self,m,off):
        self.m=m;z=np.load(HERE/'candidates/native_constants.npz')
        self.h0=m.tensor(z['h0']);self.constant=float(z['constant']);self.off=m.tensor(off)
    def vals(self,u,w,off):
        m=self.m;P=m.feats(u);packed=torch.einsum('rbc,rc,rbp->p',w,off,P)/m.fac
        hc=torch.zeros((76,76),**m.d);hc[m.ii,m.jj]=packed;hc=hc+hc.T-torch.diag(hc.diag())
        ev=torch.linalg.eigvalsh(self.h0-hc);l1=ev.abs().sum()
        gap=self.constant+l1-ev.sum()+.5*off.square().sum()
        lam=l1+.25*(w.abs().sum(1)+(off+w.sum(1)).abs()).square().sum()
        return gap,lam,gap*(2*lam-gap)
    def evaluate(self,u,w,steps=12):
        # Offsets are continuously refitted for every candidate as requested.
        # The smooth LCU objective is used only to screen; final exact values
        # use abs and the native gap, with a strict-floor penalty.
        with torch.enable_grad():
            off=self.off.clone().requires_grad_(True)
            opt=torch.optim.LBFGS([off],lr=.5,max_iter=steps,history_size=8,line_search_fn='strong_wolfe',tolerance_change=1e-12)
            def closure():
                opt.zero_grad();g,l,s=self.vals(u.detach(),w.detach(),off)
                loss=s/3000+100*torch.relu(3.9153-g).square()
                loss.backward();return loss
            opt.step(closure)
        with torch.no_grad():
            g,l,s=self.vals(u,w,off)
            score=float(s)+300000*max(0,3.91505-float(g))**2
            return score,off.detach(),dict(E_gap=float(g),Lambda=float(l),lambda_eff=float(s.clamp(min=0).sqrt()))
