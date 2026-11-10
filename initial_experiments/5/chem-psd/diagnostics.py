"""Independent span/rank and envelope-gradient controls."""
import argparse
from scipy.linalg import eigh
from core import *

@torch.no_grad()
def joint_span(m):
    # This explicitly permits cross-block coefficients and is only a lower
    # bound. SVD handles intentionally shared, duplicate cross-block atoms.
    U,s,_=torch.linalg.svd(m.Q.flatten(0,1).T,full_matrices=False)
    cutoff=float(s.max())*1e-10;keep=s>cutoff;Q=U[:,keep]
    T=Q.T@m.H@Q;ev,V=torch.linalg.eigh((T+T.T)/2)
    approximation=Q@((V*ev.clamp(min=0))@V.T)@Q.T
    return dict(union_span_dimension=int(keep.sum()),union_span_PSD_lower_bound=float((approximation-m.H).norm()),
                union_span_singular_cutoff=cutoff,smallest_retained_singular_value=float(s[keep][-1]))

def control():
    z=np.load(HERE/'parameters/prior_r4b180_psd.npz');u=z['U'];m=Covariance(u);S=tensor(z['S'])
    with torch.no_grad():
        E=m.direct(S)[0]-m.H;KP=torch.linalg.solve_triangular(m.Rm,S@m.Q,upper=True);gp=4*KP@E
        M=torch.zeros((*gp.shape[:2],76,76),**D);M[:,:,ii,jj]=gp/tensor(fac);M[:,:,jj,ii]=gp/tensor(fac)
        gu=2*(M@tensor(u)[:,:,:,None]).squeeze(-1);gu-=tensor(u)*(gu*tensor(u)).sum(2,keepdim=True)
        w=m.weights(S);K=w@w.transpose(1,2)
        rng=np.random.default_rng(227);checks=[]
        for _ in range(3):
            v=rng.normal(size=u.shape);v-=u*(v*u).sum(2,keepdims=True);v/=np.linalg.norm(v)
            def f(eps):
                un=u+eps*v;un/=np.linalg.norm(un,axis=2,keepdims=True);P=tensor(features(un));H=(P.transpose(1,2)@K@P).sum(0)
                return float((H-m.H).square().sum())
            eps=1e-5;fd=(f(eps)-f(-eps))/(2*eps);analytic=float((gu*tensor(v)).sum());checks.append(dict(finite_difference=fd,analytic=analytic,absolute_error=abs(fd-analytic)))
        assert max(c['absolute_error'] for c in checks)<1e-7
        report(HERE/'results/envelope_gradient_control.json',dict(checks=checks))

def main():
    p=argparse.ArgumentParser();p.add_argument('inputs',nargs='*');p.add_argument('--control',action='store_true');a=p.parse_args()
    if a.control:control()
    records=[]
    H=target();ev=eigh(H,eigvals_only=True)[::-1];negative=float(np.square(ev[ev<0]).sum());positive=ev[ev>0]
    for path in a.inputs:
        z=np.load(path);m=Covariance(z['U']);rec=dict(input=path,**joint_span(m))
        if 'S' in z:
            S=tensor(z['S']);e=torch.linalg.eigvalsh(S);rec['PSD_ranks_at_1e8_relative']=(e>e[:,-1:]*1e-8).sum(1).tolist()
        records.append(rec)
    ranks=[88,128,160,175,192,224,256,289,290,320,350,352,512,640,720]
    report(HERE/'results/span_diagnostics.json',dict(unrestricted_PSD_rank_bounds=[dict(rank=k,error_2=float(np.sqrt(negative+np.square(positive[k:]).sum()))) for k in ranks],dictionaries=records))
if __name__=='__main__':main()
