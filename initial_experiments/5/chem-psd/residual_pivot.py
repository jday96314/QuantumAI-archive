"""One resource-screened pivot: append exact integer residual projector atoms.

Positive residual modes propose atoms only; the optimized objective always
retains the complete signed native target. Each block is conditioned on its
existing projector span. Four high-leverage additions are shared, with the
remaining additions selected from ensemble-specific residual modes.
"""
from scipy.linalg import eigh
from core import *
from angle_gauge_probe import encode
from solve import solve
from dictionaries import pivot

@torch.no_grad()
def main():
    tag='residual_r2b375';z=dict(np.load(HERE/'parameters/residual_pivot_seed.npz'));m=Covariance(z['U']);tree=z['tree_edges']
    S,h=solve(m,tensor(z['S']),800,None,'fista',200,quiet=True);approx,initial=m.direct(S)
    residual=(m.H-approx).cpu().numpy();ev,V=eigh(residual,subset_by_index=[2926-64,2925]);positive=ev>0
    F=(V[:,positive]*np.sqrt(ev[positive])[None]).T[::-1].copy()
    _,vec=np.linalg.eigh(unpack(F));base=np.load(INC)
    incumbent=reconstruct_u_vectors(base['angle_words'][76:]/128,tree)
    pool=np.r_[incumbent,vec.transpose(0,2,1).reshape(-1,76)]
    poolwords=encode(pool,tree,8,True);pool=reconstruct_u_vectors(poolwords/128,tree)
    P=tensor(features(pool));out=[];indices=[];common=[]
    # Global leverage pivots supply the shared high-leverage residual vectors.
    _,common=pivot(pool,F,4)
    for r in range(2):
        target=tensor(F[r::2]);Q=m.Q[r];cor=P@Q.T
        rhs=P@target.T-cor@(Q@target.T);diag=(P.square().sum(1)-cor.square().sum(1)).clamp(min=0)
        residualP=P-cor@Q;chosen=[];extraQ=[]
        for k in range(24):
            gain=rhs.square().sum(1)/diag.clamp(min=1e-20);gain[diag<1e-10]=-1
            if chosen:gain[chosen]=-1
            j=common[k] if k<4 and float(diag[common[k]])>1e-10 else int(gain.argmax())
            chosen.append(j);v=residualP[j].clone()
            if extraQ:
                qq=torch.stack(extraQ);v-=qq.T@(qq@v);v-=qq.T@(qq@v)
            v/=v.norm();extraQ.append(v);c=P@v;rhs-=c[:,None]*(target@v)[None,:];diag=(diag-c*c).clamp(min=0)
        out.append(np.r_[z['U'][r],pool[chosen]]);indices.append(chosen)
    u=np.stack(out);nm=Covariance(u);w=m.weights(S);warm=torch.cat([w,torch.zeros((2,24,w.shape[2]),**D)],dim=1)
    Sn,hn=solve(nm,nm.from_weights(warm),2400,None,'fista',200,quiet=True)
    words=encode(u.reshape(-1,76),tree,8,True)
    assert np.max(abs(features(reconstruct_u_vectors(words/128,tree).reshape(u.shape))-features(u)))<1e-12
    np.savez(HERE/f'parameters/{tag}_psd.npz',U=u,w=nm.weights(Sn).cpu().numpy(),S=Sn.cpu().numpy(),angle_words=words,tree_edges=tree)
    report(HERE/f'results/{tag}_psd.json',dict(initial_error_2=initial,final=nm.diagnostics(Sn,True),source='residual_pivot_seed.npz',
         strategy='Exact eight-bit residual-eigenatom augmentation, four shared plus twenty ensemble-specific atoms per block',pool_size=len(pool),pool_indices=indices,initial_history=h,history=hn))
if __name__=='__main__':main()
