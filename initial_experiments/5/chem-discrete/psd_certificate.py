"""Dual lower bound for fixed-dictionary PSD block reconstruction, any ranks.

For any D with Q_r^T D Q_r negative semidefinite, the squared physical
residual is at least 2<G,D>-||D||². Sequentially subtracting each positive
block compression preserves all previously satisfied dual inequalities.
"""
import argparse,json
from core import *
p=argparse.ArgumentParser();p.add_argument('fit');a=p.parse_args()
z,u,w,tree,words=load(a.fit);R,B,C=w.shape
u=reconstruct_u_vectors(words/128,tree).reshape(R,B,76);sp=np.load(SP);ii,jj,fac=sp['ii'],sp['jj'],sp['fac']
G=np.load(HERE/'candidates/native_pair_target.npz')['G'];P=u[:,:,ii]*u[:,:,jj]*fac
X=np.matmul(w.transpose(0,2,1),P).reshape(-1,len(ii));M=X.T@X;D0=G-M;qs=[np.linalg.qr(pp.T)[0] for pp in P]
spectral=(sp['eigenvectors']*sp['eigenvalues'])@sp['eigenvectors'].T
proposals=[];best=-float('inf');bestD=None
for name,target in [('native_residual',G),('spectral_residual',spectral)]:
    D=target-M
    for Q in qs:
        h=Q.T@D@Q;ev,V=np.linalg.eigh((h+h.T)/2);F=Q@(V*np.sqrt(np.maximum(ev,0)));D-=F@F.T
    maxev=[float(np.linalg.eigvalsh(Q.T@D@Q)[-1]) for Q in qs]
    shift=max(0,max(maxev))+1e-11;D-=shift*np.eye(len(D))
    # Both proposals are certified against the direct native target. The
    # spectral residual is merely an alternative feasible dual multiplier.
    bound2=float(2*np.sum(G*D)-np.square(D).sum());finalev=[float(np.linalg.eigvalsh(Q.T@D@Q)[-1]) for Q in qs]
    proposals.append(dict(proposal=name,lower_bound_squared=bound2,block_max_eigenvalues=finalev,shift=shift))
    if max(finalev)<0 and bound2>best:best=bound2;bestD=D.copy();winner=proposals[-1]
D=bestD;bound2=best;finalev=winner['block_max_eigenvalues'];shift=winner['shift']
row=dict(fit=a.fit,shape=[R,B,C],target='direct_native_reference',primal_error=float(np.linalg.norm(D0)),dual_lower_bound_squared=bound2,dual_lower_bound=float(np.sqrt(max(0,bound2))),block_max_eigenvalues=finalev,numerical_safety_shift=shift,excludes_20mHa=bool(bound2>.020**2 and max(finalev)<0),dual_proposals=proposals)
np.savez_compressed(HERE/f'candidates/{Path(a.fit).stem}_dual.npz',D=D,tree_edges=tree,angle_words=words)
(HERE/f'results/{Path(a.fit).stem}_dual.json').write_text(json.dumps(row,indent=2));print(json.dumps(row))
