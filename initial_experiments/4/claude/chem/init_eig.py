"""Initial block vectors from eigenvectors of the target's leading pair-space directions."""
import argparse, numpy as np
p=argparse.ArgumentParser();p.add_argument('--R',type=int);p.add_argument('--B',type=int);p.add_argument('--C',type=int)
p.add_argument('--out');p.add_argument('--interleave',action='store_true');a=p.parse_args()
sp=np.load('/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/parameters/reference_spectrum.npz')
ev,V,ii,jj,fac=sp['eigenvalues'],sp['eigenvectors'],sp['ii'],sp['jj'],sp['fac']
order=np.argsort(ev)[::-1];n=76;R,B,C=a.R,a.B,a.C
K=R*C
dirs=order[:K]
groups=[dirs[r::R] for r in range(R)] if a.interleave else [dirs[r*C:(r+1)*C] for r in range(R)]
U=np.zeros((R,B,n))
for r,g in enumerate(groups):
    cands=[];scores=[]
    for c in g:
        m=np.zeros((n,n));m[ii,jj]=V[:,c]/fac;m=(m+m.T)-np.diag(np.diag(m))
        mu,vec=np.linalg.eigh(m)
        for i in range(n):cands.append(vec[:,i]);scores.append(np.sqrt(max(ev[c],0))*abs(mu[i]))
    idx=np.argsort(scores)[::-1][:B];U[r]=np.array(cands)[idx]
np.savez(a.out,U=U,w=np.zeros((R,B,C)))
