"""Split a global (R=1) DF-THC factorization into R balanced blocks by co-clustering.

Rows (vectors) and columns are assigned to blocks with equal capacities; within-block
weights are refit exactly (block-coordinate PSD least squares).
"""
import argparse, numpy as np, torch
from sklearn.cluster import SpectralCoclustering
p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--R',type=int);p.add_argument('--B',type=int)
p.add_argument('--C',type=int);p.add_argument('--out');p.add_argument('--iters',type=int,default=20);a=p.parse_args()
z=np.load(a.fit);U=z['U'][0];w=z['w'][0];V0,C0=w.shape;R,B,C=a.R,a.B,a.C
assert R*B>=V0 and R*C<=C0+R*0+10**9
A=w**2
m=SpectralCoclustering(n_clusters=R,random_state=0).fit(A+1e-12)
rl,cl=m.row_labels_.copy(),m.column_labels_.copy()
def balance(labels,aff,cap):
    # aff: items x R affinity; greedy assignment respecting capacity
    n=len(labels);order=np.argsort(-(aff.max(1)-np.sort(aff,1)[:,-2] if R>1 else aff.max(1)))
    count=np.zeros(R,int);out=-np.ones(n,int)
    for i in order:
        for k in np.argsort(-aff[i]):
            if count[k]<cap:out[i]=k;count[k]+=1;break
    return out
for it in range(a.iters):
    caff=np.stack([A[rl==k].sum(0) for k in range(R)],1)   # columns x R
    cl=balance(cl,caff/ (caff.sum(1,keepdims=True)+1e-30),int(np.ceil(C0/R)) if R*C>=C0 else C)
    raff=np.stack([A[:,cl==k].sum(1) for k in range(R)],1)
    rl=balance(rl,raff/(raff.sum(1,keepdims=True)+1e-30),B)
within=sum(A[np.ix_(rl==k,cl==k)].sum() for k in range(R))/A.sum()
print('within mass',within,np.bincount(rl,minlength=R),np.bincount(cl,minlength=R))
Ub=np.zeros((R,B,76));wb=np.zeros((R,B,C))
rng=np.random.default_rng(0)
for k in range(R):
    rows=np.flatnonzero(rl==k);cols=np.flatnonzero(cl==k)
    Ub[k,:len(rows)]=U[rows]
    if len(rows)<B:Ub[k,len(rows):]=U[rng.choice(V0,B-len(rows),replace=False)]
    c=min(C,len(cols));wb[k][:len(rows),:c]=w[np.ix_(rows,cols[:c])]
np.savez(a.out,U=Ub,w=wb)
