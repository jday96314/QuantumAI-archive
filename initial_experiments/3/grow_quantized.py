"""Enrich a coarse dictionary with neighboring quantized rank-one features.

Selection maximizes the exact two-sided Frobenius projection gain, including
both the new diagonal component and cross terms with the current subspace.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
import torch
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--target');p.add_argument('--pool')
p.add_argument('--basis',type=int,default=600);p.add_argument('--out',required=True);a=p.parse_args()
torch.set_num_threads(3);HERE=Path(__file__).resolve().parent;z=np.load(a.fit);u=z['U'][0];bits=int(z['brot']);tree=z['tree_edges']
if a.pool:
 pool=np.load(a.pool);assert int(pool['brot'])==bits and np.array_equal(pool['tree_edges'],tree)
 words=np.unique(pool['angle_words'].reshape(-1,75),axis=0)
else:
 assert a.target
 t=np.load(a.target)['U'][0].copy();th=np.zeros((len(t),75))
 for j in reversed(range(75)):
  pa,ch=tree[j];th[:,j]=np.arctan2(t[:,ch],t[:,pa]);t[:,pa]=np.hypot(t[:,pa],t[:,ch]);t[:,ch]=0
 units=th/np.pi*(1<<(bits-1));lower=np.floor(units);fraction=units-lower;rng=np.random.default_rng(719)
 tables=[np.rint(units)]
 for threshold in [.125,.25,.375,.625,.75,.875]:tables.append(lower+(fraction>threshold))
 for k in range(3):tables.append(lower+(fraction>rng.random(fraction.shape)))
 words=np.unique(np.concatenate(tables).astype(np.int64)%(1<<bits),axis=0)
uc=reconstruct_u_vectors(words/(1<<(bits-1)),tree);sp=np.load(HERE/'parameters/reference_spectrum.npz')
ii,jj,fac=sp['ii'],sp['jj'],sp['fac'];ev,V=sp['eigenvalues'],sp['eigenvectors']
d=dict(device='cuda',dtype=torch.float64)
P=torch.tensor(u[:,ii]*u[:,jj]*fac,**d);Pc=torch.tensor(uc[:,ii]*uc[:,jj]*fac,**d)
G=torch.tensor((V*np.maximum(ev,0))@V.T,**d);s=P@P.T;si=torch.linalg.inv(s)
ee,vv=torch.linalg.eigh(s);ir=(vv/ee.sqrt())@vv.T
c=P@Pc.T;t=si@c;h=P@G@P.T;hc=P@G@Pc.T;hcc=Pc@G@Pc.T
gr=Pc@Pc.T-c.T@t;hr=hcc-t.T@hc-hc.T@t+t.T@h@t
cross=ir@(hc-h@t);h0=ir@h@ir;const=float(np.square(ev).sum());error2=const-float(h0.square().sum())
chosen=[];start=time.time()
for step in range(a.basis-len(u)):
 diag=gr.diag().clamp(min=1e-15);hd=hr.diag()
 score=2*cross.square().sum(0)/diag+(hd/diag).square()
 score[diag<1e-10]=-1
 if chosen:score[torch.tensor(chosen,device='cuda')]=-1
 j=int(score.argmax());gain=float(score[j]);chosen.append(j);error2-=gain
 v=gr[:,j].clone();dj=float(diag[j]);beta=v/dj;hj=hr[:,j].clone();hjj=float(hr[j,j])
 newrow=(hj-hjj*beta)/np.sqrt(dj)
 cross=cross-cross[:,j:j+1]*beta[None,:];cross=torch.cat([cross,newrow[None,:]],0)
 hr=hr-beta[:,None]*hj[None,:]-hj[:,None]*beta[None,:]+hjj*beta[:,None]*beta[None,:]
 gr=gr-v[:,None]*v[None,:]/dj
 if step%10==0:print(json.dumps(dict(B=len(u)+step+1,error_2=float(np.sqrt(max(0,error2))),gain=gain,seconds=time.time()-start)),flush=True)
out=np.r_[u,uc[chosen]];basewords=z['angle_words'].reshape(len(u),75)
np.savez(a.out,U=out[None,:,:],angle_words=np.r_[basewords,words[chosen]],tree_edges=tree,brot=bits)
