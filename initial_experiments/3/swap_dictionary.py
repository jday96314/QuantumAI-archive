"""Monotone single-atom swaps from a residual-directed lattice candidate pool."""
import argparse,json,time
from pathlib import Path
import numpy as np
import torch
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--pool',required=True);p.add_argument('--out',required=True)
p.add_argument('--swaps',type=int,default=64);a=p.parse_args();torch.set_num_threads(1)
root=Path(__file__).resolve().parent;z=dict(np.load(a.fit));pool=np.load(a.pool);tree=z['tree_edges'];bits=int(z['brot'])
assert np.array_equal(tree,pool['tree_edges']) and bits==int(pool['brot'])
words=z['angle_words'].reshape(-1,75);cw=np.unique(np.r_[words,pool['angle_words'].reshape(-1,75)],axis=0)
uc=reconstruct_u_vectors(cw/(1<<(bits-1)),tree);lookup={tuple(w):j for j,w in enumerate(cw)};active=[lookup[tuple(w)] for w in words]
sp=np.load(root/'parameters/reference_spectrum.npz');ii,jj,fac=sp['ii'],sp['jj'],sp['fac'];d=dict(device='cuda',dtype=torch.float64)
pc=torch.tensor(uc[:,ii]*uc[:,jj]*fac,**d);e,V=sp['eigenvalues'],sp['eigenvectors'];G=torch.tensor((V*np.maximum(e,0))@V.T,**d);const=G.square().sum()
CC=pc@pc.T;HC=pc@G@pc.T;idx=torch.tensor(active,device='cuda');S=CC[idx][:,idx];H=HC[idx][:,idx]
si=torch.linalg.inv(S+1e-10*torch.eye(len(S),**d));tried=set();accepted=0;began=time.time()
for step in range(a.swaps):
 cap=si@H;K=cap@si;dd=si.diag();loss=2*(cap*K.T).sum(1)/dd-(K.diag()/dd).square()
 previous=(si,H,S,active);KS=K@S;prior_error=float((const-2*(K*H).sum()+(KS*KS.T).sum()).clamp(min=0).sqrt())
 for j,atom in enumerate(active):
  if atom in tried:loss[j]=float('inf')
 remove=int(loss.argmin());old=active[remove];tried.add(old)
 keep=torch.arange(len(active),device='cuda')!=remove;sib=(si-torch.outer(si[:,remove],si[remove,:])/si[remove,remove])[keep][:,keep]
 hb=H[keep][:,keep];base=[atom for j,atom in enumerate(active) if j!=remove];ix=torch.tensor(base,device='cuda')
 cross=CC[ix];t=sib@cross;hc=HC[ix];diag=(CC.diag()-(cross*t).sum(0)).clamp(min=1e-18)
 nr=hc-hb@t;hr=HC.diag()-2*(t*hc).sum(0)+(t*(hb@t)).sum(0)
 gain=2*(nr*(sib@nr)).sum(0)/diag+(hr/diag).square();gain[diag<1e-9]=-1
 selected=int(gain.argmax());oldgain=float(gain[old]);newgain=float(gain[selected])
 if newgain<=oldgain+1e-12:selected=old
 else:accepted+=1;tried.add(selected)
 vec=t[:,selected];gap=float(diag[selected]);newsi=torch.empty_like(si)
 newsi[:-1,:-1]=sib+torch.outer(vec,vec)/gap;newsi[:-1,-1]=-vec/gap;newsi[-1,:-1]=-vec/gap;newsi[-1,-1]=1/gap;si=newsi
 active=base+[selected];idx=torch.tensor(active,device='cuda');H=HC[idx][:,idx];S=CC[idx][:,idx]
 cap=si@H;K=cap@si;KS=K@S;error=float((const-2*(K*H).sum()+(KS*KS.T).sum()).clamp(min=0).sqrt())
 if error>prior_error+1e-10:
  if selected!=old:accepted-=1
  si,H,S,active=previous;error=prior_error
  cap=si@H;K=cap@si
 out=dict(U=uc[active][None,:,:],angle_words=cw[active],tree_edges=tree,brot=bits,err2=error)
 np.savez(a.out,**out)
 print(json.dumps(dict(step=step,error_2=error,accepted=accepted,gain=newgain-oldgain,gram_bound=float(.25*K.abs().sum()),seconds=time.time()-began)),flush=True)
