"""Add residual-directed quantized atoms and reproject to a physical rank.

Only writes below the current experiment. The native full verifier is still
required after normalization/offset and one-body angle optimization.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
import torch
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors

p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--pool',default='candidates/selectswap7_anf/submission.npz')
p.add_argument('--out',required=True);p.add_argument('--basis',type=int,default=710);p.add_argument('--columns',type=int,default=350)
p.add_argument('--ridge',type=float,default=1e-6);a=p.parse_args();torch.set_num_threads(2)
root=Path(__file__).resolve().parent;legacy=(root/'../../3').resolve();z=dict(np.load(a.fit));pool=dict(np.load(a.pool))
tree=z['tree_edges'];assert np.array_equal(tree,pool['tree_edges']);bits=int(z['brot']);words=z['angle_words'].reshape(-1,len(tree))
pw=pool['angle_words'];pw=pw[-pool['U'].shape[0]*pool['U'].shape[1]:]
cw=np.unique(np.r_[words,pw],axis=0);lookup={tuple(w):i for i,w in enumerate(cw)};active=[lookup[tuple(w)] for w in words]
U=reconstruct_u_vectors(cw/(1<<(bits-1)),tree);sp=np.load(legacy/'parameters/reference_spectrum.npz');ev,V=sp['eigenvalues'],sp['eigenvectors']
d=dict(device='cuda',dtype=torch.float64);P=torch.tensor(U[:,sp['ii']]*U[:,sp['jj']]*sp['fac'],**d)
G=torch.tensor((V*ev)@V.T,**d);Gpos=torch.tensor((V*np.maximum(ev,0))@V.T,**d)
S=P@P.T;H=P@Gpos@P.T;Hactual=P@G@P.T;const=G.square().sum();began=time.time()
ids=torch.tensor(active,device='cuda');ss=S[ids][:,ids];hh=H[ids][:,ids];si=torch.linalg.inv(ss+a.ridge*torch.eye(len(ids),**d))
out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True)
while True:
 ids=torch.tensor(active,device='cuda');ss=S[ids][:,ids];hh=H[ids][:,ids];ha=Hactual[ids][:,ids]
 L=torch.linalg.cholesky(ss+a.ridge*torch.eye(len(ids),**d));white=torch.linalg.solve_triangular(L,hh,upper=False)
 white=torch.linalg.solve_triangular(L,white.T,upper=False).T;e,v=torch.linalg.eigh(.5*(white+white.T))
 w=torch.linalg.solve_triangular(L.T,v[:,-a.columns:]*e[-a.columns:].clamp(min=0).sqrt(),upper=True)
 K=w@w.T;KS=K@ss;err=float((const-2*(K*ha).sum()+(KS*KS.T).sum()).clamp(min=0).sqrt())
 name=out.with_name(out.stem+f'_b{len(ids)}.npz')
 np.savez(name,U=U[active][None,:,:],w=w.cpu().numpy()[None,:,:],W=np.zeros((1,a.columns)),tree_edges=tree,brot=bits,angle_words=cw[active],err2=err)
 print(json.dumps(dict(B=len(ids),error_2=err,gram_bound=float(.25*K.abs().sum()),path=str(name),seconds=time.time()-began)),flush=True)
 if len(active)>=a.basis:break
 cross=S[ids];t=si@cross;hc=H[ids];diag=(S.diag()+a.ridge-(cross*t).sum(0)).clamp(min=1e-18)
 nr=hc-hh@t;hr=H.diag()-2*(t*hc).sum(0)+(t*(hh@t)).sum(0)
 gain=2*(nr*(si@nr)).sum(0)/diag+(hr/diag).square();gain[diag<1e-9]=-1;gain[ids]=-1
 j=int(gain.argmax());gap=float(diag[j]);vec=t[:,j];N=len(active)
 nextsi=torch.empty((N+1,N+1),**d);nextsi[:N,:N]=si+torch.outer(vec,vec)/gap
 nextsi[N,:N]=-vec/gap;nextsi[:N,N]=-vec/gap;nextsi[N,N]=1/gap;si=nextsi;active.append(j)
