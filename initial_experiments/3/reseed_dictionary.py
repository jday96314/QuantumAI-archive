"""Replace redundant lattice atoms with quantized residual eigendirections."""
import argparse,io,json,time,zipfile
from pathlib import Path
import numpy as np
import torch
from compress_chem import encode
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--out',required=True)
p.add_argument('--modes',type=int,default=32);p.add_argument('--per-mode',type=int,default=8)
p.add_argument('--replace',type=int,default=16,help='also replace the least useful existing atoms')
a=p.parse_args();torch.set_num_threads(1);root=Path(__file__).resolve().parent
for _ in range(60):
 try:
  with np.load(io.BytesIO(Path(a.fit).read_bytes())) as z:z={k:z[k] for k in z.files}
  break
 except (OSError,ValueError,zipfile.BadZipFile,EOFError):time.sleep(.05)
u=z['U'][0];target=len(u);bits=int(z['brot']);tree=z['tree_edges'];words=z['angle_words'].reshape(-1,75)
sp=np.load(root/'parameters/reference_spectrum.npz');ii,jj,fac=sp['ii'],sp['jj'],sp['fac']
P=u[:,ii]*u[:,jj]*fac;_,ix=np.unique(np.round(P,10),axis=0,return_index=True);ix=np.sort(ix)
u,words,P=u[ix],words[ix],P[ix];print(json.dumps(dict(original=target,unique=len(u))),flush=True)
d=dict(device='cuda',dtype=torch.float64);pt=torch.tensor(P,**d);ev,V=sp['eigenvalues'],sp['eigenvectors']
G=torch.tensor((V*np.maximum(ev,0))@V.T,**d);S=pt@pt.T;si=torch.linalg.inv(S+1e-10*torch.eye(len(S),**d));H=pt@G@pt.T;K=si@H@si
chosen=torch.arange(len(u),device='cuda')
while len(chosen)>target-a.replace:
 cap=si@H;kk=cap@si;diag=si.diag();loss=2*(cap*kk.T).sum(1)/diag-(kk.diag()/diag).square()
 remove=int(loss.argmin());keep=torch.arange(len(chosen),device='cuda')!=remove
 si=(si-torch.outer(si[:,remove],si[remove,:])/si[remove,remove])[keep][:,keep]
 H=H[keep][:,keep];chosen=chosen[keep]
ix=chosen.cpu().numpy();u,words,P=u[ix],words[ix],P[ix];pt=pt[chosen];S=pt@pt.T;K=si@H@si
print(json.dumps(dict(retained=len(u),replacing=target-len(u))),flush=True)
residual=G-pt.T@K@pt;e,v=torch.linalg.eigh(.5*(residual+residual.T));vf=v[:,-a.modes:].T.cpu().numpy()
mat=np.zeros((len(vf),76,76));mat[:,ii,jj]=vf/fac;mat=mat+mat.transpose(0,2,1)-np.eye(76)[None,:,:]*mat
ee,vv=np.linalg.eigh(mat);order=np.argsort(abs(ee),axis=1)[:,-a.per_mode:]
vectors=np.stack([vv[r,:,order[r,j]] for r in range(len(vv)) for j in range(a.per_mode)])
work=vectors.copy();theta=np.zeros((len(work),75))
for j in reversed(range(75)):
 pa,ch=tree[j];theta[:,j]=np.arctan2(work[:,ch],work[:,pa]);work[:,pa]=np.hypot(work[:,pa],work[:,ch]);work[:,ch]=0
units=theta/np.pi*(1<<(bits-1));low=np.floor(units);fraction=units-low;rng=np.random.default_rng(711029)
tables=[np.rint(units)]+[low+(fraction>t) for t in [.125,.25,.375,.625,.75,.875]]
tables.extend(low+(fraction>rng.random(fraction.shape)) for _ in range(3))
cw=np.unique(np.concatenate(tables).astype(np.int64)%(1<<bits),axis=0)
uc=reconstruct_u_vectors(cw/(1<<(bits-1)),tree);pc=torch.tensor(uc[:,ii]*uc[:,jj]*fac,**d)
np.savez(str(Path(a.out).with_suffix(''))+'_pool.npz',U=uc[None,:,:],angle_words=cw,tree_edges=tree,brot=bits)
cross=pt@pc.T;t=si@cross;hc=pt@G@pc.T;hcc=pc@G@pc.T;gr=pc@pc.T-cross.T@t
hr=hcc-t.T@hc-hc.T@t+t.T@H@t;es,vs=torch.linalg.eigh(S);ir=(vs/es.clamp(min=1e-12).sqrt())@vs.T
cg=ir@(hc-H@t);h0=ir@H@ir;error2=float(G.square().sum()-h0.square().sum());chosen=[]
for step in range(target-len(u)):
 diag=gr.diag().clamp(min=1e-15);hd=hr.diag();score=2*cg.square().sum(0)/diag+(hd/diag).square();score[diag<1e-10]=-1
 if chosen:score[torch.tensor(chosen,device='cuda')]=-1
 j=int(score.argmax());chosen.append(j);gain=float(score[j]);error2-=gain
 column=gr[:,j].clone();dj=float(diag[j]);beta=column/dj;hj=hr[:,j].clone();hjj=float(hr[j,j])
 newrow=(hj-hjj*beta)/np.sqrt(dj);cg=cg-cg[:,j:j+1]*beta[None,:];cg=torch.cat([cg,newrow[None,:]],0)
 hr=hr-beta[:,None]*hj[None,:]-hj[:,None]*beta[None,:]+hjj*beta[:,None]*beta[None,:];gr=gr-column[:,None]*column[None,:]/dj
 print(json.dumps(dict(B=len(u)+step+1,error_2=float(np.sqrt(max(0,error2))),gain=gain)),flush=True)
np.savez(a.out,U=np.r_[u,uc[chosen]][None,:,:],angle_words=np.r_[words,cw[chosen]],tree_edges=tree,brot=bits)
