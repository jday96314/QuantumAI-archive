"""Select a compact shared dictionary from complementary quantized fits."""
import argparse,io,json,time,zipfile
from pathlib import Path
import numpy as np
import torch
p=argparse.ArgumentParser();p.add_argument('--first',required=True);p.add_argument('--second',required=True)
p.add_argument('--out',required=True);p.add_argument('--basis',type=int,default=700)
p.add_argument('--bias',type=float,default=0.,help='removal preference for second-source rows')
p.add_argument('--target-rank',type=int,help='select for the leading target modes')
a=p.parse_args();torch.set_num_threads(1);root=Path(__file__).resolve().parent
def snapshot(path):
 for _ in range(60):
  try:
   with np.load(io.BytesIO(Path(path).read_bytes())) as z:return {k:z[k] for k in z.files}
  except (OSError,ValueError,zipfile.BadZipFile,EOFError):time.sleep(.05)
 raise RuntimeError('Could not obtain a coherent snapshot')
z=snapshot(a.first);v=snapshot(a.second);assert int(z['brot'])==int(v['brot'])
tree=z['tree_edges'];assert np.array_equal(tree,v['tree_edges'])
u=np.concatenate([z['U'][0],v['U'][0]]);words=np.concatenate([z['angle_words'].reshape(-1,75),v['angle_words'].reshape(-1,75)])
source=np.r_[np.zeros(len(z['U'][0])),np.ones(len(v['U'][0]))]
sp=np.load(root/'parameters/reference_spectrum.npz');ii,jj,fac=sp['ii'],sp['jj'],sp['fac']
P=u[:,ii]*u[:,jj]*fac
_,unique=np.unique(np.round(P,10),axis=0,return_index=True);unique=np.sort(unique)
u,words,source,P=u[unique],words[unique],source[unique],P[unique]
d=dict(device='cuda',dtype=torch.float64);pt=torch.tensor(P,**d)
e,V=sp['eigenvalues'],sp['eigenvectors'];positive=np.maximum(e,0)
Gfull=torch.tensor((V*positive)@V.T,**d)
if a.target_rank:
 selected=np.argsort(e)[-a.target_rank:];G=torch.tensor((V[:,selected]*positive[selected])@V[:,selected].T,**d)
else:G=Gfull
const=G.square().sum()
S=pt@pt.T;H=pt@G@pt.T;si=torch.linalg.inv(S+1e-10*torch.eye(len(S),**d))
chosen=torch.arange(len(u),device='cuda');src=torch.tensor(source,**d);began=time.time()
while len(chosen)>a.basis:
 cap=si@H;K=cap@si;diag=si.diag();ksk=(cap*K.T).sum(1)
 loss=2*ksk/diag-(K.diag()/diag).square();score=loss-a.bias*src
 ix=int(score.argmin());keep=torch.arange(len(chosen),device='cuda')!=ix
 si=(si-torch.outer(si[:,ix],si[ix,:])/si[ix,ix])[keep][:,keep]
 H=H[keep][:,keep];chosen=chosen[keep];src=src[keep]
 if len(chosen)%50==0:
  captured=si@H;err=float((const-(captured*captured.T).sum()).clamp(min=0).sqrt())
  print(json.dumps(dict(B=len(chosen),error_2=err,second_rows=int(src.sum()),last_loss=float(loss[ix]),seconds=time.time()-began)),flush=True)
ids=chosen.cpu().numpy();pp=pt[chosen];ss=pp@pp.T;cap=si@H;K=cap@si;KS=K@ss
Hfull=pp@Gfull@pp.T
err=float((Gfull.square().sum()-2*(K*Hfull).sum()+(KS*KS.T).sum()).clamp(min=0).sqrt())
np.savez(a.out,U=u[ids][None,:,:],angle_words=words[ids],tree_edges=tree,brot=int(z['brot']),err2=err)
print(json.dumps(dict(final_error=err,gram_bound=float(.25*K.abs().sum()),second_rows=int(src.sum()))),flush=True)
