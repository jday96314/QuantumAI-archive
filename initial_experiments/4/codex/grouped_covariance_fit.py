"""Alternate exact covariance blocks with dictionary-vector gradient steps.

This numerical family optimizes PSD covariance matrices directly rather than
factoring them through fixed-rank weights during gradient updates. It is a
research fit; angle quantization and native acceptance remain separate gates.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
import torch


def main():
 p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True);p.add_argument('--steps',type=int,default=4000);p.add_argument('--lr',type=float,default=.0003);a=p.parse_args();torch.set_num_threads(2);begin=time.time()
 z=dict(np.load(a.input));sp=np.load('../../3/parameters/reference_spectrum.npz');d=dict(device='cuda',dtype=torch.float64)
 u=torch.tensor(z['U'],requires_grad=True,**d);R,B,N=u.shape;ii,jj=sp['ii'],sp['jj'];fac=torch.tensor(sp['fac'],**d);eye=torch.eye(B,**d)
 G=torch.tensor((sp['eigenvectors']*sp['eigenvalues'])@sp['eigenvectors'].T,**d);constant=G.square().sum();diag=torch.tensor((ii==jj).astype(float),**d);trace=diag@G@diag
 def features():
  U=u/u.norm(dim=2,keepdim=True);P=U[:,:,ii]*U[:,:,jj]*fac;L=torch.linalg.cholesky(P@P.transpose(1,2)+1e-12*eye);Q=torch.linalg.solve_triangular(L,P,upper=False);return U,P,L,Q
 with torch.no_grad():
  _,_,L,_=features();F=L.transpose(1,2)@torch.tensor(z['w'],**d);H=F@F.transpose(1,2)
 def contracted(Q,H):
  A=(Q.reshape(-1,Q.shape[-1])@Q.reshape(-1,Q.shape[-1]).T).reshape(R,B,R,B).permute(0,2,1,3).contiguous();T=Q@G@Q.transpose(1,2)
  cross=(A@H[None])@A.transpose(2,3)
  err=constant-2*(H*T).sum()+(H[:,None]*cross.transpose(2,3)).sum()
  return err,A,T
 def trace_scale(Q,H):
  t=Q@diag;current=(t*(H@t[:,:,None]).squeeze(-1)).sum();return H*(trace/current)
 opt=torch.optim.Adam([u],lr=a.lr);best=float('inf')
 for step in range(a.steps):
  with torch.no_grad():
   U,P,L,Q=features();_,A,T=contracted(Q,H)
   for r in range(R):
    cross=(A[r]@H)@A[r].transpose(1,2);target=T[r]-cross.sum(0)+cross[r];e,v=torch.linalg.eigh((target+target.T)/2);H[r]=(v*e.clamp(min=0))@v.T
   H=trace_scale(Q,H)
  opt.zero_grad();U,P,L,Q=features();scaled=trace_scale(Q,H);err,_,_=contracted(Q,scaled);err.backward();opt.step();opt.param_groups[0]['lr']=a.lr*(.1+.9*.5*(1+np.cos(np.pi*step/a.steps)))
  if step%20==0 or step==a.steps-1:
   with torch.no_grad():
    U,P,L,Q=features();HH=trace_scale(Q,H);er,_,_=contracted(Q,HH);error=float(er.clamp(min=0).sqrt())
    if error<best:
     best=error;e,v=torch.linalg.eigh(HH);F=v*e.clamp(min=0).sqrt()[:,None,:];w=torch.linalg.solve_triangular(L.transpose(1,2),F,upper=True)
     if step==0:
      X=(w.transpose(1,2)@P).reshape(R*B,-1);direct=(X.T@X-G).square().sum();assert abs(float(direct-er))<1e-8
     path=Path(a.output);tmp=path.with_name(path.stem+'.tmp.npz');np.savez(tmp,U=U.cpu().numpy(),w=w.cpu().numpy(),W=np.zeros((R,B)),tree_edges=z['tree_edges'],err2=error);tmp.replace(path)
    print(json.dumps(dict(step=step,R=R,B=B,C=B,error_2=error,best=best,seconds=time.time()-begin)),flush=True)
 print(json.dumps(dict(status='FINISHED',scope='unquantized alternating covariance/dictionary fit; no native acceptance',best_error_2=best,seconds=time.time()-begin)),flush=True)
if __name__=='__main__':main()
