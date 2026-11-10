"""Exact PSD covariance block updates for a fixed grouped THC dictionary.

Whiten each rank-one dictionary. The least-squares update of one covariance
block is a PSD projection of its residual, optionally truncated to C columns.
This is a numerical initializer, not a native circuit or acceptance check.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
import torch


def main():
 p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True);p.add_argument('--columns',type=int,default=64);p.add_argument('--cycles',type=int,default=200);a=p.parse_args();torch.set_num_threads(2);begin=time.time()
 z=dict(np.load(a.input));sp=np.load('../../3/parameters/reference_spectrum.npz');d=dict(device='cuda',dtype=torch.float64)
 U=torch.tensor(z['U'],**d);R,B,N=U.shape;C=min(a.columns,B);ii,jj=sp['ii'],sp['jj'];fac=torch.tensor(sp['fac'],**d)
 P=U[:,:,ii]*U[:,:,jj]*fac;L=torch.linalg.cholesky(P@P.transpose(1,2)+1e-12*torch.eye(B,**d));Q=torch.linalg.solve_triangular(L,P,upper=False)
 G=torch.tensor((sp['eigenvectors']*sp['eigenvalues'])@sp['eigenvectors'].T,**d);constant=G.square().sum()
 A=(Q.reshape(-1,Q.shape[-1])@Q.reshape(-1,Q.shape[-1]).T).reshape(R,B,R,B).permute(0,2,1,3).contiguous()
 T=Q@G@Q.transpose(1,2);F=L.transpose(1,2)@torch.tensor(z['w'],**d);H=F@F.transpose(1,2)
 def error():
  term=torch.einsum('rsij,sjk,rslk->rsil',A,H,A)
  return (constant-2*(H*T).sum()+(H[:,None]*term.transpose(2,3)).sum()).clamp(min=0).sqrt()
 best=float(error());bestH=H.clone();initial=best;print(json.dumps(dict(cycle=-1,error_2=best,R=R,B=B,C=C)),flush=True)
 for cycle in range(a.cycles):
  for r in range(R):
   cross=(A[r]@H)@A[r].transpose(1,2)
   residual=T[r]-cross.sum(0)+cross[r];residual=(residual+residual.T)/2
   e,v=torch.linalg.eigh(residual);e=e[-C:].clamp(min=0);v=v[:,-C:]
   H[r]=(v*e)@v.T
  if cycle%10==0 or cycle==a.cycles-1:
   er=float(error())
   if er<best:best=er;bestH=H.clone()
   print(json.dumps(dict(cycle=cycle,error_2=er,best=best,seconds=time.time()-begin)),flush=True)
 e,v=torch.linalg.eigh(bestH);F=v[:,:,-C:]*e[:,-C:].clamp(min=0).sqrt()[:,None,:];w=torch.linalg.solve_triangular(L.transpose(1,2),F,upper=True)
 X=(w.transpose(1,2)@P).reshape(R*C,-1);direct=float((X.T@X-G).square().sum().sqrt());assert abs(direct-best)<1e-7,(direct,best)
 np.savez(a.output,U=z['U'],w=w.cpu().numpy(),W=np.zeros((R,C)),tree_edges=z['tree_edges'],err2=direct)
 report=dict(status='FINISHED',scope='fixed-dictionary PSD block fit; no native acceptance',R=R,B=B,C=C,starting_error_2=initial,best_error_2=direct,seconds=time.time()-begin)
 Path('results/'+Path(a.output).stem+'.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if __name__=='__main__':main()
