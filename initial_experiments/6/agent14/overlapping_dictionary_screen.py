"""Shared-atom R2 dictionary screen using exact incumbent angle words.

Possibly indefinite unrestricted block-covariance projection is a lower
bound only for each saved dictionary. Shared intersections use the exact
Moore-Penrose normal-equation chart, with checked residuals.
"""
import sys,json,time
from pathlib import Path
ROOT=Path.cwd().parents[2]
for path in ['initial_experiments/5/chem-qfit','initial_experiments/5/chem-audit','initial_experiments/5/chem-census']:sys.path.insert(0,str(ROOT/path))
import numpy as np
import torch
from census_common import initial
start=time.time();z,m,u,w,off,words=initial();torch.set_num_threads(2);rows=[]
with torch.no_grad():
 P=m.feats(u)[0];G=m.G;K=w[0]@w[0].T;N=len(P);graph=K.square();graph.fill_diagonal_(0);hub=graph.sum(1);normw=K.diag();generator=torch.Generator(device=P.device);generator.manual_seed(142057)
 for B in [384,400,416,448]:
  common_count=2*B-N
  for method in ['hub','weight','random']:
   rank=torch.argsort(hub if method=='hub' else normw,descending=True) if method!='random' else torch.randperm(N,generator=generator,device=P.device);common=rank[:common_count];rest=rank[common_count:];gg=graph[rest][:,rest];dd=gg.sum(1).clamp(min=1e-30).rsqrt();gg=gg*dd[:,None]*dd[None,:];ev,V=torch.linalg.eigh(gg);rest=rest[torch.argsort(V[:,-2])];half=len(rest)//2;ids=[torch.cat([common,rest[:half]]),torch.cat([common,rest[half:]])];Q=[torch.linalg.qr(P[ii].T)[0] for ii in ids];C=Q[0].T@Q[1];G1=Q[0].T@G@Q[0];G2=Q[1].T@G@Q[1];lam,V=torch.linalg.eigh(C@C.T);den=1-lam[:,None]*lam[None,:];A=G1-C@G2@C.T;rhs=V.T@A@V;sol=torch.where(den.abs()>1e-10,rhs/den,torch.zeros_like(rhs));k1=V@sol@V.T;k2=G2-C.T@k1@C;fitted=Q[0]@k1@Q[0].T+Q[1]@k2@Q[1].T;res=G-fitted;error=float(res.norm());normal=max(float((Q[i].T@res@Q[i]).norm()) for i in range(2));assert normal<1e-8,normal;r=dict(B=B,common_atoms=common_count,method=method,unrestricted_error=error,normal_equation_residual=normal,null_chart_dimension=int((den.abs()<=1e-10).sum()),below_20mHa=error<.02,seconds=time.time()-start);rows.append(r);print(json.dumps(r),flush=True);np.savez(f'parameters/shared_B{B}_{method}.npz',ids=np.array([ii.cpu().numpy() for ii in ids]),scope='fixed exact-word dictionary, unrestricted covariance screen');Path('results/overlapping_dictionary_screen.json').write_text(json.dumps(rows,indent=2)+'\n')
