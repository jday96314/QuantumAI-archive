"""Native signed-capacity checks for new positive-factor dictionaries.

The starting nonnegative fit is not a PSD optimum. These exact QR/eigen
block updates test both original rank and a full-rank relaxation, while a
union-span projection is a necessary lower bound. No native acceptance.
"""
import argparse,time,json,hashlib
from pathlib import Path
import numpy as np
import torch
from audit_core import Metric,load,save
HERE=Path(__file__).resolve().parent
@torch.no_grad()
def main():
 p=argparse.ArgumentParser();p.add_argument('--init',type=Path,required=True);p.add_argument('--tag',required=True);p.add_argument('--cycles',type=int,default=20);a=p.parse_args();torch.set_num_threads(1);start=time.time()
 z,u,w,tree,words=load(a.init);m=Metric(tree);u=m.decode(m.tensor(words.reshape(*w.shape[:2],75)));w=m.tensor(w);P=m.feats(u);Q,_=torch.linalg.qr(P.flatten(0,1).T);H=Q.T@m.G@Q
 floor=float((m.G.square().sum()-H.square().sum()).clamp(min=0).sqrt());rows=[]
 for relaxed in (False,True):
  ww=w if not relaxed else torch.cat([w,torch.zeros((*w.shape[:2],w.shape[1]-w.shape[2]),**m.d)],-1)
  for cycles in (1,4,a.cycles):
   ww,e=m.project(u,ww,cycles);rec=dict(relaxed=relaxed,C=ww.shape[-1],added_cycles=cycles,error2=e,negative_weights=int((ww<0).sum()),weight_norm=float(ww.norm()),seconds=time.time()-start)
   rows.append(rec);print(json.dumps(rec),flush=True)
  save(HERE/'parameters'/f'{a.tag}_{"relaxed" if relaxed else "ranked"}.npz',u,ww,tree,words,W=(-ww.sum(1)).cpu().numpy())
 out=dict(status='FINITE_POSITIVE_DICTIONARY_SIGNED_PSD_REFIT',arguments={k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()},union_span_lower_bound=floor,rows=rows,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='fixed exact dictionary only; full rank is a relaxation, finite block fits are upper bounds')
 (HERE/'results'/f'{a.tag}.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
