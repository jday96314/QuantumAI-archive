"""Check realistic rank increases on exact integer dictionaries; no float fit."""
import json
from core import *
for tag in ['legacy_exact_B700','anneal_B671']:
    z,u,w,tree,words=load(HERE/f'candidates/{tag}.npz');m=Metric(tree);ut=m.tensor(u);wt=m.tensor(w)
    for C in [350,360,365,375,400]:
        ww=torch.zeros((1,u.shape[1],C),**m.d);ww[:,:,:350]=wt
        ww,err=m.project(ut,ww,1)
        save(HERE/f'candidates/{tag}_C{C}.npz',ut,ww,tree,words,err2=err)
        print(json.dumps(dict(tag=tag,C=C,error=err)),flush=True)
