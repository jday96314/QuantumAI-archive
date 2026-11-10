"""Align rank-compressed physical generators to original column subsets."""
import json
from pathlib import Path
import numpy as np
import torch
from audit_core import Metric
from chem_rank_screen import quantities,BANK,HERE


if __name__=='__main__':
    torch.set_num_threads(1);z=dict(np.load(BANK));m=Metric(z['tree_edges'])
    words=z['angle_words'][76:];u=m.decode(m.tensor(words)[None]);w0=m.tensor(z['W'])
    P=m.feats(u)[0];X0=w0[0].T@P
    rankings=dict(first=torch.arange(350,device=u.device),
                  physical=X0.square().sum(1).argsort(descending=True),
                  weight=w0[0].square().sum(0).argsort(descending=True))
    reports=[]
    with torch.no_grad():
        for C in [348,344,342,340]:
            old=dict(np.load(HERE/'parameters'/f'chem_rank_C{C}.npz'))
            w=m.tensor(old['w']);off=m.tensor(old['W']);best=float('inf')
            for name,order in rankings.items():
                warm=w0[0][:,order[:C]]
                a,_,b=torch.linalg.svd(w[0].T@warm,full_matrices=False);rotation=a@b
                neww=w@rotation;newoff=off@rotation
                row=dict(C=C,gauge=name,**quantities(m,u,neww,newoff))
                reports.append(row);print(json.dumps(row),flush=True)
                if row['lambda_eff']<best:
                    best=row['lambda_eff']
                    np.savez(HERE/'parameters'/f'chem_rank_C{C}_aligned.npz',
                             U=u.cpu().numpy(),w=neww.cpu().numpy(),W=newoff.cpu().numpy(),
                             tree_edges=z['tree_edges'],angle_words=words,brot=8)
            (HERE/'results/chem_rank_gauges.json').write_text(json.dumps(reports,indent=2)+'\n')
