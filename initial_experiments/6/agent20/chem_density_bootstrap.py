"""Native signed-target density/Schur bootstrap in an exact angle dictionary.

The first76 square-root generators reproduce diagonal-density blocks and
all density/off-diagonal couplings exactly. Remaining generators come from
the signed Schur complement, which is truncated only with a reported error.
Each block then selects quantized atoms before any PSD fit, avoiding a final
rounding shock. This is a numerical experiment, not a native submission.
"""
import argparse,json,time,hashlib
from pathlib import Path
import numpy as np
import torch
from audit_core import Metric,encode,save
from chem_rank_screen import BANK,quantities
HERE=Path(__file__).resolve().parent


def select(m,pool,words,target,B,ridge):
    P=m.feats(pool);S=P@P.T;original=S.clone();cross=P@target.T;chosen=[]
    diag_floor=1e-8
    for _ in range(B):
        diag=S.diag().clamp(min=1e-13);score=cross.square().sum(1)/(diag+ridge)
        score[diag<diag_floor]=-1
        if chosen:score[chosen]=-1
        j=int(score.argmax());chosen.append(j)
        v=S[:,j].clone()/diag[j].sqrt();h=cross[j].clone()/diag[j].sqrt()
        S-=v[:,None]*v[None,:];cross-=v[:,None]*h[None,:]
    ix=torch.tensor(chosen,device=P.device);selected=P[ix]
    w=torch.linalg.lstsq(selected.T,target.T).solution
    return pool[ix],words[np.asarray(chosen)],w,float((w.T@selected-target).norm())


@torch.no_grad()
def main():
    p=argparse.ArgumentParser();p.add_argument('--B',type=int,default=240);p.add_argument('--C',type=int,default=88)
    p.add_argument('--ridge',type=float,default=.001);p.add_argument('--cycles',type=int,default=10);p.add_argument('--tag',default='density_R4B240C88');a=p.parse_args()
    torch.set_num_threads(1);torch.manual_seed(2026102510);start=time.time();z=dict(np.load(BANK));m=Metric(z['tree_edges'])
    ii=m.ii;jj=m.jj;diag=torch.nonzero(ii==jj).flatten();off=torch.nonzero(ii!=jj).flatten();G=m.G
    D=G[diag[:,None],diag[None,:]];C=G[diag[:,None],off[None,:]]
    de,V=torch.linalg.eigh(D);assert float(de.min())>0
    root=(V*de.sqrt())@V.T;invroot=(V/de.sqrt())@V.T
    Xd=torch.zeros((76,len(G)),**m.d);Xd[:,diag]=root;Xd[:,off]=invroot@C
    S=G[off[:,None],off[None,:]]-Xd[:,off].T@Xd[:,off]
    ev,E=torch.linalg.eigh((S+S.T)/2);K=4*a.C-76;take=torch.arange(len(ev)-K,len(ev),device=ev.device)
    Xr=torch.zeros((K,len(G)),**m.d);Xr[:,off]=(E[:,take]*ev[take].clamp(min=0).sqrt()).T
    X=torch.cat([Xd,Xr]);target=X.reshape(4,a.C,-1)
    floor=float((X.T@X-G).norm());records=[];urows=[];wrows=[];arows=[]
    inc=m.decode(m.tensor(z['angle_words'][76:]));coordinate=torch.eye(76,**m.d)
    for r in range(4):
        mats=m.unpack(target[r]);strength=target[r].norm(dim=1)
        order=strength.argsort(descending=True)[:16]
        _,vec=torch.linalg.eigh(mats[order]);vectors=vec.transpose(1,2).reshape(-1,76)
        # Physical Schur generators provide new directions; exact incumbent
        # atoms and coordinates are controls in the same quantized pool.
        raw=torch.cat([vectors,coordinate]);new_words=encode(raw.cpu().numpy(),z['tree_edges'],8,True)
        allwords=np.concatenate([new_words,z['angle_words'][76:]],0)
        decoded=m.decode(m.tensor(allwords))
        u,aw,w,error=select(m,decoded,allwords,target[r],a.B,a.ridge)
        urows.append(u);arows.append(aw);wrows.append(w)
        row=dict(block=r,factor_projection_error=error,weight_norm=float(w.norm()),seconds=time.time()-start)
        records.append(row);print(json.dumps(row),flush=True)
    u=torch.stack(urows);w=torch.stack(wrows);aw=np.stack(arows);offset=-w.sum(1)
    before=quantities(m,u,w,offset)
    w,_=m.project(u,w,a.cycles);m.project_gap(u,w,offset)
    after=quantities(m,u,w,offset)
    save(HERE/'parameters'/f'{a.tag}.npz',u,w,z['tree_edges'],aw,W=offset.cpu().numpy(),target=target.cpu().numpy())
    result=dict(status='FINITE_DENSITY_SCHUR_BOOTSTRAP',arguments=vars(a),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
       density_eigenvalue_range=[float(de.min()),float(de.max())],schur_norm=float(S.norm()),
       signed_schur_min_eigenvalue=float(ev.min()),square_root_truncation_error=floor,
       density_root_negative_entries=int((root<0).sum()),blocks=records,before=before,after=after,seconds=time.time()-start,
       scope='new exact-word initializer and signed-target PSD refit; no circuit/native acceptance')
    (HERE/'results'/f'{a.tag}.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)

if __name__=='__main__':main()
