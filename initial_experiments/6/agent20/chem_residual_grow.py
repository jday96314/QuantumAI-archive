"""Grow exact-word multiblock dictionaries using missing generator gradients.

Unlike continuous multiblock rounding, every proposed atom is first encoded
and decoded on the actual 8-bit tree lattice. Missing feature directions are
the leading singular vectors of X_r*(G-X^T X), projected outside the current
dictionary span. Their matrix eigenvectors propose new rank-one atoms.
This is a finite capacity/conditioning screen, not native acceptance.
"""
import argparse,json,time,hashlib
from pathlib import Path
import numpy as np
import torch
from audit_core import Metric,load,save,encode
from chem_rank_screen import BANK,quantities

HERE=Path(__file__).resolve().parent

@torch.no_grad()
def main():
    p=argparse.ArgumentParser();p.add_argument('--init',type=Path,required=True)
    p.add_argument('--B',type=int,default=257);p.add_argument('--batch',type=int,default=4)
    p.add_argument('--cycles',type=int,default=6);p.add_argument('--tag',default='residual_grow_R4B257')
    a=p.parse_args();torch.set_num_threads(1);torch.manual_seed(2026110702)
    z,_,w,tree,words=load(a.init);R,B,C=w.shape;assert a.B>=B
    m=Metric(tree);aw=np.asarray(words).reshape(R,B,75).copy();u=m.decode(m.tensor(aw));wt=m.tensor(w)
    off=m.tensor(z['H2_I']-w.sum(1) if 'H2_I' in z else z.get('offset',z.get('W')))
    base=dict(np.load(BANK));iw=base['angle_words'][76:]
    if not np.array_equal(tree,base['tree_edges']):iw=encode(base['U'][0],tree,8,True)
    incumbent=m.decode(m.tensor(iw));start=time.time();history=[]
    def record(stage):
        stats=quantities(m,u,wt,off);row=dict(stage=stage,shape=list(wt.shape),seconds=time.time()-start,**stats)
        history.append(row);print(json.dumps(row),flush=True)
        save(HERE/'parameters'/f'{a.tag}_B{u.shape[1]}.npz',u,wt,tree,aw,W=off.cpu().numpy())
        (HERE/'results'/f'{a.tag}.json').write_text(json.dumps(dict(arguments={k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()},source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),initializer_sha256=hashlib.sha256(a.init.read_bytes()).hexdigest(),history=history),indent=2)+'\n')
    record('initial_actual_words')
    wt,_=m.project(u,wt,a.cycles);m.project_gap(u,wt,off);record('initial_PSD_refit')
    while u.shape[1]<a.B:
        P=m.feats(u);X=wt.transpose(1,2)@P;flat=X.flatten(0,1);residual=m.G-flat.T@flat
        add=min(a.batch,a.B-u.shape[1]);chosen_u=[];chosen_w=[]
        for r in range(R):
            Q,_=torch.linalg.qr(P[r].T);Z=X[r]@residual;Z-=Z@Q@Q.T
            _,_,V=torch.linalg.svd(Z,full_matrices=False)
            matrices=m.unpack(V[:6]);eig,vec=torch.linalg.eigh(matrices)
            order=eig.abs().argsort(-1,descending=True)[:,:4]
            proposals=torch.gather(vec,2,order[:,None,:].expand(-1,76,-1)).transpose(1,2).reshape(-1,76)
            pw=encode(proposals.cpu().numpy(),tree,8,True)
            pool_u=torch.cat([m.decode(m.tensor(pw)),incumbent],0)
            pool_words=np.concatenate([pw,iw],0)
            features=m.feats(pool_u);outside=features-features@Q@Q.T
            cross=outside@Z.T;selected=[]
            for _ in range(add):
                norm=outside.square().sum(1)
                score=cross.square().sum(1)/(norm+.05)
                score[norm<.01]=-1
                if selected:score[selected]=-1
                index=int(score.argmax());assert float(score[index])>=0
                selected.append(index);q=outside[index]/norm[index].sqrt()
                gain=outside@q;direction=(q@Z.T)
                outside-=gain[:,None]*q[None,:];cross-=gain[:,None]*direction[None,:]
            chosen_u.append(pool_u[selected]);chosen_w.append(pool_words[selected])
        u=torch.cat([u,torch.stack(chosen_u)],1);aw=np.concatenate([aw,np.stack(chosen_w)],1)
        wt=torch.cat([wt,torch.zeros((R,add,C),**m.d)],1)
        wt,_=m.project(u,wt,a.cycles);m.project_gap(u,wt,off);record('residual_atoms_and_PSD_refit')
    print(json.dumps(dict(status='FINITE_EXACT_WORD_GROWTH_COMPLETE',final=history[-1])),flush=True)

if __name__=='__main__':main()
