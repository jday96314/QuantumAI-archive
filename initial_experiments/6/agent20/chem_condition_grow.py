"""Grow an exact dictionary to repair cancellation, allowing in-span directions.

Ridge generator residuals price atoms that replace expensive coefficient
cancellations. Unlike outside-span-only growth, useful redundant directions
are eligible; exact duplicate projectors are excluded. Original target
factors and words stay fixed; numerical native-formula evidence is recorded
for every ridge and PSD checkpoint. This is not a submission acceptance.
"""
import argparse,json,time,hashlib
from pathlib import Path
import numpy as np
import torch
from audit_core import Metric,load,encode,save
from chem_rank_screen import BANK,quantities
HERE=Path(__file__).resolve().parent


@torch.no_grad()
def main():
    p=argparse.ArgumentParser();p.add_argument('--init',type=Path,required=True);p.add_argument('--B',type=int,default=511)
    p.add_argument('--batch',type=int,default=4);p.add_argument('--ridge',type=float,default=.0001)
    p.add_argument('--tag',default='condition_R2B511');a=p.parse_args();torch.set_num_threads(1)
    torch.manual_seed(2026102511);start=time.time();z,raw,ww,tree,words=load(a.init);R,B,C=ww.shape
    m=Metric(tree);aw=np.asarray(words).reshape(R,B,75);u=m.decode(m.tensor(aw));w=m.tensor(ww)
    off=m.tensor(z['H2_I']-ww.sum(1) if 'H2_I' in z else z['W']);P=m.feats(u);target=w.transpose(1,2)@P
    iz=dict(np.load(BANK));iw=iz['angle_words'][76:]
    if not np.array_equal(tree,iz['tree_edges']):iw=encode(iz['U'][0],tree,8,True)
    incumbent=m.decode(m.tensor(iw));history=[]
    def record(stage,wt,ot):
        row=dict(stage=stage,B=u.shape[1],seconds=time.time()-start,**quantities(m,u,wt,ot));history.append(row);print(json.dumps(row),flush=True)
        (HERE/'results'/f'{a.tag}.json').write_text(json.dumps(dict(arguments={k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()},source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),initializer_sha256=hashlib.sha256(a.init.read_bytes()).hexdigest(),history=history),indent=2)+'\n')
    record('original_exact_words',w,off)
    while u.shape[1]<a.B:
        P=m.feats(u);newu=[];newaw=[];add=min(a.batch,a.B-u.shape[1])
        for r in range(R):
            S=P[r]@P[r].T;ridgew=torch.linalg.solve(S+a.ridge*torch.eye(len(S),**m.d),P[r]@target[r].T)
            residual=target[r]-ridgew.T@P[r]
            _,_,V=torch.linalg.svd(residual,full_matrices=False)
            eig,vec=torch.linalg.eigh(m.unpack(V[:8]));ix=eig.abs().argsort(-1,descending=True)[:,:8]
            proposed=torch.gather(vec,2,ix[:,None,:].expand(-1,76,-1)).transpose(1,2).reshape(-1,76)
            pw=encode(proposed.cpu().numpy(),tree,8,True);pu=m.decode(m.tensor(pw));pool=torch.cat([pu,incumbent]);pwords=np.concatenate([pw,iw]);features=m.feats(pool)
            near=(pool@u[r].T).square().amax(1)>.9999
            # Ridge Schur complement: redundant but nonidentical atoms can
            # lower coefficient cost even if their exact-span residual is zero.
            cross=features@P[r].T;solve=torch.linalg.solve(S+a.ridge*torch.eye(len(S),**m.d),cross.T)
            A=features@features.T-cross@solve;rhs=features@residual.T;chosen=[]
            for _ in range(add):
                diag=A.diag().clamp(min=1e-12);score=rhs.square().sum(1)/(diag+a.ridge)
                score[near]=-1
                if chosen:score[chosen]=-1
                j=int(score.argmax());assert float(score[j])>=0;chosen.append(j)
                col=A[:,j].clone()/(diag[j]+a.ridge).sqrt();h=rhs[j].clone()/(diag[j]+a.ridge).sqrt()
                A-=col[:,None]*col[None,:];rhs-=col[:,None]*h[None,:]
            newu.append(pool[chosen]);newaw.append(pwords[chosen])
        u=torch.cat([u,torch.stack(newu)],1);aw=np.concatenate([aw,np.stack(newaw)],1);w=torch.cat([w,torch.zeros((R,add,C),**m.d)],1)
        record('added_zero_weight_atoms',w,off)
    P=m.feats(u)
    for ridge in [0.,1e-8,1e-7,1e-6,1e-5,1e-4,.001]:
        fitted=[]
        for r in range(R):
            if ridge:wr=torch.linalg.solve(P[r]@P[r].T+ridge*torch.eye(a.B,**m.d),P[r]@target[r].T)
            else:wr=torch.linalg.lstsq(P[r].T,target[r].T).solution
            fitted.append(wr)
        wt=torch.stack(fitted);ot=off.clone();m.project_gap(u,wt,ot);record(f'fixed_generators_ridge_{ridge:g}',wt,ot)
        save(HERE/'parameters'/f'{a.tag}_ridge_{ridge:g}.npz',u,wt,tree,aw,W=ot.cpu().numpy())
        for cycles in [4,16]:
            wt,_=m.project(u,wt,cycles,ridge=ridge);m.project_gap(u,wt,ot);record(f'PSD_ridge_{ridge:g}_cycles_{cycles}',wt,ot)
            save(HERE/'parameters'/f'{a.tag}_PSD_{ridge:g}_{cycles}.npz',u,wt,tree,aw,W=ot.cpu().numpy())
    print(json.dumps(dict(status='FINITE_CONDITIONING_GROWTH_COMPLETE',seconds=time.time()-start)),flush=True)

if __name__=='__main__':main()
