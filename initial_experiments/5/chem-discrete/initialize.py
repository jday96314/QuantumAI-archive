"""Combinatorial block support and factor-gauge initialization on exact words.

Alternate balanced factor subspaces and dictionary row subset selection.
Orbital vectors are never trained: all atoms are exact incumbent integer rows.
The factor-space surrogate proposes dictionaries; full rank-C tensor refits
measure their actual quality and retain every improved initialization.
"""
import argparse,json,time
from scipy.cluster.hierarchy import linkage,leaves_list
from core import *
p=argparse.ArgumentParser();p.add_argument('--B',type=int,default=479);p.add_argument('--rounds',type=int,default=16)
p.add_argument('--seed',type=int,default=0);p.add_argument('--tensor',action='store_true');a=p.parse_args()
torch.manual_seed(a.seed);z,un,wn,tree,words=load(INC/'submission.npz');m=Metric(tree)
u=m.decode(m.tensor(words));P=m.feats(u);w=m.tensor(wn[0]);F=w.T@P;C=175;N=len(P)
support=abs(wn[0]);support/=np.linalg.norm(support,axis=0,keepdims=True)
ids=leaves_list(linkage(support.T,method='ward')) if a.seed==0 else np.random.default_rng(a.seed).permutation(350)
rot=torch.eye(350,**m.d)[:,torch.tensor(ids.copy(),device='cuda')];start=time.time();best=float('inf');traj=[]
Sn=(P@P.T).cpu().numpy();ivn=np.linalg.inv(Sn)
@torch.no_grad()
def subset(wg):
    ix=np.arange(N);iv=ivn.copy();ww=wg.cpu().numpy().copy();S=Sn.copy()
    while len(ix)>a.B:
        if a.tensor:
            K=ww@ww.T;score=2*np.einsum('ij,ij->i',K@S,K)/iv.diagonal()-(K.diagonal()/iv.diagonal())**2
        else:score=np.square(ww).sum(-1)/iv.diagonal()
        b=int(score.argmin());keep=np.arange(len(ix))!=b
        col=iv[keep,b]/iv[b,b];ww=ww[keep]-col[:,None]*ww[b]
        iv=iv[keep][:,keep]-iv[keep,b,None]*iv[b,keep][None,:]/iv[b,b]
        S=S[keep][:,keep];ix=ix[keep]
    return torch.tensor(ix,device='cuda'),m.tensor(ww)
with torch.no_grad():
 for step in range(a.rounds):
    weights=w@rot;keep=[];ws=[];qs=[]
    for r in range(2):
        ix,ww=subset(weights[:,r*C:(r+1)*C]);keep.append(ix);ws.append(ww);qs.append(torch.linalg.qr(P[ix].T).Q)
    ut=torch.stack([u[ix] for ix in keep]);wt=torch.stack(ws);aw=torch.stack([m.tensor(words)[ix] for ix in keep])
    wt,er=m.project(ut,wt,20)
    tag=f'init_B{a.B}_seed{a.seed}_'+('tensor' if a.tensor else 'factor')
    save(HERE/f'trajectories/{tag}_{step}.npz',ut,wt,tree,aw,err2=er)
    if er<best:best=er;save(HERE/f'candidates/{tag}.npz',ut,wt,tree,aw,err2=er)
    Hs=[(F@q)@(F@q).T for q in qs];ev,V=torch.linalg.eigh(Hs[0]-Hs[1]);rot=torch.cat([V[:,-C:],V[:,:C]],1)
    row=dict(step=step,B=a.B,seed=a.seed,tensor=a.tensor,error_2=er,best=best,surrogate=float((F.square().sum()-((rot[:,:C].T@F@qs[0]).square().sum()+(rot[:,C:].T@F@qs[1]).square().sum()))),seconds=time.time()-start)
    print(json.dumps(row),flush=True);traj.append(row);(HERE/f'results/{tag}.json').write_text(json.dumps(traj,indent=2))
