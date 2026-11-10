"""Optimize integer one-body rows jointly against their weighted matrix."""
import argparse,json,time
from pathlib import Path
import numpy as np
import torch
from compress_chem import encode
p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--out',required=True)
p.add_argument('--brot',type=int,default=8);p.add_argument('--steps',type=int,default=3000);p.add_argument('--target',type=float,default=0.);a=p.parse_args()
torch.set_num_threads(2);root=Path(__file__).resolve().parent;z=dict(np.load(a.fit));ref=np.load(Path('/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/parameters')/'chem_reference.npz')
tree=z.get('tree_edges',np.load(Path('/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/candidates')/'chem35c45_rot.npz')['tree_edges'])
u,w,W=z['U'],z['w'],z['W'];X=np.einsum('rbp,rbq,rbc->rcpq',u,u,w,optimize=True)
h=ref['h1_thc']-.5*int(ref['electrons'])*ref['S']+float(ref['beta_1'])*np.eye(76)-np.einsum('rc,rcpq->pq',W,X)
ev,v=np.linalg.eigh(.5*(h+h.T));words=encode(v.T,tree,a.brot)
d=dict(device='cuda',dtype=torch.float64);theta=torch.tensor(words,**d);e=torch.tensor(ev,**d);H=torch.tensor(h,**d)
alpha=np.pi/(1<<(a.brot-1));start=torch.zeros((76,76),**d);start[:,0]=1
def state():
    v=start.clone();J=torch.zeros((76,75,76),**d);Q=J.clone()
    for t,(pa,ch) in enumerate(tree):
        pa,ch=int(pa),int(ch);ang=theta[:,t]*alpha;c,s=ang.cos(),ang.sin()
        for D in [J,Q]:
            jp,jc=D[:,:,pa].clone(),D[:,:,ch].clone();D[:,:,pa]=c[:,None]*jp-s[:,None]*jc;D[:,:,ch]=s[:,None]*jp+c[:,None]*jc
        vp,vc=v[:,pa].clone(),v[:,ch].clone();v[:,pa]=c*vp-s*vc;v[:,ch]=s*vp+c*vc
        J[:,t,pa]=-v[:,ch];J[:,t,ch]=v[:,pa];Q[:,t,pa]=-v[:,pa];Q[:,t,ch]=-v[:,ch]
    residual=(v.T*e)@v-H
    return v,J,Q,residual
began=time.time();accepted=0;best=float('inf')
with torch.no_grad():
 for step in range(a.steps):
    u,J,Q,R=state();loss=R.square().sum();err=float(loss.sqrt())
    if err<best:
        best=err;out=dict(z);out['onebody_angle_words']=theta.cpu().numpy().astype(np.int64)%(1<<a.brot);out['onebody_brot']=a.brot;out['error_1_optimized']=err;np.savez(a.out,**out)
    if err<=a.target:
        print(json.dumps(dict(step=step,error_1=err,best=best,target=a.target,seconds=time.time()-began)),flush=True);break
    candidate=u[:,None,None,:]+np.sin(alpha)*torch.stack([J,-J],dim=2)+(1-np.cos(alpha))*Q[:,:,None,:]
    old=(u*(u@R)).sum(1);quad=(candidate*(candidate@R)).sum(3);dot=(candidate*u[:,None,None,:]).sum(3)
    changes=2*e[:,None,None]*(quad-old[:,None,None])+2*e[:,None,None].square()*(1-dot.square())
    values,ix=changes.reshape(76,-1).min(1);order=values.argsort();order=order[values[order]<-1e-14]
    moved=False
    for count in sorted(set([len(order),32,16,8,4,2,1]),reverse=True):
        if not count or count>len(order):continue
        rows=order[:count];cols=ix[rows]//2;signs=ix[rows]%2
        new=u.clone();new[rows]=candidate[rows,cols,signs];newloss=((new.T*e)@new-H).square().sum()
        if float(newloss)<float(loss)-1e-14:
            theta[rows,cols]+=torch.where(signs==0,1.,-1.);accepted+=count;moved=True;break
    if step%10==0 or not moved:print(json.dumps(dict(step=step,error_1=err,best=best,accepted=accepted,moved=moved,seconds=time.time()-began)),flush=True)
    if not moved:break
