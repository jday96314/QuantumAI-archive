"""Joint angle/weight continuation with an explicit LCU one-norm penalty."""
import argparse,json,time
from pathlib import Path
import numpy as np
import torch
p=argparse.ArgumentParser();p.add_argument('--init',required=True);p.add_argument('--out',required=True)
p.add_argument('--brot',type=int,default=8);p.add_argument('--steps',type=int,default=16000)
p.add_argument('--lr',type=float,default=.02);p.add_argument('--weight-lr',type=float,default=.0002)
p.add_argument('--gamma',type=float,default=1e-5);p.add_argument('--strength',type=float,default=30.)
p.add_argument('--start-strength',type=float,default=1e-5);p.add_argument('--whitened',action='store_true')
p.add_argument('--project-every',type=int,default=0);p.add_argument('--project-cycles',type=int,default=3)
p.add_argument('--positive-weights',action='store_true')
p.add_argument('--hard-forward',action='store_true')
a=p.parse_args();torch.set_num_threads(3);root=Path(__file__).resolve().parent;z=np.load(a.init)
R,bases,n=z['U'].shape;u0=z['U'].reshape(-1,n);w0=z['w'];B=len(u0);tree=z.get('tree_edges',np.load(root/'candidates/chem35c45_rot.npz')['tree_edges'])
if a.positive_weights:
    assert not a.whitened
    w0=np.maximum(w0,0)
v=u0.copy();th=np.zeros((B,75))
for j in reversed(range(75)):
    pa,ch=tree[j];th[:,j]=np.arctan2(v[:,ch],v[:,pa]);v[:,pa]=np.hypot(v[:,pa],v[:,ch]);v[:,ch]=0
alpha=np.pi/(1<<(a.brot-1));d=dict(device='cuda',dtype=torch.float32)
theta=torch.tensor(th/alpha,requires_grad=True,**d);w=torch.tensor(w0,requires_grad=True,**d)
sp=np.load(root/'parameters/reference_spectrum.npz');ev,V=sp['eigenvalues'],sp['eigenvectors'];G=torch.tensor((V*np.maximum(ev,0))@V.T,**d)
ii=torch.tensor(sp['ii'],device='cuda');jj=torch.tensor(sp['jj'],device='cuda');fac=torch.tensor(sp['fac'],**d)
traces=torch.tensor(w0.sum(1),**d)
iv=(sp['ii']==sp['jj']).astype(float);vt=V.T@iv;trace_norm=float(np.maximum(ev,0)@(vt*vt))
traces*=np.sqrt(trace_norm/float(traces.square().sum()))
def weights(t=None,u=None):
    ww=w
    if a.whitened:
        if u is None:u=decode(theta if t is None else t)
        P=(u[:,ii]*u[:,jj]*fac).reshape(R,bases,-1).double()
        L=torch.linalg.cholesky(P@P.transpose(1,2)+1e-10*torch.eye(bases,device='cuda'))
        ww=torch.linalg.solve_triangular(L.transpose(1,2),w.double(),upper=True).float()
    return ww*torch.sqrt(traces.square().sum()/ww.sum(1).square().sum().clamp(min=1e-15))
start=torch.zeros((B,n),**d);start[:,0]=1;stages=[]
for j,(pa,ch) in enumerate(tree):
    pa,ch=int(pa),int(ch)
    if not stages or any(pa in (x,y) or ch in (x,y) for _,x,y in stages[-1]):stages.append([])
    stages[-1].append((j,pa,ch))
stages=[tuple(torch.tensor(x,device='cuda') for x in zip(*s)) for s in stages]
def decode(t):
    if a.hard_forward:t=t+(t.round()-t).detach()
    u=start
    for ts,ps,cs in stages:
        ang=t[:,ts]*alpha;c,s=ang.cos(),ang.sin();up,uc=u[:,ps],u[:,cs]
        new=u.clone();new[:,ps]=c*up-s*uc;new[:,cs]=s*up+c*uc;u=new
    return u
if a.whitened:
    with torch.no_grad():
        pu=(decode(theta)[:,ii]*decode(theta)[:,jj]*fac).reshape(R,bases,-1).double()
        L=torch.linalg.cholesky(pu@pu.transpose(1,2)+1e-10*torch.eye(bases,device='cuda'))
        w.copy_((L.transpose(1,2)@w.double()).float())
def physical(t):
    u=decode(t);P=(u[:,ii]*u[:,jj]*fac).reshape(R,bases,-1);X=torch.bmm(weights(t,u).transpose(1,2),P).reshape(-1,len(sp['ii']));return (X.T@X-G).square().sum(),u
with torch.no_grad():
    u=decode(theta);ur=u.reshape(R,bases,n);S=torch.bmm(ur,ur.transpose(1,2)).square();ww=weights();K=torch.bmm(ww,ww.transpose(1,2));yn=(torch.bmm(K,S)*K).sum(2).reshape(-1).clamp(min=1e-14)
    v=start.clone();J=torch.zeros((B,75,n),**d)
    for j,(pa,ch) in enumerate(tree):
        pa,ch=int(pa),int(ch);ang=theta[:,j]*alpha;c,s=ang.cos(),ang.sin();jp,jc=J[:,:,pa].clone(),J[:,:,ch].clone()
        J[:,:,pa]=c[:,None]*jp-s[:,None]*jc;J[:,:,ch]=s[:,None]*jp+c[:,None]*jc
        up,uc=v[:,pa].clone(),v[:,ch].clone();v[:,pa]=c*up-s*uc;v[:,ch]=s*up+c*uc
        J[:,j,pa]=-alpha*v[:,ch];J[:,j,ch]=alpha*v[:,pa]
    weight=(8*J.square().sum(2)*yn[:,None]).clamp(min=1e-14)
opt=torch.optim.Adam([{'params':[theta],'lr':a.lr},{'params':[w],'lr':a.weight_lr}]);best=float('inf');started=time.time()
@torch.no_grad()
def project_weights():
    assert a.whitened
    u=decode(theta);P=(u[:,ii]*u[:,jj]*fac).reshape(R,bases,-1).double()
    L=torch.linalg.cholesky(P@P.transpose(1,2)+1e-10*torch.eye(bases,device='cuda'))
    Q=torch.linalg.solve_triangular(L,P,upper=False);flat=Q.reshape(R*bases,-1)
    ov=(flat@flat.T).reshape(R,bases,R,bases).permute(0,2,1,3).contiguous()
    cross=(flat@G.double()@flat.T).reshape(R,bases,R,bases)
    target=torch.stack([cross[r,:,r,:] for r in range(R)])
    aa=L.transpose(1,2)@weights(theta,u).double();H=aa@aa.transpose(1,2)
    residual=target-torch.einsum('rsab,sbc,rsdc->rad',ov,H,ov)
    C=w.shape[2]
    for _ in range(a.project_cycles):
        for r in range(R):
            hh=.5*(residual[r]+residual[r].T)+H[r];ee,vv=torch.linalg.eigh(hh)
            ff=vv[:,-C:]*ee[-C:].clamp(min=0).sqrt();new=ff@ff.T;change=new-H[r];H[r]=new
            A=ov[:,r];residual-=A@change@A.transpose(1,2)
    ee,vv=torch.linalg.eigh(H);ff=vv[:,:,-C:]*ee[:,-C:].clamp(min=0).sqrt()[:,None,:]
    left,_,right=torch.linalg.svd(ff.transpose(1,2)@w.double(),full_matrices=False)
    w.copy_((ff@(left@right)).float())
    if w in opt.state:
        opt.state[w]['exp_avg'].zero_();opt.state[w]['exp_avg_sq'].zero_()
for step in range(a.steps):
    opt.zero_grad();error2,u=physical(theta);ww=weights();l1=.25*torch.sqrt(ww.square()+1e-12).sum(1).square().sum()
    strength=10**(np.log10(a.start_strength)+(np.log10(a.strength)-np.log10(a.start_strength))*step/max(1,a.steps-1))
    penalty=(weight*(1-torch.cos(2*np.pi*theta))/(4*np.pi*np.pi)).sum() if not a.hard_forward else torch.zeros((),**d)
    loss=error2+a.gamma*l1+strength*penalty;loss.backward();opt.step()
    if a.positive_weights:
        with torch.no_grad():w.clamp_(min=0)
    if a.project_every and step%a.project_every==0:project_weights()
    decay=.05+.95*.5*(1+np.cos(np.pi*step/a.steps));opt.param_groups[0]['lr']=a.lr*decay;opt.param_groups[1]['lr']=a.weight_lr*decay
    if step%100==0 or step==a.steps-1:
        with torch.no_grad():
            words=theta.round();hard,uq=physical(words);value=float(hard.sqrt())
            out=dict(U=uq.cpu().numpy().astype(float).reshape(R,bases,n),w=weights(words,uq).cpu().numpy().astype(float),W=np.zeros((R,w.shape[2])),tree_edges=tree,brot=a.brot,angle_words=words.cpu().numpy().astype(np.int64)%(1<<a.brot),err2=value)
            if value<best:best=value;np.savez(a.out,**out)
            soft_error,soft_u=physical(theta)
            soft=dict(out,U=soft_u.cpu().numpy().astype(float).reshape(R,bases,n),w=weights(theta,soft_u).cpu().numpy().astype(float),err2=float(soft_error.sqrt()));soft.pop('angle_words');soft.pop('brot');np.savez(str(Path(a.out).with_suffix(''))+'_soft.npz',**soft)
        print(json.dumps(dict(step=step,soft_error=float(error2.detach().sqrt()),hard_error=value,best=best,lambda_h2_lower=float(l1.detach()),strength=strength,seconds=time.time()-started)),flush=True)
