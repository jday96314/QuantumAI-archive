"""Explore nonnegative shared factors, whose LCU norm avoids cancellations."""
import argparse,json,time
from pathlib import Path
import numpy as np
import torch
p=argparse.ArgumentParser();p.add_argument('--init',required=True);p.add_argument('--out',required=True)
p.add_argument('--steps',type=int,default=25000);p.add_argument('--lr',type=float,default=.001)
p.add_argument('--optimizer',choices=['adam','lbfgs'],default='adam')
p.add_argument('--fixed-traces',action='store_true')
p.add_argument('--fixed-vectors',action='store_true')
a=p.parse_args();torch.set_num_threads(3);root=Path(__file__).resolve().parent;z=np.load(a.init)
sp=np.load(root/'parameters/reference_spectrum.npz');ev,V=sp['eigenvalues'],sp['eigenvectors']
d=dict(device='cuda',dtype=torch.float32);G=torch.tensor((V*np.maximum(ev,0))@V.T,**d)
ii=torch.tensor(sp['ii'],device='cuda');jj=torch.tensor(sp['jj'],device='cuda');fac=torch.tensor(sp['fac'],**d)
u=torch.tensor(z['U'][0],requires_grad=not a.fixed_vectors,**d);w=torch.tensor(abs(z['w'][0]),requires_grad=True,**d)
traces=torch.tensor(abs(z['w'][0]).sum(0),**d)
if a.fixed_traces:
    iv=(sp['ii']==sp['jj']).astype(float);vt=V.T@iv;target=float(np.maximum(ev,0)@(vt*vt))
    traces*=np.sqrt(target/float(traces.square().sum()))
def vals():
    ww=w.square() if a.optimizer=='lbfgs' else w
    if a.fixed_traces:ww=ww*torch.sqrt(traces.square().sum()/ww.sum(0).square().sum().clamp(min=1e-15))
    un=u/u.norm(dim=1,keepdim=True);P=un[:,ii]*un[:,jj]*fac;X=ww.T@P;pred=X.T@X
    return (pred-G).square().sum(),un,pred,ww
with torch.no_grad():
    if a.optimizer=='lbfgs':w.copy_(w.clamp(min=1e-9).sqrt())
    _,_,pred,_=vals();scale=float(((pred*G).sum()/pred.square().sum()).clamp(min=1e-12).pow(.25))
    if not a.fixed_traces:w.mul_(np.sqrt(scale) if a.optimizer=='lbfgs' else scale)
opt=(torch.optim.LBFGS(([w] if a.fixed_vectors else [u,w]),lr=a.lr,max_iter=12,history_size=30,line_search_fn='strong_wolfe',tolerance_grad=1e-10,tolerance_change=1e-13) if a.optimizer=='lbfgs' else torch.optim.Adam([{'params':[u],'lr':0 if a.fixed_vectors else a.lr},{'params':[w],'lr':a.lr*.3}]))
best=float('inf');start=time.time();current={}
def closure():
    opt.zero_grad();loss,un,pred,ww=vals();loss.backward();current.update(un=un.detach(),ww=ww.detach(),error=float(loss.detach().sqrt()));return loss
for step in range(a.steps):
    if a.optimizer=='lbfgs':opt.step(closure)
    else:closure()
    un,ww,err=current['un'],current['ww'],current['error']
    if step%(1 if a.optimizer=='lbfgs' else 25)==0 and err<best:
        best=err;out=dict(U=un.cpu().numpy().astype(float)[None,:,:],w=ww.cpu().numpy().astype(float)[None,:,:],W=np.zeros((1,w.shape[1])),err2=err)
        if 'tree_edges' in z:out['tree_edges']=z['tree_edges']
        np.savez(a.out,**out)
    if step%(5 if a.optimizer=='lbfgs' else 100)==0:print(json.dumps(dict(step=step,error_2=err,best=best,lambda_h2_no_offset=float(.25*ww.sum(0).square().sum()),seconds=time.time()-start)),flush=True)
    if a.optimizer=='adam':
        opt.step()
        with torch.no_grad():w.clamp_(min=0)
        decay=.05+.95*.5*(1+np.cos(np.pi*step/a.steps));opt.param_groups[0]['lr']=0 if a.fixed_vectors else a.lr*decay;opt.param_groups[1]['lr']=a.lr*.3*decay
