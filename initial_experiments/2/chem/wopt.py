"""Re-optimize the identity offsets W (error_2 and the two-body angle words do not depend on W).

h1' changes with W, so the one-body rows are re-encoded from its new eigenvectors.
"""
import sys,numpy as np,torch
src,ref,dst=sys.argv[1],sys.argv[2],sys.argv[3];steps=int(sys.argv[4]);floor=float(sys.argv[5]) if len(sys.argv)>5 else 3.9165
z=dict(np.load(src));r=np.load(ref);tree=z['tree_edges'];N=76
w=z['W'];Wu=z['H2_I']-w.sum(1);R,B,C=w.shape
d=dict(dtype=torch.float64,device='cuda')
words=z['angle_words'];ang=words[N:]/2048.
v=np.zeros((len(ang),N));v[:,0]=1
for t,(p,c) in enumerate(tree):
    th=np.pi*ang[:,t];cs,sn=np.cos(th),np.sin(th);vp,vc=v[:,p].copy(),v[:,c].copy();v[:,p]=cs*vp-sn*vc;v[:,c]=sn*vp+cs*vc
U=torch.tensor(v.reshape(R,B,N),**d);wt=torch.tensor(w,**d)
X=torch.einsum('rbp,rbq,rbc->rcpq',U,U,wt)
h1=torch.tensor(r['h1_thc'],**d);S=torch.tensor(r['S'],**d);b1=float(r['beta_1']);eta=int(r['electrons']);ec=float(r['ecore'])
A0=h1-0.5*eta*S+b1*torch.eye(N,**d)
Wt=torch.tensor(Wu,**d,requires_grad=True)
def met(Wt):
    h=A0-torch.einsum('rc,rcpq->pq',Wt,X);h=0.5*(h+h.T)
    ev=torch.linalg.eigvalsh(h.cpu()).to('cuda');lh1=ev.abs().sum()
    lh2=0.25*((wt.abs().sum(1)+(Wt+wt.sum(1)).abs())**2).sum();L=lh1+lh2
    gap=-22140.27335-(ec-eta*b1-(lh1-ev.sum())-0.5*(Wt**2).sum())
    return L,gap,torch.sqrt(gap*(2*L-gap)),h
L,g,le,_=met(Wt);print('start',le.item(),g.item(),L.item(),flush=True)
best=(le.item(),Wt.detach().clone())
for lr in [3e-4,1e-4,3e-5]:
    opt=torch.optim.Adam([Wt],lr=lr)
    for s in range(steps):
        opt.zero_grad();L,g,le,_=met(Wt);(le+1e4*torch.relu(floor-g)**2).backward();opt.step()
        if g.item()>=floor and le.item()<best[0]:best=(le.item(),Wt.detach().clone())
    Wt=best[1].clone().requires_grad_(True);print(lr,'best',best[0],flush=True)
L,g,le,h=met(best[1]);print('final',le.item(),g.item(),L.item())
Wn=best[1].cpu().numpy();ev,vec=np.linalg.eigh(h.detach().cpu().numpy())
vv=vec.T.copy();th=np.zeros((N,len(tree)))
for t in reversed(range(len(tree))):
    p,c=tree[t];th[:,t]=np.arctan2(vv[:,c],vv[:,p]);vv[:,p]=np.hypot(vv[:,p],vv[:,c]);vv[:,c]=0
w1=np.mod(np.rint(th/np.pi*2048).astype(np.int64),4096)
z['angle_words']=np.concatenate([w1,words[N:]]);z['H2_I']=Wn+w.sum(1)
np.savez(dst,**z)
