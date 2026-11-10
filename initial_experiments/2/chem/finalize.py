"""Quantize a fitted factorization to 12-bit Givens words, refit w and W, and write submission.npz.

Stages: encode two-body vectors on the tree; refit w with the decoded vectors fixed (error_2);
optimize W and w for lambda_eff under error_2 and E_gap constraints; encode the one-body
eigenvectors of h1'; evaluate with the host metric code.
"""
import argparse,json,numpy as np,torch
p=argparse.ArgumentParser()
p.add_argument('--ref',required=True);p.add_argument('--fit',required=True);p.add_argument('--base',required=True)
p.add_argument('--out',required=True);p.add_argument('--steps',type=int,default=1500)
p.add_argument('--lamsteps',type=int,default=1500);p.add_argument('--err2cap',type=float,default=0.0192)
p.add_argument('--gapfloor',type=float,default=3.9160);p.add_argument('--lr',type=float,default=1e-3)
a=p.parse_args()
torch.backends.cuda.matmul.allow_tf32=False;dev='cuda'
ref=np.load(a.ref);base=dict(np.load(a.base));tree=base['tree_edges']
fit=np.load(a.fit);U=fit['U'];w=fit['w'];W=fit['W'];R,B,C=w.shape;N=76;BROT=12

def encode(v):
    """Words for unit vectors v (n, N) on `tree` (inverse of reconstruct_u_vectors)."""
    v=v.copy();th=np.zeros((len(v),len(tree)))
    for t in reversed(range(len(tree))):
        pa,ch=tree[t];th[:,t]=np.arctan2(v[:,ch],v[:,pa]);v[:,pa]=np.hypot(v[:,pa],v[:,ch]);v[:,ch]=0
    return np.mod(np.rint(th/np.pi*2**(BROT-1)).astype(np.int64),1<<BROT)
def decode(words):
    ang=words/2**(BROT-1);v=np.zeros((len(words),N));v[:,0]=1
    for t,(pa,ch) in enumerate(tree):
        th=np.pi*ang[:,t];c,s=np.cos(th),np.sin(th);vp,vc=v[:,pa].copy(),v[:,ch].copy()
        v[:,pa]=c*vp-s*vc;v[:,ch]=s*vp+c*vc
    return v
w2=encode(U.reshape(-1,N));Uq=decode(w2).reshape(R,B,N)
print('quantization max |dU|',np.abs(np.abs((Uq*U).sum(-1))-1).max(),flush=True)

G=torch.tensor(ref['G'],device=dev)
uq=torch.tensor(Uq,device=dev,dtype=torch.float32)
P=torch.einsum('rbp,rbq->rbpq',uq,uq).reshape(R,B,N*N)
def err2(wt):
    X=torch.einsum('rbk,rbc->rck',P,wt).reshape(-1,N*N);return (X.T@X-G).norm()
wt=torch.tensor(w,device=dev,dtype=torch.float32,requires_grad=True)
opt=torch.optim.Adam([wt],lr=3e-4);best=(1e9,None)
for s in range(a.steps):
    opt.zero_grad();e=err2(wt);(e*e).backward();opt.step()
    if e.item()<best[0]:best=(e.item(),wt.detach().clone())
    if s%300==0:print('refit',s,e.item(),flush=True)
print('refit err2',best[0],flush=True)
wt=best[1].clone().requires_grad_(True)

# lambda_eff stage (float64 for metric pieces).
d64=dict(dtype=torch.float64,device=dev)
h1=torch.tensor(ref['h1_thc'],**d64);S=torch.tensor(ref['S'],**d64);beta1=float(ref['beta_1'])
eta=int(ref['electrons']);ecore=float(ref['ecore']);EUP=-22140.27335
P64=P.to(torch.float64)
Wt=torch.tensor(W,**d64).requires_grad_(True)
def metrics(wt,Wt):
    w64=wt.to(torch.float64)
    X=torch.einsum('rbk,rbc->rck',P64,w64).reshape(R,C,N,N)
    h1p=h1-0.5*eta*S+beta1*torch.eye(N,**d64)-torch.einsum('rc,rcpq->pq',Wt,X)
    h1p=0.5*(h1p+h1p.T)
    ev=torch.linalg.eigvalsh(h1p.cpu()).to(dev)
    lh1=ev.abs().sum();wB=Wt+w64.sum(1)
    lh2=0.25*((w64.abs().sum(1)+wB.abs())**2).sum()
    Lam=lh1+lh2
    esos=ecore-eta*beta1-(lh1-ev.sum())-0.5*(Wt**2).sum()
    gap=EUP-esos
    return Lam,gap,torch.sqrt(gap*(2*Lam-gap)),h1p
wt=wt.detach()
L0,g0,le0,_=metrics(wt,Wt);print('start Lambda %.4f gap %.4f lam_eff %.4f'%(L0.item(),g0.item(),le0.item()),flush=True)
# error_2 does not depend on W: optimize W alone.
opt=torch.optim.Adam([Wt],lr=a.lr);best=None
for s_ in range(a.lamsteps):
    opt.zero_grad()
    Lam,gap,le,_=metrics(wt,Wt)
    loss=le+1e3*torch.relu(a.gapfloor-gap)**2
    if not torch.isfinite(loss):print('nonfinite',s_);break
    loss.backward();opt.step()
    if gap.item()>=a.gapfloor-2e-4 and (best is None or le.item()<best[0]):best=(le.item(),Wt.detach().clone())
    if s_%500==0:print('lam',s_,'lam_eff %.4f Lambda %.3f gap %.4f'%(le.item(),Lam.item(),gap.item()),flush=True)
if best is not None:Wt=best[1]
Lam,gap,le,h1p=metrics(wt,Wt);print('chosen lam_eff',le.item(),'err2',err2(wt).item(),flush=True)
wf=wt.cpu().numpy().astype(np.float64);Wf=Wt.detach().cpu().numpy()
ev,vec=np.linalg.eigh(h1p.detach().cpu().numpy())
w1=encode(vec.T)
words=np.concatenate([w1,w2]).astype(np.int64)
out=dict(base);out.update(U=Uq,W=wf,H2_I=Wf+wf.sum(1),angle_words=words)
np.savez(a.out,**out)
try:
    from dfthcsos.verify import load_inputs,evaluate_words
    inp=load_inputs(a.out);m=evaluate_words(inp,words,tree,BROT)
    print('FINAL',json.dumps(m),flush=True)
except ImportError as exc:
    print('evaluate separately:',exc)
