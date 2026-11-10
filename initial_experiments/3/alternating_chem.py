"""Alternate exact tensor-factor projections and local dictionary fitting.

The spectral projection preserves the best rank-K reference tensor. For each
local dictionary, eliminate its linear weights by least squares before taking
an optimizer step on the unit vectors. This avoids optimizing poorly scaled
weights and vectors together.
"""
import argparse,time,json
from pathlib import Path
import numpy as np
import torch
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
p=argparse.ArgumentParser();p.add_argument('--R',type=int,default=80);p.add_argument('--B',type=int,default=25)
p.add_argument('--C',type=int,default=4);p.add_argument('--rounds',type=int,default=200)
p.add_argument('--inner',type=int,default=100);p.add_argument('--init');p.add_argument('--out',required=True)
p.add_argument('--lr',type=float,default=.003)
p.add_argument('--weight',type=float,default=0.,help='ridge in the G-weighted feature metric; zero uses Frobenius')
a=p.parse_args();torch.set_num_threads(3);HERE=Path(__file__).resolve().parent
dev='cuda';dt=torch.float64;d=dict(device=dev,dtype=dt);n=76;R,B,C=a.R,a.B,a.C;K=R*C
ref=np.load(HERE/'parameters/reference_spectrum.npz');ev=ref['eigenvalues'];V=ref['eigenvectors']
ii,jj=ref['ii'],ref['jj'];fac=ref['fac'];D=len(ii)
F=torch.tensor((V[:,:K]*np.sqrt(np.maximum(ev[:K],0))).T,**d)
G=torch.tensor((V*ev)@V.T,**d)
ti=torch.tensor(ii,device=dev);tj=torch.tensor(jj,device=dev);ff=torch.tensor(fac,**d)
def asmat(Y):
    m=torch.zeros((*Y.shape[:-1],n,n),**d)
    m[...,ti,tj]=Y/ff;m[...,tj,ti]=Y/ff
    return m
target=asmat(F.reshape(R,C,D))
e,q=torch.linalg.eigh((target@target).sum(1));U=q[:,:,-B:].transpose(1,2).contiguous()
if a.init:
    z=np.load(a.init);old=z['U'];cnt=min(R,len(old));rng=np.random.default_rng(19)
    for r in range(cnt):
        score=np.linalg.norm(z['w'][r],axis=1);ix=np.argsort(score)[-min(B,len(score)):]
        U[r,:len(ix)]=torch.tensor(old[r,ix],**d)
        if len(ix)<B:U[r,len(ix):]=torch.tensor(rng.standard_normal((B-len(ix),n)),**d)
U=U.clone().requires_grad_(True)
eye=torch.eye(B,**d)
def project_weights(y):
    un=U/U.norm(dim=2,keepdim=True)
    gram=(un@un.transpose(1,2)).square()
    if a.weight:
        P=un[:,:,ti]*un[:,:,tj]*ff
        Pf=P.reshape(R*B,D)@F.T
        Pf=Pf.reshape(R,B,K)
        gram=Pf@Pf.transpose(1,2)+a.weight*gram
    gram=gram+1e-10*eye
    cross=torch.einsum('rbp,rcpq,rbq->rbc',un,y,un)
    w=torch.linalg.solve(gram,cross)
    return un,w,gram,cross
def xvals(un,w):
    P=un[:,:,ti]*un[:,:,tj]*ff
    return torch.bmm(w.transpose(1,2),P).reshape(K,D)
opt=torch.optim.Adam([U],lr=a.lr);best=1e9;start=time.time()
# A passed factorization supplies the one-body correction to project into
# the new auxiliary span after tensor fitting.
base=np.load(HERE/'candidates/chem35c45_rot.npz')
bu=reconstruct_u_vectors(base['angle_words'][76:]/2048,base['tree_edges']).reshape(35,60,n)
bp=bu[:,:,ii]*bu[:,:,jj]*fac
bx=np.einsum('rbc,rbp->rcp',base['W'],bp).reshape(-1,D)
bo=(base['H2_I']-base['W'].sum(1)).reshape(-1)
offset_target=torch.tensor(bo@bx,**d)
if a.init:
    z=np.load(a.init);initialw=torch.zeros((R,B,C),**d)
    for r in range(min(R,len(z['U']))):
        score=np.linalg.norm(z['w'][r],axis=1);ix=np.argsort(score)[-min(B,len(score)):]
        initialw[r,:len(ix),:min(C,z['w'].shape[2])]=torch.tensor(z['w'][r,ix,:C],**d)
    un=U/U.norm(dim=2,keepdim=True);X=xvals(un,initialw).detach()
else:
    with torch.no_grad():un,w,_,_=project_weights(target);X=xvals(un,w)
for round in range(a.rounds):
    with torch.no_grad():
        mat=X@F.T
        if a.weight:mat=mat*torch.tensor(np.maximum(ev[:K],0)+a.weight,**d)
        left,_,right=torch.linalg.svd(mat,full_matrices=False)
        Y=(left@right@F).reshape(R,C,D);target=asmat(Y)
        if a.weight:
            Yflat=Y.reshape(K,D)
            Yw=(Yflat@F.T)@F+a.weight*Yflat
            norm=(Yflat*Yw).sum();target=asmat(Yw.reshape(R,C,D))
        else:norm=Y.square().sum()
    for inner in range(a.inner):
        opt.zero_grad();un,w,gram,cross=project_weights(target)
        # The envelope theorem allows fixed optimal weights in the local
        # derivative; the solve is still differentiable for this formulation.
        loss=norm-(w*cross).sum()
        loss.backward();opt.step()
    with torch.no_grad():
        un,w,_,_=project_weights(target);X=xvals(un,w)
        error=float((X.T@X-G).norm())
        if error<best:
            best=error
            gram=X@X.T;v=X@offset_target
            # Small Tikhonov term controls almost-null modes.
            W=torch.linalg.solve(gram+1e-8*torch.eye(K,**d),v).reshape(R,C)
            np.savez(a.out,U=un.cpu().numpy(),w=w.cpu().numpy(),W=W.cpu().numpy(),err2=error)
    if round%5==0:
        print(json.dumps(dict(round=round,error_2=error,best=best,local_loss=float(loss),seconds=time.time()-start)),flush=True)
    opt.param_groups[0]['lr']=a.lr*(.1+.9*(1-round/a.rounds))
