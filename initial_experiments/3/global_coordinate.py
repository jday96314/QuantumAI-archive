"""Monotone integer-angle descent with variable projection and local curvature."""
import argparse,json,time
from pathlib import Path
import numpy as np
import torch
from compress_chem import encode
p=argparse.ArgumentParser();p.add_argument('--init',required=True);p.add_argument('--out',required=True)
p.add_argument('--brot',type=int,default=12);p.add_argument('--steps',type=int,default=500)
p.add_argument('--ridge',type=float,default=1e-5)
p.add_argument('--tree',help='optional numpy array of Givens tree edges')
p.add_argument('--columns',type=int,help='rank constraint during variable projection')
p.add_argument('--moves',type=int,default=16)
p.add_argument('--pair-candidates',type=int,default=0);p.add_argument('--pair-trials',type=int,default=12)
p.add_argument('--candidate-energy',action='store_true');p.add_argument('--curvature-padding',type=float,default=.03)
a=p.parse_args();torch.set_num_threads(3)
HERE=Path(__file__).resolve().parent;z=np.load(a.init);u0=z['U'][0];B=len(u0);n=76
tree=(np.load(a.tree) if a.tree else z['tree_edges'] if 'tree_edges' in z else
      np.load(HERE/'candidates/chem35c45_rot.npz')['tree_edges']);T=len(tree)
sp=np.load(HERE/'parameters/reference_spectrum.npz');ev,V=sp['eigenvalues'],sp['eigenvectors'];ii,jj=sp['ii'],sp['jj'];fac=sp['fac']
d=dict(device='cuda',dtype=torch.float64);G=torch.tensor((V*np.maximum(ev,0))@V.T,**d);const=G.square().sum()
theta=torch.tensor(encode(u0,tree,a.brot),**d,requires_grad=True)
ti=torch.tensor(ii,device='cuda');tj=torch.tensor(jj,device='cuda');ff=torch.tensor(fac,**d)
eye=torch.eye(B,**d);alpha=np.pi/(1<<(a.brot-1))
stages=[]
for t,(pa,ch) in enumerate(tree):
    pa,ch=int(pa),int(ch)
    if not stages or any(pa in (x,y) or ch in (x,y) for _,x,y in stages[-1]):stages.append([])
    stages[-1].append((t,pa,ch))
stages=[tuple(torch.tensor(x,device='cuda') for x in zip(*s)) for s in stages]
startvec=torch.zeros((B,n),**d);startvec[:,0]=1
def decode():
    angles=theta*alpha;u=startvec
    for ts,ps,cs in stages:
        ang=angles[:,ts];co,sn=ang.cos(),ang.sin();up,uc=u[:,ps],u[:,cs]
        new=u.clone();new[:,ps]=co*up-sn*uc;new[:,cs]=sn*up+co*uc;u=new
    return u
def state(grad=False):
    u=decode();P=u[:,ti]*u[:,tj]*ff;S=(u@u.T).square();H=(P@G)@P.T;Sr=S+a.ridge*eye
    with torch.no_grad():
        L=torch.linalg.cholesky(Sr)
        if a.columns and a.columns<B:
            hh=torch.linalg.solve_triangular(L,H,upper=False)
            hh=torch.linalg.solve_triangular(L,hh.T,upper=False).T
            ee,vv=torch.linalg.eigh(.5*(hh+hh.T))
            factor=torch.linalg.solve_triangular(L.T,vv[:,-a.columns:]*ee[-a.columns:].clamp(min=0).sqrt(),upper=True)
            K=factor@factor.T
        else:K=torch.cholesky_solve(torch.cholesky_solve(H,L).T,L).T
        K=.5*(K+K.T)
    KS=K@Sr;loss=const-2*(K*H).sum()+(KS*KS.T).sum()
    phys=K@S;error=(const-2*(K*H).sum()+(phys*phys.T).sum()).clamp(min=0).sqrt()
    return loss,error,u,K,S,L
def curvature(u,K,S,L,g):
    v=startvec.clone();J=torch.zeros((B,T,n),**d)
    for t,(pp,cc) in enumerate(tree):
        p,c=int(pp),int(cc);ang=theta[:,t]*alpha;co,sn=ang.cos(),ang.sin()
        jp,jc=J[:,:,p].clone(),J[:,:,c].clone()
        J[:,:,p]=co[:,None]*jp-sn[:,None]*jc;J[:,:,c]=sn[:,None]*jp+co[:,None]*jc
        vp,vc=v[:,p].clone(),v[:,c].clone();v[:,p]=co*vp-sn*vc;v[:,c]=sn*vp+co*vc
        J[:,t,p]=-alpha*v[:,c];J[:,t,c]=alpha*v[:,p]
    raw=2*J.square().sum(2)
    cross=2*(J@u.T)*(u@u.T)[:,None,:]
    flat=cross.reshape(B*T,B)
    projected=(flat*torch.cholesky_solve(flat.T,L).T).sum(1).reshape(B,T)
    effective=(raw-projected).clamp(min=0)+a.curvature_padding*raw
    ynorm=(K@S*K).sum(1)
    h=(4*effective*ynorm[:,None]).clamp(min=1e-18)
    pairs=[]
    if a.pair_candidates:
        barrier=.5*h-g.abs();ratio=barrier/h.clamp(min=1e-18)
        score=barrier.clone();score[h<1e-14]=float('inf');ratio[h<1e-14]=float('inf')
        candidates=[score.flatten().argsort()[:a.pair_candidates],ratio.flatten().argsort()[:a.pair_candidates//2]]
        if a.candidate_energy:candidates.append((g.square()/h.clamp(min=1e-18)).flatten().argsort(descending=True)[:a.pair_candidates//2])
        ids=torch.unique(torch.cat(candidates))
        rows=ids//T;deriv=J.reshape(B*T,76)[ids];uv=u[rows]
        ju=deriv@uv.T
        gram=2*(deriv@deriv.T)*(uv@uv.T)+2*ju*ju.T
        cr=flat[ids];gram-=cr@torch.cholesky_solve(cr.T,L)
        yy=K@S@K;hh=4*gram*yy[rows[:,None],rows[None,:]]
        signs=-g.flatten()[ids].sign();bb=.5*h.flatten()[ids]-g.flatten()[ids].abs()
        predicted=bb[:,None]+bb[None,:]+hh*signs[:,None]*signs[None,:]
        predicted[torch.tril(torch.ones_like(predicted,dtype=torch.bool))]=float('inf')
        choices=predicted.flatten().argsort()[:a.pair_trials]
        pairs=[[int(ids[int(ix)//len(ids)]),int(ids[int(ix)%len(ids)])] for ix in choices.cpu().tolist()]
        # Optimize a small binary quadratic model of simultaneous moves.
        # Its off-diagonal terms include the weight refit's tangent projection.
        interaction=(hh*signs[:,None]*signs[None,:]).cpu().numpy()
        np.fill_diagonal(interaction,0)
        field=bb.cpu().numpy().copy();selected=np.zeros(len(ids),dtype=bool)
        if float(predicted.flatten()[choices[0]])<0:
            first=int(choices[0]);seeds=[first//len(ids),first%len(ids)]
            for j in seeds:selected[j]=True;field+=interaction[:,j]
        for _ in range(max(512,2*a.moves)):
            gain=field*(1-2*selected.astype(int))
            if selected.sum()>=max(32,a.moves):gain[~selected]=np.inf
            j=int(gain.argmin())
            if gain[j]>=-1e-16:break
            direction=-1 if selected[j] else 1;selected[j]=not selected[j];field+=direction*interaction[:,j]
        together=ids[torch.tensor(np.flatnonzero(selected),device='cuda')].cpu().tolist()
        if together:pairs.insert(0,together)
    return h,pairs
best=float('inf');accepted=0;start=time.time();stale=0
for step in range(a.steps):
    theta.grad=None;loss,error,u,K,S,L=state();value=float(loss.detach());loss.backward();g=theta.grad.detach()
    with torch.no_grad():
        err=float(error)
        if err<best:
            best=err
            np.savez(a.out,U=u.cpu().numpy()[None,:,:],K=K.cpu().numpy(),err2=err,
                     angle_words=theta.cpu().numpy().astype(np.int64)%(1<<a.brot),tree_edges=tree,brot=a.brot)
        h,pairs=curvature(u.detach(),K,S.detach(),L,g)
        scores=-g.abs()+.5*h
        order=torch.argsort(scores.flatten())[:max(256,8*a.moves)]
        choices=[];seen=set()
        for x in order.cpu().tolist():
            row=x//T
            if row not in seen:choices.append(x);seen.add(row)
            if len(choices)>=a.moves:break
        changed=False
        sizes=sorted(set([a.moves>>k for k in range(a.moves.bit_length())]+[32,16,8,4,2,1]),reverse=True)
        coupled=([pairs[0][:k] for k in sorted(set([len(pairs[0])]+sizes),reverse=True) if k<=len(pairs[0])] if pairs else [])
        groups=coupled+[choices[:k] for k in sizes if len(choices)>=k]
        groups.extend([[int(x)] for x in order[1:33].cpu().tolist()])
        groups.extend(pairs[1:])
        for ids in groups:
            index=torch.tensor(ids,device='cuda');flat=theta.view(-1);direction=-g.flatten()[index].sign()
            flat[index]+=direction
            new=state()[0];nv=float(new)
            if nv<value-1e-13:
                changed=True;accepted+=len(ids);break
            flat[index]-=direction
        stale=0 if changed else stale+1
    if step%5==0:print(json.dumps(dict(step=step,error_2=err,best=best,accepted=accepted,changed=changed,
                                      lambda_gram_bound=float(.25*K.abs().sum()),seconds=time.time()-start)),flush=True)
    if stale>=3:break
z=dict(np.load(a.out));e,v=np.linalg.eigh(z.pop('K'));ix=np.argsort(e)[::-1];e=np.maximum(e[ix],0);v=v[:,ix]
C=a.columns or B
z.update(w=(v[:,:C]*np.sqrt(e[:C]))[None,:,:],W=np.zeros((1,C)));np.savez(a.out,**z)
