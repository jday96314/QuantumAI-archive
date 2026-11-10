"""Exact-word, multi-vector local repair on any legal shared tree.

Proposal construction adapts the completed chem-qfit/discrete.py geometry.
Native target, cached compact refits, and matched fixed-word controls are local.
"""
import argparse
from core import *

class Geometry:
    def __init__(self,tree):
        self.tree=tree;self.i=torch.tensor(ii,device=DEVICE);self.j=torch.tensor(jj,device=DEVICE);self.f=tensor(fac)
    def decode(self,t):
        u=torch.zeros((*t.shape[:-1],N),**D);u[...,0]=1
        for j,(p,c) in enumerate(self.tree):
            co=(t[...,j]*np.pi/128).cos();si=(t[...,j]*np.pi/128).sin()
            a=u[...,p].clone();b=u[...,c].clone();u[...,p]=co*a-si*b;u[...,c]=si*a+co*b
        return u
    def feats(self,u):return u[...,self.i]*u[...,self.j]*self.f
    def derivatives(self,t):
        tf=t.reshape(-1,75);u=torch.zeros((len(tf),N),**D);u[:,0]=1
        J=torch.zeros((len(tf),75,N),**D)
        for j,(p,c) in enumerate(self.tree):
            co=(tf[:,j]*np.pi/128).cos();si=(tf[:,j]*np.pi/128).sin()
            a=J[:,:,p].clone();b=J[:,:,c].clone();J[:,:,p]=co[:,None]*a-si[:,None]*b;J[:,:,c]=si[:,None]*a+co[:,None]*b
            a=u[:,p].clone();b=u[:,c].clone();u[:,p]=co*a-si*b;u[:,c]=si*a+co*b
            J[:,j,p]=-np.pi/128*u[:,c];J[:,j,c]=np.pi/128*u[:,p]
        return u.reshape(*t.shape[:-1],N),J.reshape(*t.shape[:-1],75,N)
    @torch.no_grad()
    def proposals(self,t,w,H,step):
        u,J=self.derivatives(t);R,B,C=w.shape;P=self.feats(u);F=(w.transpose(1,2)@P).reshape(-1,len(ii))
        E=F.T@F-H;Y=w@(w.transpose(1,2)@P);YD=Y@E
        M=torch.zeros((R,B,N,N),**D);M[:,:,self.i,self.j]=YD/self.f
        M=M+M.transpose(-1,-2);M.diagonal(dim1=-1,dim2=-2).mul_(.5)
        g=8*(J*(M@u[...,None]).squeeze(-1)[:,:,None,:]).sum(-1)
        raw=2*J.square().sum(-1);h=raw.clone();Qs=[]
        for r in range(R):
            Q,Rm=torch.linalg.qr(P[r].T);Qs.append(Q)
            cross=2*(J[r]@u[r].T)*(u[r]@u[r].T)[:,None,:]
            cr=cross.reshape(B*75,B)
            coeff=torch.linalg.solve_triangular(Rm.T,cr.T,upper=False).T
            perp=(raw[r]-coeff.square().sum(-1).reshape(B,75)).clamp(min=0)
            h[r]=4*(perp+.001*raw[r])*Y[r].square().sum(-1)[:,None]
        score=.5*h-g.abs();score[h<1e-18]=torch.inf;order=score.argsort(-1)
        cp=[]
        for k in [1,2,4,8]:
            inds=order[:,:,:k];delta=torch.zeros_like(t);delta.scatter_(-1,inds,-g.gather(-1,inds).sign())
            if step%9==8:delta*=2
            cp.append(t+delta)
        cp=torch.stack(cp,2);dp=self.feats(self.decode(cp))-P[:,:,None,:]
        score=4*(dp*(Y@E)[:,:,None,:]).sum(-1)
        for r in range(R):
            perp=dp[r]-((dp[r]@Qs[r])@Qs[r].T)
            score[r]+=2*Y[r].square().sum(-1)[:,None]*perp.square().sum(-1)
        score[dp.square().sum(-1)<1e-22]=torch.inf
        best,ki=score.min(-1);flatcp=cp.reshape(R*B,4,75)
        choices=flatcp[torch.arange(R*B,device=DEVICE),ki.reshape(-1)]
        return choices,best.reshape(-1).argsort()

def main():
    p=argparse.ArgumentParser();p.add_argument('--init',required=True);p.add_argument('--tag',required=True)
    p.add_argument('--steps',type=int,default=100);p.add_argument('--cycles',type=int,default=6)
    p.add_argument('--seconds',type=float,default=600);p.add_argument('--refit-budget',type=int,default=0);a=p.parse_args()
    z,u,w,tree,words=load(a.init);R,B,C=w.shape
    if words is None:words=quantize(u,tree)[0]
    words=words.reshape(R,B,75);geom=Geometry(tree);H=target();HT=tensor(H)
    w,err,_=Refit(u,H).fit(w,100);initial=err;history=[];accepted=0;refits=100*R;start=time.time()
    out=HERE/f'candidates/{a.tag}.npz';save(out,u,w,tree,words,err2=err)
    stale=0
    for step in range(a.steps):
        w,baseline,_=Refit(u,H).fit(w,a.cycles);refits+=a.cycles*R;err=baseline
        choices,order=geom.proposals(tensor(words),tensor(w),HT,step);changed=False;attempts=0
        for count in [128,32,8,2,1]:
            tr=tensor(words).reshape(-1,75).clone();ids=order[:count];tr[ids]=choices[ids]
            nw=tr.cpu().numpy().astype(np.int64).reshape(R,B,75)%256
            nu=decode(nw*np.pi/128,tree);ww,ne,_=Refit(nu,H).fit(w,a.cycles);refits+=a.cycles*R;attempts+=1
            if ne<baseline-1e-10:
                u,w,words,err=nu,ww,nw,ne;changed=True;accepted+=len(ids);break
            if a.refit_budget and refits>=a.refit_budget:break
        save(out,u,w,tree,words,err2=err)
        stale=0 if changed else stale+1
        rec=dict(step=step,error_2=err,fixed_word_baseline=baseline,changed=changed,accepted_rows=accepted,
                 attempts=attempts,block_eigensolves=refits,seconds=time.time()-start)
        history.append(rec);report(HERE/f'results/{a.tag}.json',history,False)
        if step%10==0 or not changed:print(json.dumps(rec),flush=True)
        if step%25==0:save(HERE/f'trajectories/{a.tag}_{step}.npz',u,w,tree,words,err2=err)
        if stale>=4 or time.time()-start>a.seconds or (a.refit_budget and refits>=a.refit_budget):break
    w,err,hist=Refit(u,H).fit(w,300);words=signed_words(words,tree)
    save(out,decode(words*np.pi/128,tree),w,tree,words,err2=err)
    report(HERE/f'results/{a.tag}_final.json',dict(initial=initial,final=err,steps=len(history),
           block_eigensolves=refits,polish_block_eigensolves=300*R,accepted_rows=accepted,polish_history=hist,seconds=time.time()-start,
           stop_reason='refit_budget' if a.refit_budget and refits>=a.refit_budget else 'stale_or_step_or_time_limit'))
if __name__=='__main__':main()
