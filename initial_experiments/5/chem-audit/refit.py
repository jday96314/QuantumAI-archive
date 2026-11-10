"""Fixed integer dictionary weight/offset refit, round-4 Codex AL method.

Uses a QR-compressed *exact native* residual, gauge-preserving weight updates,
and explicit gap constraints. Only numerical candidates; native acceptance is
performed separately. Adapted from round-3 constrained_factor.py and qfit.
"""
import argparse,time,json
from audit_core import *

p=argparse.ArgumentParser()
p.add_argument('--fit',required=True);p.add_argument('--out',required=True)
p.add_argument('--steps',type=int,default=8000);p.add_argument('--lr',type=float,default=2e-6)
p.add_argument('--cap',type=float,default=.019990);p.add_argument('--gap',type=float,default=3.91515)
p.add_argument('--mue',type=float,default=.8);p.add_argument('--mug',type=float,default=.6)
p.add_argument('--offset-factor',type=float,default=3.)
p.add_argument('--project-gap',action='store_true')
a=p.parse_args();z,u,w0,tree,words=load(a.fit);m=Metric(tree)
ut=m.tensor(u);P=m.feats(ut);R,B,C=w0.shape
assert R==1,'compressed fixed-dictionary refit currently supports R1'
with torch.no_grad():
    Q,Rm=torch.linalg.qr(P[0].T);H=Q.T@m.G@Q
    outside=(m.G-Q@H@Q.T).square().sum()
w=m.tensor(w0).requires_grad_();off=m.tensor(z['W']).requires_grad_()
opt=torch.optim.Adam([{'params':[w],'lr':a.lr},{'params':[off],'lr':a.lr*a.offset_factor}])
rho=2.;mue=a.mue;mug=a.mug;best=float('inf');began=time.time();history=[]

def vals():
    A=Rm@w[0]
    e2=outside+(A@A.T-H).square().sum()
    hc=m.unpack(((w@off[...,None]).squeeze(-1)[:,None,:]@P).sum((0,1)))
    eig=torch.linalg.eigvalsh(m.h0-hc);ab=torch.sqrt(eig.square()+1e-14).sum()
    gap=m.constant+ab-eig.sum()+.5*off.square().sum()
    l1=torch.sqrt(w.square()+1e-16).sum(1)
    wb=torch.sqrt((off+w.sum(1)).square()+1e-16)
    lam=ab+.25*(l1+wb).square().sum()
    return e2,gap,lam,gap*(2*lam-gap)

def store(er):
    save(a.out,ut,w,tree,words,W=off.detach().cpu().numpy(),err2=er,lagrange_error=mue,lagrange_gap=mug)

def feasible_probe():
    global best
    # Exact PSD/rank-C projection in the current column gauge supplies a
    # feasible endpoint; interpolate only if needed to restore error margin.
    with torch.no_grad():
        ww=w.detach().clone();oo=off.detach().clone()
        if float(m.loss(ut,ww).sqrt())>a.cap:
            wp,_=m.project(ut,ww,1)
            if float(m.loss(ut,wp).sqrt())>a.cap:return
            lo,hi=0.,1.
            for _ in range(35):
                mid=(lo+hi)/2;trial=ww+mid*(wp-ww)
                if float(m.loss(ut,trial).sqrt())>a.cap:lo=mid
                else:hi=mid
            ww=ww+hi*(wp-ww)
        m.project_gap(ut,ww,oo,a.gap)
        g,l,s=m.physical(ut,ww,oo);sv=float(s.clamp(min=0).sqrt())
        if float(g)>3.91495 and sv<best:
            er=float(m.loss(ut,ww).sqrt());best=sv
            save(a.out,ut,ww,tree,words,W=oo.cpu().numpy(),err2=er,lagrange_error=mue,lagrange_gap=mug)

for step in range(a.steps+1):
    if step%20==0:
        with torch.no_grad():
            ee,g,l,s=vals();er=float(ee.clamp(min=0).sqrt());gv=float(g);sv=float(s.clamp(min=0).sqrt())
            if er<=a.cap and gv>=3.91495 and gv<2*float(l) and sv<best:
                best=sv;store(er)
            mue=max(0,mue+rho*(float(ee)/a.cap**2-1))
            mug=max(0,mug+rho*(a.gap-gv)/a.gap)
            if step%200==0 or step==a.steps:
                feasible_probe()
                row=dict(step=step,error_2=er,E_gap=gv,Lambda=float(l),lambda_eff=sv,best=best,
                    mue=mue,mug=mug,seconds=time.time()-began)
                history.append(row);print(json.dumps(row),flush=True)
    if step==a.steps:break
    opt.zero_grad();ee,g,l,s=vals()
    def pen(v,mu):return (torch.relu(mu+rho*v).square()-mu*mu)/(2*rho)
    loss=s/3000+pen(ee/a.cap**2-1,mue)+pen((a.gap-g)/a.gap,mug)
    loss.backward();opt.step()
    if a.project_gap and step%20==19:m.project_gap(ut,w,off,a.gap)
    factor=.1+.9*(1-step/a.steps)
    opt.param_groups[0]['lr']=a.lr*factor
    opt.param_groups[1]['lr']=a.offset_factor*a.lr*factor

save(str(Path(a.out).with_suffix(''))+'_last.npz',ut,w,tree,words,W=off.detach().cpu().numpy(),lagrange_error=mue,lagrange_gap=mug)
if np.isfinite(best):
    zz,uu,ww,tt,aw=load(a.out)
    cp,_=checkpoint(m,m.tensor(uu),m.tensor(ww),m.tensor(zz['W']),tt,2)
else:cp=None
write_json(HERE/'results'/f'{Path(a.out).stem}_refit.json',dict(arguments=vars(a),history=history,best=best if np.isfinite(best) else None,checkpoint=cp,seconds=time.time()-began))
