"""Globalized PSD variable projection with feasible geometry line searches.

An inexact inner solve can fool ordinary reduced-objective L-BFGS. Here a
geometry trial retains the current atom covariance, providing a feasible
upper bound on the variable-projected objective. Armijo acceptance of this
bound guarantees descent before the next convex PSD optimization. Curvature
history and a physical Gauss-Newton diagonal precondition the unit vectors.
"""
import argparse,time,signal,fcntl
from core import *
from solve import solve

@torch.no_grad()
def gradient(m,S):
    approx,err=m.direct(S);E=approx-m.H
    KP=torch.linalg.solve_triangular(m.Rm,S@m.Q,upper=True);gp=4*KP@E
    M=torch.zeros((*gp.shape[:2],76,76),**D);M[:,:,ii,jj]=gp/tensor(fac);M[:,:,jj,ii]=gp/tensor(fac)
    u=tensor(m.u);g=2*(M@u[:,:,:,None]).squeeze(-1);g-=u*(g*u).sum(2,keepdim=True)
    diagonal=8*KP.square().sum(2,keepdim=True)
    return g,diagonal,err

@torch.no_grad()
def fit(a):
    z=dict(np.load(a.input));m=Covariance(z['U']);u=tensor(z['U']);start=time.time();hist=[];memory=[]
    S0=m.from_weights(z['w']) if 'w' in z else None
    S,initialhist=solve(m,S0,a.initial,None,'fista',200,quiet=True)
    initial=m.diagnostics(S,True);g,diag,err=gradient(m,S);best=err
    def save():
        np.savez(HERE/f'parameters/{a.tag}.npz',U=u.cpu().numpy(),S=S.cpu().numpy(),w=m.weights(S).cpu().numpy(),tree_edges=z['tree_edges'])
    save();stopping=False
    def stop(signum,frame):
        nonlocal stopping
        stopping=True
    signal.signal(signal.SIGTERM,stop)
    for k in range(a.steps):
        if stopping:
            reason='Graceful checkpoint stop';break
        # Limited-memory curvature update of the reduced physical objective.
        q=g.clone();alphas=[]
        for ss,yy,rho in reversed(memory):
            alpha=rho*(ss*q).sum();alphas.append(alpha);q-=alpha*yy
        direction=q/diag.clamp(min=1e-6)
        for (ss,yy,rho),alpha in zip(memory,reversed(alphas)):
            beta=rho*(yy*direction).sum();direction+=ss*(alpha-beta)
        direction=-direction;direction-=u*(direction*u).sum(2,keepdim=True)
        if float((g*direction).sum())>=0:
            memory=[];direction=-g/diag.clamp(min=1e-6)
        direction*=min(1.,a.radius/float(direction.norm(dim=2).max()))
        gd=float((g*direction).sum());w=m.weights(S);old_u=u.clone();old_g=g.clone();old_err=err
        accepted=False;step=1.;trials=[]
        for ls in range(24):
            un=u+step*direction;un/=un.norm(dim=2,keepdim=True)
            P=un[:,:,ii]*un[:,:,jj]*tensor(fac);X=(w.transpose(1,2)@P).flatten(0,1)
            # Direct physical pair-space norm is used for acceptance; cached
            # contraction cancellation cannot invent an improvement.
            trial=float((X.T@X-m.H).square().sum());trials.append(trial)
            if trial<err*err+1e-4*step*gd:
                accepted=True;u=un;break
            step*=.5
        if not accepted:
            if memory:memory=[];continue
            reason='No feasible Armijo geometry step';break
        nm=Covariance(u.cpu().numpy(),m.H)
        warm=nm.from_weights(w)
        inner=a.inner if k%25 else max(400,a.inner)
        S,innerhist=solve(nm,warm,inner,None,'fista',inner,quiet=True);m=nm
        g,diag,err=gradient(m,S)
        # The covariance inner step must preserve the accepted feasible bound.
        assert err*err<=trial+1e-8,(err*err,trial)
        ss=u-old_u;yy=g-old_g;sy=float((ss*yy).sum())
        if sy>1e-14*float(ss.norm()*yy.norm()):
            memory.append((ss,yy,1/sy));memory=memory[-20:]
        if err<best:best=err;save()
        rec=dict(step=k,error_2=err,previous_error_2=old_err,fixed_covariance_trial_error_2=trial**.5,
             line_search_scale=step,line_search_trials=len(trials),seconds=time.time()-start,
             inner_projected_gradient=innerhist[-1]['projected_gradient_norm'],weight_squared_norm=float(w.square().sum()),
             maximum_vector_change=float((u-old_u).norm(dim=2).max()))
        hist.append(rec)
        if k%5==0:
            print(json.dumps(rec),flush=True)
            (HERE/f'results/{a.tag}_progress.json').write_text(json.dumps(dict(initial=initial,history=hist,device=D['device']),indent=2)+'\n')
    else:reason='Finite outer iteration budget'
    # An independent long convex solve certifies the final fixed geometry.
    S,finalhist=solve(m,S,a.final,None,'fista',200,quiet=True);g,diag,err=gradient(m,S);save()
    report(HERE/f'results/{a.tag}.json',dict(input=a.input,method='PSD variable projection with feasible trust line search',
         initial=initial,final=m.diagnostics(S,True),termination=reason,history=hist,initial_history=initialhist,final_history=finalhist,
         outer_steps=len(hist),device=D['device'],seconds=time.time()-start))

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--tag',required=True)
    p.add_argument('--steps',type=int,default=250);p.add_argument('--inner',type=int,default=80)
    p.add_argument('--initial',type=int,default=1200);p.add_argument('--final',type=int,default=2400)
    p.add_argument('--radius',type=float,default=.02);a=p.parse_args()
    with (HERE/f'results/{a.tag}.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if (HERE/f'results/{a.tag}.json').exists():
            print('Existing completed result: '+a.tag,flush=True)
        else:fit(a)
if __name__=='__main__':main()
