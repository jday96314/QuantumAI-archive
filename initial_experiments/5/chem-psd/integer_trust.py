"""Coupled integer-vector geometry moves with exact PSD variable projection.

All words of a row may change together. Physical unit-vector gradients propose
several row moves; exact tensor increments screen each proposal. Combined
moves require direct physical descent, then full block-PSD covariances refit.
This never uses floating-point rounding loss as an acceptance criterion.
The exact row-increment identity follows completed chem-blocks/quantized_screen.
"""
import argparse,time
from core import *
from solve import solve
from geometry_trust import gradient
from angle_gauge_probe import encode

@torch.no_grad()
def main():
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--tag',required=True)
    p.add_argument('--steps',type=int,default=100);p.add_argument('--inner',type=int,default=30);p.add_argument('--initial',type=int,default=400)
    p.add_argument('--newton',action='store_true')
    p.add_argument('--projected-trials',action='store_true')
    a=p.parse_args();z=dict(np.load(a.input));tree=z.get('tree_edges',np.load(INC)['tree_edges']);shape=z['U'].shape;R,B,_=shape
    words=z.get('angle_words',encode(z['U'].reshape(-1,76),tree,8,True)).reshape(-1,75)[-R*B:]
    initialwords=words.copy()
    u=reconstruct_u_vectors(words/128,tree).reshape(shape);m=Covariance(u)
    S0=m.from_weights(z['w']) if 'w' in z else None
    S,hi=solve(m,S0,a.initial,None,'fista',200,quiet=True);initial=m.diagnostics(S,True);start=time.time();records=[]
    def save():np.savez(HERE/f'parameters/{a.tag}.npz',U=u,w=m.weights(S).cpu().numpy(),S=S.cpu().numpy(),tree_edges=tree,angle_words=words)
    save()
    for k in range(a.steps):
        g,diagonal,err=gradient(m,S);newton_info=None
        if a.newton:
            from geometry_newton import direction as newton_direction
            direction,newton_info=newton_direction(m,S)
        else:direction=-g/diagonal.clamp(min=1e-6)
        direction*=torch.minimum(torch.ones_like(diagonal),.02/direction.norm(dim=2,keepdim=True).clamp(min=1e-20))
        approximation,_=m.direct(S);E=approximation-m.H;w=m.weights(S)
        Y=torch.linalg.solve_triangular(m.Rm,S@m.Q,upper=True).flatten(0,1)
        alpha=w.square().sum(2).flatten();EY=Y@E;yy=Y.square().sum(1)
        P=m.P.flatten(0,1);bestgain=torch.zeros(R*B,**D);bestwords=words.copy();row_checks=[]
        for scale in [.5,1.,2.,4.]:
            proposal=tensor(u)+scale*direction;proposal/=proposal.norm(dim=2,keepdim=True)
            nw=encode(proposal.cpu().numpy().reshape(-1,76),tree,8,True)
            nu=reconstruct_u_vectors(nw/128,tree);d=tensor(features(nu))-P
            dd=d.square().sum(1);yd=(Y*d).sum(1)
            inc=4*(d*EY).sum(1)+2*alpha*(d*(d@E)).sum(1)
            inc+=2*yy*dd+2*yd.square()+4*alpha*yd*dd+alpha.square()*dd.square()
            improve=inc<bestgain;bestgain[improve]=inc[improve];mask=improve.cpu().numpy();bestwords[mask]=nw[mask]
        order=bestgain.argsort();order=order[bestgain[order]<-1e-13];accepted=False;newwords=None;projectedS=None;projected_model=None;projected_attempts=0
        if k==0 and len(order):
            j=int(order[0]);control=words.copy();control[j]=bestwords[j]
            cu=reconstruct_u_vectors(control/128,tree).reshape(shape);cx=(w.transpose(1,2)@tensor(features(cu))).flatten(0,1)
            direct=float((cx.T@cx-m.H).square().sum());predicted=err*err+float(bestgain[j])
            assert abs(direct-predicted)<1e-8,(direct,predicted)
        for count in sorted(set([len(order),256,128,64,32,16,8,4,2,1]),reverse=True):
            if not 0<count<=len(order):continue
            selected=order[:count].cpu().numpy();testwords=words.copy();testwords[selected]=bestwords[selected]
            testu=reconstruct_u_vectors(testwords/128,tree).reshape(shape);tp=tensor(features(testu));X=(w.transpose(1,2)@tp).flatten(0,1)
            trial=float((X.T@X-m.H).square().sum())
            if count==1:
                predicted=err*err+float(bestgain[order[0]])
                assert abs(trial-predicted)<1e-8,(trial,predicted)
            if trial<err*err-1e-13:
                accepted=True;u=testu;newwords=testwords;break
            if a.projected_trials and count>=16 and projected_attempts<5:
                projected_attempts+=1;tm=Covariance(testu,m.H)
                ts,th=solve(tm,tm.from_weights(w),50,None,'fista',50,quiet=True)
                projected_error=tm.direct(ts)[1]
                if projected_error<err-1e-9:
                    accepted=True;u=testu;newwords=testwords;projectedS=ts;projected_model=tm;break
        if not accepted:reason='No improving decoded vector proposal';break
        words=newwords;nm=projected_model if projected_model is not None else Covariance(u,m.H)
        warm=projectedS if projectedS is not None else nm.from_weights(w)
        S,h=solve(nm,warm,a.inner,None,'fista',a.inner,quiet=True);m=nm
        direct=m.direct(S)[1];assert direct*direct<=trial+1e-8;save()
        rec=dict(step=k,error_2=direct,previous_error_2=err,accepted_rows=count,proposed_improving_rows=len(order),
                 changed_integer_words_from_initial=int(np.count_nonzero(testwords!=initialwords)),
                 fixed_covariance_trial_error_2=trial**.5,projected_trial_acceptance=projectedS is not None,projected_attempts=projected_attempts,
                 newton=newton_info,inner_projected_gradient=h[-1]['projected_gradient_norm'],seconds=time.time()-start)
        records.append(rec)
        if k%5==0:print(json.dumps(rec),flush=True)
    else:reason='Finite integer geometry budget'
    S,hf=solve(m,S,2000,None,'fista',200,quiet=True);save()
    report(HERE/f'results/{a.tag}.json',dict(input=a.input,initial=initial,final=m.diagnostics(S,True),termination=reason,
          history=records,initial_history=hi,final_history=hf,device=D['device'],seconds=time.time()-start))
if __name__=='__main__':main()
