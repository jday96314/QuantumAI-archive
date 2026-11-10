"""Quantized physical controls, metric-aware covariance and offset refits.

Native evaluate_words is used for final physical checks. Circuit simulation is
reserved for promising fits. This script cannot claim native oracle acceptance.
"""
import argparse,time
import numpy as np
import torch
from core import *
from angle_gauge_probe import encode
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors,build_h1_prime,reconstruct_x
from dfthcsos.verify import load_inputs,evaluate_words,passes_accuracy
from resources import row as resource_row

def covariance_refit(U,w,sp,cycles,C):
    d=dict(device='cuda',dtype=torch.float64);R,B,N=U.shape
    P=torch.tensor(feature(U,sp),**d);I=torch.eye(B,**d)
    L=torch.linalg.cholesky(P@P.transpose(1,2)+1e-12*I)
    Q=torch.linalg.solve_triangular(L,P,upper=False)
    F=L.transpose(1,2)@torch.tensor(w,**d);H=F@F.transpose(1,2)
    G=torch.tensor(gram(sp),**d)
    A=(Q.reshape(R*B,-1)@Q.reshape(R*B,-1).T).reshape(R,B,R,B).permute(0,2,1,3).contiguous()
    T=Q@G@Q.transpose(1,2);constant=G.square().sum();best=float('inf');bestH=None;records=[]
    def err(H):
        cross=(A@H[None])@A.transpose(2,3)
        return (constant-2*(H*T).sum()+(H[:,None]*cross.transpose(2,3)).sum()).clamp(min=0).sqrt()
    for cycle in range(cycles+1):
        if cycle:
          for r in range(R):
            cross=(A[r]@H)@A[r].transpose(1,2)
            V=T[r]-cross.sum(0)+cross[r];e,v=torch.linalg.eigh((V+V.T)/2)
            H[r]=(v[:,-C:]*e[-C:].clamp(min=0))@v[:,-C:].T
        er=float(err(H))
        if er<best:best=er;bestH=H.clone()
        if cycle%10==0 or cycle==cycles:records.append(dict(cycle=cycle,error_2=er))
    e,v=torch.linalg.eigh(bestH);F=v[:,:,-C:]*e[:,-C:].clamp(min=0).sqrt()[:,None,:]
    ww=torch.linalg.solve_triangular(L.transpose(1,2),F,upper=True).cpu().numpy()
    X=np.einsum('rbc,rbd->rcd',ww,feature(U,sp)).reshape(R*C,-1)
    direct=float(np.linalg.norm(X.T@X-gram(sp)));assert abs(direct-best)<1e-6
    return ww,direct,records

def quantization_aware(U,w,tree,bits,sp,steps,seed):
    # Differentiate through variable projection of physical generator targets.
    # The forward pass always uses actual rounded integer words. This is a
    # factorization-aware heuristic, not a proof of low-bit capacity.
    d=dict(device='cuda',dtype=torch.float64);R,B,N=U.shape
    target=torch.tensor(np.einsum('rbc,rbd->rcd',w,feature(U,sp)),**d)
    words=encode(U.reshape(-1,76),tree,bits,True)
    init=words.astype(float);init[init>=1<<(bits-1)]-=1<<bits
    rng=np.random.default_rng(seed);theta=torch.tensor(init+rng.uniform(-.15,.15,init.shape),requires_grad=True,**d)
    fac=torch.tensor(sp['fac'],**d);I=torch.eye(B,**d);G=torch.tensor(gram(sp),**d);constant=G.square().sum()
    stages=[]
    for t,(p,c) in enumerate(tree):
        p,c=int(p),int(c)
        if not stages or any(p in (pp,cc) or c in (pp,cc) for _,pp,cc in stages[-1]):stages.append([])
        stages[-1].append((t,p,c))
    stages=[tuple(torch.tensor(x,device='cuda') for x in zip(*stage)) for stage in stages]
    start=torch.zeros((R*B,76),**d);start[:,0]=1
    opt=torch.optim.Adam([theta],lr=.04);best=float('inf');beststate=None;records=[]
    for step in range(steps+1):
        opt.zero_grad();rounded=theta+(theta.round()-theta).detach();u=start
        for ts,ps,cs in stages:
            ang=rounded[:,ts]*(np.pi/(1<<(bits-1)));co,sn=ang.cos(),ang.sin();up,uc=u[:,ps],u[:,cs]
            new=u.clone();new[:,ps]=co*up-sn*uc;new[:,cs]=sn*up+co*uc;u=new
        Uq=u.reshape(R,B,76);P=Uq[:,:,sp['ii']]*Uq[:,:,sp['jj']]*fac
        ww=torch.linalg.solve(P@P.transpose(1,2)+1e-10*I,P@target.transpose(1,2))
        X=(ww.transpose(1,2)@P).reshape(R*w.shape[2],-1)
        error=constant-2*(X*(X@G)).sum()+(X@X.T).square().sum()
        er=float(error.detach().clamp(min=0).sqrt())
        if er<best:
            best=er;beststate=(Uq.detach().cpu().numpy(),ww.detach().cpu().numpy(),theta.detach().round().cpu().numpy().astype(np.int64)%(1<<bits))
        if step%25==0 or step==steps:records.append(dict(step=step,error_2=er,best=best))
        if step<steps:error.backward();opt.step()
    return beststate,best,records

def offset_fit(U,w,inputs,sp,incumbent):
    # Constrained scalar offset refit: least-squares to the incumbent one-body
    # correction, followed by minimum-gap convex optimization and a tiny
    # outward adjustment if needed to satisfy the strict lower gap.
    X=np.einsum('rbc,rbd->rcd',w,feature(U,sp)).reshape(-1,len(sp['ii']))
    old=feature(incumbent['U'][0],sp);Xold=incumbent['w'][0].T@old
    correction=incumbent['offset'].ravel()@Xold
    off=np.linalg.solve(X@X.T+1e-8*np.eye(len(X)),X@correction)
    d=dict(device='cuda',dtype=torch.float64);m=torch.tensor(unpack(X,sp),**d)
    h0=torch.tensor(inputs.ref['h1_thc']-.5*inputs.ref['electrons']*inputs.S+inputs.beta_1*np.eye(76),**d)
    const=inputs.energy_upper-inputs.ref['ecore']+inputs.ref['electrons']*inputs.beta_1
    v=torch.tensor(off,requires_grad=True,**d)
    def gap():
        eig=torch.linalg.eigvalsh(h0-torch.einsum('c,cpq->pq',v,m))
        return const+2*torch.relu(-eig).sum()+.5*v.square().sum()
    initial=float(gap().detach());opt=torch.optim.LBFGS([v],max_iter=120,history_size=20,line_search_fn='strong_wolfe',tolerance_grad=1e-11,tolerance_change=1e-14)
    def closure():opt.zero_grad();g=gap();g.backward();return g
    opt.step(closure);minimum=float(gap().detach())
    # If the attainable gap falls below the floor, choose a nearby valid point.
    if minimum<3.9153:
        with torch.no_grad():
            direction=torch.tensor(np.random.default_rng(846).standard_normal(len(X)),**d);direction/=direction.norm()
            origin=v.clone();lo,hi=0.,1.
            while True:
                v.copy_(origin+hi*direction)
                if float(gap())>3.9153:break
                hi*=2
            for _ in range(45):
                mid=(lo+hi)/2;v.copy_(origin+mid*direction)
                if float(gap())<3.9153:lo=mid
                else:hi=mid
            v.copy_(origin+hi*direction)
    return v.detach().cpu().numpy().reshape(w.shape[0],w.shape[2]),dict(initial_gap=initial,minimum_gap=minimum,refitted_gap=float(gap().detach()))

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--bits',type=int,nargs='+',default=[7,8])
    p.add_argument('--cycles',type=int,default=60);p.add_argument('--qat-steps',type=int,default=100);a=p.parse_args();torch.set_num_threads(2);beg=time.time()
    z=dict(np.load(a.input));inc,sp=load();R,B,N=z['U'].shape;C=z['w'].shape[2]
    inputs=load_inputs(SOURCE);rows=[]
    for bits in a.bits:
        words=encode(z['U'].reshape(-1,76),z['tree_edges'],bits,True)
        U=reconstruct_u_vectors(words/(1<<(bits-1)),z['tree_edges']).reshape(R,B,N)
        X=np.einsum('rbc,rbd->rcd',z['w'],feature(U,sp)).reshape(R*C,-1)
        direct=float(np.linalg.norm(X.T@X-gram(sp)))
        ww,er,records=covariance_refit(U,z['w'],sp,a.cycles,C)
        options=[('nearest_covariance',U,ww,words,er)]
        qat_records=[]
        if a.qat_steps:
            state,qa,qat_records=quantization_aware(z['U'],z['w'],z['tree_edges'],bits,sp,a.qat_steps,178+bits)
            qu,qw,qwords=state;qw,qe,qrec=covariance_refit(qu,qw,sp,a.cycles,C)
            options.append(('factor_aware_covariance',qu,qw,qwords,qe))
        checks=[]
        for name,qu,qw,qwords,qe in options:
            # Restore the exact signed gauge if a discrete step crossed a
            # quadrant boundary; only the vector's global sign may change.
            qwords=encode(qu.reshape(-1,N),z['tree_edges'],bits,True)
            decoded=reconstruct_u_vectors(qwords/(1<<(bits-1)),z['tree_edges']).reshape(R,B,N)
            assert np.max(np.abs(feature(decoded,sp)-feature(qu,sp)))<1e-12
            offset,gapinfo=offset_fit(qu,qw,inputs,sp,inc)
            XX=reconstruct_x(qu,qw)
            hp=build_h1_prime(inputs.ref['h1_thc'],inputs.S,inputs.beta_1,inputs.ref['electrons'],offset,XX)
            _,v=np.linalg.eigh(hp);one=encode(v.T,z['tree_edges'],bits,True)
            full=np.concatenate([one,qwords]);params=dict(U=qu,W=qw,H2_I=offset+qw.sum(1),H1_I=inc['H1_I'],Hsym_lower=inc['Hsym_lower'],tree_edges=z['tree_edges'],angle_words=full)
            native=load_inputs(params);metrics=evaluate_words(native,full,z['tree_edges'],bits)
            assert abs(metrics['error_2']-qe)<2e-6,(metrics['error_2'],qe)
            entry=dict(method=name,metrics=metrics,passes_physical=passes_accuracy(metrics),margins=dict(error_1=.02-metrics['error_1'],error_2=.02-metrics['error_2'],Lambda=metrics['Lambda']-49,E_gap=metrics['E_gap']-3.914898),offset_refit=gapinfo)
            choices=[resource_row(R,B,C,bits,J) for J in range(4,16)]
            choices=[c for c in choices if c]
            rr=min(choices,key=lambda c:c['cost_conditional_upper'])
            entry['resource_bound']={k:rr[k] for k in ['J','k2','oracle_kind','T_upper','Q_upper','oracle_Q_upper']}
            entry['resource_bound']['cost_at_measured_lambda_eff']=metrics['lambda_eff']*rr['T_upper']*rr['Q_upper']
            entry['accepted_task_score']=None
            checks.append(entry);emit(dict(bits=bits,**entry))
            np.savez(HERE/'candidates'/f'{Path(a.input).stem}_q{bits}_{name}.npz',**params)
        rows.append(dict(bits=bits,direct_rounded_error_2=direct,covariance_refit_error_2=er,covariance_records=records,quantization_aware_records=qat_records,physical_checks=checks))
    data=dict(status='NATIVE_PHYSICAL_WORD_EVALUATION_ONLY',oracle_simulation='NOT_RUN; no passing physical candidate asserted',input=a.input,shape=[R,B,C],rows=rows,seconds=time.time()-beg)
    report(Path(a.input).stem+'_quantized.json',data);emit(data)

if __name__=='__main__':main()
