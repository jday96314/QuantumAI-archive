"""Actual decoded-angle screens with metric-correct constrained weight refits.

The PSD updates whiten the projector Gram, and retain at most C columns.
A bounded integer angle search accepts only true full-tensor improvements;
it re-fits ensemble covariances again afterward. No rounded-fit impossibility
claim is supported by this local screen.
"""
import argparse
import time
from pathlib import Path
import numpy as np
import torch
from common import HERE,load,SOURCE,feature,unpack,report
from angle_gauge_probe import encode
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors,evaluate_dfthc_angle_words
from dfthcsos.verify import load_inputs

def refit(U,w,G,C,ii,jj,fac,cycles):
    R,B,N=U.shape;d=dict(device=U.device,dtype=U.dtype)
    P=U[:,:,ii]*U[:,:,jj]*fac
    L=torch.linalg.cholesky(P@P.transpose(1,2)+1e-12*torch.eye(B,**d))
    Q=torch.linalg.solve_triangular(L,P,upper=False)
    F=L.transpose(1,2)@w;H=F@F.transpose(1,2)
    flat=Q.reshape(R*B,-1)
    A=(flat@flat.T).reshape(R,B,R,B).permute(0,2,1,3).contiguous()
    T=Q@G@Q.transpose(1,2)
    def error(H):
        cross=(A@H[None])@A.transpose(2,3)
        return G.square().sum()-2*(H*T).sum()+(H[:,None]*cross.transpose(2,3)).sum()
    best=float(error(H));bestH=H.clone();start=max(0,best)**.5
    for k in range(cycles):
        for r in range(R):
            cross=(A[r]@H)@A[r].transpose(1,2)
            residual=T[r]-cross.sum(0)+cross[r]
            e,v=torch.linalg.eigh((residual+residual.T)/2)
            e=e[-C:].clamp(min=0);v=v[:,-C:]
            H[r]=(v*e)@v.T
        err=float(error(H))
        if err<best:best=err;bestH=H.clone()
    e,v=torch.linalg.eigh(bestH)
    F=v[:,:,-C:]*e[:,-C:].clamp(min=0).sqrt()[:,None,:]
    wn=torch.linalg.solve_triangular(L.transpose(1,2),F,upper=True)
    X=(wn.transpose(1,2)@P).reshape(R*C,-1)
    direct=float((X.T@X-G).norm())
    assert abs(direct-max(0,best)**.5)<1e-7
    return wn,dict(cycles=cycles,initial_error_2=start,best_error_2=direct)

def decode_torch(words,tree,bits):
    u=torch.zeros((len(words),76),device=words.device,dtype=words.dtype);u[:,0]=1
    for t,(p,c) in enumerate(tree):
        angle=words[:,t]*(np.pi/(1<<(bits-1)))
        co,sn=angle.cos(),angle.sin();up,uc=u[:,p],u[:,c]
        v=u.clone();v[:,p]=co*up-sn*uc;v[:,c]=sn*up+co*uc;u=v
    return u

def integer_search(words,w,G,tree,bits,ii,jj,fac,moves):
    """Greedy coordinate proposals selected by full physical tensor gradients.

    Gradients rank candidates only; all accepted +/-1 integer changes must
    reduce the exact full tensor norm after redecoding the Givens tree.
    """
    R,B,C=w.shape;d=dict(device=G.device,dtype=G.dtype)
    theta=torch.tensor(words,dtype=G.dtype,device=G.device).requires_grad_()
    history=[];initial=None;accepted=0
    for move in range(moves):
        u=decode_torch(theta,tree,bits)
        P=(u[:,ii]*u[:,jj]*fac).reshape(R,B,-1)
        X=(w.transpose(1,2)@P).reshape(R*C,-1)
        delta=X.T@X-G;loss=delta.square().sum()
        if initial is None:initial=float(loss.detach().sqrt())
        loss.backward()
        with torch.no_grad():
            grad=theta.grad.detach();idx=grad.abs().flatten().topk(16).indices
            best=float(loss);bestchange=None
            for flatidx in idx.tolist():
                row,t=divmod(flatidx,75);r,b=divmod(row,B)
                step=-1 if grad[row,t]>0 else 1
                trial=theta[row:row+1].clone();trial[0,t]+=step
                unew=decode_torch(trial,tree,bits)[0]
                pnew=unew[ii]*unew[jj]*fac
                change=pnew-P[r,b]
                y=w[r,b]@X[r*C:(r+1)*C]
                alpha=w[r,b].square().sum()
                # Exact norm of E + y*d^T + d*y^T + alpha*d*d^T.
                dd=change.square().sum();yd=y@change
                increment=4*change@(delta@y)+2*alpha*change@(delta@change)
                increment+=2*y.square().sum()*dd+2*yd.square()+4*alpha*yd*dd+alpha.square()*dd.square()
                val=float(loss+increment)
                if val<best-1e-13:best=val;bestchange=(row,t,step)
            theta.grad=None
            if bestchange is None:break
            row,t,step=bestchange;theta[row,t]+=step;accepted+=1
            checku=decode_torch(theta,tree,bits)
            checkP=(checku[:,ii]*checku[:,jj]*fac).reshape(R,B,-1)
            checkX=(w.transpose(1,2)@checkP).reshape(R*C,-1)
            direct=float((checkX.T@checkX-G).square().sum())
            assert abs(direct-best)<1e-8,(direct,best)
            assert direct<float(loss),('nonimproving integer move',direct,float(loss))
            history.append(dict(move=move,row=row,angle=t,step=step,error_2=best**.5))
    out=theta.detach().cpu().numpy().round().astype(np.int64)%(1<<bits)
    # Restore the signed gauge exactly, preserving all row projectors.
    uu=reconstruct_u_vectors(out/(1<<(bits-1)),tree)
    out=encode(uu,tree,bits,True)
    return out,dict(proposed_rounds=moves,accepted_moves=accepted,initial_error_2=initial,history=history)

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--bits',type=int,nargs='+',default=[7,8])
    p.add_argument('--cycles',type=int,default=40);p.add_argument('--moves',type=int,default=32);a=p.parse_args()
    torch.set_num_threads(2);begin=time.time();z,sp=load();fit=dict(np.load(a.input));tree=fit['tree_edges'];R,B,N=fit['U'].shape;C=fit['w'].shape[2]
    d=dict(device='cuda',dtype=torch.float64);ii,jj=sp['ii'],sp['jj'];fac=torch.tensor(sp['fac'],**d)
    G=torch.tensor((sp['eigenvectors']*sp['eigenvalues'])@sp['eigenvectors'].T,**d)
    native=load_inputs(SOURCE)
    # Independently bind cached spectrum to the actual FCIDUMP and fixed shift.
    tei=native.ref['tei'];S=native.S
    ident=(ii==jj).astype(float);svec=S[ii,jj]*sp['fac']
    target=tei[ii[:,None],jj[:,None],ii[None,:],jj[None,:]]*sp['fac'][:,None]*sp['fac'][None,:]
    target+=.5*(svec[:,None]*ident[None,:]+ident[:,None]*svec[None,:])
    discrepancy=float(np.max(np.abs(target-G.cpu().numpy())))
    frob_discrepancy=float(np.linalg.norm(target-G.cpu().numpy()))
    print(__import__('json').dumps(dict(cached_target_vs_host_max_discrepancy=discrepancy,cached_target_vs_host_frobenius_discrepancy=frob_discrepancy)),flush=True)
    # The inherited target is extremely close but not bit-identical. Use the
    # actual host integrals and fixed shift for every quantized objective.
    assert frob_discrepancy<1e-6,frob_discrepancy
    G=torch.tensor(target,**d)
    baseF=np.einsum('rbc,rbd->rcd',z['w'],feature(z['U'],sp)).reshape(350,-1)
    offset_target=z['offset'].reshape(-1)@baseF
    records=[]
    for bits in a.bits:
        words=encode(fit['U'].reshape(-1,76),tree,bits,True)
        un=reconstruct_u_vectors(words/(1<<(bits-1)),tree).reshape(R,B,76)
        U=torch.tensor(un,**d);w=torch.tensor(fit['w'],**d)
        P=U[:,:,ii]*U[:,:,jj]*fac;X=(w.transpose(1,2)@P).reshape(R*C,-1)
        rounded=float((X.T@X-G).norm())
        with torch.no_grad():w,first=refit(U,w,G,C,ii,jj,fac,a.cycles)
        words,integer=integer_search(words,w,G,tree,bits,ii,jj,fac,a.moves)
        un=reconstruct_u_vectors(words/(1<<(bits-1)),tree).reshape(R,B,76)
        U=torch.tensor(un,**d)
        with torch.no_grad():w,second=refit(U,w,G,C,ii,jj,fac,a.cycles)
        wn=w.cpu().numpy();Pf=feature(un,sp);F=np.einsum('rbc,rbd->rcd',wn,Pf).reshape(R*C,-1)
        # Offset refit uses the associated full symmetric contraction, with a
        # small ridge only on nearly null modes, never componentwise clipping.
        off=np.linalg.solve(F@F.T+1e-8*np.eye(R*C),F@offset_target).reshape(R,C)
        Xmat=unpack(F,sp).reshape(R,C,76,76)
        h0=native.ref['h1_thc']-.5*native.ref['electrons']*S+native.beta_1*np.eye(76)
        h=h0-np.einsum('rc,rcpq->pq',off,Xmat)
        _,v=np.linalg.eigh((h+h.T)/2)
        one=encode(v.T,tree,bits,True);allwords=np.r_[one,words]
        metrics=evaluate_dfthc_angle_words(allwords,wn,off,S,native.beta_1,native.ref,brot=bits,tree_edges=tree,energy_upper=native.energy_upper)
        assert abs(metrics['error_2']-second['best_error_2'])<1e-7
        out=dict(U=un,w=wn,W=off,angle_words=words,onebody_angle_words=one,onebody_brot=bits,brot=bits,tree_edges=tree,err2=metrics['error_2'])
        stem=Path(a.input).stem+f'_q{bits}';np.savez(HERE/f'parameters/{stem}.npz',**out)
        margins=dict(error_1=.02-metrics['error_1'],error_2=.02-metrics['error_2'],Lambda=metrics['Lambda']-49,E_gap=metrics['E_gap']-3.914898)
        valid=(margins['error_1']>=0 and margins['error_2']>=0 and margins['Lambda']>0 and margins['E_gap']>0 and bool(metrics['gap_in_domain']))
        record=dict(bits=bits,rounded_error_2=rounded,first_covariance_refit=first,integer_search=integer,second_covariance_refit=second,
                    native_physical_metrics=metrics,margins=margins,passes_physical_screen=valid,
                    offset_contraction_refit_residual=float(np.linalg.norm(off.reshape(-1)@F-offset_target)))
        records.append(record)
        print(__import__('json').dumps(record),flush=True)
    report(HERE/f'results/{Path(a.input).stem}_quantized.json',dict(status='QUANTIZED_PHYSICAL_SCREEN_ONLY',input=a.input,
        cached_target_vs_host_max_discrepancy=discrepancy,cached_target_vs_host_frobenius_discrepancy=frob_discrepancy,
        quantized_objective='exact native FCIDUMP plus fixed-shift target',R=R,B=B,C=C,records=records,seconds=time.time()-begin,
        limitation='Bounded local quantizer/refits; failure is not a precision impossibility proof. No candidate oracle/native task acceptance is claimed.'))

if __name__=='__main__':main()
