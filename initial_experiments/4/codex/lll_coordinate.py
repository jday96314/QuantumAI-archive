"""Coordinated integer-angle jumps using a reduced local curvature lattice.

This replaces unit-step/ternary proposals with potentially large moves across
all angles of an atom. Each move must lower the actual rank-constrained projected
loss. It writes research fits only; native chemistry verification is separate.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
import torch
from lattice_babai import proposals

p=argparse.ArgumentParser();p.add_argument('--input',default='candidates/wide700_q8_c350.npz')
p.add_argument('--output',default='candidates/lll700_q8_c350.npz');p.add_argument('--steps',type=int,default=100)
p.add_argument('--rows',type=int,default=12);p.add_argument('--trials',type=int,default=6)
p.add_argument('--padding',type=float,default=.005);p.add_argument('--max-step',type=int,default=24)
p.add_argument('--ridge',type=float,default=1e-6);p.add_argument('--columns',type=int,default=350)
p.add_argument('--global-dim',type=int,default=0)
p.add_argument('--paired-rows',type=int,default=0)
p.add_argument('--moves-per-step',type=int,default=1)
a=p.parse_args();torch.set_num_threads(3);root=Path(__file__).resolve().parent
z=dict(np.load(a.input));tree=z['tree_edges'];bits=int(z['brot']);words=z['angle_words'].reshape(-1,75)
B,T=words.shape;n=76;d=dict(device='cuda',dtype=torch.float64)
theta=torch.tensor(words,requires_grad=True,**d);alpha=np.pi/(1<<(bits-1))
sp=np.load(root/'../../3/parameters/reference_spectrum.npz');ev,V=sp['eigenvalues'],sp['eigenvectors']
G=torch.tensor((V*np.maximum(ev,0))@V.T,**d);Gactual=torch.tensor((V*ev)@V.T,**d)
const=G.square().sum();actual_const=Gactual.square().sum()
ii=torch.tensor(sp['ii'],device='cuda');jj=torch.tensor(sp['jj'],device='cuda');fac=torch.tensor(sp['fac'],**d)
start=torch.zeros((B,n),**d);start[:,0]=1;eye=torch.eye(B,**d);stages=[]
for t,(pa,ch) in enumerate(tree):
    pa,ch=int(pa),int(ch)
    if not stages or any(pa in (x,y) or ch in (x,y) for _,x,y in stages[-1]):stages.append([])
    stages[-1].append((t,pa,ch))
stages=[tuple(torch.tensor(x,device='cuda') for x in zip(*s)) for s in stages]


def decode():
    u=start
    for ts,ps,cs in stages:
        ang=theta[:,ts]*alpha;co,sn=ang.cos(),ang.sin();up,uc=u[:,ps],u[:,cs]
        new=u.clone();new[:,ps]=co*up-sn*uc;new[:,cs]=sn*up+co*uc;u=new
    return u


def state():
    u=decode();P=u[:,ii]*u[:,jj]*fac;S=(u@u.T).square();H=P@G@P.T;Sr=S+a.ridge*eye
    with torch.no_grad():
        L=torch.linalg.cholesky(Sr);hh=torch.linalg.solve_triangular(L,H,upper=False)
        hh=torch.linalg.solve_triangular(L,hh.T,upper=False).T;e,v=torch.linalg.eigh(.5*(hh+hh.T))
        w=torch.linalg.solve_triangular(L.T,v[:,-a.columns:]*e[-a.columns:].clamp(min=0).sqrt(),upper=True)
        K=w@w.T
    KS=K@Sr;loss=const-2*(K*H).sum()+(KS*KS.T).sum()
    with torch.no_grad():
        Hactual=P@Gactual@P.T;phys=K@S
        error=(actual_const-2*(K*Hactual).sum()+(phys*phys.T).sum()).clamp(min=0).sqrt()
    return loss,error,u,K,S,L,w


def jacobian():
    u=start.clone();J=torch.zeros((B,T,n),**d)
    for t,(pa,ch) in enumerate(tree):
        pa,ch=int(pa),int(ch);ang=theta[:,t]*alpha;co,sn=ang.cos(),ang.sin()
        up,uc=u[:,pa].clone(),u[:,ch].clone();jp,jc=J[:,:,pa].clone(),J[:,:,ch].clone()
        J[:,:,pa]=co[:,None]*jp-sn[:,None]*jc;J[:,:,ch]=sn[:,None]*jp+co[:,None]*jc
        u[:,pa]=co*up-sn*uc;u[:,ch]=sn*up+co*uc
        J[:,t,pa]=-alpha*u[:,ch];J[:,t,ch]=alpha*u[:,pa]
    return J


def save_fit(er,u,w):
    # Keep readers from observing a partially written NPZ.
    target=root/a.output
    temporary=target.with_name(target.stem+'.tmp.npz')
    np.savez(temporary,U=u.cpu().numpy()[None,:,:],w=w.cpu().numpy()[None,:,:],W=np.zeros((1,a.columns)),
             angle_words=theta.detach().cpu().numpy().astype(np.int64)%(1<<bits),tree_edges=tree,brot=bits,err2=er)
    temporary.replace(target)


best=float('inf');began=time.time();total_moves=0
for step in range(a.steps):
    theta.grad=None;loss,error,u,K,S,L,w=state();loss.backward();g=theta.grad.detach();value=float(loss.detach())
    with torch.no_grad():
        er=float(error);record=dict(step=step,error_2=er,loss=value,seconds=time.time()-began)
        if er<best:
            best=er
            save_fit(er,u,w)
        J=jacobian();cross=2*(J@u.T)*(u@u.T)[:,None,:]
        flat=cross.reshape(B*T,B);solved=torch.cholesky_solve(flat.T,L).T.reshape(B,T,B)
        ynorm=(K@S*K).sum(1);raw=2*J.square().sum(2)
        diagonal=(4*(raw-(cross*solved).sum(2)).clamp(min=0)*ynorm[:,None]).clamp(min=1e-25)
        score=(g.square()/diagonal).sum(1);ids=score.argsort(descending=True)[:a.rows]
        choices=[]
        for row in ids.cpu().tolist():
            rawH=2*J[row]@J[row].T
            gram=rawH-cross[row]@solved[row].T
            H=4*ynorm[row]*(gram+a.padding*rawH)
            for predicted,move in proposals(H.cpu().numpy(),g[row].cpu().numpy(),max_step=a.max_step):
                choices.append((predicted,row,move))
        def global_proposals(selected):
            index=torch.tensor(selected,device='cuda');r=index//T
            deriv=J.reshape(B*T,n)[index];uv=u[r];ju=deriv@uv.T
            rawgram=2*(deriv@deriv.T)*(uv@uv.T)+2*ju*ju.T
            crossflat=cross.reshape(B*T,B)[index];solvedflat=solved.reshape(B*T,B)[index]
            gram=rawgram-crossflat@solvedflat.T;YY=K@S@K
            H=4*gram*YY[r[:,None],r[None,:]]
            H+=a.padding*torch.diag(4*ynorm[r]*raw.flatten()[index])
            return [(predicted,tuple(selected),move) for predicted,move in
                    proposals(H.cpu().numpy(),g.flatten()[index].cpu().numpy(),max_step=a.max_step)]
        if a.global_dim:
            count=a.global_dim
            energy=(g.square()/diagonal).flatten();ratio=(g.abs()/diagonal).flatten()
            selected=list(dict.fromkeys(energy.argsort(descending=True)[:count//2].cpu().tolist()+
                                         ratio.argsort(descending=True)[:count//2].cpu().tolist()))
            for ix in energy.argsort(descending=True)[:2*count].cpu().tolist():
                if len(selected)>=count:break
                if ix not in selected:selected.append(ix)
            choices.extend(global_proposals(selected))
        if a.paired_rows:
            anchors=score.argsort(descending=True)[:a.paired_rows].cpu().tolist();done=set()
            overlap=(u@u.T).abs()
            for anchor in anchors:
                neighbors=overlap[anchor].argsort(descending=True).cpu().tolist()
                neighbor=next(j for j in neighbors if j!=anchor and tuple(sorted((j,anchor))) not in done)
                done.add(tuple(sorted((neighbor,anchor))))
                selected=list(range(anchor*T,(anchor+1)*T))+list(range(neighbor*T,(neighbor+1)*T))
                choices.extend(global_proposals(selected))
        choices.sort(key=lambda x:x[0]);accepted=False;seen={};trials=0;accepted_moves=0
        original_value=value;changed_total=0;largest_step=0
        for predicted,row,move in choices:
            # Avoid spending all trials on a single atom and near-identical jumps.
            if seen.get(row,0)>=(2 if isinstance(row,int) else 4):continue
            seen[row]=seen.get(row,0)+1;trials+=1;delta=torch.tensor(move,**d)
            if isinstance(row,int):theta[row]+=delta
            else:theta.view(-1)[list(row)]+=delta
            actual_state=state();actual_value=float(actual_state[0])
            if actual_value<value-1e-13:
                accepted=True;accepted_moves+=1;changed=int(np.count_nonzero(move));total_moves+=changed;changed_total+=changed
                largest_step=max(largest_step,int(abs(move).max()));value=actual_value
                actual_error=float(actual_state[1])
                if actual_error<best:
                    best=actual_error;save_fit(best,actual_state[2],actual_state[6])
                record.update(row=row if isinstance(row,int) else 'global',predicted=predicted,actual_gain=original_value-value,
                              changed_angles=changed_total,largest_step=largest_step,final_error_2=actual_error)
            else:
                if isinstance(row,int):theta[row]-=delta
                else:theta.view(-1)[list(row)]-=delta
            if trials>=a.trials or accepted_moves>=a.moves_per_step:break
        record.update(accepted=accepted,accepted_moves=accepted_moves,trials=trials,best=best,total_moves=total_moves,candidate_count=len(choices))
        print(json.dumps(record),flush=True)
        if not accepted:break
with torch.no_grad():
    _,error,u,_,_,_,w=state()
    if float(error)<best:
        best=float(error);save_fit(best,u,w)
print(json.dumps(dict(status='FINISHED',scope='classical integer-angle research, not native verification',best_error_2=best,seconds=time.time()-began)),flush=True)
