"""Multi-block native-normalization optimization in exact integer dictionaries.

An augmented stage can approach the physical feasible set; --feasible uses
backtracked projected directions and accepts only physical feasible iterates.
Fixed native shifts, exact SOS gap, one-body integer rows and strict gap margin
are included. The checkpoint always retains the best feasible normalization.
"""
import argparse,json,time,subprocess
from scipy.optimize import minimize,linear_sum_assignment
from compact import *
from physical import offset_fun,onebody
from dfthcsos.verify import load_inputs
from dfthcsos.scoreinfra.dfthc_math import reconstruct_x,build_h1_prime
p=argparse.ArgumentParser();p.add_argument('--init',required=True);p.add_argument('--out',required=True)
p.add_argument('--steps',type=int,default=3000);p.add_argument('--lr',type=float,default=1e-5)
p.add_argument('--cap',type=float,default=.01999);p.add_argument('--gap',type=float,default=3.9153)
p.add_argument('--jacobi',type=int,default=6);p.add_argument('--feasible',action='store_true')
p.add_argument('--weight',type=float,default=1.);p.add_argument('--seconds',type=float,default=1500)
a=p.parse_args();z,un,wn,tree,words=load(a.init);R,B,C=wn.shape;m=Metric(tree)
t=m.tensor(words.reshape(R,B,75));u=m.decode(t);un=u.cpu().numpy()
if not (HERE/'candidates/native_constants.npz').exists():
    inp=load_inputs(INC/'submission.npz');base=dict(np.load(INC/'submission.npz'))
    h0=build_h1_prime(inp.ref['h1_thc'],inp.S,inp.beta_1,inp.ref['electrons'],np.zeros_like(inp.W),reconstruct_x(base['U'],inp.w))
    constant=inp.energy_upper-inp.ref['ecore']+inp.ref['electrons']*inp.beta_1
    np.savez(HERE/'candidates/native_constants.npz',h0=h0,constant=constant)
const=np.load(HERE/'candidates/native_constants.npz');h0=const['h0'];constant=float(const['constant'])
X=reconstruct_x(un,wn);xf=X.reshape(-1,76*76);tr=np.trace(X,axis1=-2,axis2=-1).reshape(-1)
fun=lambda off:offset_fun(off,xf,tr,h0,constant,wn)
if 'W' in z and z['W'].shape==(R,C) and np.any(z['W']):oo=z['W'].reshape(-1).copy()
else:oo=minimize(fun,np.zeros(R*C),jac=True,method='L-BFGS-B',options=dict(maxiter=600,ftol=1e-14,gtol=1e-9)).x
if fun(oo)[0]<a.gap:
    origin=oo.copy();direction=-wn.sum(1).reshape(-1)-origin
    if fun(origin+direction)[0]<a.gap:direction=-origin
    lo,hi=0.,1.
    for _ in range(70):
        mid=(lo+hi)/2
        if fun(origin+mid*direction)[0]<a.gap:lo=mid
        else:hi=mid
    oo=origin+hi*direction
oo=oo.reshape(R,C)
if a.jacobi:
    from jacobi_l1 import encode as bencode,decode as bdecode
    for r in range(R):
        mat=np.r_[wn[r],(oo[r]+wn[r].sum(0))[None]]
        src=HERE/f'candidates/{Path(a.out).stem}_jacobi{r}_in.bin';dst=src.with_name(src.name.replace('_in','_out'))
        bencode(mat,src)
        with (HERE/f'results/{Path(a.out).stem}_jacobi{r}.log').open('w') as f:
            subprocess.run([str(ROOT/'initial_experiments/4/codex/jacobi_l1'),str(src),str(dst),str(a.jacobi),'12000','913104'],stdout=f,check=True)
        rot=bdecode(dst);assert abs(rot@rot.T-mat@mat.T).max()<1e-9
        wn[r]=rot[:-1];oo[r]=rot[-1]-wn[r].sum(0)
h=h0-np.einsum('rc,rcpq->pq',oo,reconstruct_x(un,wn))
one,oc=onebody(h,tree,target=.0145);vone=m.decode(m.tensor(one));_,vv=np.linalg.eigh(h)
rr,cc=linear_sum_assignment(-abs(vone.cpu().numpy()@vv));match=np.empty(76,dtype=int);match[rr]=cc;match=torch.tensor(match,device='cuda')
c=Compact(m,u);f=c.factors(m.tensor(wn)).detach().requires_grad_(True)
# w_B parametrization keeps the normalization's identity coefficient explicit.
wb=m.tensor(oo+wn.sum(1)).requires_grad_(True);h0t=m.tensor(h0)
def vals():
    w=c.weights(f);off=wb-w.sum(1);packed=torch.einsum('rbc,rc,rbp->p',w,off,c.P)/m.fac
    hc=torch.zeros((76,76),**m.d);hc[m.ii,m.jj]=packed;hc=hc+hc.T-torch.diag(hc.diag())
    h=h0t-hc;e,evec=torch.linalg.eigh(h);sm=(e*e+1e-16).sqrt()
    gap=constant+sm.sum()-e.sum()+.5*off.square().sum()
    lam=sm.sum()+.25*((w*w+1e-16).sqrt().sum(1)+(wb*wb+1e-16).sqrt()).square().sum()
    ar,ac=linear_sum_assignment(-abs((vone@evec).detach().cpu().numpy()))
    assignment=np.empty(76,dtype=int);assignment[ar]=ac;assignment=torch.tensor(assignment,device='cuda')
    e2=c.loss(f);e1=((vone.T*e[assignment])@vone-h).square().sum()
    score=gap*(2*lam-gap)
    return e2,e1,gap,lam,score,w,off
start=time.time();best=float('inf');besterr=float('inf');traj=[];mue=0.;mug=0.;rho=4.;params=[f,wb]
opt=torch.optim.Adam(params,lr=a.lr)
def record(step):
    global best,besterr
    with torch.no_grad():e2,e1,gap,lam,score,w,off=vals()
    er=float(e2.sqrt());eone=float(e1.sqrt());gv=float(gap);lv=float(lam);sv=float(score.clamp(min=0).sqrt())
    feasible=er<=.019999 and eone<=.019999 and gv>=3.91505 and lv>49 and gv<2*lv
    kw=dict(W=off.detach().cpu().numpy(),onebody_angle_words=one,err2=er)
    if feasible and sv<best:best=sv;save(a.out,u,w,tree,t,**kw)
    if er<besterr:besterr=er;save(str(Path(a.out).with_suffix(''))+'_minerror.npz',u,w,tree,t,**kw)
    row=dict(step=step,error_2=er,error_1=eone,E_gap=gv,Lambda=lv,lambda_eff=sv,feasible=feasible,best_feasible=best,seconds=time.time()-start)
    traj.append(row);print(json.dumps(row),flush=True);Path(str(a.out)+'.json').write_text(json.dumps(traj,indent=2))
    if step%200==0:save(HERE/f'trajectories/{Path(a.out).stem}_{step}.npz',u,w,tree,t,**kw)
    return row
record(-1)
for step in range(a.steps):
    if not a.feasible:
        opt.zero_grad();e2,e1,gap,lam,score,w,off=vals()
        ge=e2/a.cap**2-1;gg=(a.gap-gap)/a.gap
        def pen(g,mu):return (torch.relu(mu+rho*g).square()-mu*mu)/(2*rho)
        loss=a.weight*score/3000+pen(ge,mue)+pen(gg,mug)
        loss.backward();torch.nn.utils.clip_grad_norm_(params,1e6);opt.step()
        if step%20==0:mue=max(0,mue+rho*float(ge.detach()));mug=max(0,mug+rho*float(gg.detach()))
    else:
        e2,e1,gap,lam,score,w,off=vals();objective=score/3000
        grads=torch.autograd.grad(objective,params,retain_graph=True);g=torch.cat([x.flatten() for x in grads]);direction=-g
        constraints=[e2/a.cap**2-1,e1/.01995**2-1,(a.gap-gap)/a.gap]
        cg=[]
        for con in constraints:
            grad=torch.autograd.grad(con,params,retain_graph=True,allow_unused=True)
            cg.append(torch.cat([(x if x is not None else torch.zeros_like(pp)).flatten() for x,pp in zip(grad,params)]))
        # Active-set half-space projection in the whitened factor coordinates.
        for _ in range(8):
            for con,grad in zip(constraints,cg):
                bound=max(0,-float(con.detach()))*.25/max(a.lr,1e-20)
                excess=grad@direction-bound
                if float(excess)>0:direction-=excess/(grad@grad).clamp(min=1e-30)*grad
        old=[p.detach().clone() for p in params];oldscore=float(objective.detach());scale=a.lr;success=False
        for ls in range(24):
            with torch.no_grad():
                k=0
                for pp,base in zip(params,old):pp.copy_(base+scale*direction[k:k+pp.numel()].reshape(pp.shape));k+=pp.numel()
                ne,ne1,ng,nl,ns,_,_=vals()
            if float(ne)<=.019999**2 and float(ne1)<=.019999**2 and float(ng)>=3.91505 and float(nl)>49 and float(ns/3000)<oldscore:
                success=True;break
            scale*=.5
        if not success:
            with torch.no_grad():
                for pp,base in zip(params,old):pp.copy_(base)
    if step%20==0 or step==a.steps-1:record(step)
    if time.time()-start>=a.seconds:break
