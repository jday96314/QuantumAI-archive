"""Convex accelerated projected gradient and exact block coordinate solves.

Rank-C coordinate updates are nonconvex; dual bounds always refer to the
unrestricted block-PSD problem, never certify the rank-constrained optimum.
"""
import argparse,time
from core import *

@torch.no_grad()
def solve(model,S=None,steps=2000,C=None,method='fista',every=100,tol=1e-9,quiet=False):
    R,B=model.R,model.B
    S=project(model.A,C) / R if S is None else project(S,C)
    Y=S.clone();t=1.;best=float(model.objective(S));bestS=S.clone();records=[];start=time.time()
    # ||A* A|| <= R since each Q_r is an isometry. Power iteration sharpens
    # this estimate; retain R as the safe step bound for convex FISTA.
    for k in range(steps):
        if method=='coordinate':
            for r in range(R):
                rhs=model.A[r]-(model.cross[r]@S@model.cross[r].transpose(1,2)).sum(0)+S[r]
                S[r]=project(rhs[None],C)[0]
        else:
            new=project(Y-(model.apply(Y)-model.A)/R,C)
            nt=(1+np.sqrt(1+4*t*t))/2
            # Adaptive gradient restart handles the heavily overlapping cone.
            if float(((Y-new)*(new-S)).sum())>0:nt=1.
            Y=new+((t-1)/nt)*(new-S);S=new;t=nt
        if k%every==0 or k==steps-1:
            val=float(model.objective(S))
            if val<best:best=val;bestS=S.clone()
            rec=dict(step=k,seconds=time.time()-start,**model.diagnostics(bestS));records.append(rec)
            if not quiet:print(json.dumps(rec),flush=True)
            if C is None and rec['primal_dual_gap_squared']<tol:break
    return bestS,records

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--out',required=True)
    p.add_argument('--steps',type=int,default=2000);p.add_argument('--rank',type=int)
    p.add_argument('--method',default='fista',choices=['fista','coordinate']);p.add_argument('--every',type=int,default=100)
    a=p.parse_args();z=dict(np.load(a.input));m=Covariance(z['U']);S=None
    if 'S' in z:S=tensor(z['S'])
    elif 'w' in z:S=m.from_weights(z['w'])
    S,hist=solve(m,S,a.steps,a.rank,a.method,a.every)
    w=m.weights(S,a.rank).cpu().numpy()
    np.savez(a.out,U=z['U'],S=S.cpu().numpy(),w=w,tree_edges=np.load(INC)['tree_edges'])
    report(HERE/f'results/{Path(a.out).stem}.json',dict(input=a.input,rank=a.rank,method=a.method,
           diagnostics=m.diagnostics(S,True,a.rank),history=hist))
if __name__=='__main__':main()
