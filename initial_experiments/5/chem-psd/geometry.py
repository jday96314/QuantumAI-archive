"""Unit-vector variable projection through a full convex PSD inner solve.

Outer L-BFGS improves projector geometry without a shared eigenbasis. Inner
covariance solves replace short alternating factor-gradient fits. The envelope
gradient is exact at a converged inner optimum; recorded inner stationarity
quantifies its approximation during geometry search. Final fixed-dictionary
solves and dual certificates are independent of outer convergence claims.
"""
import argparse,time
from scipy.optimize import minimize
from core import *
from solve import solve

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--tag',required=True)
    p.add_argument('--steps',type=int,default=300);p.add_argument('--inner',type=int,default=80);p.add_argument('--rank',type=int)
    a=p.parse_args();z=dict(np.load(a.input));u0=z['U'];shape=u0.shape;H=tensor(target());start=time.time()
    warm=z.get('w');best=float('inf');records=[];evals=0;latest={}
    def evaluate(raw):
        nonlocal warm,best,evals,latest
        u=raw.reshape(shape);u=u/np.linalg.norm(u,axis=2,keepdims=True)
        m=Covariance(u,H);S0=m.from_weights(warm) if warm is not None else None
        with torch.no_grad():
            S,hist=solve(m,S0,a.inner,a.rank,'fista',a.inner,1e-10,quiet=True)
            w=m.weights(S,a.rank);warm=w.cpu().numpy()
            approx,err=m.direct(S);E=approx-H
            # Differentiate the physical objective at fixed atom covariance.
            # d f / d P = 4 K P E, K P = Rm^-1 S Q.
            KP=torch.linalg.solve_triangular(m.Rm,S@m.Q,upper=True)
            gp=4*(KP@E)
            M=torch.zeros((*gp.shape[:2],76,76),**D)
            ft=tensor(fac);M[:,:,ii,jj]=gp/ft;M[:,:,jj,ii]=gp/ft
            ut=tensor(u);gu=2*(M@ut[:,:,:,None]).squeeze(-1)
            gu-=ut*(gu*ut).sum(2,keepdim=True)
            gu=gu/torch.linalg.vector_norm(tensor(raw.reshape(shape)),dim=2,keepdim=True)
            grad=gu.cpu().numpy();evals+=1
            latest=dict(evaluation=evals,error_2=err,seconds=time.time()-start,
                        weight_squared_norm=float(w.square().sum()),inner_projected_gradient=hist[-1]['projected_gradient_norm'])
            if err<best:
                best=err
                np.savez(HERE/f'parameters/{a.tag}.npz',U=u,w=warm,S=S.cpu().numpy(),tree_edges=z['tree_edges'],error_2=err)
            if evals==1 or evals%10==0:
                records.append(latest);print(json.dumps(latest),flush=True)
            return err*err*1e4,grad.ravel()*1e4
    result=minimize(evaluate,u0.ravel(),jac=True,method='L-BFGS-B',options=dict(maxiter=a.steps,maxls=20,maxcor=20,ftol=1e-13,gtol=1e-8))
    bestz=np.load(HERE/f'parameters/{a.tag}.npz');m=Covariance(bestz['U'])
    S,hist=solve(m,tensor(bestz['S']),2400,a.rank,'fista',200,quiet=True)
    final=m.diagnostics(S,True)
    np.savez(HERE/f'parameters/{a.tag}.npz',U=bestz['U'],w=m.weights(S,a.rank).cpu().numpy(),S=S.cpu().numpy(),tree_edges=z['tree_edges'])
    report(HERE/f'results/{a.tag}.json',dict(input=a.input,outer_method='L-BFGS variable projection',inner_steps=a.inner,
        rank=a.rank,outer_iterations=result.nit,evaluations=evals,termination=str(result.message),best_geometry_error_2=best,
        final_fixed_dictionary=final,final_inner_history=hist,records=records,seconds=time.time()-start))
if __name__=='__main__':main()
