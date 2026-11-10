"""Native physical evaluation and complete ANGLE construction for research fits.

Build invalid research submissions too, so full native recovery can document
which physical gate rejects them. No validity is inferred from static counts.
Offset derivative/one-body solver attributed to chem-qfit/physical.py (AGY).
"""
import argparse
from dataclasses import replace
from scipy.optimize import minimize
from census_common import *
import physical as prior
from oracles import build
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas
p=argparse.ArgumentParser();p.add_argument('fit');p.add_argument('--tag',required=True);a=p.parse_args()
z,u,w,tree,words=load(a.fit);m=Metric(tree)
u=reconstruct_u_vectors(words/128,tree).reshape(u.shape)
words=encode(u.reshape(-1,76),tree,8,True);u=reconstruct_u_vectors(words/128,tree).reshape(u.shape)
X=reconstruct_x(u,w);xf=X.reshape(-1,76*76);tr=np.trace(X,axis1=-2,axis2=-1).reshape(-1)
h0=m.h0.cpu().numpy();off=np.asarray(z['W'],dtype=float).reshape(-1)
f=lambda o:prior.offset_fun(o,xf,tr,h0,m.constant,w)
res=minimize(f,off,jac=True,method='L-BFGS-B',options=dict(maxiter=500,ftol=1e-14,gtol=1e-9,maxls=40))
off=res.x
# Convex gap minimization is a concrete attainable normalization screen. If
# below the strict native floor, move on a specified line to a strict margin.
if f(off)[0]<3.91515:
    origin=off.copy();direction=-w.sum(1).reshape(-1)-origin
    if f(origin+direction)[0]<3.91515:direction=-origin
    assert f(origin+direction)[0]>3.91515
    lo,hi=0.,1.
    for _ in range(65):
        mid=(lo+hi)/2
        if f(origin+mid*direction)[0]<3.91515:lo=mid
        else:hi=mid
    off=origin+hi*direction
off=off.reshape(w.shape[0],w.shape[-1])
h=build_h1_prime(m.inp.ref['h1_thc'],m.inp.S,m.inp.beta_1,m.inp.ref['electrons'],off,X)
prior.Metric=lambda _:m
one,oc=prior.onebody(h,tree,steps=3000,target=.019)
allwords=np.concatenate([one,words]);base=m.inp.params
params=dict(U=u,W=w,H2_I=off+w.sum(1),H1_I=base['H1_I'],Hsym_lower=base['Hsym_lower'])
inp=replace(m.inp,params=params,w=w,W=off)
mt=evaluate_words(inp,allwords,tree,8)
out=HERE/'candidates'/a.tag;out.mkdir(exist_ok=True)
np.savez(out/'submission.npz',**params,tree_edges=tree,angle_words=allwords)
cs=[]
for k,d in enumerate(batch_deltas(allwords,7)):
    c=build(d);assert c.num_qubits<=128
    c.to_file(str(out/f'angle_k_{k}.kmb'),format='kmb');cs.append(c)
rc=resource_costs(inp,cs,8,2*sum(c.operation_statistics().num_toffoli_like for c in cs))
record=dict(fit=a.fit,metrics=mt,passes_physical=passes_accuracy(mt),resources=rc,
    conditional_objective=mt['lambda_eff']*rc['toffoli']*rc['qubits'],onebody=oc,
    offset=dict(success=bool(res.success),iterations=res.nit,message=str(res.message)),
    files={p.name:sha(p) for p in sorted(out.iterdir())})
write(f'{a.tag}_physical.json',record);print(json.dumps(record),flush=True)
