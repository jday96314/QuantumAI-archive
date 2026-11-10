"""Native physical preflight, optional exact offset refinement, signed ANGLE.

The fast offset objective/gradient follows AGY as generalized in qfit physical.
SLSQP here minimizes lambda_eff squared subject to the strict-gap margin.
Only a separate locked full verifier can accept the resulting directory.
"""
import argparse,time,json
from dataclasses import replace
from scipy.linalg import eigh
from scipy.optimize import minimize
from audit_core import *
import physical as prior_physical
# Keep the reused one-body solver on the audit's selected device and exact
# native target. This changes only this process's helper binding, no source.
prior_physical.Metric=Metric
onebody=prior_physical.onebody
from oracles import build
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas
from dfthcsos.verify import resource_costs

p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--tag',required=True)
p.add_argument('--offset-steps',type=int,default=0);p.add_argument('--build',action='store_true');a=p.parse_args()
torch.set_num_threads(2)
z,u,w,tree,words=load(a.fit);m=Metric(tree);base=dict(np.load(INC/'submission.npz'))
u=np.asarray(u,dtype=np.float64);w=np.asarray(w,dtype=np.float64)
if words is None:words=encode(u.reshape(-1,76),tree,8,True)
u=reconstruct_u_vectors(words/128,tree).reshape(u.shape)
words=encode(u.reshape(-1,76),tree,8,True)
canonical_u=reconstruct_u_vectors(words/128,tree).reshape(u.shape)
assert np.max(np.abs(np.einsum('rbp,rbq->rbpq',u,u)-np.einsum('rbp,rbq->rbpq',canonical_u,canonical_u)))<1e-12
u=canonical_u
X=reconstruct_x(u,w);xf=X.reshape(-1,76*76);tr=np.trace(X,axis1=-2,axis2=-1).reshape(-1)
h0=m.h0.cpu().numpy();l1=np.abs(w).sum(1).reshape(-1);t=w.sum(1).reshape(-1)
off=np.asarray(z['W'],dtype=np.float64).reshape(-1).copy()

def objective(o):
    h=h0-(o@xf).reshape(76,76);ev,V=eigh((h+h.T)/2)
    norm=np.abs(ev).sum();g=m.constant+norm-ev.sum()+.5*o@o
    l=norm+.25*np.square(l1+np.abs(o+t)).sum()
    ax=xf@((V*np.sign(ev))@V.T).reshape(-1)
    gg=o+tr-ax;lg=-ax+.5*(l1+np.abs(o+t))*np.sign(o+t)
    return g,l,g*(2*l-g),gg,2*(l-g)*gg+2*g*lg

before=objective(off);offset_record=dict(before=dict(E_gap=before[0],Lambda=before[1],lambda_eff=np.sqrt(max(0,before[2]))))
if before[0]<=3.914898:
    repaired=m.tensor(off.reshape(w.shape[0],w.shape[-1]))
    m.project_gap(m.tensor(u),m.tensor(w),repaired,3.91515)
    off=repaired.cpu().numpy().reshape(-1);before=objective(off)
    offset_record['gap_repair']=dict(E_gap=before[0],Lambda=before[1],lambda_eff=np.sqrt(max(0,before[2])))
checks=[];rng=np.random.default_rng(836121)
# The native objective has genuine |w_B| kinks. Check differentiable points
# away from those kinks rather than falsely rejecting a valid subgradient.
probe=off+rng.normal(scale=.001,size=len(off));probe_values=objective(probe)
for _ in range(4):
    v=rng.normal(size=len(off));v/=np.linalg.norm(v);eps=1e-6
    plus=objective(probe+eps*v);minus=objective(probe-eps*v)
    gap_fd=(plus[0]-minus[0])/(2*eps);product_fd=(plus[2]-minus[2])/(2*eps)
    checks.append(dict(gap_difference=abs(gap_fd-probe_values[3]@v),product_difference=abs(product_fd-probe_values[4]@v)))
assert max(x['gap_difference'] for x in checks)<1e-5
assert max(x['product_difference'] for x in checks)<.02
offset_record['directional_gradient_checks']=checks
if a.offset_steps:
    fun=lambda o:(objective(o)[2]/3000,objective(o)[4]/3000)
    con=dict(type='ineq',fun=lambda o:objective(o)[0]-3.91515,jac=lambda o:objective(o)[3])
    start=time.time();res=minimize(fun,off,jac=True,method='SLSQP',constraints=[con],options=dict(maxiter=a.offset_steps,ftol=1e-12))
    after=objective(res.x)
    offset_record.update(success=bool(res.success),message=res.message,iterations=res.nit,seconds=time.time()-start,
        proposed=dict(E_gap=after[0],Lambda=after[1],lambda_eff=np.sqrt(max(0,after[2]))))
    if after[0]>=3.91514 and after[2]>0 and after[2]<before[2]:off=res.x
off=off.reshape(w.shape[0],w.shape[-1])
ut=m.tensor(u);wt=m.tensor(w);ot=m.tensor(off)
cp,_=checkpoint(m,ut,wt,ot,tree,2)
h=build_h1_prime(m.inp.ref['h1_thc'],m.inp.S,m.inp.beta_1,m.inp.ref['electrons'],off,X)
one,oc=onebody(h,tree,steps=2000,target=.012)
allwords=np.concatenate([one,words]);params=dict(U=u,W=w,H2_I=off+w.sum(1),H1_I=base['H1_I'],Hsym_lower=base['Hsym_lower'])
inp=replace(m.inp,params=params,w=w,W=off)
metrics=evaluate_words(inp,allwords,tree,8)
record=dict(arguments=vars(a),checkpoint=cp,offset=offset_record,onebody=oc,metrics=metrics,physical_pass=passes_accuracy(metrics))
if a.build and record['physical_pass']:
    out=HERE/'candidates'/a.tag;out.mkdir(exist_ok=True)
    np.savez(out/'submission.npz',**params,tree_edges=tree,angle_words=allwords)
    circuits=[]
    for k,delta in enumerate(batch_deltas(allwords,7)):
        circuit=build(delta);assert circuit.num_qubits<=128
        circuit.to_file(str(out/f'angle_k_{k}.kmb'),format='kmb');circuits.append(circuit)
    T=2*sum(c.operation_statistics().num_toffoli_like for c in circuits)
    costs=resource_costs(inp,circuits,8,T)
    record.update(resources_static=costs,estimated_score=metrics['lambda_eff']*costs['toffoli']*costs['qubits'],circuit_directory=str(out))
write_json(HERE/f'results/{a.tag}_physical.json',record);print(json.dumps(record),flush=True)
