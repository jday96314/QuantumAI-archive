import hashlib,json,time
from audit_core import *

start=time.time();base=dict(np.load(INC/'submission.npz'));tree=base['tree_edges'];m=Metric(tree,compare_cached=True)
u=m.decode(m.tensor(base['angle_words'][76:]).reshape(1,719,75));w=m.tensor(m.inp.w);off=m.tensor(m.inp.W)
native=evaluate_words(m.inp,base['angle_words'],tree,8)
local=checkpoint(m,u,w,off,tree)[0]
assert abs(local['float']['error_2']-native['error_2'])<1e-9
g,l,s=m.physical(u,w,off)
assert abs(float(g)-native['E_gap'])<1e-8
assert abs(float(l)-native['Lambda'])<1e-8
bound=(3.914898*(-m.constant))**.5
records=[]
for name,path in [('incumbent',INC/'submission.npz'),('discrete',ROOT/'initial_experiments/5/chem-qfit/candidates/discrete_control.npz')]:
    z,uu,ww,tt,words=load(path);ut=m.tensor(uu);wt=m.tensor(ww)
    if words is not None:ut=m.decode(m.tensor(words).reshape(*ut.shape[:-1],75))
    # Closest physical generator gauge to original columns; rotate offsets
    # only if supplied. discrete_control uses the incumbent warm gauge.
    if name=='discrete':
        F=w.transpose(1,2)@m.feats(u);F1=wt.transpose(1,2)@m.feats(ut)
        left,_,right=torch.linalg.svd(F1[0]@F[0].T)
        wt=wt@(left@right)
    oo=off.clone();m.project_gap(ut,wt,oo)
    cp,(uq,wq,oq,aw)=checkpoint(m,ut,wt,oo,tt,2)
    save(HERE/f'parameters/{name}_seed.npz',ut,wt,tt,words,W=oo.cpu().numpy())
    save(HERE/f'parameters/{name}_projected.npz',uq,wq,tt,aw,W=oq.cpu().numpy())
    records.append(dict(name=name,metrics=cp))
    print(json.dumps(records[-1]),flush=True)

# Independent algebraic identity and nonnegative brackets on native factors.
t=m.inp.w.sum(1);l1=np.abs(m.inp.w).sum(1);O=m.inp.W
X=reconstruct_x(u.cpu().numpy(),m.inp.w)
ev=np.linalg.eigvalsh(build_h1_prime(m.inp.ref['h1_thc'],m.inp.S,m.inp.beta_1,m.inp.ref['electrons'],O,X))
bracket=.5*(t*t+l1*l1)+O*t+l1*np.abs(O+t)
rhs=-m.constant+2*np.maximum(ev,0).sum()+bracket.sum()
lhs=2*native['Lambda']-native['E_gap']
assert abs(rhs-lhs)<1e-8
rec=dict(native=native,exact_pair=local,cached_target_difference=m.cached_target_difference,
    bound=dict(c=m.constant,lambda_eff_lower_bound=bound,scope='all native factorizations at fixed shifts with the strict positive gap floor',
        assumptions=['normalized unit projector rows','real finite weights and identity offsets','native LCU 1-norm and SOS formulas','native fixed S, beta_1, ecore, electron count and energy_upper','E_gap > 3.914898 and gap in physical domain'],
        identity_lhs=lhs,identity_rhs=rhs,minimum_bracket=float(bracket.min()),
        maximum_normalization_gain=native['lambda_eff']/bound),
    seeds=records,source_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in INC.iterdir() if p.is_file()},seconds=time.time()-start)
write_json(HERE/'results/calibration.json',rec)
print(json.dumps({k:v for k,v in rec.items() if k not in ['source_hashes','seeds']}),flush=True)
