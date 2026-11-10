"""Column-gauge basin escapes preserving every physical Gram invariant.

Part of normalization stage, not an orbital-vector fitting method. Perturb
several column pairs together, then perform exact pairwise L1 minimization.
All weights, offsets, and integer words are retained in each checkpoint.
"""
import argparse,json,subprocess,time
from jacobi_l1 import encode as bencode,decode as bdecode
from core import *
p=argparse.ArgumentParser();p.add_argument('--init',required=True);p.add_argument('--out',required=True)
p.add_argument('--rounds',type=int,default=10);p.add_argument('--seed',type=int,default=711);a=p.parse_args()
z,u,w,tree,words=load(a.init);R,B,C=w.shape;off=z['W'];rng=np.random.default_rng(a.seed)
A=np.concatenate([w,(off+w.sum(1))[:,None,:]],axis=1);base=A@A.transpose(0,2,1)
def val(A):return float(.25*np.square(abs(A).sum(1)).sum())
best=val(A);current=best;bestA=A.copy();trajectory=[];start=time.time()
save(a.out,u,w,tree,words,W=off,onebody_angle_words=z['onebody_angle_words'],err2=z['err2'])
for step in range(a.rounds):
    trial=A.copy()
    if step:
        for r in range(R):
            for _ in range(64):
                j,k=rng.choice(C,2,replace=False);angle=rng.normal(scale=[.01,.03,.1,.3][step%4]);co,si=np.cos(angle),np.sin(angle)
                x,y=trial[r,:,j].copy(),trial[r,:,k].copy();trial[r,:,j]=co*x+si*y;trial[r,:,k]=-si*x+co*y
    for r in range(R):
        src=HERE/f'candidates/{Path(a.out).stem}_{step}_{r}_in.bin';dst=src.with_name(src.name.replace('_in','_out'))
        bencode(trial[r],src)
        with (HERE/f'results/{Path(a.out).stem}_{step}_{r}.log').open('w') as log:
            subprocess.run([str(ROOT/'initial_experiments/4/codex/jacobi_l1'),str(src),str(dst),'8','12000',str(a.seed+step*31+r)],stdout=log,check=True)
        trial[r]=bdecode(dst)
    nv=val(trial);delta=float(abs(trial@trial.transpose(0,2,1)-base).max());assert delta<1e-9
    accepted=bool(nv<current or rng.random()<np.exp(min(0,(current-nv)/.5)))
    if accepted:A=trial;current=nv
    if nv<best:
        best=nv;bestA=trial.copy();ww=trial[:,:-1];oo=trial[:,-1]-ww.sum(1)
        save(a.out,u,ww,tree,words,W=oo,onebody_angle_words=z['onebody_angle_words'],err2=z['err2'])
        save(HERE/f'trajectories/{Path(a.out).stem}_{step}.npz',u,ww,tree,words,W=oo,onebody_angle_words=z['onebody_angle_words'],err2=z['err2'])
    row=dict(step=step,lambda_h2=nv,best_lambda_h2=best,accepted=accepted,gram_max_error=delta,seconds=time.time()-start)
    trajectory.append(row);print(json.dumps(row),flush=True);Path(str(a.out)+'.json').write_text(json.dumps(trajectory,indent=2))
    if step%4==3:A=bestA.copy();current=best
