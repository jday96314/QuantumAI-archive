"""Gauge-only L1 Jacobi sweep, reusing round-4 Codex's executable."""
import argparse,subprocess,json
from audit_core import *
from jacobi_l1 import encode as binary_encode,decode as binary_decode
p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--out',required=True);p.add_argument('--sweeps',type=int,default=30);a=p.parse_args()
z,u,w,tree,words=load(a.fit);u=np.asarray(u,dtype=np.float64);w=np.asarray(w,dtype=np.float64);off=np.asarray(z['W'],dtype=np.float64);m=Metric(tree)
before=checkpoint(m,m.tensor(u),m.tensor(w),m.tensor(off),tree,1)[0]
changes=[]
for r in range(len(w)):
    A=np.concatenate([w[r],(off[r]+w[r].sum(0))[None]],axis=0)
    src=HERE/f'parameters/{Path(a.out).stem}_{r}_input.bin';dst=src.with_name(src.name.replace('input','output'))
    binary_encode(A,src)
    with (HERE/f'results/{Path(a.out).stem}_{r}_jacobi.log').open('w') as log:
        subprocess.run([str(ROOT/'initial_experiments/4/codex/jacobi_l1'),str(src),str(dst),str(a.sweeps),'12000','193145'],stdout=log,check=True)
    B=binary_decode(dst);delta=np.max(np.abs(A@A.T-B@B.T));assert delta<1e-9
    changes.append(dict(gram_max_error=float(delta),lambda_h2_before=float(.25*np.abs(A).sum(0).dot(np.abs(A).sum(0))),lambda_h2_after=float(.25*np.abs(B).sum(0).dot(np.abs(B).sum(0)))))
    w[r]=B[:-1];off[r]=B[-1]-w[r].sum(0)
save(a.out,u,w,tree,words,W=off)
after=checkpoint(m,m.tensor(u),m.tensor(w),m.tensor(off),tree,1)[0]
rec=dict(arguments=vars(a),before=before,after=after,gauge_changes=changes)
write_json(HERE/f'results/{Path(a.out).stem}_gauge.json',rec);print(json.dumps(rec),flush=True)
