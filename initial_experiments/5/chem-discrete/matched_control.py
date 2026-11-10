"""Fixed-word weight refit at the same 652 block eigensolves as search."""
import json,time,argparse
from core import *
p=argparse.ArgumentParser();p.add_argument('--prefix',default='');a=p.parse_args()
z,u,w,tree,words=load(HERE/'candidates/refit100_legacy_exact_split_B479.npz');m=Metric(tree);t=m.tensor(words.reshape(*u.shape[:2],75));u=m.decode(t);w=m.tensor(w)
start=time.time();records=[]
for sweep in range(327):
    if sweep%20==0 or sweep==326:
        er=float(m.loss(u,w).sqrt());r=dict(sweep=sweep,PSD_block_sweeps=2*sweep,error_2=er,seconds=time.time()-start)
        records.append(r);print(json.dumps(r),flush=True)
    if sweep<326:w,_=m.project(u,w,1)
save(HERE/f'candidates/{a.prefix}matched_fixed_B479.npz',u,w,tree,t,err2=er)
(HERE/f'results/{a.prefix}matched_fixed_B479.json').write_text(json.dumps(records,indent=2))
