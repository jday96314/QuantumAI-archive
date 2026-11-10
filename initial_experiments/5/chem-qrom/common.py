import hashlib,json
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
PROJECT=HERE.parents[2]
INC=PROJECT/'initial_experiments/4/codex/candidates/norm_cycle2_signed'
LAMBDA=48.59417497035283
BASE_COST=708547388.9553173

def save(name,data):
    path=HERE/'results'/name
    path.write_text(json.dumps(data,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x))+'\n')

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def predicted(batch,rows):
    t=19900+2*sum(r['T'] for r in rows);q=422+8*batch
    valid=max(r['Q'] for r in rows)<=128
    return dict(T=t,Q=q,score=LAMBDA*t*q if valid else None,within_cap=valid,
                improvement=BASE_COST/(LAMBDA*t*q) if valid else None)
