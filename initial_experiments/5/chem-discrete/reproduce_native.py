"""Final exact native-target reproduction of the completed qfit word controls."""
import json,hashlib
from core import *
tags=['incumbent_control','discrete_control','legacy_exact_split_B479','discrete_claude_B479','qd_b447_C175','hardfreeze_B447','anneal_B671','pruned_B639','constrained_B711']
records=[];m=None
for tag in tags:
    source=ROOT/f'initial_experiments/5/chem-qfit/candidates/{tag}.npz';z,u,w,tree,words=load(source)
    if m is None:m=Metric(tree)
    actual=m.decode(m.tensor(words.reshape(*u.shape[:2],75)));err=float(m.loss(actual,m.tensor(w)).sqrt())
    row=dict(tag=tag,source=str(source),shape=list(w.shape),brot=8,native_pair_error_2=err,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
    records.append(row);print(json.dumps(row),flush=True)
(HERE/'results/native_reproduction.json').write_text(json.dumps(records,indent=2))
