"""Recompute all reported numerical residuals from saved exact words in float64."""
import json
from core import *
paths=['incumbent_control','incumbent_projected','legacy_block_control','discrete_control',
       'legacy_ste_control','pruned_B511','pruned_B639','pruned_B671','pruned_B700',
       'freeze_B479','freeze_B447','freeze_B671','freeze_B700',
       'anglefreeze_B479','anglefreeze_B511','anglefreeze_B639',
       'randomfreeze_B479','hardfreeze_B447','discrete_B700','legacy_exact_B700',
       'anneal_B671','exchange_B700','exchange_B671','exchange_qd_B479',
       'legacy_exact_split_B479','discrete_claude_B479','constrained_B711',
       'legacy_exact_B700_C360','legacy_exact_B700_C365','greedy_B711',
       'ablation_align_B447','ablation_noalign_B447']
rows=[];m=None
for tag in paths:
    path=HERE/f'candidates/{tag}.npz'
    if not path.exists():continue
    z,u,w,tree,words=load(path)
    if m is None:m=Metric(tree)
    if words is not None:
        uq=reconstruct_u_vectors(words/128,tree).reshape(u.shape)
        signs=np.where((uq*u).sum(-1)<0,-1.,1.)
        decode_error=float(np.max(abs(uq-signs[...,None]*u)))
    else:uq=u;decode_error=None
    row=dict(tag=tag,shape=list(w.shape),error_2=float(m.loss(m.tensor(uq),m.tensor(w)).sqrt()),
             decoded_vector_error_up_to_sign=decode_error,saved_error=float(z.get('err2',np.nan)))
    rows.append(row);print(json.dumps(row),flush=True)
(HERE/'results/numerical_summary.json').write_text(json.dumps(rows,indent=2)+'\n')
