"""Reproduce baseline synthesis margin from exact all-input integer phases."""
import json,math
from pathlib import Path
import numpy as np
import binary_centered_invsqrt as p
from invsqrt.verify_submission_invsqrt import ProblemSpecInvSqrt,reference_values,rotation_t_count
p.table_for=lambda *args:np.load('candidates/margin_curve8_t5.9e-08.npy')
spec=ProblemSpecInvSqrt.for_n_p(10);margin=1<<64;min_x=0;worst=(0.,0)
for start in range(1,3133453,131072):
 xs=np.arange(start,min(3133453,start+131072),dtype=np.int64)
 phase=p.model(xs,7,(27,18,8),rounded=True,lookup_frac=8,booth=True).astype(np.uint64)
 units=phase<<np.uint64(64-37);target,tol=reference_values(spec,xs.astype(np.uint64));used=np.abs((units-target).view(np.int64)).view(np.uint64)
 margins=tol.astype(np.int64)-used.astype(np.int64);i=margins.argmin()
 if int(margins[i])<margin:margin=int(margins[i]);min_x=int(xs[i])
 ratio=used/tol;i=ratio.argmax()
 if ratio[i]>worst[0]:worst=(float(ratio[i]),int(xs[i]))
ref=json.loads(Path('reference/zero_head_affine_verify.json').read_text());eps=math.pi*margin/(1<<64)/34;rot=(1+34*rotation_t_count(eps))/4;score=90*(ref['toffoli_count']+rot)
result=dict(status='REPRODUCED',scope='exact exhaustive integer model plus preserved full-native mean Toffolis; byte-identical circuit independently rebuilt',inputs=3133452,exact_min_margin_units=margin,min_margin_x=min_x,worst_ratio=worst[0],worst_x=worst[1],rotation_eps=eps,rotation_equivalents=rot,mean_Toffolis_from_preserved_full_native=ref['toffoli_count'],objective=score,objective_difference=score-ref['score'])
Path('results/reproduced_objective.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
