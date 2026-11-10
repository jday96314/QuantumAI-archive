import json
from pathlib import Path
import kickmix as km
import affine_invsqrt as affine
import fixed_top_affine_invsqrt as fixed
import binary_centered_invsqrt as poly
import zero_head_affine_invsqrt as zero
from optimize_workspace import compact
affine.restore_affine=fixed.restore_affine;poly.cshift=zero.cshift
rows=[];fail=[]
for limit in [78,80,82,84,86,88,90]:
 for keep in [0,1,2,3,4]:
  for delay in [3,4,5]:
   for product in [12,17,None]:
    try:
     c=compact(affine.build(limit=limit,keep=keep,delay=delay,product=product,deferred=3),Path('candidates/delay_tmp.kmb'));T=c.operation_statistics().num_toffoli_like;cost=c.num_qubits*(T-11+250.850210910699)
     row=dict(limit=limit,keep=keep,delay=delay,product=product,Q=c.num_qubits,T=T,estimate=cost);rows.append(row)
     if cost<min([r['estimate'] for r in rows[:-1]] or [float('inf')]):c.to_file('candidates/delay_best.kmb',format='kmb');print(json.dumps(row),flush=True)
    except (ValueError,AssertionError) as e:fail.append(dict(limit=limit,keep=keep,delay=delay,product=product,error=str(e)))
Path('results/delay_probe.json').write_text(json.dumps(dict(rows=sorted(rows,key=lambda r:r['estimate']),failures=fail),indent=2)+'\n')
print(json.dumps(dict(valid=len(rows),failures=len(fail),best=sorted(rows,key=lambda r:r['estimate'])[:3])),flush=True)
