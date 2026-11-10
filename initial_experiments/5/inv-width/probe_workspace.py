"""Resource screen exact incumbent arithmetic schedules; no model changes."""
import argparse, json
from pathlib import Path
import kickmix as km
import affine_invsqrt as affine
import fixed_top_affine_invsqrt as fixed
import binary_centered_invsqrt as poly
import zero_head_affine_invsqrt as zero
from optimize_workspace import compact

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',default='results/width_probe.json');a=p.parse_args()
 affine.restore_affine=fixed.restore_affine;poly.cshift=zero.cshift
 rows=[]
 for limit in [70,74,78,80,82,84,86,88,90]:
  for keep in [0,2,4]:
   for product in [4,8,12,17]:
    try:
     c=affine.build(limit=limit,keep=keep,delay=2,product=product,deferred=3)
     c=compact(c,Path('candidates/probe_tmp.kmb'));T=c.operation_statistics().num_toffoli_like;cost=c.num_qubits*(T-11+250.850210910699)
     row=dict(limit=limit,keep=keep,product=product,Q=c.num_qubits,static_Toffolis=T,projected_objective=cost)
     if not rows or cost<min(r['projected_objective'] for r in rows):c.to_file('candidates/probe_best.kmb',format='kmb');print(json.dumps(row),flush=True)
     rows.append(row)
    except (ValueError,AssertionError) as e:rows.append(dict(limit=limit,keep=keep,product=product,error=str(e))) if False else None
 valid=sorted(rows,key=lambda r:r['projected_objective']);Path(a.output).write_text(json.dumps(valid,indent=2)+'\n');print(json.dumps(dict(valid=len(valid),best=valid[:3])),flush=True)
