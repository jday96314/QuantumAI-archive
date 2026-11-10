"""Re-run schedules after the adjacent-product arithmetic change."""
import inspect,json,itertools,time
from pathlib import Path
import affine_invsqrt as affine
import binary_centered_invsqrt as poly
import fixed_top_affine_invsqrt as fixed
import zero_head_affine_invsqrt as zero
from optimize_workspace import compact
source=inspect.getsource(affine.build).replace("    (ROOT/'parameters/affine_build_source.py').write_text(source)","    source=source.replace('from booth_adder import accumulate','from booth_pair_borrow_adder import accumulate')")
source=source.replace('    namespace=dict(p.__dict__)',"    source=source.replace('from booth_adder import accumulate','from booth_pair_borrow_adder import accumulate')\n    namespace=dict(p.__dict__)")
namespace=dict(vars(affine));namespace.update(restore_affine=fixed.restore_affine);exec(source,namespace);poly.cshift=zero.cshift
build=namespace['build'];rows=[];start=time.time();best=175049.66869
with Path('results/pair_booth_frontier.jsonl').open('w') as f:
 for delay,keep,limit,product in itertools.product([2,3,4,5,6],[2,3,4,5],[84,86,88,89,90],[14,16,17,18]):
  try:
   path=Path('candidates/pair_frontier_temp.kmb');c=compact(build(limit=limit,keep=keep,delay=delay,product=product,deferred=3),path);T=c.operation_statistics().num_toffoli_like;Q=c.num_qubits;cost=Q*(T-11+250.85021);r=dict(delay=delay,keep=keep,limit=limit,product=product,Q=Q,T=T,cost_estimate=cost)
   if cost<best:
    best=cost;dest=Path(f'candidates/pair_frontier_q{Q}_t{T}.kmb');c.to_file(str(dest),format='kmb');r['path']=str(dest);print(json.dumps(dict(**r,seconds=time.time()-start)),flush=True);Path('results/pair_booth_frontier_best.json').write_text(json.dumps(r,indent=2)+'\n')
  except (ValueError,AssertionError) as e:r=dict(delay=delay,keep=keep,limit=limit,product=product,error=str(e))
  f.write(json.dumps(r)+'\n');f.flush()
print(json.dumps(dict(status='FINISHED',best_estimate=best,seconds=time.time()-start)),flush=True)
