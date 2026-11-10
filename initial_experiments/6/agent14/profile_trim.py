import json,time
from pathlib import Path
import profile_schedule as ps
from trimmed_middle import TrimmedArithmetic
from mobius_series import BoundedArithmetic
from dqi.verify_submission_dqi import clopper_pearson_upper_bound,repetition_overhead
# All new artifacts live here; prior experiments are imports only.
Path('results').mkdir(exist_ok=True)
start=time.time();rows=[]
for cls in [BoundedArithmetic,TrimmedArithmetic]:
 ps.BoundedArithmetic=cls
 for k,c in [(1,3),(4,3),(7,3),(10,3),(13,2),(15,2),(17,2)]+[(j,1) for j in range(19,35)]:
  row=ps.stage((k,c));row['arithmetic']=cls.__name__;rows.append(row)
  print(json.dumps(row),flush=True)
 Path('results/trim_stage_profiles.json').write_text(json.dumps(rows,indent=2)+'\n')
a=sum(r['T'] for r in rows if r['arithmetic']=='BoundedArithmetic');b=sum(r['T'] for r in rows if r['arithmetic']=='TrimmedArithmetic')
print(json.dumps(dict(original=a,trimmed=b,saving=a-b,seconds=time.time()-start)),flush=True)
