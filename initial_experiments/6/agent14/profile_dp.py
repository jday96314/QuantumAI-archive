import json,time
from pathlib import Path
import profile_schedule as ps
from triangular_dp import DPArithmetic
from trimmed_middle import TrimmedArithmetic
ps.BoundedArithmetic=DPArithmetic
start=time.time();rows=[]
for k,c in [(1,3),(4,3),(7,3),(10,3),(13,2),(15,2),(17,2)]+[(j,1) for j in range(19,35)]:
 row=ps.stage((k,c));rows.append(row);print(json.dumps(row),flush=True)
 Path('results/dp_stage_profiles.json').write_text(json.dumps(rows,indent=2)+'\n')
print(json.dumps(dict(T=sum(r['T'] for r in rows),seconds=time.time()-start)),flush=True)
