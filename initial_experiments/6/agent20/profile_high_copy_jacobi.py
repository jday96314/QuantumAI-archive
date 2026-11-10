"""Exact count-emitter profiles compared with unchanged general blocks."""
import json
from pathlib import Path
import profile_continuant_frame as profile
import high_copy_jacobi as high
HERE=Path(__file__).resolve().parent
old=profile.general.block;rows=[]
for b in (2,3,4,5,6):
 for kind in ('prior','high_copy'):
  profile.general.block=old if kind=='prior' else high.block
  r=profile.run(70,b,'prior');r['kind']=kind;rows.append(r);print(json.dumps(r),flush=True);(HERE/'results/high_copy_jacobi_block_profiles.json').write_text(json.dumps(rows,indent=2)+'\n')
