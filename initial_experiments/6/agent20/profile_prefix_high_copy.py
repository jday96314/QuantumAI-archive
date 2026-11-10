"""Fixed early-block resources, with the exact arithmetic count emitter."""
from pathlib import Path
import json
import profile_continuant_frame as profile
import prefix_high_copy_jacobi as prefix
HERE=Path(__file__).resolve().parent
profile.general.block=prefix.block;rows=[]
for b in (3,4,5,6,8,10):
 r=profile.run(70,b,'prior');r['kind']='prefix_inplace';rows.append(r);print(json.dumps(r),flush=True);(HERE/'results/prefix_high_copy_block_profiles.json').write_text(json.dumps(rows,indent=2)+'\n')
