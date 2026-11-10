"""Check eager padding clearance on all saved full-circuit failure syndromes."""
from pathlib import Path
import json
import numpy as np
from raw_cf_model import extract,synthesize
from quadratic_schur_model import F
rows=[]
for p in Path('results/raw_cf_seedless_failures').glob('*.npz'):
 z=np.load(p)
 for i,s in zip(z['shot_indices'],z['syndromes']):
  h,_=extract(s.tolist());A,B,e,_=synthesize(h);rows.append(dict(shot=int(i),exception=e,degree=len(B)-1,constant_nonzero=bool(B and B[0])))
Path('results/eager_padding_failure_cases.json').write_text(json.dumps(rows,indent=2)+'\n');print(rows)
