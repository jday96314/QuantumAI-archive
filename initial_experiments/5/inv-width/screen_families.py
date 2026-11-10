"""Tighten independent model precisions; cross-pollinate split phase schedules."""
import json
from pathlib import Path
import numpy as np
from family_seeds import fit,model_check,resource_estimate
rows=[]
for kind,configs in [('direct',[(7,(27,18,8)),(7,(28,19,9)),(8,(27,18,7))]),('serial',[(7,(27,18,8)),(7,(28,19,9)),(8,(27,18,7))]),('linear',[(10,(28,18)),(10,(27,17)),(11,(27,17))])]:
 for frac,widths in configs:
  name=f'{kind}_f{frac}_w'+ '_'.join(map(str,widths));table=fit(frac,widths,kind);check=model_check(table,frac,widths,kind);r=resource_estimate(table,frac,widths,kind,check)
  r['rotations']=(widths[0]+10-3) if kind!='direct' else sum([widths[0]+10-3,19+10,14+10])
  np.save('parameters/'+name+'.npy',table)
  row=dict(name=name,kind=kind,frac=frac,widths=widths,model=check,resources=r);rows.append(row);print(json.dumps(row),flush=True)
  Path('results/family_precision_screen.json').write_text(json.dumps(rows,indent=2)+'\n')
