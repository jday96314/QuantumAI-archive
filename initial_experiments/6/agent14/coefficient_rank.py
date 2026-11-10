import json
from pathlib import Path
import numpy as np
P=Path('../../4/codex')
t=np.load(P/'candidates/margin_curve8_t5.9e-08.npy').copy()
i=np.arange(len(t));p=i>>8;coarse=(i&255)>>1
t[:,0]+=(1<<25)+(coarse<<18)+(1<<17)-np.array([961<<17,1217<<17])[p]
columns=[(np.ones(len(t),dtype=int),'one')]+[((i>>j)&1,f'a{j}') for j in range(9)]
for col in range(3):
 for j in range(int(t[:,col].max()).bit_length()):columns.append(((t[:,col]>>j)&1,f'c{col}.{j}'))
basis={};relations=[]
for bits,name in columns:
 v=int.from_bytes(np.packbits(bits,bitorder='little').tobytes(),'little');expr={name}
 while v:
  b=v.bit_length()-1
  if b not in basis:basis[b]=(v,expr);break
  bv,be=basis[b];v^=bv;expr^=be
 if not v:relations.append(sorted(expr))
res=dict(columns=len(columns),rank=len(basis),relations=relations,coefficient_ranges=[(int(x.min()),int(x.max())) for x in t.T])
Path('results/coefficient_rank.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
