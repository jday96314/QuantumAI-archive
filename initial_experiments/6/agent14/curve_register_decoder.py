"""Exact quadratic interpolation of curvature from live coefficient/address."""
import sys,json,random,time
from pathlib import Path
sys.path.append(str(Path.cwd().parents[2]/'initial_experiments/5/chem-qrom'))
import numpy as np
import affine_invsqrt as affine
from quadratic import solve
T=affine.table();start=time.time();records=[]
# C0 is unchanged at curvature release and replay. Include the lookup address.
for keep in [0,1,2,3]:
 width=23+9+keep;mask=(1<<keep)-1;bank=[int(a)|((i&511)<<23)|((int(c)&mask)<<32) for i,(a,b,c) in enumerate(T)];targets=[int(c)>>keep for a,b,c in T];mapping={}
 for x,y in zip(bank,targets):
  assert x not in mapping or mapping[x]==y
  mapping[x]=y
 bank=list(mapping);targets=[mapping[x] for x in bank];d=7-keep;n=width+1;aug=[x|(1<<width) for x in bank];rng=random.Random(142141+keep);found=None
 for count in range(15,66):
  for trial in range(4):
   left=[rng.getrandbits(n) for _ in range(count)];output=[rng.randrange(1,1<<d) for _ in range(count)];rows=[];values=[]
   for x,y in zip(aug,targets):
    active=[k for k,l in enumerate(left) if (l&x).bit_count()%2]
    for j in range(d):
     row=x<<(j*n)
     for k in active:
      if output[k]>>j&1:row^=x<<((d+k)*n)
     rows.append(row);values.append(y>>j&1)
   sol,rank=solve(rows,values)
   if sol is not None:
    masks=[(sol>>(j*n))&((1<<n)-1) for j in range(d+count)];found=dict(width=width,keep=keep,linear=masks[:d],left=left,right=masks[d:],output=output,products=count,trial=trial);break
  if found:break
 if found:Path(f'parameters/curve_register_decoder_k{keep}.json').write_text(json.dumps(found,indent=2)+'\n')
 row=dict(keep=keep,width=width,rows=len(bank),products=found['products'] if found else None,seconds=time.time()-start);records.append(row);print(json.dumps(row),flush=True)
Path('results/curve_register_decoder.json').write_text(json.dumps(records,indent=2)+'\n')
