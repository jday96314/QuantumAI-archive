"""Exact curvature decoding from live C0 low bits, slope high bits, address."""
import sys,json,random,time
from pathlib import Path
sys.path.append(str(Path.cwd().parents[2]/'initial_experiments/5/chem-qrom'))
import numpy as np
import affine_invsqrt as affine
from quadratic import solve
T=affine.table();start=time.time();records=[]
for cut in [14,13,12]:
 width=21+8+17-cut;mapping={};bad=0
 for i,(c0,c1,c2) in enumerate(T):
  zz=(np.arange(8192,dtype=np.int64)+((i&1)<<13))^8192;updated=np.full(len(zz),c1,dtype=np.int64)
  for g in range(7):
   bit=lambda k:0 if k<0 else (zz>>min(k,13))&1
   digit=bit(2*g-1)+bit(2*g)-2*bit(2*g+1)
   for mag in [1,2]:
    ee=-11+2*g+(mag==2);v=c2<<ee if ee>=0 else (c2>>-ee)+((c2>>(-ee-1))&1);updated-=np.where(abs(digit)==mag,np.sign(digit)*v,0)
  for hi in np.unique(updated>>cut):
   x=(int(c0)&((1<<21)-1))|((i>>1)<<21)|(int(hi)<<29);y=int(c2) if i!=69 else int(T[68,2])
   assert x not in mapping or mapping[x]==y
   mapping[x]=y
 bank=list(mapping);targets=[mapping[x] for x in bank];d=7;n=width+1;aug=[x|(1<<width) for x in bank];rng=random.Random(142143+cut);found=None
 for count in range(30,81):
  for trial in range(8):
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
    masks=[(sol>>(j*n))&((1<<n)-1) for j in range(d+count)];found=dict(width=width,cut=cut,linear=masks[:d],left=left,right=masks[d:],output=output,products=count,trial=trial);break
  if found:break
 if found:Path(f'parameters/curve_live_cut{cut}.json').write_text(json.dumps(found,indent=2)+'\n')
 row=dict(cut=cut,width=width,rows=len(bank),products=found['products'] if found else None,seconds=time.time()-start);records.append(row);print(json.dumps(row),flush=True)
Path('results/curve_live_fit.json').write_text(json.dumps(records,indent=2)+'\n')
