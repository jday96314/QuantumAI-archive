"""Full-domain numeric screen of once-rounded fused curvature/linear terms."""
import json,time
from pathlib import Path
import numpy as np
TABLE=np.load('../../5/inv-census/candidates/margin_curve8_t5.9e-08.npy');N=3133452
with np.load('../../5/inv-phase/reference/exact_target.npz') as d:target=d['target'];tol=d['tol']
rows=[];start=time.time()
for extra in [0,1,2,3]:
 for mode in ['separate','joint']:
  mm=1<<64;mx=0;where=0
  for lo in range(1,N+1,131072):
   xs=np.arange(lo,min(N+1,lo+131072),dtype=np.int64);k=21-np.floor(np.log2(xs)).astype(np.int64);mant=xs<<k;s=k//2;c0,c1,c2=TABLE[((mant>>13)&255)+((k&1)<<8)].T;z=(mant&16383)-8192;u=(z+16)>>5;h=(u*u+32)>>6
   if mode=='separate':p=(c0<<extra)-(((c1*z)<<extra)+(1<<11))//4096+(((c2*h)<<extra)+(1<<6))//128
   else:
    p=(((c0<<12)-(c1*z)+(c2*h<<5))*(1<<extra)+(1<<11))//4096
   phase=p.astype(np.uint64)<<(s+27-extra).astype(np.uint64);err=np.abs((phase-target[lo-1:lo-1+len(xs)]).view(np.int64)).astype(np.uint64);to=tol[lo-1:lo-1+len(xs)];ratio=err/to;j=int(ratio.argmax())
   if float(ratio[j])>mx:mx=float(ratio[j]);where=int(xs[j])
   mm=min(mm,int((to.astype(np.int64)-err.astype(np.int64)).min()))
  r=dict(extra_fraction_bits=extra,mode=mode,status='MODEL_PASS' if mm>0 else 'MODEL_FAIL',worst_ratio=mx,worst_x=where,min_margin_units=mm,seconds=time.time()-start);rows.append(r);print(json.dumps(r),flush=True)
Path('results/fused_round_screen.json').write_text(json.dumps(rows,indent=2)+'\n')
