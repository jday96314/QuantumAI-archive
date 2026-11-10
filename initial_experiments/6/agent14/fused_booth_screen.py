"""Price-free numeric screen for actual streamed partial-product rounding.

Only c0 is refitted, one integer bias per existing prefix leaf. Duplicate
coefficient rows retain identical biases, so the lookup leaf structure stays.
"""
import json,time
from pathlib import Path
import numpy as np
N=3133452;TABLE=np.load('../../5/inv-census/candidates/margin_curve8_t5.9e-08.npy');_,leaf=np.unique(TABLE,axis=0,return_inverse=True);leaves=int(leaf.max())+1
with np.load('../../5/inv-phase/reference/exact_target.npz') as d:target=d['target'];tol=d['tol']

def booth(a,z,L,exponent):
 out=np.zeros_like(z)
 for g in range((L+1)//2):
  bit=lambda i:0 if i<0 else (z>>min(i,L-1))&1
  digit=bit(2*g-1)+bit(2*g)-2*bit(2*g+1)
  for mag in [1,2]:
   e=exponent+2*g+int(mag==2);v=a<<e if e>=0 else ((a>>-e)+((a>>(-e-1))&1));out+=np.where(abs(digit)==mag,np.sign(digit)*v,0)
 return out

def calc(lo,q,extra,cm):
 xs=np.arange(lo,min(N+1,lo+131072),dtype=np.int64);k=21-np.floor(np.log2(xs)).astype(np.int64);mant=xs<<k;s=k//2;addr=((mant>>13)&255)+((k&1)<<8);c0,c1,c2=TABLE[addr].T;z=(mant&16383)^8192;signed=(mant&16383)-8192;u=(signed+(1<<(q-1)))>>q;h=((u*u)+(1<<(15-2*q)))>>(16-2*q)
 linear=booth(c1,z,14,-12+extra)
 if cm=='booth':curve=booth(c2,h,12,-7+extra)
 else:
  curve=np.zeros_like(h)
  for bit in range(11):
   ee=bit-7+extra;v=c2<<ee if ee>=0 else ((c2>>-ee)+((c2>>(-ee-1))&1));curve+=((h>>bit)&1)*v
 p=(c0<<extra)-linear+curve;scale=np.left_shift(np.int64(1),(s+27-extra));phase=p.astype(np.uint64)<<((s+27-extra).astype(np.uint64));err=(phase-target[lo-1:lo-1+len(xs)]).view(np.int64).astype(float);to=tol[lo-1:lo-1+len(xs)].astype(float);return p,addr,scale,err,to,xs

rows=[];start=time.time()
for q in [4,5]:
 for extra in [0,1,2]:
  for cm in ['booth','binary']:
   lower=np.full(leaves,-np.inf);upper=np.full(leaves,np.inf)
   for lo in range(1,N+1,131072):
    p,a,scale,e,t,x=calc(lo,q,extra,cm);which=leaf[a];unit=scale*(1<<extra);np.maximum.at(lower,which,(-t-e)/unit);np.minimum.at(upper,which,(t-e)/unit)
   feasible=np.ceil(lower)<=np.floor(upper);bias=np.clip(np.rint((lower+upper)/2),np.ceil(lower),np.floor(upper)).astype(np.int64);mx=0;mm=1<<64;where=0
   for lo in range(1,N+1,131072):
    p,a,scale,e,t,x=calc(lo,q,extra,cm);e=e+bias[leaf[a]]*scale*(1<<extra);ratio=np.abs(e)/t;j=int(ratio.argmax())
    if float(ratio[j])>mx:mx=float(ratio[j]);where=int(x[j])
    mm=min(mm,int((t-np.abs(e)).min()))
   r=dict(coordinate_round_bits=q,extra=extra,curvature_rows=cm,status='MODEL_PASS' if mm>0 else 'MODEL_FAIL',worst_ratio=mx,worst_x=where,min_margin_units=mm,infeasible_leaf_count=int((~feasible).sum()),leaf_count=leaves,bias_min=int(bias.min()),bias_max=int(bias.max()),seconds=time.time()-start);rows.append(r);print(json.dumps(r),flush=True)
   if mm>0:
    tab=TABLE.copy();tab[:,0]+=bias[leaf];np.save(f'parameters/fused_q{q}_e{extra}_{cm}.npy',tab)
Path('results/fused_booth_screen.json').write_text(json.dumps(rows,indent=2)+'\n')
