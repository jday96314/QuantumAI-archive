"""Refit streamed fused arithmetic without extra output fraction bits.

Each original prefix leaf gets integer c0/c1/c2 adjustments. The exact
rounded Booth-row model is screened on every original integer input.
"""
import json,time
from pathlib import Path
import numpy as np
TABLE=np.load('../../5/inv-census/candidates/margin_curve8_t5.9e-08.npy');N=3133452;orig,leaf=np.unique(TABLE,axis=0,return_inverse=True)
with np.load('../../5/inv-phase/reference/exact_target.npz') as d:target=d['target'];tol=d['tol'].astype(np.int64)
xs=np.arange(1,N+1,dtype=np.int64);k=21-np.floor(np.log2(xs)).astype(np.int64);mant=xs<<k;addr=((mant>>13)&255)+((k&1)<<8);which=leaf[addr];z=(mant&16383)^8192;u=(((mant&16383)-8192+16)>>5);h=(u*u+32)>>6;s=k//2;scale=np.left_shift(np.int64(1),s+27)

def booth(a,z,L,exponent):
 out=np.zeros_like(z)
 for g in range((L+1)//2):
  bit=lambda i:0 if i<0 else (z>>min(i,L-1))&1
  digit=bit(2*g-1)+bit(2*g)-2*bit(2*g+1)
  for mag in [1,2]:
   e=exponent+2*g+int(mag==2);v=a<<e if e>=0 else ((a>>-e)+((a>>(-e-1))&1));out+=np.where(abs(digit)==mag,np.sign(digit)*v,0)
 return out

start=time.time();rows=[]
for mode in ['booth','binary']:
 result=orig.copy();bad=[];changed=[];margins=[];maxratio=0.;worstx=0
 for L,coef in enumerate(orig):
  ix=np.flatnonzero(which==L);zz=z[ix];hh=h[ix];ss=s[ix];sc=scale[ix];tt=target[ix];to=tol[ix];best=(-float('inf'),None,None);order=sorted([(a,b) for a in range(-5,6) for b in range(-5,6)],key=lambda v:abs(v[0])+abs(v[1]));trials=0
  for ds,dc in order:
   c0,c1,c2=map(int,coef);c1+=ds;c2+=dc
   if not 0<c1<(1<<17) or not 0<c2<(1<<7):continue
   linear=booth(np.full_like(zz,c1),zz,14,-12)
   if mode=='booth':curve=booth(np.full_like(hh,c2),hh,12,-7)
   else:
    curve=np.zeros_like(hh)
    for bit in range(11):
     ee=bit-7;v=c2<<ee if ee>=0 else ((c2>>-ee)+((c2>>(-ee-1))&1));curve+=((hh>>bit)&1)*v
   p=c0-linear+curve;phase=p.astype(np.uint64)<<((ss+27).astype(np.uint64));err=(phase-tt).view(np.int64);lo=int(np.ceil(np.max((-to.astype(float)-err)/sc)));hi=int(np.floor(np.min((to.astype(float)-err)/sc)));trials+=1
   if lo>hi:continue
   for bias in range(lo,hi+1):
    used=err+bias*sc;margin=int((to-np.abs(used)).min())
    if margin>best[0]:best=(margin,[c0+bias,c1,c2],used)
   if best[0]>=40000000:break
  if best[1] is None:bad.append(L);margins.append(-1);continue
  margin,c,err=best;result[L]=c;margins.append(int(margin));ratio=np.abs(err)/to;j=int(ratio.argmax())
  if float(ratio[j])>maxratio:maxratio=float(ratio[j]);worstx=int(xs[ix[j]])
  if c!=coef.tolist():changed.append(dict(leaf=L,old=coef.tolist(),new=c,margin=int(margin),trials=trials))
  if L%32==0:print(json.dumps(dict(mode=mode,leaf=L,failed=len(bad),min_margin=min(margins),seconds=time.time()-start)),flush=True)
 r=dict(mode=mode,status='MODEL_PASS' if not bad and min(margins)>0 else 'MODEL_FAIL',inputs=N,leaf_count=len(orig),failed_leaves=bad,min_margin_units=min(margins),worst_ratio=maxratio,worst_x=worstx,changed_leaf_count=len(changed),changes=changed,seconds=time.time()-start);rows.append(r);Path(f'results/fused_leaf_repair_{mode}.json').write_text(json.dumps(r,indent=2)+'\n');np.save(f'parameters/fused_q5_e0_{mode}_repaired.npy',result[leaf]);print(json.dumps({key:value for key,value in r.items() if key!='changes'}),flush=True)
Path('results/fused_leaf_repair.json').write_text(json.dumps(rows,indent=2)+'\n')
