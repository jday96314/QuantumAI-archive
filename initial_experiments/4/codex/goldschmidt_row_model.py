"""Seed/residual family using the actual signed radix-four row rounding.

Both products round individual Booth rows. The small quadratic correction
uses an independently rounded high part of R, making a small integer square
sufficient. This is an exhaustive numerical screen, not a quantum circuit.
"""
import argparse,json
from pathlib import Path
import numpy as np
from goldschmidt_model import table,rounded


def booth(coef,z,zbits,shift):
 out=np.zeros_like(z)
 for g in range((zbits+1)//2):
  bit=lambda k:0 if k<0 else (z>>min(k,zbits-1))&1
  digit=bit(2*g-1)+bit(2*g)-2*bit(2*g+1)
  for mag in [1,2]:
   row=rounded(coef,shift-2*g-(mag==2))
   out+=np.where(abs(digit)==mag,np.sign(digit)*row,0)
 return out


def monomial_square(t,bits,shift):
 out=np.zeros_like(t)
 for i in range(bits):
  out+=((t>>i)&1)*int(rounded(np.int64(3),shift-2*i))
  for j in range(i+1,bits):
   value=int(rounded(np.int64(3),shift-i-j-1))
   out+=((t>>i)&(t>>j)&1)*(-value if j==bits-1 else value)
 return out


def model(xs,frac,seedbits,resbits,outbits,squarebits,monomial=False,floor=False):
 tab=table(frac,seedbits,resbits)
 k=21-np.floor(np.log2(xs)).astype(np.int64);x=xs<<k;L=21-frac
 q=(x>>L)&((1<<frac)-1);p=k&1;Y,V,R0=tab[q+(p<<frac)].T
 z=(x&((1<<L)-1))-(1<<(L-1))
 R=R0-booth(V,z,L,2*seedbits+21-resbits)
 t=R>>(resbits-squarebits) if floor else rounded(R,resbits-squarebits)
 if monomial:
  tb=squarebits-frac+1
  assert t.min()>=-(1<<(tb-1)) and t.max()<(1<<(tb-1))
  S=R+monomial_square(t,tb,2+2*squarebits-resbits)
 else:S=R+rounded(3*t*t,2+2*squarebits-resbits)
 sbits=resbits-frac+1
 assert S.min()>=-(1<<(sbits-1)) and S.max()<(1<<(sbits-1))
 corr=booth(Y,S,sbits,seedbits+resbits+1-outbits)
 y=(Y<<(outbits-seedbits))+corr
 return y<<(k//2),R,S,t


def main():
 p=argparse.ArgumentParser();p.add_argument('--frac',type=int,default=7);p.add_argument('--seedbits',type=int,default=10);p.add_argument('--resbits',type=int,default=28);p.add_argument('--outbits',type=int,default=28);p.add_argument('--squarebits',type=int,default=18);p.add_argument('--monomial',action='store_true');p.add_argument('--floor',action='store_true');a=p.parse_args()
 worst=0.;where=0;rmin=0;rmax=0;smin=0;smax=0;tmin=0;tmax=0
 for start in range(1,3133453,131072):
  xs=np.arange(start,min(3133453,start+131072));y,R,S,t=model(xs,**vars(a));error=y.astype(float)*np.sqrt(xs)/2**(a.outbits+10)-1;j=int(abs(error).argmax())
  if abs(error[j])>worst:worst=float(abs(error[j]));where=int(xs[j])
  rmin=min(rmin,int(R.min()));rmax=max(rmax,int(R.max()));smin=min(smin,int(S.min()));smax=max(smax,int(S.max()));tmin=min(tmin,int(t.min()));tmax=max(tmax,int(t.max()))
 record=dict(status='MODEL_PASS' if worst<2**-24 else 'MODEL_FAIL',scope='all-x Booth-row approximation with a small square; no quantum circuit',parameters=vars(a),worst_relative_error=worst,worst_x=where,residual_range=[rmin,rmax],corrected_residual_range=[smin,smax],square_input_range=[tmin,tmax])
 tag=f'gold_row_f{a.frac}_s{a.seedbits}_r{a.resbits}_o{a.outbits}_sq{a.squarebits}'+('_monomial' if a.monomial else '')+('_floor' if a.floor else '');np.save(f'parameters/{tag}_table.npy',table(a.frac,a.seedbits,a.resbits));Path(f'results/{tag}_model.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record),flush=True)
if __name__=='__main__':main()
