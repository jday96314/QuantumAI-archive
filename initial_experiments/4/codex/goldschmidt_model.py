"""Independent fixed-point seed/residual inverse-square-root approximation.

A short seed determines its exact square. The residual r=1-t*y0^2/a^2 is
small; y0*(1+r/2+3*r^2/8) is a different implementation family from direct
piecewise quadratic interpolation. This file establishes numerical feasibility
only; no quantum resource saving or native acceptance is claimed.
"""
import argparse,json
from pathlib import Path
import numpy as np


def rounded(v,shift):return (v+(1<<(shift-1)))>>shift if shift>0 else v<<(-shift)


def table(frac=8,seedbits=10,resbits=26):
 q=np.tile(np.arange(1<<frac),2);p=np.repeat(np.arange(2),1<<frac);center=1+(q+.5)/2**frac
 seed=np.rint(np.where(p,1.,1/np.sqrt(2))/np.sqrt(center)*2**seedbits).astype(np.int64)
 seed=np.minimum(seed,(1<<seedbits)-1)
 square=seed*seed*(2-p)
 # Exact dyadic residual center, rounded once to its signed accumulator.
 shift=2*seedbits+frac+1-resbits
 R0=(1<<resbits)-rounded(square*((1<<(frac+1))+2*q+1),shift)
 return np.stack([seed,square,R0],1)


def model(xs,tab,frac,seedbits,resbits,outbits,combined=False):
 k=21-np.floor(np.log2(xs)).astype(np.int64);x=xs<<k;L=21-frac;q=(x>>L)&((1<<frac)-1);p=k&1
 Y,V,R0=tab[q+(p<<frac)].T;z=(x&((1<<L)-1))-(1<<(L-1))
 R=R0-rounded(V*z,2*seedbits+21-resbits)
 if combined:
  correction=rounded(Y*(R*(1<<(resbits+2))+3*R*R),seedbits+2*resbits+3-outbits)
 else:correction=rounded(Y*R,seedbits+resbits+1-outbits)+rounded(3*Y*R*R,seedbits+2*resbits+3-outbits)
 y=(Y<<(outbits-seedbits))+correction
 return y<<(k//2),R


def main():
 p=argparse.ArgumentParser();p.add_argument('--frac',type=int,default=8);p.add_argument('--seedbits',type=int,default=10);p.add_argument('--resbits',type=int,default=26);p.add_argument('--outbits',type=int,default=28);p.add_argument('--combined',action='store_true');a=p.parse_args()
 tab=table(a.frac,a.seedbits,a.resbits);worst=0.;where=0;rmin=0;rmax=0
 for start in range(1,3133453,131072):
  xs=np.arange(start,min(3133453,start+131072));y,R=model(xs,tab,a.frac,a.seedbits,a.resbits,a.outbits,a.combined)
  error=y.astype(float)*np.sqrt(xs)/2**(a.outbits+10)-1;j=int(abs(error).argmax())
  if abs(error[j])>worst:worst=float(abs(error[j]));where=int(xs[j])
  rmin=min(rmin,int(R.min()));rmax=max(rmax,int(R.max()))
 record=dict(status='MODEL_PASS' if worst<2**-24 else 'MODEL_FAIL',scope='all-x integer approximation only; no circuit',parameters=vars(a),worst_relative_error=worst,worst_x=where,residual_min=rmin,residual_max=rmax,table_bits=[int(abs(tab[:,i]).max()).bit_length()+(i==2) for i in range(3)])
 tag=f'goldschmidt_f{a.frac}_s{a.seedbits}_r{a.resbits}_o{a.outbits}'+('_combined' if a.combined else '');np.save(f'parameters/{tag}_table.npy',tab);Path(f'results/{tag}_model.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record),flush=True)
if __name__=='__main__':main()
