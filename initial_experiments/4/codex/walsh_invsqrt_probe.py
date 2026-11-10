"""Screen direct Walsh phase synthesis of the normalized mantissa.

This favorable front-end screen omits parity and exponent modulation. It does
not construct a quantum circuit or establish a bound on other approximations.
"""
import argparse,json,time
from pathlib import Path
import numpy as np


def fwht(v):
 out=v.copy();n=len(out);h=1
 while h<n:
  view=out.reshape(-1,2*h);left=view[:,:h].copy();right=view[:,h:].copy();view[:,:h]=left+right;view[:,h:]=left-right;h*=2
 return out


def main():
 p=argparse.ArgumentParser();p.add_argument('--bits',type=int,default=21);a=p.parse_args();begin=time.time();n=1<<a.bits;x=1+np.arange(n)/n;target=1/np.sqrt(x);coef=fwht(target)/n;order=np.argsort(-abs(coef));records=[]
 for count in [32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536]:
  if count>n:continue
  kept=np.zeros(n);kept[order[:count]]=coef[order[:count]];approx=fwht(kept);error=abs(approx/target-1);record=dict(terms=count,rotations=count-1,max_relative_error=float(error.max()),worst_index=int(error.argmax()),omitted_coefficient_l1=float(abs(coef[order[count:]]).sum()));records.append(record);print(json.dumps(record),flush=True)
  if error.max()<2**-24:break
 result=dict(scope='normalized mantissa only, largest Walsh coefficients, floating-point approximation; no quantum resources or native acceptance',bits=a.bits,negative_coefficients=int(np.sum(coef < -1e-15)),records=records,seconds=time.time()-begin);Path(f'results/walsh_mantissa_{a.bits}_probe.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
