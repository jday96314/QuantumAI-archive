"""Independent full rational-series block model for the half-head copy."""
import argparse,json
from pathlib import Path
import check_blocked_cfraction as check
import check_jacobi
from check_continuant_head_frame import model as cofactors
import high_copy_jacobi as high
from compact_blocked_cfraction import Arithmetic
HERE=Path(__file__).resolve().parent
PAIRS=3

def head(m,H,inverse=False):
 H=list(H)
 for i in (reversed(range(0,len(H),2)) if inverse else range(0,len(H),2)):
  if inverse:
   a,*tail=H[i:];v=check.divide(m,[1],[1,*tail],len(tail)+1);H[i:]=[a,*[m.gmul(a or 1,x) for x in v[1:]]]
  else:H[i:]=check_jacobi.pair(m,H[i:])
 return H

def model(m,words,inverse=False):
 b=PAIRS;n=len(words)-2*b;H,R=list(words[:2*b]),list(words[2*b:]);raw=head(m,H,True) if inverse else H.copy()
 if not inverse:H=head(m,H)
 info=cofactors(m,H);A=[1,*info[:b]];D=[1,*info[b:2*b-1]];p=info[2*b-1] or 1
 P=check.product(m,A[1:],raw[b:],2*b-1)[b-1:]
 if n:
  if inverse:
   R=check.divide(m,R,A,n);R=check.product(m,D,R,n);R=check.divide(m,R,[1,0,*R],n);R=[m.gmul(p,v) for v in R]
   R=check.divide(m,R,D,n)
   for i,v in enumerate(P[:n]):R[i]^=v
   R=check.divide(m,R,A,n)
  else:
   R=check.product(m,A,R,n)
   for i,v in enumerate(P[:n]):R[i]^=v
   R=check.divide(m,check.product(m,A,R,n),[p,0,*check.product(m,D,R,max(0,n-2))],n)
 if inverse:H=raw
 return H+R

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--pairs',type=int,default=3);p.add_argument('--n',type=int,nargs='+',default=[8,14,34]);a=p.parse_args();PAIRS=a.pairs
 check.MODE=f'high_copy_jacobi_b{PAIRS}';check.PAIRS=PAIRS;check.block.block=lambda m,E,math,W,inverse=False:high.block(m,E,math,W,PAIRS,inverse);check.block.CubeArithmetic=Arithmetic
 check.model=model;check.pair=check_jacobi.pair;rows=[]
 for n in a.n:
  rows.append(check.check(n));(HERE/'results'/f'high_copy_jacobi_b{PAIRS}_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
