"""Native exact totalized block model, including armed diagonal exceptions."""
import argparse,json
from pathlib import Path
import check_blocked_cfraction as check
from check_high_copy_jacobi import head
import check_jacobi
import prefix_high_copy_jacobi as prefix
from compact_blocked_cfraction import Arithmetic
HERE=Path(__file__).resolve().parent
PAIRS=3

def frame_model(m,H):
 b=len(H)//2;p=H[0];P=[];B=[H[1]]
 for beta,alpha in zip(H[2::2],H[3::2]):
  effective=beta or 1;out=[];p=m.gmul(p,effective)
  for i,v in enumerate(P):out.append(m.gmul(effective,v)^m.gmul(alpha,B[i])^(B[i-1] if i else 0))
  out.extend([beta^m.gmul(alpha,B[-1])^(B[-2] if len(B)>1 else 0),alpha^B[-1]]);P,B=B,out
 return [1,*reversed(B)],[1,*reversed(P)],p

def model(m,words,inverse=False):
 b=PAIRS;n=len(words)-2*b;H,R=list(words[:2*b]),list(words[2*b:]);raw=head(m,H,True) if inverse else H.copy()
 if not inverse:H=head(m,H)
 A,D,p=frame_model(m,H);p=p or 1
 P=[m.gmul(A[-1] or 1,raw[b+i]) for i in range(b)]
 for i in range(b):
  for k in range(i+1,b):P[i]^=m.gmul(A[b+i-k],raw[b+k])
 if n:
  if inverse:
   R=check.divide(m,R,A,n);R=check.product(m,D,R,n);R=check.divide(m,R,[1,0,*R],n);R=[m.gmul(p,v) for v in R];R=check.divide(m,R,D,n)
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
 check.MODE=f'prefix_high_copy_b{PAIRS}';check.PAIRS=PAIRS;check.block.block=lambda m,E,math,W,inverse=False:prefix.block(m,E,math,W,PAIRS,inverse);check.block.CubeArithmetic=Arithmetic
 check.model=model;check.pair=check_jacobi.pair;rows=[]
 for n in a.n:
  rows.append(check.check(n));(HERE/'results'/f'prefix_high_copy_b{PAIRS}_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
