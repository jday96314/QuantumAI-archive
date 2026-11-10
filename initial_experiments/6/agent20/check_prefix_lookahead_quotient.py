"""Totalized scalar reference, independent head model and native block checks."""
from pathlib import Path
import argparse,json
import check_blocked_cfraction as check
from check_triple_gap_head import look_model
import prefix_lookahead_quotient_block as block
from compact_blocked_cfraction import Arithmetic
import check_jacobi
from projective_polynomial_inverse import bind
HERE=Path(__file__).resolve().parent
PAIRS=6

def head(m,H,inverse=False):
 H=list(H)
 for k in (reversed(range(0,len(H),4)) if inverse else range(0,len(H),4)):H[k:]=look_model(m,H[k:],inverse)
 return H

def parameters(m,H):
 n=len(H)//2;W=list(H);prev=1;flags=[]
 for k in range(0,n,2):
  a,b,h,e=map(int,W[2*k:2*k+4]);gap=not a;flags.append(gap)
  if gap:
   v=m.gmul(b,prev or 1);q=m.ginv(v);W[2*k:2*k+4]=[q,m.gmul(e,q or 1),0,m.gmul(h,q or 1)];prev=q
  else:
   d=e^b;c=h^m.gmul(b,d);v=m.gmul(a,prev or 1);q=m.ginv(v);v=m.gmul(c,q or 1);r=m.ginv(v);W[2*k:2*k+4]=[q,m.gmul(b,q or 1),r,m.gmul(d,r or 1)];prev=r
 pairs=list(reversed(list(zip(W[::2],W[1::2]))));revflags=[False]*n
 for i,gap in enumerate(flags):
  j=n-2-2*i;revflags[j]=gap
  if gap:pairs[j],pairs[j+1]=pairs[j+1],pairs[j]
 def one(q0,q1,f,g):
  out=[0]*(len(f)+1)
  for i,v in enumerate(g):out[i]=v^m.gmul(q0,f[i])^(m.gmul(q1,f[i-1]) if i else 0)
  out[-2]=m.gmul(q0,f[-1] or 1)^(m.gmul(q1,f[-2]) if len(f)>1 else 0);out[-1]=m.gmul(q1,f[-1] or 1);return out,f
 f,g=[1],[]
 for k in reversed(range(0,n,2)):
  if not revflags[k]:f,g=one(pairs[k+1][1],pairs[k+1][0],f,g);f,g=one(pairs[k][1],pairs[k][0],f,g)
  else:
   q2,q1=pairs[k];q0=pairs[k+1][1];gg=[*f,0];out=[*g,q0,q1,q2];nn=len(gg)
   for i in range(nn-2):
    for j,q in enumerate((q0,q1,q2)):
     if i>=j:out[i]^=m.gmul(q,gg[i-j])
   lead=gg[-2] or 1;out[-3]=m.gmul(q0,lead)^(m.gmul(q1,gg[-3]) if nn>=3 else 0)^(m.gmul(q2,gg[-4]) if nn>=4 else 0);out[-2]=m.gmul(q1,lead)^(m.gmul(q2,gg[-3]) if nn>=3 else 0);out[-1]=m.gmul(q2,lead);f,g=out,gg
 lam=f[-1];d=g[-2] if flags[-1] else g[-1];ff=[m.gmul(v,m.ginv(lam or 1)) for v in f[:-1]];gg=[m.gmul(v,m.ginv(d or 1)) for v in g];D=[1,*reversed(gg[:-2]),0] if flags[-1] else [1,*reversed(gg[:-1])];p=m.ginv(m.gmul(lam,d or 1));return [1,*reversed(ff)],D,p,flags[-1]

def model(m,words,inverse=False):
 n=len(words)-2*PAIRS;H,R=list(words[:2*PAIRS]),list(words[2*PAIRS:]);raw=head(m,H,True) if inverse else H.copy()
 if not inverse:H=head(m,H)
 A,D,p,last=parameters(m,H);p=p or 1;b=PAIRS;P=[m.gmul(A[-1] or 1,raw[b+i]) for i in range(b)]
 for i in range(b):
  for k in range(i+1,b):P[i]^=m.gmul(A[b+i-k],raw[b+k])
 shift=2+last
 if inverse:
  R=check.divide(m,R,A,n);R=check.product(m,D,R,n);R=check.divide(m,R,[1,*([0]*(shift-1)),*R],n);R=[m.gmul(p,v) for v in R];R=check.divide(m,R,D,n)
  for i,v in enumerate(P[:n]):R[i]^=v
  R=check.divide(m,R,A,n)
 else:
  R=check.product(m,A,R,n)
  for i,v in enumerate(P[:n]):R[i]^=v
  R=check.divide(m,check.product(m,A,R,n),[p,*([0]*(shift-1)),*check.product(m,D,R,max(0,n-shift))],n)
 return (raw if inverse else H)+R

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--pairs',type=int,default=6);p.add_argument('--n',type=int,nargs='+',default=[14,34]);a=p.parse_args();PAIRS=a.pairs;check.MODE=f'prefix_lookahead_quotient_b{PAIRS}';check.PAIRS=PAIRS;check.block.block=lambda m,E,math,W,inverse=False:block.block(m,E,math,W,PAIRS,inverse);check.block.CubeArithmetic=Arithmetic;check.block.bind=lambda m:bind(m,limit=1024);check.model=model;check.pair=check_jacobi.pair;rows=[]
 for n in a.n:rows.append(check.check(n));(HERE/'results'/f'prefix_lookahead_quotient_b{PAIRS}_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
