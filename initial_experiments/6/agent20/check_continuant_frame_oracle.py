"""Complete SMALL oracle integration of the compact continuant parameter frame."""
from pathlib import Path
import argparse,json
import check_jacobi_oracle as wrapper
import general_jacobi as general
import continuant_head_frame as compact
import check_general_jacobi as ref
import check_jacobi
from compact_blocked_cfraction import Arithmetic
from check_continuant_head_frame import model as frame_model
HERE=Path(__file__).resolve().parent
PAIRS=2

def model(F,s):
 class Adapter:
  gmul=staticmethod(F.mul);gpow=staticmethod(F.pow);ginv=staticmethod(lambda a:F.inv(a) if a else 0)
 m=Adapter;W=list(s)
 for k in range(0,len(W),2*PAIRS):
  b=min(PAIRS,(len(W)-k)//2)
  if b==1:W[k:]=check_jacobi.pair(m,W[k:]);continue
  H,R=W[k:k+2*b],W[k+2*b:];n=len(R)
  for j in range(0,len(H),2):H[j:]=check_jacobi.pair(m,H[j:])
  info=frame_model(m,H);A=[1,*info[:b]];D=[1,*info[b:2*b-1]];p=info[2*b-1] or 1
  P=ref.parameters(m,H)[3]
  if n:
   R=ref.check.product(m,A,R,n)
   for j,v in enumerate(P[:n]):R[j]^=v
   R=ref.check.divide(m,ref.check.product(m,A,R,n),[p,0,*ref.check.product(m,D,R,max(n-2,0))],n)
  W[k:]=H+R
 def fold(w):
  if len(w)==2:return [w[0]],[w[1]]
  a,b=w[:2];P,Q=fold(w[2:]);nn=len(P)
  return [a,*[F.mul(a or 1,q) for q in Q]],[b^Q[0],*[P[i]^F.mul(b,Q[i])^(Q[i+1] if i+1<nn else 0) for i in range(nn)]]
 P,Q=fold(W)
 from algebra import trim
 return trim([1,*Q]),trim(P)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--pairs',type=int,default=2);p.add_argument('--n',type=int,default=14);p.add_argument('--m',type=int,default=127);p.add_argument('--shots',type=int,default=4096);a=p.parse_args();PAIRS=a.pairs
 general.parameters=compact.parameters
 wrapper.EMIT_EXTRACT=lambda m,E,math,W,inverse=False:general.extract(m,E,math,W,inverse,pairs=PAIRS)
 wrapper.EMIT_ARITHMETIC=Arithmetic
 workspace=HERE/f'continuant_frame_b{PAIRS}'
 for sub in ('parameters','candidates','results'):(workspace/sub).mkdir(parents=True,exist_ok=True)
 cf=wrapper.cf;cf.HERE=workspace;cf.TRANSLATOR=HERE/'../../4/codex/translate_prefix';cf.emit=wrapper.emit;cf.model=model
 cf.test(a.n,a.m,a.shots,'agent20-continuant-frame-small-a')
