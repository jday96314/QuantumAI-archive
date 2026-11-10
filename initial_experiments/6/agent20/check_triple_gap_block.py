"""Independent scalar block/inverse and native phase/scratch checks."""
from pathlib import Path
import argparse,json
import check_blocked_cfraction as check
import triple_gap_block as block
from check_triple_gap_head import model as head_model
from check_continuant_head_frame import model as cofactor_model
from compact_blocked_cfraction import Arithmetic
import check_jacobi
HERE=Path(__file__).resolve().parent

def parameters(m,H):
 if not H[0]:
  z,u,h,e,a,b=H;return [1,e^b,h^m.gmul(b,e),a^m.gmul(b,h)],[1,e,h],m.gmul(u,a or 1)
 if not H[2]:
  a,b,z,u,h,e=H;return [1,b^e,h^m.gmul(b,e),u^m.gmul(b,h)],[1,b,0],m.gmul(a,u or 1)
 info=cofactor_model(m,H);return [1,*info[:3]],[1,*info[3:5]],info[5]

def model(m,words,inverse=False):
 H,R=list(words[:6]),list(words[6:]);raw=head_model(m,H,True) if inverse else H.copy()
 if not inverse:H=head_model(m,H)
 A,D,p=parameters(m,H);p=p or 1;n=len(R);extra=int(bool(H[0]) and not H[2]);shift=2+extra
 P=check.product(m,A[1:],raw[3:],5)[2:]
 if n:
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
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,nargs='+',default=[8,14,34]);a=p.parse_args();check.MODE='triple_gap_block';check.PAIRS=3;check.block.block=block.block;check.block.CubeArithmetic=Arithmetic;check.model=model;check.pair=check_jacobi.pair;rows=[]
 for n in a.n:
  rows.append(check.check(n));(HERE/'results/triple_gap_block_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
