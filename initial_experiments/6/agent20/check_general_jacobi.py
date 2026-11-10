"""Independent polynomial transfer model for larger emitted Jacobi blocks."""
import argparse,json
from pathlib import Path
import check_blocked_cfraction as check
import check_jacobi
import general_jacobi as general

PAIRS=3


def fraction(m,H):
    if len(H)==2:return [H[0]],[1,H[1]]
    P,Q=fraction(m,H[2:]);a,b=H[:2]
    out=[m.gmul(a,v) for v in Q]
    denominator=[0]*(max(len(Q)+1,len(P)+2))
    for i,v in enumerate(Q):denominator[i]^=v;denominator[i+1]^=m.gmul(b,v)
    for i,v in enumerate(P):denominator[i+2]^=v
    return out,denominator


def parameters(m,H):
    B,A=fraction(m,H);D=fraction(m,H[:-2])[1];length=len(H)
    series=check.divide(m,B,A,length);p=1
    for a in H[::2]:p=m.gmul(p,a)
    affine=[]
    for j in range(len(H)//2):
        v=0
        for i in range(j+1,len(A)):v^=m.gmul(A[i],series[length+j-i])
        affine.append(v)
    return A,D,p,affine


def model(m,words,inverse=False):
    H,R=list(words[:2*PAIRS]),list(words[2*PAIRS:]);n=len(R)
    if not inverse:
        for i in range(PAIRS):H=H[:2*i]+check_jacobi.pair(m,H[2*i:])
    A,D,p,P=parameters(m,H);p=p or 1
    if n:
        if inverse:
            R=check.divide(m,R,A,n);R=check.product(m,D,R,n)
            R=check.divide(m,R,[1,0,*R],n);R=[m.gmul(p,v) for v in R]
            R=check.divide(m,R,D,n)
            for i,v in enumerate(P[:n]):R[i]^=v
            R=check.divide(m,R,A,n)
        else:
            R=check.product(m,A,R,n)
            for i,v in enumerate(P[:n]):R[i]^=v
            R=check.divide(m,check.product(m,A,R,n),[p,0,*check.product(m,D,R,max(n-2,0))],n)
    if inverse:
        for i in reversed(range(PAIRS)):
            a,*tail=H[2*i:];reciprocal=check.divide(m,[1],[1,*tail],len(tail)+1)
            H[2*i:]=[a,*[m.gmul(v,a or 1) for v in reciprocal[1:]]]
    return H+R


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--pairs',type=int,default=3);p.add_argument('--n',type=int,nargs='+',default=[8,14,34]);args=p.parse_args();PAIRS=args.pairs
    check.MODE=f'general_jacobi_b{PAIRS}';check.PAIRS=PAIRS
    check.block.block=lambda m,E,math,W,inverse=False:general.block(m,E,math,W,PAIRS,inverse)
    check.model=model;check.pair=check_jacobi.pair;check.parameters=parameters
    rows=[]
    for n in args.n:
        rows.append(check.check(n));Path(f'results/{check.MODE}_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
