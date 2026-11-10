"""Independent rational-series model and native checks for Jacobi blocks."""
import json
from pathlib import Path
import check_blocked_cfraction as check
import jacobi_decoder as jacobi


def pair(m,words):
    a,*tail=words
    normalized=[1,*[m.gmul(v,m.ginv(a or 1)) for v in tail]]
    return [a,*check.divide(m,[1],normalized,len(words))[1:]]


def parameters(m,H):
    a,b,c,d=H;e=b^d;h=m.gmul(b,d)^c;p=m.gmul(a,c)
    p0=m.gmul(a,m.gpow(b,4)^m.gmul(c,c)^m.gmul(c,m.gmul(e,e)))
    p1=m.gmul(h,m.gmul(a,m.gpow(b,3)^m.gmul(c,d)))
    return [1,e,h],[1,b],p,[p0,p1]


def model(m,words,inverse=False):
    H,R=list(words[:4]),list(words[4:]);n=len(R)
    if not inverse:H=pair(m,H);H=H[:2]+pair(m,H[2:])
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
            denominator=[p,0,*check.product(m,D,R,max(n-2,0))]
            R=check.divide(m,check.product(m,A,R,n),denominator,n)
    if inverse:
        def unpair(words):
            a,*tail=words
            reciprocal=check.divide(m,[1],[1,*tail],len(words))
            return [a,*[m.gmul(v,a or 1) for v in reciprocal[1:]]]
        H=H[:2]+unpair(H[2:]);H=unpair(H)
    return H+R


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--compact',action='store_true');args=p.parse_args()
    check.MODE='jacobi_block';check.block.block=jacobi.block
    if args.compact:
        import compact_jacobi
        check.MODE='compact_jacobi_block';check.block.block=compact_jacobi.block;check.block.CubeArithmetic=compact_jacobi.Arithmetic
    check.model=model;check.pair=pair;check.parameters=parameters
    rows=[]
    for n in (6,14,34):
        rows.append(check.check(n))
        Path(f'results/{check.MODE}_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
