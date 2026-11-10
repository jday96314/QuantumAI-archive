"""Independent mixed-chart model, stressing a zero first syndrome."""
import json
from pathlib import Path
import check_blocked_cfraction as check
import check_jacobi
import lookahead_jacobi as lookahead


def model(m,words,inverse=False):
    H,R=list(words[:4]),list(words[4:]);n=len(R);flag=not H[0]
    if not inverse:
        if flag:
            a,u,v,w=H;e=m.gmul(v,m.ginv(u or 1));h=m.gmul(w,m.ginv(u or 1))^m.gmul(e,e);H=[a,u,h,e]
        else:
            H=check_jacobi.pair(m,H);H=H[:2]+check_jacobi.pair(m,H[2:]);a,b,c,d=H;H=[a,b,c^m.gmul(b,d),b^d]
    a,b,h,e=H;bg=b if flag else 0;c=h^m.gmul(b,e)^m.gmul(b,b);d=e^b
    A=[1,e,h];D=[1,b^bg];p=(m.gmul(a,c)^bg) or 1
    p0=m.gmul(a,m.gpow(b,4)^m.gmul(c,c)^m.gmul(c,m.gmul(e,e)))^m.gmul(bg,m.gpow(e,3))
    p1=m.gmul(h,m.gmul(a,m.gpow(b,3)^m.gmul(c,d))^m.gmul(bg,h^m.gmul(e,e)))
    P=[p0,p1];shift=2+int(flag)
    if n:
        if inverse:
            R=check.divide(m,R,A,n);R=check.product(m,D,R,n)
            R=check.divide(m,R,[1,*([0]*(shift-1)),*R],n);R=[m.gmul(p,v) for v in R]
            R=check.divide(m,R,D,n)
            for i,v in enumerate(P[:n]):R[i]^=v
            R=check.divide(m,R,A,n)
        else:
            R=check.product(m,A,R,n)
            for i,v in enumerate(P[:n]):R[i]^=v
            R=check.divide(m,check.product(m,A,R,n),[p,*([0]*(shift-1)),*check.product(m,D,R,max(n-shift,0))],n)
    if inverse:
        if flag:H=[a,b,m.gmul(b or 1,e),m.gmul(b or 1,h^m.gmul(e,e))]
        else:
            d=e^b;c=h^m.gmul(b,d);H=[a,b,c,d]
            for i in (2,0):
                aa,*tail=H[i:];reciprocal=check.divide(m,[1],[1,*tail],len(tail)+1)
                H[i:]=[aa,*[m.gmul(v,aa or 1) for v in reciprocal[1:]]]
    return H+R


if __name__=='__main__':
    check.MODE='lookahead_jacobi';check.FORCE_ZERO_HEAD=True
    check.block.block=lookahead.block;check.block.CubeArithmetic=lookahead.Arithmetic
    check.model=model;check.pair=check_jacobi.pair
    rows=[]
    for n in (6,14,34):
        row=check.check(n);row['scope']='mixed leading-zero chart; half of all fixtures have a zero first syndrome; native forward/inverse value, phase and scratch; encoded head differs from serial'
        rows.append(row);Path('results/lookahead_jacobi_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
