"""Larger Jacobi blocks with short-series preparation of the affine shift.

The transfer polynomials A and B describe the short-head Padé fraction B/A.
Compute its first 2b Taylor coefficients in scratch, and take the high band
of A times that truncated series to prepare P. This avoids expanding the
large P polynomials into an expensive repeated expression tree.
"""
from contextlib import contextmanager,ExitStack
from functools import lru_cache
from pathlib import Path
import argparse,json
import jacobi_decoder as jacobi
from continuant_polynomials import derive_jacobi,Compiler
from field_dag import TreeEvaluator


@lru_cache(None)
def formulas(pairs):
    polys=derive_jacobi(pairs);c=Compiler(2*pairs)
    roots={name:[c.compile(p) for p in (value if isinstance(value,list) else [value])] for name,value in polys.items() if name in ('A','B','D','p')}
    return c,roots


@contextmanager
def parameters(m,E,math,H):
    pairs=len(H)//2;compiler,roots=formulas(pairs)
    def formula(name,out):
        ev=TreeEvaluator(m,E,math,compiler,H)
        for root,word in zip(roots[name][1:] if name in ('A','D') else roots[name],out):ev.xor(root,word)
    with ExitStack() as stack:
        P=[stack.enter_context(m.opened(E,f'general_jacobi_P{i}',12)) for i in range(pairs)]
        A=[(),*[stack.enter_context(m.opened(E,f'general_jacobi_A{i}',12)) for i in range(1,pairs+1)]]
        formula('A',A[1:])
        def affine(clear=False):
            with ExitStack() as scratch:
                low=[scratch.enter_context(m.opened(E,f'general_jacobi_H{i}',12)) for i in range(pairs)]
                formula('B',low)
                with ExitStack() as high_scratch:
                    high=[high_scratch.enter_context(m.opened(E,f'general_jacobi_H{i}',12)) for i in range(pairs,2*pairs)]
                    math.triangular(A,1,[*low,*high],True)
                    targets=[[] for _ in range(2*pairs-1)]
                    for i,out in enumerate(P):targets[pairs-1+i]=[out]
                    rt=m.polyprod_rt(A[1:],high,targets)
                    (E.free if clear else E.mac)(rt)
                    math.triangular(A,1,[*low,*high])
                formula('B',low)
        affine()
        with ExitStack() as extra:
            D=[(),*[extra.enter_context(m.opened(E,f'general_jacobi_D{i}',12)) for i in range(1,pairs)]]
            p=extra.enter_context(m.opened(E,'general_jacobi_determinant',12))
            formula('D',D[1:]);formula('p',[p])
            yield A,D,p,P
            formula('p',[p]);formula('D',D[1:])
        affine(True);formula('A',A[1:])


def block(m,E,math,W,pairs=3,inverse=False):
    H,R=W[:2*pairs],W[2*pairs:]
    if not inverse:jacobi.extract(m,E,math,H)
    if R:
        with parameters(m,E,math,H) as (A,D,p,P):
            def shift():
                for src,dst in zip(P,R):
                    for x,y in zip(src,dst):E.cx(x,y)
            def reciprocal():
                math.known_zero.add(());math.inverse([(),(),*R]);math.known_zero.remove(())
            if inverse:
                math.triangular(A,1,R,True);math.triangular(D,1,R);reciprocal()
                jacobi.cf.scale(m,E,R,p);math.triangular(D,1,R,True);shift();math.triangular(A,1,R,True)
            else:
                math.triangular(A,1,R);shift();math.triangular(D,1,R)
                jacobi.cf.scale(m,E,R,p,True);reciprocal();math.triangular(D,1,R,True);math.triangular(A,1,R)
    if inverse:jacobi.extract(m,E,math,H,True)


def extract(m,E,math,W,inverse=False,pairs=3):
    for i in (reversed(range(0,len(W),2*pairs)) if inverse else range(0,len(W),2*pairs)):
        count=min(pairs,(len(W)-i)//2)
        if count==1:jacobi.pair(m,E,math,W[i:],inverse)
        else:block(m,E,math,W[i:],count,inverse)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--pairs',type=int,nargs='+',default=[3,4]);args=parser.parse_args();rows=[]
    jacobi.cfblock.cf.synthesize=jacobi.synthesize
    for pairs in args.pairs:
        jacobi.cfblock.extract=lambda m,E,math,W,inverse=False,pairs=pairs:extract(m,E,math,W,inverse,pairs)
        r=jacobi.cfblock.profile();r['pairs']=pairs;rows.append(r);print(json.dumps(r),flush=True)
        Path('results/general_jacobi_profiles.json').write_text(json.dumps(rows,indent=2)+'\n')
