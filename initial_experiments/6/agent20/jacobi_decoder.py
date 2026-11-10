"""In-register Jacobi fractions: f=a/(1+b*x+x^2*g).

One peel consumes two moments but uses only the inverse of a. Coefficients
remain in the syndrome slots. Zero a is armed, so every input restores;
noncompleted zero-pivot charts still require separate decoding analysis.
"""
from contextlib import contextmanager
from pathlib import Path
import json
import blocked_cfraction as cfblock

cf=cfblock.cf


def pair(m,E,math,W,inverse=False):
    a,*rest=W
    if inverse:
        math.inverse([a,*rest]);cf.scale(m,E,rest,a)
    else:
        cf.scale(m,E,rest,a,True);math.inverse([a,*rest])


def extract(m,E,math,W,inverse=False):
    for i in (reversed(range(0,len(W),2)) if inverse else range(0,len(W),2)):
        pair(m,E,math,W[i:],inverse)


def synthesize(m,E,W,inverse=False):
    if len(W)==2:return W[::2],W[1::2]
    a,b=W[:2];P,Q=W[2::2],W[3::2]
    def fold(back=False):
        def shear():
            # Prior Q has implicit constant one; b stays unchanged until
            # all products have used it.
            rt=[]
            for i,(p,q) in enumerate(zip(P,Q)):
                rt+=m.prod_rt(p,b,q)
                if i+1<len(Q):rt += [('L',(x,),0,(y,)) for x,y in zip(Q[i+1],p)]
            E.xor(rt)
        if back:
            cf.scale(m,E,P,a,True)
            for x,y in zip(b,P[0]):E.cx(y,x)
            E.swap_fn(None,list(zip(sum(P,[]),sum(Q,[]))))()
            shear()
        else:
            shear();E.swap_fn(None,list(zip(sum(P,[]),sum(Q,[]))))()
            for x,y in zip(b,P[0]):E.cx(y,x)
            cf.scale(m,E,P,a)
    if inverse:fold(True);synthesize(m,E,W[2:],True)
    else:synthesize(m,E,W[2:]);fold()
    return W[::2],W[1::2]


@contextmanager
def parameters(m,E,math,H):
    a,b,c,d=H
    # A=1+(b+d)x+(bd+c)x^2, D=1+bx, p=ac.
    # P0=a*(b^4+c^2+c*(b+d)^2), P1=(bd+c)*a*(b^3+c*d).
    def prod(left,right,out):return m.prod_rt(out,left,right)
    with m.opened(E,'jacobi_e',12) as e,m.opened(E,'jacobi_h',12) as h,m.opened(E,'jacobi_p',12) as p:
        er=[('L',(x,),0,(y,)) for w in (b,d) for x,y in zip(w,e)];E.mac(er)
        hr=prod(b,d,h)+[('L',(x,),0,(y,)) for x,y in zip(c,h)];E.mac(hr)
        pr=prod(a,c,p);E.mac(pr)
        with m.opened(E,'jacobi_P0',12) as p0,m.opened(E,'jacobi_P1',12) as p1:
            def prepare(clear=False):
                with m.opened(E,'jacobi_aux',12) as aux:
                    rr=[];m.bil_terms(rr,m.KT,math.spec([(1,c,1)]),math.spec([(1,e,2)]),tgt=aux)
                    for co,word,q in [(1,b,4),(1,c,2)]:m.lin_rows_terms(rr,math.rows(co,q),word,aux)
                    E.mac(rr)
                    (E.free if clear else E.mac)(prod(a,aux,p0));E.free(rr)
                    rr=prod(c,d,aux)+math.cube_rts(m.spec_plain(b),lambda bits:[aux[j] for j in bits]);E.mac(rr)
                    with m.opened(E,'jacobi_scaled_aux',12) as scaled:
                        sr=prod(a,aux,scaled);E.mac(sr)
                        (E.free if clear else E.mac)(prod(h,scaled,p1));E.free(sr)
                    E.free(rr)
            prepare();yield [(),e,h],[(),b],p,(p0,p1);prepare(True)
        E.free(pr);E.free(hr);E.free(er)


def block(m,E,math,W,inverse=False):
    H,R=W[:4],W[4:]
    def head(back=False):
        for i in ((2,0) if back else (0,2)):pair(m,E,math,H[i:],back)
    if not inverse:head()
    if R:
        with parameters(m,E,math,H) as (A,D,p,P):
            def shift():
                for src,dst in zip(P,R):
                    for x,y in zip(src,dst):E.cx(x,y)
            def reciprocal():
                math.known_zero.add(())
                math.inverse([(),(),*R])
                math.known_zero.remove(())
            if inverse:
                math.triangular(A,1,R,True);math.triangular(D,1,R)
                reciprocal();cf.scale(m,E,R,p);math.triangular(D,1,R,True)
                shift();math.triangular(A,1,R,True)
            else:
                math.triangular(A,1,R);shift();math.triangular(D,1,R)
                cf.scale(m,E,R,p,True);reciprocal()
                math.triangular(D,1,R,True);math.triangular(A,1,R)
    if inverse:head(True)


def blocked_extract(m,E,math,W,inverse=False):
    for i in (reversed(range(0,len(W),4)) if inverse else range(0,len(W),4)):
        if len(W)-i==2:pair(m,E,math,W[i:],inverse)
        else:block(m,E,math,W[i:],inverse)


def profile(blocked=False,N=70):
    cfblock.extract=blocked_extract if blocked else extract
    cfblock.cf.synthesize=synthesize
    return cfblock.profile(N)


if __name__=='__main__':
    reports=[]
    for blocked in (False,True):
        report=profile(blocked);report['blocked']=blocked;reports.append(report);print(json.dumps(report),flush=True)
    Path('results/jacobi_profiles.json').write_text(json.dumps(reports,indent=2)+'\n')
