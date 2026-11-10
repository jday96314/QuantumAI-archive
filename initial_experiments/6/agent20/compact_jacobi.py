"""Jacobi block with its denominator stored in the consumed head slots.

The head change (a,b,c,d)->(a,b,c+bd,b+d) is a pair of shears, so it needs
neither a retained branch bit nor a field division, including at zero heads.
"""
from contextlib import contextmanager
from pathlib import Path
import json
import jacobi_decoder as jacobi
from compact_blocked_cfraction import Arithmetic,XorWord


@contextmanager
def frame(m,E,H):
    _,b,c,d=H
    rt=m.prod_rt(c,b,d);E.xor(rt)
    for x,y in zip(b,d):E.cx(x,y)
    yield
    for x,y in zip(b,d):E.cx(x,y)
    E.xor(rt)


@contextmanager
def parameters(m,E,math,H):
    a,b,h,e=H;d=XorWord(e,b)
    def product(left,right,target):
        rt=[];m.bil_terms(rt,m.KT,math.spec([(1,left,1)]),math.spec([(1,right,1)]),tgt=target);return rt
    with m.opened(E,'compact_jacobi_p',12) as p,m.opened(E,'compact_jacobi_P0',12) as p0,m.opened(E,'compact_jacobi_P1',12) as p1:
        def prepare(clear=False):
            with m.opened(E,'compact_jacobi_original_c',12) as c:
                cr=product(b,e,c)+[('L',(x,),0,(y,)) for x,y in zip(h,c)]
                m.lin_rows_terms(cr,math.rows(1,2),b,c);E.mac(cr)
                (E.free if clear else E.mac)(product(a,c,p))
                with m.opened(E,'compact_jacobi_aux',12) as aux:
                    rr=[];m.bil_terms(rr,m.KT,math.spec([(1,c,1)]),math.spec([(1,e,2)]),tgt=aux)
                    for word,q in [(b,4),(c,2)]:m.lin_rows_terms(rr,math.rows(1,q),word,aux)
                    E.mac(rr);(E.free if clear else E.mac)(product(a,aux,p0));E.free(rr)
                    rr=product(c,d,aux)+math.cube_rts(m.spec_plain(b),lambda bits:[aux[j] for j in bits]);E.mac(rr)
                    with m.opened(E,'compact_jacobi_scaled_aux',12) as scaled:
                        sr=product(a,aux,scaled);E.mac(sr)
                        (E.free if clear else E.mac)(product(h,scaled,p1));E.free(sr)
                    E.free(rr)
                E.free(cr)
        prepare();yield [(),e,h],[(),b],p,(p0,p1);prepare(True)


def block(m,E,math,W,inverse=False):
    H,R=W[:4],W[4:]
    def head(back=False):
        for i in ((2,0) if back else (0,2)):jacobi.pair(m,E,math,H[i:],back)
    if not inverse:head()
    if R:
        with frame(m,E,H),parameters(m,E,math,H) as (A,D,p,P):
            def shift():
                for src,dst in zip(P,R):
                    for x,y in zip(src,dst):E.cx(x,y)
            def reciprocal():
                math.known_zero.add(());math.inverse([(),(),*R]);math.known_zero.remove(())
            if inverse:
                math.triangular(A,1,R,True);math.triangular(D,1,R)
                reciprocal();jacobi.cf.scale(m,E,R,p);math.triangular(D,1,R,True)
                shift();math.triangular(A,1,R,True)
            else:
                math.triangular(A,1,R);shift();math.triangular(D,1,R)
                jacobi.cf.scale(m,E,R,p,True);reciprocal()
                math.triangular(D,1,R,True);math.triangular(A,1,R)
    if inverse:head(True)


def extract(m,E,math,W,inverse=False):
    for i in (reversed(range(0,len(W),4)) if inverse else range(0,len(W),4)):
        if len(W)-i==2:jacobi.pair(m,E,math,W[i:],inverse)
        else:block(m,E,math,W[i:],inverse)


if __name__=='__main__':
    jacobi.cfblock.extract=extract;jacobi.cfblock.CubeArithmetic=Arithmetic
    jacobi.cfblock.cf.synthesize=jacobi.synthesize
    report=jacobi.cfblock.profile()
    Path('results/compact_jacobi_profile.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
