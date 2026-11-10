"""Reuse the consumed head as Möbius parameters during blocked extraction.

The temporary head frame stores u=c0*armed(c1), h=c1*c3,
d=c1+c2, e=c3 (or c2 on the c3=0 chart), plus one zero-chart bit.
Only p,P0,P1 remain outside the head during the long tail transform.
This is exact on generic and completed charts and a permutation everywhere.
"""
from contextlib import contextmanager
import json
from pathlib import Path

import blocked_cfraction as reference


class XorWord(tuple):
    def __new__(cls,a,b):return tuple.__new__(cls,[frozenset((x,y)) if x!=y else frozenset() for x,y in zip(a,b)])


class Arithmetic(reference.CubeArithmetic):
    def spec(self,terms):
        if not any(isinstance(w,XorWord) for _,w,_ in terms):return super().spec(terms)
        out=[set() for _ in range(12)]
        for c,word,q in terms:
            if not c or tuple(word) in self.known_zero:continue
            for bit,row in enumerate(self.rows(int(c),int(q))):
                for j in self.m.bitlist(row,12):
                    out[bit].symmetric_difference_update(word[j] if isinstance(word,XorWord) else (word[j],))
        return list(map(frozenset,out))


@contextmanager
def original_c1(m,E,H,zero):
    _,h,d,e=H
    with m.opened(E,'compact_original_c1',12) as out:
        def terms(iv):
            rt=m.prod_rt(out,h,iv)
            rt += [('P',(zero,),0,tuple((a,b)) if a!=b else (),0,(target,)) for a,b,target in zip(d,e,out)]
            return rt
        with m.inv_word(E,e) as iv:E.mac(terms(iv))
        yield out
        with m.inv_word(E,e) as iv:E.free(terms(iv))


@contextmanager
def frame(m,E,H):
    a,b,c,d=H
    with m.opened(E,'compact_head_zero',1) as flag:
        zero=flag[0]
        with m.eq_test(E,d,0) as test:E.cx(test,zero)
        reference.cf.scale(m,E,[a],b)
        for x,y in zip(b,c):E.cx(x,y)
        reference.cf.scale(m,E,[b],d)
        E.swaps(zero,list(zip(b,d)))
        for x,y in zip(c,d):E.b.ccx(zero,x,y)
        yield
        for x,y in zip(c,d):E.b.ccx(zero,x,y)
        E.swaps(zero,list(zip(b,d)))
        reference.cf.scale(m,E,[b],d,True)
        for x,y in zip(b,c):E.cx(x,y)
        reference.cf.scale(m,E,[a],b,True)
        with m.eq_test(E,d,0) as test:E.cx(test,zero)


@contextmanager
def parameters(m,E,math,H):
    u,h,d,e=H;epsilon=XorWord(d,e)
    def product(a,b,out,qa=1,qb=1):
        rt=[];m.bil_terms(rt,m.KT,math.spec([(1,a,qa)]),math.spec([(1,b,qb)]),tgt=out);return rt
    def linear(a,b):return [('L',(x,),0,(y,)) for x,y in zip(a,b)]
    def qterms(q):return product(d,e,q)+linear(h,q)
    def vw(v,w,p):return product(u,d,v),product(u,d,w,qb=2)+linear(p,w)
    def outputs(v,w,p0,p1):return product(epsilon,w,p0)+product(h,v,p0),product(h,w,p1)
    with m.opened(E,'compact_block_p',12) as p:
        with m.opened(E,'compact_block_q',12) as q:
            qr=qterms(q);E.mac(qr);E.mac(product(u,q,p));E.free(qr)
        with m.opened(E,'compact_block_P0',12) as p0,m.opened(E,'compact_block_P1',12) as p1:
            with m.opened(E,'compact_block_v',12) as v,m.opened(E,'compact_block_w',12) as w:
                vr,wr=vw(v,w,p);E.mac(vr);E.mac(wr)
                r0,r1=outputs(v,w,p0,p1);E.mac(r0);E.mac(r1)
                E.free(wr);E.free(vr)
            yield [(),epsilon,h],[(),d],p,(p0,p1)
            with m.opened(E,'compact_block_v',12) as v,m.opened(E,'compact_block_w',12) as w:
                vr,wr=vw(v,w,p);E.mac(vr);E.mac(wr)
                r0,r1=outputs(v,w,p0,p1);E.free(r1);E.free(r0)
                E.free(wr);E.free(vr)
        with m.opened(E,'compact_block_q',12) as q:
            qr=qterms(q);E.mac(qr);E.free(product(u,q,p));E.free(qr)


def block(m,E,math,W,inverse=False):
    H,R=W[:4],W[4:]
    def head(back=False):
        if back:reference.cf.peel_pair(m,E,math,H[2:],True);reference.cf.peel_pair(m,E,math,H,True)
        else:reference.cf.peel_pair(m,E,math,H);reference.cf.peel_pair(m,E,math,H[2:])
    if not inverse:head()
    if R:
        with frame(m,E,H),parameters(m,E,math,H) as (A,D,p,P):
            def shift():
                for src,dst in zip(P,R):
                    for a,b in zip(src,dst):E.cx(a,b)
            if inverse:
                math.triangular(A,1,R,True);math.triangular(D,1,R)
                math.inverse([(),*R]);reference.cf.scale(m,E,R,p)
                math.triangular(D,1,R,True);shift();math.triangular(A,1,R,True)
            else:
                math.triangular(A,1,R);shift();math.triangular(D,1,R)
                reference.cf.scale(m,E,R,p,True);math.inverse([(),*R])
                math.triangular(D,1,R,True);math.triangular(A,1,R)
    if inverse:head(True)


def extract(m,E,math,W,inverse=False):
    steps=list(range(0,len(W),4))
    for k in (reversed(steps) if inverse else steps):
        if len(W)-k==2:reference.cf.peel_pair(m,E,math,W[k:],inverse)
        else:block(m,E,math,W[k:],inverse)


if __name__=='__main__':
    reference.CubeArithmetic=Arithmetic;reference.extract=extract
    report=reference.profile()
    Path('results/compact_blocked_cfraction_profile.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)
