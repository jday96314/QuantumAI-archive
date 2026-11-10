"""Jacobi rational synthesis through unnormalized polynomial quotients.

For R(t)=a0/(t+b0+a1/(t+b1+...)), let k0=1/a0,
ki=1/(ai*kprev), qi=ki*(t+bi). Then R=1/(q0+1/(q1+...)).
Reconstruct consecutive continuants with one degree-one convolution per fold,
rather than scaling an entire polynomial bank at every Jacobi fold.
Zero scalar pivots are armed; generic input correctness and total inverse
are distinct contracts. Completion/gap charts are not handled by this module.
"""
import polynomial_inverse_bilinear as poly
from reciprocal_inplace import apply as invert
from jacobi_decoder import cf


def convert(m,E,W,inverse=False):
    for i in (reversed(range(len(W)//2)) if inverse else range(len(W)//2)):
        a,b=W[2*i:2*i+2];prev=W[2*i-2] if i else None
        def multiply(back=False):
            if prev is None:return
            with cf.armed(m,E,prev),m.inv_word(E,prev) as pi,m.opened(E,'quotient_jacobi_product',12) as temp:
                m.mon_slot(E,a,pi if back else prev,prev if back else pi,temp)
        def reciprocal():
            with cf.armed(m,E,a),m.opened(E,'quotient_jacobi_scale',12) as temp:
                def scale(src,iv):m.mon_slot(E,b,iv,src,temp)
                invert(m,E,a,before_swap=scale)
        if inverse:reciprocal();multiply(True)
        else:multiply();reciprocal()


def layout(W,c):
    F,G=[c],[]
    for a,b in reversed(list(zip(W[::2],W[1::2]))):F,G=[*G,b,a],F
    return F,G


def synthesize(m,E,W,c,inverse=False):
    F,G=layout(W,c)
    if inverse:
        poly.extract(m,E,F,G);convert(m,E,W,True)
    else:
        convert(m,E,W);poly.extract(m,E,F,G,True)
    return F,G


def model(m,words):
    prev=1;qs=[]
    for a,b in zip(words[::2],words[1::2]):
        k=m.ginv(m.gmul(a,prev or 1));qs.append((m.gmul(b,k or 1),k));prev=k
    F,G=[1],[]
    for q0,q1 in reversed(qs):
        out=[0]*(len(F)+1)
        for i,v in enumerate(G):out[i]^=v
        for i,v in enumerate(F):out[i]^=m.gmul(q0,v);out[i+1]^=m.gmul(q1,v)
        F,G=out,F
    return F,G
