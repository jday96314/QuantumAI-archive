"""Polynomial modular inversion by reversing a monic Jacobi continued fraction.

All input/output words are reused. Generic Euclid step, F monic degree n:
 u=lead(G), V=G_low/u, alpha=F[n-1]+V[n-2],
 R=F_low+x*V+alpha*V. Store (u,alpha), recurse on (monic V,R).
Then F is the continuant of alpha and inner u. Reversing alpha and the inner
u leaves F fixed. Replace u0 by 1/(product(u)) to turn G into G^-1 mod F.
Armed zero pivots define a total reversible permutation, not a claim of
modular inversion on unsupported degree-gap charts.
"""
from import_published import load
from low_workspace_inverse import bind
import compact_arming
from jacobi_decoder import cf

compact_arming.install(cf)

def peel(m,E,F,G,inverse=False):
    u=G[-1];alpha=F[-1];V=G[:-1];R=F[:-1]
    def shear():
        for i,(v,r) in enumerate(zip(V,R)):
            E.xor(m.prod_rt(r,alpha,v))
            if i:
                for x,y in zip(V[i-1],r):E.cx(x,y)
    def high():
        if V:
            for x,y in zip(V[-1],alpha):E.cx(x,y)
    if inverse:shear();high();cf.scale(m,E,V,u)
    else:cf.scale(m,E,V,u,True);high();shear()
    return V,R,(u,alpha)

def extract(m,E,F,G,inverse=False):
    if not F:return []
    if inverse:
        tail=extract(m,E,G[:-1],F[:-1],True)
        _,_,head=peel(m,E,F,G,True)
    else:
        f,g,head=peel(m,E,F,G)
        tail=extract(m,E,f,g)
    return [head,*tail]

def apply(m,E,F,G):
    heads=extract(m,E,F,G)
    u=[h[0] for h in heads];alpha=[h[1] for h in heads]
    for i in range(len(alpha)//2):
        for x,y in zip(alpha[i],alpha[-1-i]):E.b.swap(x,y)
    for i in range(1,(len(u)+1)//2):
        for x,y in zip(u[i],u[-i]):E.b.swap(x,y)
    for v in u[1:]:cf.scale(m,E,[u[0]],v)
    with m.inv_word(E,u[0]) as iv:
        for x,y in zip(u[0],iv):E.b.swap(x,y)
    extract(m,E,F,G,True)

def model(m,F,G):
    def div(a,b):return m.gmul(a,m.ginv(b or 1))
    def get(f,g):
        if not f:return []
        u=g[-1];v=[div(a,u) for a in g[:-1]]
        alpha=f[-1]^(v[-1] if v else 0)
        r=[a^m.gmul(alpha,b)^(v[i-1] if i else 0) for i,(a,b) in enumerate(zip(f[:-1],v))]
        return [(u,alpha),*get(v,r)]
    def put(h):
        if not h:return [],[]
        (u,a),*rest=h;v,r=put(rest)
        f=[rr^m.gmul(a,vv)^(v[i-1] if i else 0) for i,(rr,vv) in enumerate(zip(r,v))]
        f.append(a^(v[-1] if v else 0))
        return f,[*[m.gmul(u or 1,vv) for vv in v],u]
    h=get(F,G);us=[x for x,_ in h];al=[a for _,a in h]
    p=us[0]
    for u in us[1:]:p=m.gmul(p,u or 1)
    hh=list(zip([m.ginv(p),*reversed(us[1:])],reversed(al)))
    f,g=put(hh);return f,g,all(us)

