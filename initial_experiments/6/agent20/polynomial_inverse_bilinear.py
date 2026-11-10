"""Unnormalized quotient reversal with a bilinear tail middle product.

Unlike the monic prototype, quotient extraction never normalizes a long
polynomial. Both coefficients of q(x) are quantum data, so the whole tail
q*G is one rank-(n+1) convolution instead of two scalar products per word.
For F=c K(q), G=c K(tail q), reversal gives G*Grev=c^2 mod F.
Scale reversed q_i by c^(2*(-1)^i); invert c iff degree is odd. This keeps
F fixed and divides Grev by c^2. Zero pivots are armed, with no correctness
claim for degree gaps. The totalized transformation is still an involution.
"""
import polynomial_inverse_reversal as monic

load=monic.load;bind=monic.bind;cf=monic.cf

def peel(m,E,F,G,inverse=False):
    q1,q0=F[-1],F[-2];R=F[:-2];n=len(G)
    def quotient(back=False):
        with cf.armed(m,E,G[-1]),m.inv_word(E,G[-1]) as gi,m.opened(E,'poly_quotient_product',12) as t:
            scale=lambda w:m.mon_slot(E,w,G[-1] if back else gi,gi if back else G[-1],t)
            scale(q0 if back else q1)
            if n>1:E.xor(m.prod_rt(q0,q1,G[-2]))
            scale(q1 if back else q0)
    def remainder():
        if not R:return
        # lead(G) cannot contribute to coefficients below n-1.
        # Trimming it removes one evaluation product at every nonterminal step.
        targets=[[r] for r in R]+[[]]
        E.xor(m.polyprod_rt([q0,q1],G[:-1],targets))
    if inverse:remainder();quotient(True)
    else:quotient();remainder()
    return G,R,(q0,q1)

def extract(m,E,F,G,inverse=False):
    if not G:return [],F[0]
    if inverse:
        tail,c=extract(m,E,G,F[:-2],True)
        _,_,head=peel(m,E,F,G,True)
    else:
        f,g,head=peel(m,E,F,G)
        tail,c=extract(m,E,f,g)
    return [head,*tail],c

def apply(m,E,F,G):
    heads,c=extract(m,E,F,G)
    for i in range(len(heads)//2):
        for a,b in zip(heads[i],heads[-1-i]):
            for x,y in zip(a,b):E.b.swap(x,y)
    with m.opened(E,'quotient_reversal_c_squared',12) as c2:
        rt=[];m.lin_rows_terms(rt,tuple(sum(((m.gmul(1<<j,1<<j)>>i)&1)<<j for j in range(12)) for i in range(12)),c,c2)
        E.mac(rt)
        cf.scale(m,E,[w for pair in heads[::2] for w in pair],c2)
        cf.scale(m,E,[w for pair in heads[1::2] for w in pair],c2,True)
        E.free(rt)
    if len(heads)%2:
        with m.inv_word(E,c) as ci:
            for x,y in zip(c,ci):E.b.swap(x,y)
    extract(m,E,F,G,True)

def model(m,F,G,degree=None):
    def div(a,b):return m.gmul(a,m.ginv(b or 1))
    def get(f,g):
        if not g:return [],f[0],True
        q1=div(f[-1],g[-1]);q0=div(f[-2]^(m.gmul(q1,g[-2]) if len(g)>1 else 0),g[-1])
        r=[v^m.gmul(q0,g[i])^(m.gmul(q1,g[i-1]) if i else 0) for i,v in enumerate(f[:-2])]
        rest,c,ok=get(g,r);return [(q0,q1),*rest],c,ok and bool(g[-1])
    def put(h,c):
        if not h:return [c],[]
        (q0,q1),*rest=h;g,r=put(rest,c)
        f=[v^m.gmul(q0,g[i])^(m.gmul(q1,g[i-1]) if i else 0) for i,v in enumerate(r)]
        f.append(m.gmul(q0,g[-1] or 1)^(m.gmul(q1,g[-2]) if len(g)>1 else 0))
        f.append(m.gmul(q1,g[-1] or 1))
        return f,g
    h,c,ok=get(F,G)
    active=len(h) if degree is None else degree
    scalar=h[active][1] if active<len(h) else c
    c2=m.gmul(scalar,scalar) or 1
    hh=[tuple(m.gmul(v,m.ginv(c2) if i%2 else c2) for v in q) for i,q in enumerate(reversed(h[:active]))]
    hh+=h[active:]
    if active%2:
        if active<len(h):hh[active]=(hh[active][0],m.ginv(scalar))
        else:c=m.ginv(c)
    f,g=put(hh,c)
    return f,g,ok and bool(c)

