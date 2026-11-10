"""Six-coefficient Mobius map without conjugation by the previous cofactor.

H_K(T)=T/(1+K*T), K=z^2*D, is an involution in characteristic two.
Modulo z^6 its corrections depend only on T[0],T[1], which never change.
This works when D has zero constant coefficient as well as when it is a unit.
The default D convention, like mobius_series, has implicit constant one.
"""
from tower_inverse import map_ops


def apply(m,E,math,A,D,p,T,inverse=False,p_is_inverse=False):
    if len(T)>6:
        raise ValueError('direct Mobius supports at most six coefficients')
    def square(w,back=False):
        ops=map_ops([m.gmul(1<<i,1<<i) for i in range(12)])
        for kind,i,j in (reversed(ops) if back else ops):
            getattr(E.b,kind)(w[i],w[j])
    def sx(src,dst):
        for i,q in enumerate(src):
            v=m.gmul(1<<i,1<<i)
            for j in range(12):
                if v>>j&1:E.cx(q,dst[j])
    def mac(a,b,c):E.xor(m.prod_rt(c,a,b))
    def shear():
        n=len(T)
        if n>2:sx(T[0],T[2])
        square(T[0])
        if n>3 and len(D)>1:mac(D[1],T[0],T[3])
        if n>4:
            sx(T[1],T[4])
            if len(D)>2:mac(D[2],T[0],T[4])
            with m.opened(E,'direct_mobius_original0',12) as old:
                for q,t in zip(T[0],old):E.cx(q,t)
                square(old,True);mac(T[0],old,T[4]);square(old)
                for q,t in zip(T[0],old):E.cx(q,t)
        if n>5:
            if len(D)>3:mac(D[3],T[0],T[5])
            mac(T[0],T[1],T[5])
            if len(D)>1:
                square(T[1]);mac(D[1],T[1],T[5]);square(T[1],True)
        square(T[0],True)
    def scale(back):
        with m.inv_word(E,p) as pi,m.opened(E,'direct_mobius_product',12) as product:
            undo=bool(back)^bool(p_is_inverse)
            for w in T:m.mon_slot(E,w,p if undo else pi,pi if undo else p,product)
    if inverse:
        math.triangular(A,1,T,True);shear();scale(True)
    else:
        scale(False);shear();math.triangular(A,1,T)
