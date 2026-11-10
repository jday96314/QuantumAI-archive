"""Shared scaling and numerator product for normal and carried-gap heads."""
from quantum_control import selected


def linear(m,E,math,src,power,dst):
    rt=[];m.lin_rows_terms(rt,math.rows(1,power),src,dst);E.emit_mac(rt)


def shear6(m,E,math,T,D,difference=False):
    n=len(T)
    def lin(src,power,dst):linear(m,E,math,src,power,dst)
    def prod(a,qa,b,qb,target):
        math.product(math.spec([(1,a,qa)]),math.spec([(1,b,qb)]),lambda tb:[target[i] for i in tb])
    def coeff(j,src,power,target):
        terms=[]
        if j<len(D):terms.append((1,D[j],1))
        if difference and j>1 and j-1<len(D):terms.append((1,D[j-1],1))
        if terms:math.product(math.spec(terms),math.spec([(1,src,power)]),lambda tb:[target[i] for i in tb])
        if difference and j==1:lin(src,power,target)
    if n>2:lin(T[0],2,T[2])
    if n>3:coeff(1,T[0],2,T[3])
    if n>4:
        lin(T[1],2,T[4]);coeff(2,T[0],2,T[4]);prod(T[0],2,T[0],1,T[4])
    if n>5:
        coeff(1,T[1],2,T[5]);coeff(3,T[0],2,T[5]);prod(T[0],2,T[1],1,T[5])


def carry_shear(m,E,math,H,D):
    # H[0] is structurally zero. This is H/(1+z^2*D*H) mod z^8.
    n=len(H)
    def lin(src,power,dst):linear(m,E,math,src,power,dst)
    def prod(a,qa,b,qb,target):
        math.product(math.spec([(1,a,qa)]),math.spec([(1,b,qb)]),lambda tb:[target[i] for i in tb])
    if n>4:lin(H[1],2,H[4])
    if n>5 and len(D)>1:prod(D[1],1,H[1],2,H[5])
    if n>6:
        lin(H[2],2,H[6])
        if len(D)>2:prod(D[2],1,H[1],2,H[6])
    if n>7:
        if len(D)>1:prod(D[1],1,H[2],2,H[7])
        if len(D)>3:prod(D[3],1,H[1],2,H[7])
        prod(H[1],2,H[1],1,H[7])


def apply(m,E,math,Bp,B,p,H,carry,shift,k,inverse=False):
    A=[B[0],*reversed(B)];D=[B[0],*reversed(Bp)]
    with m.opened(E,'head_normal_control',1) as ff:
        normal=ff[0];E.x(normal);E.cx(carry,normal)
        def scale(back=False):
            with m.inv_word(E,p) as pi,m.opened(E,'gap_head_scale',12) as product:
                for w in H[1:]:m.mon_slot(E,w,pi if back else p,p if back else pi,product)
        def shears(back=False):
            def normal_shear():
                with selected(m,E,normal) as (cm,ce):shear6(cm,ce,type(math)(cm,ce),H[2:],D)
            def shifted_shear():
                with selected(m,E,shift) as (cm,ce):shear6(cm,ce,type(math)(cm,ce),H[2:],D,True)
            def carried_shear():
                if k>1:
                    with selected(m,E,carry) as (cm,ce):carry_shear(cm,ce,type(math)(cm,ce),H,D)
            if back:carried_shear();shifted_shear();normal_shear()
            else:normal_shear();shifted_shear();carried_shear()
        math.known_zero.add(tuple(H[0]))
        if inverse:math.triangular(A,1,H,True);shears(True);scale(True)
        else:scale();shears();math.triangular(A,1,H)
        math.known_zero.remove(tuple(H[0]))
        E.cx(carry,normal);E.x(normal)
