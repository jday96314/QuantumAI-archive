"""General first compact slot, including a gap beginning at the first pivot."""
from quantum_control import selected
from tower_inverse import map_ops
from gap_serial import choose_wait


def apply(m,E,S,stash,wait,complete,inverse=False):
    s,r=S[-2:]
    if not inverse:choose_wait(m,E,r,s,None,None,wait)
    with m.opened(E,'opening_normal',1) as ctrl:
        normal=ctrl[0];E.x(normal);E.cx(wait,normal)
        def completion():
            with m.eq_test(E,r,0) as rz:E.b.ccx(rz,normal,complete)
        if not inverse:completion();E.cx(complete,r[0])
        def wait_route():
            E.swaps(wait,list(zip(s,stash[1])));E.cx(wait,r[0])
        if inverse:wait_route()
        with selected(m,E,normal) as (cm,ce):
            with m.inv_word(E,r) as ri:
                if inverse:E.swaps(normal,list(zip(r,ri)))
                def scale(back=False):
                    with m.opened(E,'opening_product',12) as temp:
                        cm.mon_slot(ce,s,r if back else ri,ri if back else r,temp)
                def last():
                    ops=map_ops([m.gmul(1<<i,1<<i) for i in range(12)])
                    for op,i,j in ops:getattr(E.b,op)(s[i],s[j])
                    ce.xor(m.prod_rt(S[-3],s,r))
                    for op,i,j in reversed(ops):getattr(E.b,op)(s[i],s[j])
                if not inverse:scale()
                if inverse:last()
                for j in (reversed(range(1,len(S)-2)) if inverse else range(1,len(S)-2)):
                    ce.xor(m.prod_rt(S[j-1],s,S[j]))
                if not inverse:last()
                if inverse:scale(True)
                if not inverse:E.swaps(normal,list(zip(r,ri)))
        if not inverse:wait_route()
        if inverse:E.cx(complete,r[0]);completion()
        E.cx(wait,normal);E.x(normal)
    if inverse:choose_wait(m,E,r,s,None,None,wait)
