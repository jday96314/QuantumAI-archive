"""Single compact Euclid step with reversible monomial completion arming.

A zero pivot is replaced by inverse(pinv), and its sublead is shifted so the
cofactor quotient is beta=1, alpha=0 on completed valid words. The fake input
contributes no boundary term to the retained remainder. A common x then keeps
completion stable. This is an unselected decoder component.
"""
from compact_euclid import mul,inv
from deferred_euclid_padding import fake_undo_ops,common_padding_ops
from monic_euclid_step import apply as generic
from quantum_coprime_padding import apply as emit_ops


def apply(m,E,math,Bp,B,C,pinv,flag,inverse=False,linear_padding=False):
    s,r=C[-2:]
    def prepare(inv=False):
        with m.inv_word(E,pinv) as pi,m.opened(E,'serial_fake_sublead',12) as p:
            rt=m.prod_rt(p,B[-1],pi);E.mac(rt)
            for source,target in zip(p,s):E.b.ccx(flag,source,target)
            for source,target in zip(pi,r):E.b.ccx(flag,source,target)
            E.free(rt)
    undo=fake_undo_ops(Bp,B,r,s)
    common=common_padding_ops(Bp,B,[r,s],False)
    if not inverse:
        with m.eq_test(E,r,0) as zero:E.cx(zero,flag)
        prepare()
        generic(m,E,math,Bp,B,C,pinv,boundary_flag=flag)
        E.swaps(None,list(zip(pinv,r)))
        if linear_padding:
            from controlled_linear import apply as linear
            linear(E.b,undo+common,flag)
        else:emit_ops(E.b,undo,flag);emit_ops(E.b,common,flag)
    else:
        if linear_padding:
            from controlled_linear import apply as linear
            linear(E.b,undo+common,flag,True)
        else:emit_ops(E.b,common,flag,True);emit_ops(E.b,undo,flag,True)
        E.swaps(None,list(zip(pinv,r)))
        generic(m,E,math,Bp,B,C,pinv,True,boundary_flag=flag)
        prepare(True)
        with m.eq_test(E,r,0) as zero:E.cx(zero,flag)
    return B,Bp+[r,s],C[:-2],pinv
