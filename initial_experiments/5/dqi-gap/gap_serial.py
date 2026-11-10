"""Shared-tail one-degree compact Euclid update with a reversible gap carry.

Charts: ordinary degree one; waiting at a missing pivot; degree-two repair of
the carry. A previous cofactor one degree too short is stored multiplied by x
and divided by a conditional coefficient rotation at the next real step.
Only local head/cofactor arithmetic is selected; both large tail passes are
shared. Completion uses the incumbent's reversible common-x padding.
"""
from contextlib import contextmanager
from mobius_series import BoundedArithmetic
from direct_mobius import apply as mobius
from continued_fraction_series import apply as cf
from gap_block import head as gap_head
from quantum_control import selected,controlled_ops,permutation_ops,picked_word
from deferred_euclid_padding import fake_undo_ops,common_padding_ops
from tower_inverse import map_ops
UNIFIED=False


def xor(E,a,b):
    for q,t in zip(a,b):E.cx(q,t)


def choose_wait(m,E,r,s,carry,prevshift,out):
    with m.eq_test(E,r,0) as rz,m.eq_test(E,s,0) as sz:
        bits=[rz,sz]+([] if carry is None else [carry,prevshift])
        with m.eq_test(E,bits,1) as yes:E.cx(yes,out)


def normal_cofactor(m,E,Bp,B,beta,alpha,pinv,shift,complete,inverse=False,norm=False,routing_E=None):
    routing_E=E if routing_E is None else routing_E
    def mac(a,b,t):E.xor(m.prod_rt(t,a,b))
    def top():
        if inverse:xor(E,B[-1],alpha)
        mac(alpha,B[-1],beta)
        if len(B)>1:xor(E,B[-2],beta)
        if not inverse:xor(E,B[-1],alpha)
    def low():
        for i,w in enumerate(Bp):
            mac(alpha,B[i],w)
            if i:xor(E,B[i-1],w)
    # Scale before routing so the inverse source beta stays intact throughout
    # the measurement-based inverse-word context.
    def scale():
        if hasattr(m,'branch_flag'):
            bm,be=m.base_module,m.base_emitter
            with picked_word(bm,be,beta,m.branch_flag) as factor:
                with bm.inv_word(be,factor) as fi,bm.opened(be,'selected_cofactor_product',12) as product:
                    for w in Bp:bm.mon_slot(be,w,fi if inverse else factor,factor if inverse else fi,product)
                    if norm:bm.mon_slot(be,pinv,factor if inverse else fi,fi if inverse else factor,product)
            return
        with m.inv_word(E,beta) as bi,m.opened(E,'gap_serial_scale',12) as product:
            for w in Bp:m.mon_slot(E,w,bi if inverse else beta,beta if inverse else bi,product)
            if norm:m.mon_slot(E,pinv,beta if inverse else bi,bi if inverse else beta,product)
    def route():
        words=Bp+[beta]
        if len(words)>1:
            ops=permutation_ops(words,words[1:]+words[:1])
            # permute(source,target) sends source[i] to target[i]. To divide
            # coefficients by x, target[i] must be source[i-1].
            controlled_ops(routing_E,ops,shift,not inverse)
    undo=fake_undo_ops(Bp,B,beta,alpha)
    common=common_padding_ops(Bp,B,[beta,alpha],False)
    if inverse:
        controlled_ops(routing_E,undo+common,complete,True)
        top();low();route();scale()
    else:
        scale();route();low();top()
        controlled_ops(routing_E,undo+common,complete)


def quadratic(m,E,Bp,B,beta,q0,q1,pinv,inverse=False,norm=False):
    def mac(a,b,t):E.xor(m.prod_rt(t,a,b))
    def top():
        def beta_change():
            if len(B)>2:xor(E,B[-3],beta)
            if len(B)>1:mac(q1,B[-2],beta)
            mac(q0,B[-1],beta)
        if not inverse:beta_change()
        if inverse:xor(E,B[-1],q1)
        if len(B)>1:xor(E,B[-2],q0)
        mac(q1,B[-1],q0)
        if not inverse:xor(E,B[-1],q1)
        if inverse:beta_change()
    if inverse:top()
    with m.inv_word(E,beta) as bi,m.opened(E,'gap_serial_quadratic',12) as product:
        for i,w in enumerate(Bp):
            if not inverse:m.mon_slot(E,w,beta,bi,product)
            mac(q0,B[i],w)
            if i:mac(q1,B[i-1],w)
            if i>1:xor(E,B[i-2],w)
            if inverse:m.mon_slot(E,w,bi,beta,product)
        if norm:m.mon_slot(E,pinv,beta if inverse else bi,bi if inverse else beta,product)
    if not inverse:top()


def apply(m,E,math,Bp,B,C,pinv,stash,carry,wait,prevshift,complete,inverse=False):
    """Flags carry/wait/prevshift are supplied and retained; complete is recorded."""
    raw=list(reversed(C))[:2];tail=list(reversed(C))[2:];physical=stash+raw
    rev=[B[0],*reversed(B)];drev=[B[0],*reversed(Bp)]
    with m.opened(E,'gap_serial_controls',2) as ctrl:
        normal,divide=ctrl
        E.x(normal);E.cx(carry,normal);E.cx(wait,normal)
        # A completion preserves the padded previous-cofactor chart.
        def divide_flag():
            E.mac(m.one_rt((prevshift,),0,(complete,),1,divide))
        def divide_clear():
            E.free(m.one_rt((prevshift,),0,(complete,),1,divide))
        def local(H,back=False,record=True):
            with selected(m,E,normal) as (cm,ce):
                ar=BoundedArithmetic(cm,ce);ar.virtual_zero=True
                if back:cf(cm,ce,ar,H[2:],True,[complete],record)
                mobius(cm,ce,ar,rev,drev,pinv,H[2:],back,p_is_inverse=True)
                if not back:cf(cm,ce,ar,H[2:],False,[complete],record)
            with selected(m,E,carry) as (cm,ce):
                ar=BoundedArithmetic(cm,ce);ar.virtual_zero=True
                def gap_pair(inv=False):
                    # e,d,a,b -> e,d,u1=a/d,u0=b/d+u1^2; e=0.
                    if inv:
                        for i,q in enumerate(H[2]):
                            value=m.gmul(1<<i,1<<i)
                            for j in range(12):
                                if value>>j&1:ce.cx(q,H[3][j])
                    with cm.inv_word(ce,H[1]) as di,cm.opened(ce,'gap_pair_scale',12) as product:
                        for w in H[2:]:cm.mon_slot(ce,w,H[1] if inv else di,di if inv else H[1],product)
                    if not inv:
                        for i,q in enumerate(H[2]):
                            value=m.gmul(1<<i,1<<i)
                            for j in range(12):
                                if value>>j&1:ce.cx(q,H[3][j])
                if back:gap_pair(True)
                mobius(cm,ce,ar,rev,drev,pinv,H,back,p_is_inverse=True)
                if not back:gap_pair()
        with m.opened(E,'gap_serial_head',48) as extra:
            Q=[extra[i*12:(i+1)*12] for i in range(4)]
            new=Bp+[Q[2],Q[3]];finalrev=[new[0],*reversed(new)]
            def cof(back=False,norm=False):
                if UNIFIED and len(B)>1:
                    from unified_cofactor import apply as common_cofactor
                    common_cofactor(m,E,math,Bp,B,Q,pinv,carry,wait,prevshift,complete,back,norm)
                    return
                divide_flag()
                with selected(m,E,normal) as (cm,ce):
                    normal_cofactor(cm,ce,Bp,B,Q[2],Q[3],pinv,divide,complete,back,norm,routing_E=E)
                divide_clear()
                if len(B)>1:
                    mid=Bp[1:]+[Q[1],Q[3],Q[2]]
                    ops=permutation_ops(mid+[Q[0],Bp[0]],new+[Q[0],Q[1]])
                    if back:controlled_ops(E,ops,carry,True)
                    with selected(m,E,carry) as (cm,ce):
                        quadratic(cm,ce,Bp[1:],B[1:],Q[1],Q[3],Q[2],pinv,back,norm)
                    if not back:controlled_ops(E,ops,carry)
                else:
                    # The first gap has previous cofactor zero. Its d is
                    # absorbed into pinv; the other public operand is one.
                    ops=permutation_ops([Q[3],Q[2]],[Q[2],Q[3]])
                    def norm_first(back=False):
                        if back:E.b.cx(carry,Q[1][0])
                        if back:
                            with m.inv_word(E,pinv) as pi:E.swaps(carry,list(zip(pinv,pi)))
                        E.swaps(carry,list(zip(pinv,Q[1])))
                        if not back:
                            with m.inv_word(E,pinv) as pi:E.swaps(carry,list(zip(pinv,pi)))
                            E.b.cx(carry,Q[1][0])
                    if back and norm:norm_first(True)
                    controlled_ops(E,ops,carry,back)
                    if not back and norm:norm_first()
                wait_ops=permutation_ops([Q[0]]+Bp+[Q[1]]+B+[Q[2],Q[3]],B+new+[Q[0],Q[1]])
                # Branches are disjoint, so their mutually controlled chart
                # maps commute even though the underlying permutations overlap.
                controlled_ops(E,wait_ops,wait,back)
            def copy():
                for a,b in zip(physical,Q):xor(E,a,b)
            def swap():
                for a,b in zip(physical,Q):E.swaps(None,list(zip(a,b)))
            def boundary(poly,source,start,flag):
                # The first retained carry moment is structurally zero.
                with m.opened(E,'gap_boundary_mask',12) as masked:
                    rt=[t for a,b in zip(source,masked) for t in m.one_rt((flag,),0,(a,),0,b)]
                    E.mac(rt)
                    for i,w in enumerate(tail):
                        j=i+start-1
                        if j<len(poly):E.xor(m.prod_rt(w,poly[j],masked))
                    E.free(rt)
            if not inverse:
                copy();local(Q);cof()
                math.triangular_suffix(finalrev,1,list(reversed(C)),2)
                boundary(finalrev,stash[1],4,carry)
                cof(True)
                boundary(rev,raw[1],2,wait)
                math.triangular(rev,1,tail,True)
                local(physical,record=False)
                for a,b in zip(Q,physical):xor(E,a,b)
                cof(norm=True);swap()
            else:
                swap();cof(True,norm=True)
                for a,b in zip(Q,physical):xor(E,a,b)
                local(physical,True,record=False)
                math.triangular(rev,1,tail)
                boundary(rev,raw[1],2,wait)
                cof()
                boundary(finalrev,stash[1],4,carry)
                math.triangular_suffix(finalrev,1,list(reversed(C)),2,True)
                cof(True);local(Q,True);copy()
        E.cx(wait,normal);E.cx(carry,normal);E.x(normal)
    return B,Bp+raw,C[:-2],pinv


def selected_cofactor(m,E,math,Bp,B,Q,pinv,carry,wait,prevshift,complete,inverse=False,norm=False):
    back=inverse
    new=Bp+[Q[2],Q[3]]
    with m.opened(E,'selected_cofactor_controls',2) as cc:
        normal,divide=cc
        E.x(normal);E.cx(carry,normal);E.cx(wait,normal)
        def divide_flag():E.mac(m.one_rt((prevshift,),0,(complete,),1,divide))
        def divide_clear():E.free(m.one_rt((prevshift,),0,(complete,),1,divide))
        divide_flag()
        with selected(m,E,normal) as (cm,ce):
            normal_cofactor(cm,ce,Bp,B,Q[2],Q[3],pinv,divide,complete,back,norm,routing_E=E)
        divide_clear()
        if len(B)>1:
            mid=Bp[1:]+[Q[1],Q[3],Q[2]]
            ops=permutation_ops(mid+[Q[0],Bp[0]],new+[Q[0],Q[1]])
            if back:controlled_ops(E,ops,carry,True)
            with selected(m,E,carry) as (cm,ce):
                quadratic(cm,ce,Bp[1:],B[1:],Q[1],Q[3],Q[2],pinv,back,norm)
            if not back:controlled_ops(E,ops,carry)
        else:
            # The first gap has previous cofactor zero. Its d is
            # absorbed into pinv; the other public operand is one.
            ops=permutation_ops([Q[3],Q[2]],[Q[2],Q[3]])
            def norm_first(back=False):
                if back:E.b.cx(carry,Q[1][0])
                if back:
                    with m.inv_word(E,pinv) as pi:E.swaps(carry,list(zip(pinv,pi)))
                E.swaps(carry,list(zip(pinv,Q[1])))
                if not back:
                    with m.inv_word(E,pinv) as pi:E.swaps(carry,list(zip(pinv,pi)))
                    E.b.cx(carry,Q[1][0])
            if back and norm:norm_first(True)
            controlled_ops(E,ops,carry,back)
            if not back and norm:norm_first()
        wait_ops=permutation_ops([Q[0]]+Bp+[Q[1]]+B+[Q[2],Q[3]],B+new+[Q[0],Q[1]])
        # Branches are disjoint, so their mutually controlled chart
        # maps commute even though the underlying permutations overlap.
        controlled_ops(E,wait_ops,wait,back)
        E.cx(wait,normal);E.cx(carry,normal);E.x(normal)
