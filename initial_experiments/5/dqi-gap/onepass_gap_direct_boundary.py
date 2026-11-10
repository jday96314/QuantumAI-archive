"""One-pass gap-aware slot using raw-prefix reconstruction from output invariants.

For new inverse scale pi=1/pnew, high current coefficients H,J:
 ordinary: raw_r=pi, raw_s=pi*H;
 carry: raw_stash=pi, raw_r=pi*H, raw_s=pi*(J+H^2);
 completion: raw_r=0, raw_s=pi*Bprev[0].
Waiting prefixes remain in the carry word. These identities eliminate the
copied head and both extra cofactor passes, while retaining a literal inverse.
"""
from quantum_control import picked_word,controlled_ops
from mobius_series import BoundedArithmetic
from laurent_cofactor import apply as cofactor
from gap_serial import selected_cofactor


def apply(m,E,math,Bp,B,C,p,stash,carry,wait,shift,complete,inverse=False):
    raw=list(reversed(C))[:2];tail=list(reversed(C))[2:];H=[(),stash[1],*raw]
    new=Bp+raw;newrev=[new[0],*reversed(new)];oldrev=[B[0],*reversed(B)]
    with m.opened(E,'onepass_controls',2) as controls:
        normal,pad=controls;E.x(normal);E.cx(carry,normal);E.cx(wait,normal)
        def mobius(back=False):
            def scale():
                E.cx(carry,normal) # normal now means !wait
                with picked_word(m,E,p,normal) as factor:
                    with m.inv_word(E,factor) as fi,m.opened(E,'onepass_head_product',12) as product:
                        for w in H[1:]:m.mon_slot(E,w,fi if back else factor,factor if back else fi,product)
                E.cx(carry,normal)
            math.known_zero.add(())
            if back:math.triangular(oldrev,1,H,True);scale()
            else:scale();math.triangular(oldrev,1,H)
            math.known_zero.remove(())
        def schur(back=False,record=True):
            def normal_step():
                if not back:
                    if record:
                        with m.eq_test(E,raw[0],0) as z:E.b.ccx(normal,z,complete)
                    E.cx(complete,raw[0][0])
                with picked_word(m,E,raw[0],normal) as beta:
                    with m.inv_word(E,beta) as bi,m.opened(E,'onepass_schur_product',12) as product:
                        m.mon_slot(E,raw[1],beta if back else bi,bi if back else beta,product)
                if back:
                    E.cx(complete,raw[0][0])
                    if record:
                        with m.eq_test(E,raw[0],0) as z:E.b.ccx(normal,z,complete)
            def carry_step():
                ops=[]
                for i,q in enumerate(raw[0]):
                    value=m.gmul(1<<i,1<<i)
                    ops += [('cx',q,raw[1][j]) for j in range(12) if value>>j&1]
                if back:controlled_ops(E,ops,carry)
                with picked_word(m,E,stash[1],carry) as beta:
                    with m.inv_word(E,beta) as bi,m.opened(E,'onepass_gap_product',12) as product:
                        for w in raw:m.mon_slot(E,w,beta if back else bi,bi if back else beta,product)
                if not back:controlled_ops(E,ops,carry)
            if back:carry_step();normal_step()
            else:normal_step();carry_step()
        def local(back=False):
            if back:schur(True);mobius(True)
            else:mobius();schur()
        def cof(back=False):
            with m.opened(E,'onepass_zero_coordinate',12) as zero:
                Q=[zero,stash[1],*raw]
                fn=selected_cofactor if len(B)==1 else cofactor
                fn(m,E,math,Bp,B,Q,p,carry,wait,shift,complete,back,True)
        def extract_old():
            for a,b in zip(B,new[1:]):E.swaps(pad,list(zip(a,b)))
        def boundary():
            if not tail:return
            high,next_high=new[-1],new[-2]
            with m.inv_word(E,p) as pi:
                with m.opened(E,'reconstructed_prefix',12) as prefix,m.opened(E,'reconstruction_product',12) as product:
                    # raw_r / pi = [ordinary] + [carry]*H.
                    rt=[('L',(),1,(prefix[0],))]+[('L',(f,),0,(prefix[0],)) for f in [carry,wait,complete]]
                    rt += [t for a,b in zip(high,prefix) for t in m.one_rt((carry,),0,(a,),0,b)]
                    E.mac(rt);m.mon_slot(E,prefix,pi,p,product)
                    for i,w in enumerate(tail):
                        if i+2<len(newrev):E.xor(m.prod_rt(w,newrev[i+2],prefix))
                    m.mon_slot(E,prefix,p,pi,product);E.free(rt)
                    # raw_s / pi selects H, J+H^2, Bprev[0], or zero.
                    rt=[];rows=math.rows(1,2)
                    for bit,t in enumerate(prefix):
                        rt.append(('L',(high[bit],),0,(t,)))
                        carry_form=set([next_high[bit],high[bit]])
                        carry_form.symmetric_difference_update(high[i] for i in m.bitlist(rows[bit],12))
                        rt += m.one_rt((carry,),0,tuple(carry_form),0,t)
                        rt += m.one_rt((wait,),0,(high[bit],),0,t)
                        rt += m.one_rt((complete,),0,(B[0][bit],high[bit]),0,t)
                    E.mac(rt);m.mon_slot(E,prefix,pi,p,product)
                    for i,w in enumerate(tail):
                        if i+1<len(newrev):E.xor(m.prod_rt(w,newrev[i+1],prefix))
                    m.mon_slot(E,prefix,p,pi,product);E.free(rt)
                    # The incoming raw carry is inverse(pnew).
                    rt=[t for a,b in zip(pi,prefix) for t in m.one_rt((carry,),0,(a,),0,b)]
                    E.mac(rt)
                    for i,w in enumerate(tail):
                        if i+3<len(newrev):E.xor(m.prod_rt(w,newrev[i+3],prefix))
                    E.free(rt)
        def arm_pad():E.cx(complete,pad);E.cx(wait,pad)
        if inverse:
            arm_pad();extract_old();math.triangular(oldrev,1,tail);extract_old();arm_pad()
            boundary();math.triangular(newrev,1,tail,True);cof(True);local(True)
        else:
            local();cof();math.triangular(newrev,1,tail);boundary()
            arm_pad();extract_old();math.triangular(oldrev,1,tail,True);extract_old();arm_pad()
        E.cx(wait,normal);E.cx(carry,normal);E.x(normal)
    return B,new,C[:-2],p
