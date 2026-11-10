"""Gap-aware three-slot compact transform with shared arithmetic and tail.

Local Schur updates choose ordinary, wait, or quadratic carry coordinates.
All cofactor passes use the single short-convolution network. The two large
polynomial tail passes occur once per block. Only a terminal wait crosses a
block boundary; its raw remainder word is retained for the following block.
"""
from contextlib import contextmanager
from mobius_series import BoundedArithmetic
from quantum_control import selected,picked_word
from gap_serial import choose_wait,xor,selected_cofactor
from gap_mobius import apply as mobius
from laurent_cofactor import apply as cofactor
from blocked_monic_euclid import schedule


@contextmanager
def previous_shift(m,E,k,start,waits,flags,terminal):
    if k<2:
        yield terminal
        return
    # This history bit is not otherwise read by the current slot. Its
    # original value is restored before the next slot or flag update.
    g=waits[k-2];E.cx(terminal,g)
    for i in range(start,k):
        if i>=2:E.b.ccx(waits[i-2],flags[i],g)
    yield g
    for i in reversed(range(start,k)):
        if i>=2:E.b.ccx(waits[i-2],flags[i],g)
    E.cx(terminal,g)


def apply(m,E,math,Bp,B,C,p,stash,k,count,waits,flags,terminal,inverse=False):
    raw=list(reversed(C))[:2*count];tail=list(reversed(C))[2*count:];physical=stash+raw
    rev=[B[0],*reversed(B)]
    def local(H,back=False,record=True):
        with previous_shift(m,E,k,k,waits,flags,terminal) as g:
            if not back:mobius(m,E,math,Bp,B,p,H,waits[k-1],g,k)
        for j in (reversed(range(count)) if back else range(count)):
            kk=k+j;carry=waits[kk-1];wait=waits[kk];complete=flags[kk];F=H[2+2*j:]
            with previous_shift(m,E,kk,k,waits,flags,terminal) as g:
                if not back and record:choose_wait(m,E,F[0],F[1],carry,g,wait)
                with m.opened(E,'schur_normal',1) as ff:
                    normal=ff[0];E.x(normal);E.cx(carry,normal);E.cx(wait,normal)
                    def normal_cf():
                        if not back:
                            if record:
                                with m.eq_test(E,F[0],0) as z:E.b.ccx(normal,z,complete)
                            E.cx(complete,F[0][0])
                        if back:
                            with selected(m,E,normal) as (cm,ce):BoundedArithmetic(cm,ce).inverse(F)
                        with picked_word(m,E,F[0],normal) as pivot:
                            with m.inv_word(E,pivot) as bi,m.opened(E,'schur_product',12) as product:
                                for w in F[1:]:m.mon_slot(E,w,pivot if back else bi,bi if back else pivot,product)
                        if not back:
                            with selected(m,E,normal) as (cm,ce):BoundedArithmetic(cm,ce).inverse(F)
                        if back:
                            E.cx(complete,F[0][0])
                            if record:
                                with m.eq_test(E,F[0],0) as z:E.b.ccx(normal,z,complete)
                    def carry_cf():
                        pair=H[:2] if j==0 else H[2*j:2*j+2]
                        G=pair[1:]+F
                        if back:
                            with selected(m,E,carry) as (cm,ce):BoundedArithmetic(cm,ce).inverse(G)
                        with picked_word(m,E,G[0],carry) as pivot:
                            with m.inv_word(E,pivot) as di,m.opened(E,'gap_schur_product',12) as product:
                                for w in G[1:]:m.mon_slot(E,w,pivot if back else di,di if back else pivot,product)
                        if not back:
                            with selected(m,E,carry) as (cm,ce):BoundedArithmetic(cm,ce).inverse(G)
                    if back:carry_cf();normal_cf()
                    else:normal_cf();carry_cf()
                    E.cx(wait,normal);E.cx(carry,normal);E.x(normal)
                if back and record:choose_wait(m,E,F[0],F[1],carry,g,wait)
        with previous_shift(m,E,k,k,waits,flags,terminal) as g:
            if back:mobius(m,E,math,Bp,B,p,H,waits[k-1],g,k,True)
    with m.opened(E,'unified_block_head',12*(2*count+1)) as bank:
        H=[()]+[bank[i*12:(i+1)*12] for i in range(2*count+1)]
        ff,(np,nb)=schedule(Bp,B,H[2:]);newrev=[nb[0],*reversed(nb)]
        def last_scale(back=False):
            with picked_word(m,E,p,waits[k+count-1]) as pivot:
                with m.inv_word(E,pivot) as pi,m.opened(E,'terminal_wait_scale',12) as product:
                    m.mon_slot(E,H[1],pivot if back else pi,pi if back else pivot,product)
        def cof_inner(back=False):
            for j in (reversed(range(count)) if back else range(count)):
                kk=k+j;oldp,old,be,al=ff[j];Q=H[:2]+[be,al]
                with previous_shift(m,E,kk,k,waits,flags,terminal) as g:
                    fn=selected_cofactor if len(old)==1 else cofactor
                    fn(m,E,math,oldp,old,Q,p,waits[kk-1],waits[kk],g,flags[kk],back,True)
        def cof(back=False):
            if back:last_scale(True)
            with m.opened(E,'cofactor_zero_coordinate',12) as zero:
                H[0]=zero
                cof_inner(back)
                H[0]=()
            if not back:last_scale()
        def scalar_prefix(back=False):
            # Restore the last slot's normalization using head coefficients,
            # without replaying polynomial arithmetic.
            for j in (reversed(range(count-1)) if back else range(count-1)):
                kk=k+j;carry=waits[kk-1];wait=waits[kk]
                source=H[1] if j==0 else H[2*j+1]
                with m.opened(E,'prefix_scalar_controls',1) as gg,m.opened(E,'prefix_beta',12) as beta:
                    normal=gg[0];E.x(normal);E.cx(carry,normal);E.cx(wait,normal)
                    rt=[('L',(),1,(beta[0],)),('L',(normal,),0,(beta[0],)),('L',(carry,),0,(beta[0],))]
                    rt += [t for a,b in zip(H[2+2*j],beta) for t in m.one_rt((normal,),0,(a,),0,b)]
                    rt += [t for a,b in zip(source,beta) for t in m.one_rt((carry,),0,(a,),0,b)]
                    E.mac(rt)
                    with m.inv_word(E,beta) as bi,m.opened(E,'prefix_p_product',12) as product:
                        m.mon_slot(E,p,beta if back else bi,bi if back else beta,product)
                    E.free(rt)
                    E.cx(wait,normal);E.cx(carry,normal);E.x(normal)
        def boundary(poly,source,start,flag):
            with m.opened(E,'fused_boundary',12) as masked:
                rt=[t for a,b in zip(source,masked) for t in m.one_rt((flag,),0,(a,),0,b)]
                E.mac(rt)
                for i,w in enumerate(tail):
                    index=i+start-1
                    if index<len(poly):E.xor(m.prod_rt(w,poly[index],masked))
                E.free(rt)
        def outgoing_boundary():
            scalar_prefix()
            with selected(m,E,waits[k+count-1]) as (cm,ce):
                with cm.inv_word(ce,p) as pi,cm.opened(ce,'outgoing_raw_wait',12) as out:
                    rt=m.prod_rt(out,H[-1],pi);ce.mac(rt)
                    for i,w in enumerate(tail):
                        if i+1<len(rev):E.xor(m.prod_rt(w,rev[i+1],out))
                    ce.free(rt)
            scalar_prefix(True)
        def copy():
            for a,b in zip(physical,H):xor(E,a,b)
        def swap():
            for a,b in zip(physical,H):E.swaps(None,list(zip(a,b)))
        if not inverse:
            copy();local(H);cof()
            math.triangular_suffix(newrev,1,list(reversed(C)),2*count)
            boundary(newrev,stash[1],2*count+2,waits[k-1])
            cof(True);outgoing_boundary()
            math.triangular(rev,1,tail,True)
            local(physical,record=False)
            for a,b in zip(H,physical):xor(E,a,b)
            cof();swap()
        else:
            swap();cof(True)
            for a,b in zip(H,physical):xor(E,a,b)
            local(physical,True,False)
            math.triangular(rev,1,tail);outgoing_boundary()
            cof();boundary(newrev,stash[1],2*count+2,waits[k-1])
            math.triangular_suffix(newrev,1,list(reversed(C)),2*count,True)
            cof(True);local(H,True);copy()
    return Bp,B,C[:-2*count],p
