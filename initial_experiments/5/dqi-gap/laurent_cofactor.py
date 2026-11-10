"""Laurent-chart arithmetic network for normal, waiting and degree-two carry charts.

Both real updates use beta*P + (alpha+lambda*x)*B plus a routed monomial.
The two-term convolution has k+1 field multiplications by evaluation, rather
than 2*k independent products. Conditional coefficient rotations supply the
degree-gap layout. All ordinary arithmetic remains unconditional and charged.
"""
from quantum_control import selected,controlled_ops,permutation_ops
from deferred_euclid_padding import fake_undo_ops,common_padding_ops


def apply(m,E,math,Bp,B,Q,pinv,carry,wait,shift,complete,inverse=False,norm=False):
    if len(B)<2:raise ValueError('initial carry uses the separate opening chart')
    zero,lam,beta,alpha=Q;k=len(B);targets=Bp+[beta]
    def xor(a,b):
        for q,t in zip(a,b):E.cx(q,t)
    def sw(a,b,flag):E.swaps(flag,list(zip(a,b)))
    def wait_prep(back=False):
        if back:E.cx(wait,beta[0])
        sw(alpha,zero,wait)
        if not back:E.cx(wait,beta[0])
    def previous_route(flag,back=False):
        ops=permutation_ops(targets,targets[1:]+targets[:1])
        controlled_ops(E,ops,flag,not back)
    def scale(back=False):
        with m.inv_word(E,beta) as bi,m.opened(E,'unified_cofactor_product',12) as product:
            for w in Bp:m.mon_slot(E,w,bi if back else beta,beta if back else bi,product)
            if norm:m.mon_slot(E,pinv,beta if back else bi,bi if back else beta,product)
    def convolution():
        # A=carry*alpha. On the carry chart B[0]=0, so B/x is a
        # polynomial. The pointwise left factor is
        # alpha+lambda+x+(1+1/x)*A, with no conditional bank shift.
        with m.opened(E,'laurent_alpha',12) as A:
            rt=[t for a,b in zip(alpha,A) for t in m.one_rt((carry,),0,(a,),0,b)]
            E.mac(rt)
            pts,I=m.vinv(k)
            remainder=[0]*k
            for a,x in enumerate(pts):
                left=math.spec([(1,alpha,1),(1,lam,1),(1^m.ginv(x),A,1)])
                right=math.spec([(m.gpow(x,i),w,1) for i,w in enumerate(B)])
                def out(tb,a=a):
                    return [w[b] for i,w in enumerate(targets) for b in m.bitlist(m.gmul(I[i][a],m.tmask(tb)),12)]
                terms=[]
                for ab,bb,tb in m.KT:
                    const=sum((x>>i)&1 for i in ab)%2
                    term=m.norm_term(m.fold(left,ab),const,m.fold(right,bb),0,out(tb))
                    if term is not None:terms.append(term)
                E.emit_mac(terms)
                for i in range(k):remainder[i]^=m.gmul(I[i][a],m.gpow(x,k))
            # Correct the single high coefficient aliased by k evaluations.
            for w,c in zip(targets,remainder):
                ops=[];m.lin_rows_terms(ops,math.rows(c,1),B[-1],w);E.emit_mac(ops)
            xor(A,beta)
            E.free(rt)
    def top(back=False):
        def outshift():
            xor(B[-1],alpha)
        def cleanup():
            for a,b in zip(beta,lam):E.b.ccx(carry,a,b)
            if k>1:
                for a,b in zip(B[-2],lam):E.b.ccx(carry,a,b)
            if k>1:
                with selected(m,E,carry) as (cm,ce):ce.xor(m.prod_rt(lam,alpha,B[-1]))
        if back:outshift();cleanup();sw(alpha,lam,carry)
        else:sw(alpha,lam,carry);cleanup();outshift()
    with m.opened(E,'unified_cofactor_flags',1) as flags:
        divide=flags[0];pad=complete
        E.cx(wait,pad)
        rt=m.one_rt((shift,),0,(complete,),1,divide)
        E.mac(rt);E.cx(carry,divide)
        padding=fake_undo_ops(Bp,B,beta,alpha)+common_padding_ops(Bp,B,[beta,alpha],False)
        if inverse:
            sw(zero,lam,wait)
            controlled_ops(E,padding,pad,True)
            top(True);convolution()
            previous_route(divide,True);scale(True)
            sw(beta,lam,carry);wait_prep(True)
        else:
            wait_prep();sw(beta,lam,carry)
            scale();previous_route(divide)
            convolution();top()
            controlled_ops(E,padding,pad)
            sw(zero,lam,wait)
        E.cx(carry,divide);E.free(rt);E.cx(wait,pad)
