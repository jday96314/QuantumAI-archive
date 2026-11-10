"""One arithmetic network for normal, waiting and degree-two carry charts.

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
    def current_route(back=False):
        ops=permutation_ops(B,B[1:]+B[:1])
        if back:E.cx(carry,B[-1][0])
        controlled_ops(E,ops,carry,not back)
        if not back:E.cx(carry,B[-1][0])
    def scale(back=False):
        with m.inv_word(E,beta) as bi,m.opened(E,'unified_cofactor_product',12) as product:
            for w in Bp:m.mon_slot(E,w,bi if back else beta,beta if back else bi,product)
            if norm:m.mon_slot(E,pinv,beta if back else bi,bi if back else beta,product)
    def convolution():
        pts,I=m.vinv(k+1)
        for a,x in enumerate(pts):
            left=math.spec([(1,alpha,1),(x,lam,1)])
            right=math.spec([(m.gpow(x,i),w,1) for i,w in enumerate(B)])
            def out(tb,a=a):
                return [w[b] for i,w in enumerate(targets) for b in m.bitlist(m.gmul(I[i][a],m.tmask(tb)),12)]
            math.product(left,right,out)
    def monomial():
        ops=[]
        for i,t in enumerate(targets):
            if i:
                xor(B[i-1],t)
                for a,b in zip(B[i-1],t):ops.append(('cx',a,b))
            if i>1:
                for a,b in zip(B[i-2],t):ops.append(('cx',a,b))
        controlled_ops(E,ops,carry)
    def top(back=False):
        def outshift():
            xor(B[-1],alpha);ops=[('cx',a,b) for a,b in zip(B[-1],alpha)]
            if k>1:ops += [('cx',a,b) for a,b in zip(B[-2],alpha)]
            controlled_ops(E,ops,carry)
        def cleanup():
            for a,b in zip(beta,lam):E.b.ccx(carry,a,b)
            if k>2:
                for a,b in zip(B[-3],lam):E.b.ccx(carry,a,b)
            if k>1:
                with selected(m,E,carry) as (cm,ce):ce.xor(m.prod_rt(lam,alpha,B[-2]))
        if back:outshift();cleanup();sw(alpha,lam,carry)
        else:sw(alpha,lam,carry);cleanup();outshift()
    with m.opened(E,'unified_cofactor_flags',2) as flags:
        divide,pad=flags
        rt=m.one_rt((shift,),0,(complete,),1,divide)
        E.mac(rt);E.cx(carry,divide);E.cx(complete,pad);E.cx(wait,pad)
        padding=fake_undo_ops(Bp,B,beta,alpha)+common_padding_ops(Bp,B,[beta,alpha],False)
        if inverse:
            sw(zero,lam,wait)
            controlled_ops(E,padding,pad,True)
            current_route();top(True);monomial();convolution();current_route(True)
            previous_route(divide,True);scale(True)
            sw(beta,lam,carry);wait_prep(True)
        else:
            wait_prep();sw(beta,lam,carry)
            scale();previous_route(divide)
            current_route();convolution();monomial();top();current_route(True)
            controlled_ops(E,padding,pad)
            sw(zero,lam,wait)
        E.cx(wait,pad);E.cx(complete,pad);E.cx(carry,divide);E.free(rt)
