"""Shared polynomial action for ordinary and leading-zero Jacobi charts.

Temporary boundary words retain the original head while the tail changes.
After swapping in the new boundary, the original head is reconstructed from
the low output polynomial coefficients for exact measurement cleanup. This
uses factorized polynomial information even when an armed chart has a common
factor, rather than assuming the rational function determines its heads.
"""
from contextlib import ExitStack,contextmanager
import lookahead_jacobi as lookahead
from compact_arming import install
install(lookahead.jacobi.cf)


@contextmanager
def selected(m,E,flag):
    # The upstream Builder permits conditional CZ but not conditional CCZ.
    # Mask product operands before measurement repayment of mixed formulas.
    with lookahead.selected(m,E,flag) as (cm,ce):
        def release(rt):
            products=[term for term in rt if term[0]=='P']
            with ce.mask_terms(products) as masked:
                product_terms=products if masked is None else masked
                terms=[term for term in rt if term[0]=='L']+product_terms
                per={}
                for term in terms:
                    for target in term[-1]:per.setdefault(target,[]).append(term)
                for target,group in reversed(list(per.items())):
                    measurement=E.b.hmr(target)
                    for term in group:
                        if term[0]=='L':
                            for source in term[1]:E.b.cz(flag,source,cond=measurement)
                            if term[2]:E.b.z(flag,cond=measurement)
                        elif masked is not None:
                            E.pivot_pair(term[1],term[2],term[3],term[4],lambda a,b:E.b.cz(a,b,cond=measurement))
                        else:
                            E.fan_pair(term[1],term[2],term[3],term[4],lambda fan:E.b.cz(flag,fan,cond=measurement))
        ce.emit_free=release
        yield cm,ce


def linear(words,target):
    return [('L',(q,),0,(t,)) for word in words if word is not None for q,t in zip(word,target)]


def selected_linear(flag,polarity,words,target):
    rt=[]
    for bit,t in enumerate(target):
        source=tuple(sorted(q for word in words if word is not None for q in (word[bit],)))
        if source:rt.append(('P',(flag,),polarity,source,0,(t,)))
    return rt


def params(m,E,math,H,flag,out,clear=False,masks=None):
    a,b,h,e=H;dn,t,u,v,s=out
    with m.opened(E,'shared_synthesis_c',12) as c:
        cr=[];m.bil_terms(cr,m.KT,m.spec_plain(b),math.spec([(1,e,1),(1,b,1)]),tgt=c);cr+=linear([h],c);E.mac(cr)
        terms=selected_linear(flag,1,[e,b],dn)
        terms+=linear([c],t)+selected_linear(flag,0,[c],t)+[('L',(flag,),0,(t[0],))]
        terms+=linear([b],u)+selected_linear(flag,0,[b,h],u)
        terms += [('L',(),1,(v[0],)),('L',(flag,),0,(v[0],))]+selected_linear(flag,0,[e],v)
        terms+=linear([a],s)+selected_linear(flag,0,[b],s)
        if masks is not None:repay(E,terms,masks)
        else:(E.free if clear else E.mac)(terms)
        E.free(cr)


def boundary_terms(m,math,H,P,Q,flag,out):
    a,b,h,e=H;B1,N1,B2=out
    q0=Q[0] if Q else None;q1=Q[1] if len(Q)>1 else None;p0=P[0] if P else None
    rt=linear([e,q0],B1)+linear([h,q1],B2)
    source=math.spec([(1,w,1) for w in [e,b,q0] if w is not None])
    m.bil_terms(rt,m.KT,m.spec_plain(a),source,tgt=N1)
    rt+=selected_linear(flag,0,[b],N1)+selected_linear(flag,1,[p0],B2)
    if q0 is not None:rt+=m.prod_rt(B2,e,q0)
    return rt


def boundary(m,E,math,H,P,Q,flag,out,clear=False):
    (E.free if clear else E.mac)(boundary_terms(m,math,H,P,Q,flag,out))


def repay(E,terms,masks):
    for term in terms:
        for target in term[-1]:
            measurement=masks[target]
            if term[0]=='L':
                for source in term[1]:E.b.z(source,cond=measurement)
                if term[2]:E.b.neg(cond=measurement)
            else:E.pivot_pair(term[1],term[2],term[3],term[4],lambda a,b:E.b.cz(a,b,cond=measurement))


def original_head(m,E,math,W,flag,J,inverse=False):
    """XOR-compute (or erase) old encoded b,h,e from four low coefficients."""
    a=W[0];N1=W[2];N2=W[4] if len(W)>4 else None;N3=W[6] if len(W)>6 else None
    B1=W[1];B2=W[3];B3=W[5] if len(W)>5 else None
    jb,jh,je=J
    with ExitStack() as stack:
        q=[stack.enter_context(m.opened(E,f'shared_head_q{i}',12)) for i in range(2)]
        q2,q3=q
        def quotient(clear=False,make_head=False,clear_head=False):
            with m.opened(E,'shared_head_pivot',12) as pivot:
                pr=linear([a],pivot)+selected_linear(flag,0,[N1],pivot);E.mac(pr)
                with lookahead.jacobi.cf.armed(m,E,pivot),m.inv_word(E,pivot) as iv:
                    rt=[]
                    for source,target in zip((N2,N3),q):
                        if source is not None:rt+=m.prod_rt(target,source,iv)
                    (E.free if clear else E.mac)(rt)
                    if make_head or clear_head:
                        E.x(flag)
                        with selected(m,E,flag) as (cm,ce):
                            head_rt=linear([B1],jb)+cm.prod_rt(jb,N1,iv)
                            (ce.free if clear_head else ce.mac)(head_rt)
                        E.x(flag)
                E.free(pr)
        quotient(make_head=not inverse)
        def gap():
            # Called in forward order to compute and reverse order to erase.
            with selected(m,E,flag) as (cm,ce):
                if inverse:
                    ce.xor(linear([B2,q3],jh)+cm.prod_rt(jh,je,q2))
                    ce.xor(linear([B1,q2],je));ce.xor(linear([N1],jb))
                else:
                    ce.xor(linear([N1],jb));ce.xor(linear([B1,q2],je))
                    ce.xor(linear([B2,q3],jh)+cm.prod_rt(jh,je,q2))
        if inverse:gap()
        E.x(flag)
        with selected(m,E,flag) as (cm,ce):
            encoding=cm.prod_rt(jh,jb,je)
            m.lin_rows_terms(encoding,math.rows(1,2),jb,jh)
            if inverse:
                ce.xor(encoding);ce.xor(linear([B1],je))
            else:
                cr=linear([B2,q2],jh)+cm.prod_rt(jh,jb,B1)
                m.lin_rows_terms(cr,math.rows(1,2),jb,jh)
                ce.mac(cr)
                ce.mac(linear([B3,q3],je)+cm.prod_rt(je,jb,q2))
        E.x(flag)
        # Off the ordinary branch both the scale target and pivot are zero;
        # armed(0)=1 makes this unselected scale act as the identity there.
        lookahead.jacobi.cf.scale(m,E,[je],jh,inverse=not inverse)
        E.x(flag)
        with selected(m,E,flag) as (cm,ce):
            if inverse:
                ce.free(linear([B3,q3],je)+cm.prod_rt(je,jb,q2))
                cr=linear([B2,q2],jh)+cm.prod_rt(jh,jb,B1)
                m.lin_rows_terms(cr,math.rows(1,2),jb,jh);ce.free(cr)
            else:
                ce.xor(linear([B1],je))
                encoding=cm.prod_rt(jh,jb,je);m.lin_rows_terms(encoding,math.rows(1,2),jb,jh);ce.xor(encoding)
        E.x(flag)
        if not inverse:gap()
        quotient(True,clear_head=inverse)


def kernel(m,E,math,P,Q,flag,coefficients,inverse=False):
    dn,t,u,v,s=coefficients;length=len(P)
    swaps=list(zip(sum(P,[]),sum(Q,[])))
    def numerator():
        rt=[]
        for i,target in enumerate(P):
            rt+=m.prod_rt(target,dn,Q[i])
            if i+1<length:rt+=selected_linear(flag,1,[Q[i+1]],target)
        E.xor(rt)
    def denominator():
        rt=[]
        if length==1:rt=m.prod_rt(Q[0],u,P[0])
        elif length:
            targets=[[] for _ in range(length+1)]
            for i,target in enumerate(Q):targets[length-1-i]=[target]
            rt=m.polyprod_rt([u,v],list(reversed(P)),targets)
        for i,target in enumerate(Q):
            if i+2<length:rt+=selected_linear(flag,0,[P[i+2]],target)
        E.xor(rt)
    if inverse:
        lookahead.jacobi.cf.scale(m,E,P,s,True);denominator()
        lookahead.jacobi.cf.scale(m,E,Q,t,True);numerator();E.swaps(flag,swaps)
    else:
        E.swaps(flag,swaps);numerator();lookahead.jacobi.cf.scale(m,E,Q,t)
        denominator();lookahead.jacobi.cf.scale(m,E,P,s)


def fold(m,E,math,W,inverse=False):
    H=W[:4];P,Q=W[4::2],W[5::2]
    names=[f'shared_boundary_{i}' for i in range(3)]
    temp=[E.A.open(name,12) for name in names] if inverse else []
    def measure_boundary():
        masks={q:E.b.hmr(q) for word in temp for q in word}
        for name,word in zip(names,temp):E.A.close(name,word)
        names.clear();return masks
    with lookahead.zero_head(m,E,H[0]) as flag:
        if inverse:
            original_head(m,E,math,W,flag,temp)
            for old,new in zip(H[1:],temp):E.swap_fn(None,list(zip(old,new)))()
            masks=measure_boundary()
        with ExitStack() as stack:
            coefficients=[stack.enter_context(m.opened(E,f'shared_parameter_{i}',12)) for i in range(5)]
            params(m,E,math,H,flag,coefficients)
            if not inverse:
                temp.extend(E.A.open(name,12) for name in names)
                boundary(m,E,math,H,P,Q,flag,temp)
                for old,new in zip(H[1:],temp):E.swap_fn(None,list(zip(old,new)))()
                masks=measure_boundary()
            kernel(m,E,math,P,Q,flag,coefficients,inverse)
            if inverse:params(m,E,math,H,flag,coefficients,True)
            else:
                parameter_masks={q:E.b.hmr(q) for word in coefficients for q in word}
                stack.close()
        if inverse:
            repay(E,boundary_terms(m,math,H,P,Q,flag,temp),masks)
        else:
            shadow=[E.A.open(f'shared_shadow_{i}',12) for i in range(3)]
            original_head(m,E,math,W,flag,shadow)
            params(m,E,math,[H[0],*shadow],flag,coefficients,masks=parameter_masks)
            for old,new in zip(temp,shadow):
                for oldq,newq in zip(old,new):E.b.z(newq,cond=masks[oldq])
            original_head(m,E,math,W,flag,shadow,True)
            for i,word in enumerate(shadow):E.A.close(f'shared_shadow_{i}',word)
    for name,word in zip(names,temp):E.A.close(name,word)


def synthesize(m,E,W,inverse=False):
    if len(W)==2:return W[::2],W[1::2]
    math=lookahead.Arithmetic(m,E)
    if not inverse and len(W)>4:synthesize(m,E,W[4:])
    fold(m,E,math,W,inverse)
    if inverse and len(W)>4:synthesize(m,E,W[4:],True)
    return W[::2],W[1::2]
