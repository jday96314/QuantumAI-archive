"""Leading-zero chart for a two-pair in-register Jacobi block.

Generic head frame is (a,b,h=bd+c,e=b+d). If a=0, encode raw
(0,u,v,w) as (0,u,h=w/u+(v/u)^2,e=v/u). The gap fraction is
f=u*x/(1+e*x+h*x^2+x^3*g). A zero second Jacobi pivot is still unsupported.
"""
from contextlib import contextmanager
import compact_jacobi as compact
from quantum_control import selected
from variable_shift_inverse import apply as variable_inverse

Arithmetic=compact.Arithmetic
jacobi=compact.jacobi


@contextmanager
def zero_head(m,E,a):
    with m.opened(E,'lookahead_zero_head',1) as flag:
        with m.eq_test(E,a,0) as test:E.cx(test,flag[0])
        yield flag[0]
        with m.eq_test(E,a,0) as test:E.cx(test,flag[0])


def head(m,E,H,flag,inverse=False):
    a,b,c,d=H
    E.x(flag)
    with selected(m,E,flag) as (cm,ce):
        ar=Arithmetic(cm,ce)
        if inverse:
            for x,y in zip(b,d):ce.cx(x,y)
            ce.xor(cm.prod_rt(c,b,d));jacobi.extract(cm,ce,ar,H,True)
        else:
            jacobi.extract(cm,ce,ar,H)
            ce.xor(cm.prod_rt(c,b,d))
            for x,y in zip(b,d):ce.cx(x,y)
    E.x(flag)
    with selected(m,E,flag) as (cm,ce):
        def square():
            rt=[];m.lin_rows_terms(rt,Arithmetic(m,E).rows(1,2),c,d);ce.xor(rt)
        if inverse:ce.swap_fn(None,list(zip(c,d)))();square()
        jacobi.cf.scale(cm,ce,[c,d],b,inverse=not inverse)
        if not inverse:square();ce.swap_fn(None,list(zip(c,d)))()


@contextmanager
def parameters(m,E,math,H,flag):
    a,b,h,e=H;d=compact.XorWord(e,b)
    def product(left,right,target):
        rt=[];m.bil_terms(rt,m.KT,math.spec([(1,left,1)]),math.spec([(1,right,1)]),tgt=target);return rt
    with m.opened(E,'lookahead_gap_b',12) as bg:
        br=[m.norm_term((flag,),0,(q,),0,(t,)) for q,t in zip(b,bg)];E.mac(br)
        with m.opened(E,'lookahead_p',12) as p,m.opened(E,'lookahead_P0',12) as p0,m.opened(E,'lookahead_P1',12) as p1:
            def prepare(clear=False):
                with m.opened(E,'lookahead_original_c',12) as c:
                    cr=product(b,e,c)+[('L',(x,),0,(y,)) for x,y in zip(h,c)]
                    m.lin_rows_terms(cr,math.rows(1,2),b,c);E.mac(cr)
                    pr=product(a,c,p)+[('L',(x,),0,(y,)) for x,y in zip(bg,p)]
                    (E.free if clear else E.mac)(pr)
                    with m.opened(E,'lookahead_aux',12) as aux,m.opened(E,'lookahead_second_aux',12) as other:
                        rr=[];m.bil_terms(rr,m.KT,math.spec([(1,c,1)]),math.spec([(1,e,2)]),tgt=aux)
                        for word,q in [(b,4),(c,2)]:m.lin_rows_terms(rr,math.rows(1,q),word,aux)
                        E.mac(rr);qr=math.cube_rts(m.spec_plain(e),lambda bits:[other[j] for j in bits]);E.mac(qr)
                        r0=product(a,aux,p0)+product(bg,other,p0)
                        (E.free if clear else E.mac)(r0);E.free(qr);E.free(rr)
                        rr=product(c,d,aux)+math.cube_rts(m.spec_plain(b),lambda bits:[aux[j] for j in bits]);E.mac(rr)
                        sr=product(a,aux,other)
                        spec=[x^frozenset(y) for x,y in zip(m.spec_plain(h),math.squared_spec(m.spec_plain(e)))]
                        m.bil_terms(sr,m.KT,m.spec_plain(bg),spec,tgt=other);E.mac(sr)
                        (E.free if clear else E.mac)(product(h,other,p1));E.free(sr);E.free(rr)
                    E.free(cr)
            prepare();yield [(),e,h],[(),compact.XorWord(b,bg)],p,(p0,p1);prepare(True)
        E.free(br)


def block(m,E,math,W,inverse=False):
    H,R=W[:4],W[4:]
    with zero_head(m,E,H[0]) as flag:
        if not inverse:head(m,E,H,flag)
        if R:
            with parameters(m,E,math,H,flag) as (A,D,p,P):
                def shift():
                    for src,dst in zip(P,R):
                        for x,y in zip(src,dst):E.cx(x,y)
                if inverse:
                    math.triangular(A,1,R,True);math.triangular(D,1,R);variable_inverse(m,E,math,R,flag)
                    jacobi.cf.scale(m,E,R,p);math.triangular(D,1,R,True);shift();math.triangular(A,1,R,True)
                else:
                    math.triangular(A,1,R);shift();math.triangular(D,1,R)
                    jacobi.cf.scale(m,E,R,p,True);variable_inverse(m,E,math,R,flag)
                    math.triangular(D,1,R,True);math.triangular(A,1,R)
        if inverse:head(m,E,H,flag,True)


def extract(m,E,math,W,inverse=False):
    for i in (reversed(range(0,len(W),4)) if inverse else range(0,len(W),4)):
        if len(W)-i==2:jacobi.pair(m,E,math,W[i:],inverse)
        else:block(m,E,math,W[i:],inverse)


def fold_pair(m,E,W,inverse=False):
    if len(W)==2:return
    a,b=W[:2];P,Q=W[2::2],W[3::2]
    def shear():
        rt=[]
        for i,(p,q) in enumerate(zip(P,Q)):
            rt+=m.prod_rt(p,b,q)
            if i+1<len(Q):rt += [('L',(x,),0,(y,)) for x,y in zip(Q[i+1],p)]
        E.xor(rt)
    swaps=list(zip(sum(P,[]),sum(Q,[])))
    if inverse:
        jacobi.cf.scale(m,E,P,a,True)
        for x,y in zip(b,P[0]):E.cx(y,x)
        E.swap_fn(None,swaps)();shear()
    else:
        shear();E.swap_fn(None,swaps)()
        for x,y in zip(b,P[0]):E.cx(y,x)
        jacobi.cf.scale(m,E,P,a)


def frame(m,E,H,flag,inverse=False):
    _,b,c,d=H;E.x(flag)
    with selected(m,E,flag) as (cm,ce):
        if inverse:
            ce.xor(cm.prod_rt(c,b,d))
            for x,y in zip(b,d):ce.cx(x,y)
        else:
            for x,y in zip(b,d):ce.cx(x,y)
            ce.xor(cm.prod_rt(c,b,d))
    E.x(flag)
    with selected(m,E,flag) as (_,ce):ce.swap_fn(None,list(zip(c,d)))()


def gap_fold(m,E,W,inverse=False):
    _,u,e,h=W[:4];P,Q=W[4::2],W[5::2]
    def shear():
        rt=[]
        for i,p in enumerate(P):
            rt+=m.prod_rt(p,h,Q[i])
            if i+1<len(Q):rt+=m.prod_rt(p,e,Q[i+1])
            if i+2<len(Q):rt += [('L',(x,),0,(y,)) for x,y in zip(Q[i+2],p)]
        E.xor(rt)
    def change_e():
        if Q:
            for x,y in zip(Q[0],e):E.cx(x,y)
    def change_h():
        if Q:E.xor(m.prod_rt(h,e,Q[0]))
        if len(Q)>1:
            for x,y in zip(Q[1],h):E.cx(x,y)
    swaps=[*zip(u,e),*zip(sum(P,[]),sum(Q,[]))]
    if inverse:
        E.swap_fn(None,swaps)();jacobi.cf.scale(m,E,Q,u,True)
        change_e();change_h();shear()
    else:
        shear();change_h();change_e();jacobi.cf.scale(m,E,Q,u)
        E.swap_fn(None,swaps)()


def synthesize(m,E,W,inverse=False):
    if len(W)==2:return W[::2],W[1::2]
    if not inverse:synthesize(m,E,W[4:]) if len(W)>4 else None
    with zero_head(m,E,W[0]) as flag:
        if not inverse:frame(m,E,W[:4],flag)
        E.x(flag)
        with selected(m,E,flag) as (cm,ce):
            if inverse:fold_pair(cm,ce,W,True);fold_pair(cm,ce,W[2:],True)
            else:fold_pair(cm,ce,W[2:]);fold_pair(cm,ce,W)
        E.x(flag)
        with selected(m,E,flag) as (cm,ce):gap_fold(cm,ce,W,inverse)
        if inverse:frame(m,E,W[:4],flag,True)
    if inverse:synthesize(m,E,W[4:],True) if len(W)>4 else None
    return W[::2],W[1::2]


if __name__=='__main__':
    import json
    from pathlib import Path
    from import_published import load
    from low_workspace_inverse import bind
    import dqi_budget_probe as probe
    m=load();m.inv_word=bind(m);b=m.Builder();b._items=probe.count_module().Builder()._items
    W=[b.qubits(12) for _ in range(70)]
    class Counter(m.Em):
        extra=0
        def emit_mac(self,rt):self.extra+=sum(t[0]=='P' for t in rt)
        def emit_free(self,rt):pass
    E=Counter(b,*b.qubits(3),m.Alloc(b));math=Arithmetic(m,E)
    extract(m,E,math,W);report=dict(Q=b.peak_qubits,extract_T=b._cur_stats['toffoli']+E.extra,
        scope='leading-zero-chart extraction only; synthesis and complete decoder are not emitted')
    print(json.dumps(report),flush=True);Path('results/lookahead_jacobi_extract_profile.json').write_text(json.dumps(report,indent=2)+'\n')
