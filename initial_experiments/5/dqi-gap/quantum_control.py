"""Coherent controls for the inherited measurement-uncomputed arithmetic.

Linear operand preparation and inversion scratch are unconditional. Only the
live target action is selected. A controlled scaling uses factor 1 off-branch.
Every emitted gate is charged regardless of the branch distribution.
"""
from contextlib import contextmanager


class GateControl:
    def __init__(self,b,flag,aux):self.raw=b;self.flag=flag;self.aux=aux
    def __getattr__(self,name):return getattr(self.raw,name)
    def x(self,t,**kw):self.raw.cx(self.flag,t,**kw)
    def z(self,t,**kw):self.raw.cz(self.flag,t,**kw)
    def neg(self,**kw):self.raw.z(self.flag,**kw)
    def cx(self,c,t,**kw):self.raw.ccx(self.flag,c,t,**kw)
    def cz(self,c,t,**kw):self.raw.ccz(self.flag,c,t,**kw)
    def swap(self,a,b,**kw):self.raw.cswap(self.flag,a,b,**kw)
    def ccx(self,a,b,t,**kw):
        self.raw.ccx(self.flag,a,self.aux)
        self.raw.ccx(self.aux,b,t,**kw)
        r=self.raw.hmr(self.aux);self.raw.cz(self.flag,a,cond=r)
    def ccz(self,a,b,t,**kw):
        self.raw.ccx(self.flag,a,self.aux)
        self.raw.ccz(self.aux,b,t,**kw)
        r=self.raw.hmr(self.aux);self.raw.cz(self.flag,a,cond=r)
    def cswap(self,c,a,b,**kw):
        self.raw.cx(a,b);self.ccx(c,b,a,**kw);self.raw.cx(a,b)
    def hmr(self,t):
        self.raw.cswap(self.flag,t,self.aux)
        return self.raw.hmr(self.aux)


@contextmanager
def selected(m,E,flag):
    """Yield (field module, emitter) with an exact quantum branch control."""
    with m.opened(E,'branch_control_aux',1) as aux:
        class Emitter(m.Em):
            def __init__(self):
                super().__init__(GateControl(E.b,flag,aux[0]),E.sA,E.sB,E.sT,E.A)
            @contextmanager
            def mask_terms(self,rt):
                if rt and all(t[0]=='P' for t in rt):
                    right=list(dict.fromkeys(q for t in rt for q in t[3]))
                    if len(right)<=12:
                        with m.opened(E,'branch_mask_operand',len(right)) as bank:
                            pairs=dict(zip(right,bank))
                            terms=[t for q,w in pairs.items() for t in m.one_rt((flag,),0,(q,),0,w)]
                            E.mac(terms)
                            masked=[('P',t[1],t[2],tuple(pairs[q] for q in t[3])+((flag,) if t[4] else ()),0,t[5]) for t in rt]
                            yield masked
                            E.free(terms)
                        return
                yield None
            def fan(self,source,targets):
                if not targets:return
                first,*rest=targets
                for t in rest:E.b.cx(first,t)
                E.b.ccx(flag,source,first)
                for t in reversed(rest):E.b.cx(first,t)
            def emit_mac(self,rt):
                with self.mask_terms(rt) as masked:
                    if masked is not None:
                        E.emit_mac(masked);return
                for t in rt:
                    if t[0]=='L':
                        u,built=E._build(t[1],t[2],E.sA)
                        self.fan(u,t[-1]);E._unbuild(t[1],t[2],E.sA,built)
                    else:
                        E.fan_pair(t[1],t[2],t[3],t[4],lambda s,t=t:self.fan(s,t[-1]))
            def emit_free(self,rt):
                with self.mask_terms(rt) as masked:
                    if masked is not None:
                        E.emit_free(masked);return
                per={}
                for term in rt:
                    for target in term[-1]:per.setdefault(target,[]).append(term)
                for target,terms in reversed(list(per.items())):
                    measurement=E.b.hmr(target)
                    for term in terms:
                        if term[0]=='L':
                            for q in term[1]:E.b.cz(flag,q,cond=measurement)
                            if term[2]:E.b.z(flag,cond=measurement)
                        else:
                            E.pivot_pair(term[1],term[2],term[3],term[4],lambda a,b:E.b.ccz(flag,a,b,cond=measurement))
        CE=Emitter()
        class Module:
            branch_flag=flag
            base_module=m
            base_emitter=E
            def __getattr__(self,name):return getattr(m,name)
            def inv_word(self,ignored,source):return m.inv_word(E,source)
            def mon_slot(self,ignored,word,scale,eraser,product):
                with m.opened(E,'selected_scale',12) as ss,m.opened(E,'selected_eraser',12) as ee:
                    def copy(src,dst,back=False):
                        if not back:E.b.x(dst[0]);E.b.cx(flag,dst[0])
                        for q,t in zip(src,dst):E.b.ccx(flag,q,t)
                        if back:E.b.cx(flag,dst[0]);E.b.x(dst[0])
                    copy(scale,ss);copy(eraser,ee)
                    m.mon_slot(E,word,ss,ee,product)
                    copy(eraser,ee,True);copy(scale,ss,True)
        yield Module(),CE


def controlled_ops(E,ops,flag,inverse=False):
    """Selected affine routing, synthesized with one Toffoli per transvection."""
    from controlled_linear import apply
    apply(E.b,ops,flag,inverse)


def permutation_ops(source,target):
    from quantum_coprime_padding import permute
    ops=[];permute(ops,source,target);return ops


@contextmanager
def picked_word(m,E,source,flag):
    """A reusable field factor equal to source on-branch and one off-branch."""
    with m.opened(E,'picked_field_factor',12) as word:
        rt=[('L',(),1,(word[0],)),('L',(flag,),0,(word[0],))]
        rt += [t for q,w in zip(source,word) for t in m.one_rt((flag,),0,(q,),0,w)]
        E.mac(rt)
        yield word
        E.free(rt)
