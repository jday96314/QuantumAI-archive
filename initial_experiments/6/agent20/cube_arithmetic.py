"""Recognize Frobenius-related products in the compact series arithmetic."""
from cube_field import host_terms
from cubic_trace_field import linear
from published_dqi import isomorphism
from dqi_budget_probe import trimmed


class CubeArithmetic(trimmed.TrimmedArithmetic):
    def __init__(self,*args):
        super().__init__(*args)
        m=self.m
        forward,back,transform=isomorphism(m)
        terms,cols=host_terms()
        self.cube_terms=[]
        for a,b,out in terms:
            aa=sum(((a&transform(1<<i,back)).bit_count()%2)<<i for i in range(12))
            bb=sum(((b&transform(1<<i,back)).bit_count()%2)<<i for i in range(12))
            self.cube_terms.append((aa,bb,transform(out,forward)))
        self.cube_linear=[transform(linear(cols,transform(1<<i,back)),forward) for i in range(12)]
        self.cube_products=0

    def squared_spec(self,source):
        return [self.m.fold(source,self.m.bitlist(row,12)) for row in self.rows(1,2)]

    def product(self,left,right,targets):
        m=self.m;generic=[]
        m.bil_terms(generic,m.KT,left,right,tws_of=targets)
        if list(map(frozenset,right))==list(map(frozenset,self.squared_spec(left))):
            operand=left
        elif list(map(frozenset,left))==list(map(frozenset,self.squared_spec(right))):
            operand=right
        else:
            return self.E.emit_mac(generic)
        self.cube_products+=1
        rt=self.cube_rts(operand,targets)
        if hasattr(m,'base_emitter') and hasattr(m,'branch_flag'):
            cost=12+sum(t[0]=='P' for t in rt)
            if cost>=self.emission_cost(generic):return self.E.emit_mac(generic)
            def emit():
                base,emitter=m.base_module,m.base_emitter
                with base.opened(emitter,'selected_cube_input',12) as word:
                    masks=[]
                    for j,target in enumerate(word):
                        masks+=base.one_rt((m.branch_flag,),0,tuple(sorted(operand[j])),0,target)
                    emitter.mac(masks)
                    emitter.emit_mac(self.cube_rts(base.spec_plain(word),targets))
                    emitter.free(masks)
            if hasattr(self.E,'emit_action'):self.E.emit_action(emit,cost)
            else:emit()
        elif sum(t[0]=='P' for t in rt)<sum(t[0]=='P' for t in generic):
            self.E.emit_mac(rt)
        else:
            self.E.emit_mac(generic)

    def cube_rts(self,operand,targets):
        m=self.m;rt=[]
        for a,b,out in self.cube_terms:
            term=m.norm_term(m.fold(operand,m.bitlist(a,12)),0,
                             m.fold(operand,m.bitlist(b,12)),0,targets(m.bitlist(out,12)))
            if term:rt.append(term)
        for j,out in enumerate(self.cube_linear):
            if out and operand[j]:rt.append(('L',tuple(sorted(operand[j])),0,targets(m.bitlist(out,12))))
        return rt

    def emission_cost(self,rt):
        products=sum(t[0]=='P' for t in rt)
        if not hasattr(self.m,'branch_flag'):return products
        if rt and all(t[0]=='P' for t in rt):
            right={q for t in rt for q in t[3]}
            if len(right)<=12:return len(right)+products
        return 2*products+sum(t[0]=='L' and bool(t[-1]) for t in rt)

    def self_middle(self,A,q,B):
        if q!=2:return super().self_middle(A,q,B)
        class Collector:
            def __init__(self):self.terms=[]
            def emit_mac(self,rt):self.terms.append(rt)
            def emit_action(self,fn,cost):self.terms.append((fn,cost))
        original=self.E
        plans=[]
        for plan in (lambda:super(CubeArithmetic,self).self_middle(A,q,B),
                     lambda:self.cube_self_middle(A,q,B)):
            collected=Collector();self.E=collected
            try:plan()
            finally:self.E=original
            plans.append(collected.terms)
        terms=min(plans,key=lambda groups:sum(rt[1] if isinstance(rt,tuple) else self.emission_cost(rt) for rt in groups))
        for rt in terms:
            if isinstance(rt,tuple):rt[0]()
            else:original.emit_mac(rt)

    def cube_self_middle(self,A,q,B):
        m=self.m;start=len(A);n=start+len(B)
        pairs=[(i,j,q*j+i-start) for j in range(1,len(A)) for i in range(1,len(A))
               if start<=q*j+i<n and tuple(A[i]) not in self.known_zero
               and tuple(A[j]) not in self.known_zero]
        degree=3*(len(A)-1)
        if len(pairs)*42 <= (degree+1)*28:
            for i,j,k in pairs:
                self.product(self.spec([(1,A[i],1)]),self.spec([(1,A[j],q)]),lambda tb,k=k:[B[k][bit] for bit in tb])
        elif pairs:
            points,inverse=m.vinv(degree+1)
            for j,point in enumerate(points):
                left=self.spec([(m.gpow(point,i),A[i],1) for i in range(1,len(A))])
                right=[frozenset(v) for v in self.squared_spec(left)]
                def targets(tb,j=j):
                    return [B[k][bit] for k in range(len(B)) if start+k<=degree
                            for bit in m.bitlist(m.gmul(inverse[start+k][j],m.tmask(tb)),12)]
                self.product(left,right,targets)
        for j in range(1,len(A)):
            if tuple(A[j]) in self.known_zero:continue
            index=q*j-start
            if 0<=index<len(B):
                terms=[];m.lin_rows_terms(terms,self.rows(1,q),A[j],B[index]);self.E.emit_mac(terms)
