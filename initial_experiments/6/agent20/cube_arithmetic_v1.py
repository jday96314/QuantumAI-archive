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
        if list(map(frozenset,right))==list(map(frozenset,self.squared_spec(left))):
            operand=left
        elif list(map(frozenset,left))==list(map(frozenset,self.squared_spec(right))):
            operand=right
        else:
            return super().product(left,right,targets)
        self.cube_products+=1
        m=self.m;rt=[]
        for a,b,out in self.cube_terms:
            term=m.norm_term(m.fold(operand,m.bitlist(a,12)),0,
                             m.fold(operand,m.bitlist(b,12)),0,targets(m.bitlist(out,12)))
            if term:rt.append(term)
        for j,out in enumerate(self.cube_linear):
            if out and operand[j]:rt.append(('L',tuple(sorted(operand[j])),0,targets(m.bitlist(out,12))))
        self.E.emit_mac(rt)

    def self_middle(self,A,q,B):
        if q!=2:return super().self_middle(A,q,B)
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
