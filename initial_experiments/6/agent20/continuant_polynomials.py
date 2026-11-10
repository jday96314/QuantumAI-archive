"""Exact characteristic-two transfer polynomials for blocked C fractions.

Polynomials use exponent tuples and XOR addition. The compiled arithmetic DAG
counts only field products; Frobenius maps and additions are binary linear.
This is an algebra/resource experiment, not a decoder acceptance claim.
"""
from functools import lru_cache
import argparse,json
from pathlib import Path


def add(*values):
    out=set()
    for value in values:out.symmetric_difference_update(value)
    return frozenset(out)


def mul(a,b):
    out=set()
    for x in a:
        for y in b:
            z=tuple(i+j for i,j in zip(x,y))
            if z in out:out.remove(z)
            else:out.add(z)
    return frozenset(out)


def derive(pairs):
    n=2*pairs;zero=frozenset();one=frozenset({(0,)*n})
    variables=[frozenset({tuple(int(i==j) for j in range(n))}) for i in range(n)]
    P,Q=[variables[-1]],[one]
    for c in reversed(variables[:-1]):
        P,Q=[mul(c,v) for v in Q],[add(Q[i] if i<len(Q) else zero,P[i-1] if i else zero) for i in range(max(len(Q),len(P)+1))]
    H=[]
    for i in range(n):H.append(add(P[i] if i<len(P) else zero,*[mul(Q[j],H[i-j]) for j in range(1,min(i+1,len(Q)))]))
    tail=[add(*[mul(Q[i],H[n+j-i]) for i in range(j+1,len(Q))]) for j in range(pairs)]
    determinant=one
    for v in variables:determinant=mul(determinant,v)
    previous=[one]
    old=[one]
    for c in variables[1:-1]:
        previous,old=[add(previous[i] if i<len(previous) else zero,mul(c,old[i-1]) if i else zero) for i in range(max(len(previous),len(old)+1))],previous
    return dict(A=Q,D=previous,p=determinant,P=tail)


def derive_jacobi(pairs):
    n=2*pairs;zero=frozenset();one=frozenset({(0,)*n})
    variables=[frozenset({tuple(int(i==j) for j in range(n))}) for i in range(n)]
    def fraction(words):
        if len(words)==2:return [words[0]],[one,words[1]]
        P,Q=fraction(words[2:]);a,b=words[:2]
        numerator=[mul(a,v) for v in Q]
        denominator=[add(Q[i] if i<len(Q) else zero,mul(b,Q[i-1]) if i and i-1<len(Q) else zero,P[i-2] if i>=2 and i-2<len(P) else zero) for i in range(max(len(Q)+1,len(P)+2))]
        return numerator,denominator
    P,A=fraction(variables);D=fraction(variables[:-2])[1] if pairs>1 else [one]
    H=[]
    for i in range(n):H.append(add(P[i] if i<len(P) else zero,*[mul(A[j],H[i-j]) for j in range(1,min(i+1,len(A)))]))
    tail=[add(*[mul(A[i],H[n+j-i]) for i in range(j+1,len(A))]) for j in range(pairs)]
    determinant=one
    for v in variables[::2]:determinant=mul(determinant,v)
    return dict(A=A,B=P,D=D,p=determinant,P=tail)


def substitute(polynomial,images):
    n=len(images);one=frozenset({(0,)*n});out=frozenset()
    for exponent in polynomial:
        term=one
        for i,e in enumerate(exponent):
            shift=0
            while e:
                if e&1:term=mul(term,frozenset(tuple(v*(1<<shift) for v in monomial) for monomial in images[i]))
                shift+=1;e>>=1
        out=add(out,term)
    return out


def coordinates(polynomials,pairs,kind):
    n=2*pairs;var=[frozenset({tuple(int(i==j) for j in range(n))}) for i in range(n)]
    images=var.copy()
    if kind=='pairs':
        for i in range(2,n,2):images[i]=add(var[i],var[i-1])
    elif kind=='prefix':
        for i in range(2,n):images[i]=add(var[i],var[i-1])
    elif kind!='original':raise ValueError(kind)
    return {name:[substitute(p,images) for p in value] if isinstance(value,list) else substitute(value,images) for name,value in polynomials.items()}


class Compiler:
    def __init__(self,n):self.n=n;self.nodes=[];self.intern={};self.operations=[]
    def node(self,kind,*args):
        key=(kind,*args)
        if key in self.intern:return self.intern[key]
        i=len(self.nodes);self.nodes.append(key);self.intern[key]=i
        products={i} if kind=='mul' else set()
        for j in (args if kind in ('add','mul') else args[:1] if kind=='square' else ()):
            products.update(self.operations[j])
        self.operations.append(frozenset(products));return i
    def plus(self,ids):
        values=set()
        for i in ids:
            if self.nodes[i][0]=='zero':continue
            for j in self.nodes[i][1:] if self.nodes[i][0]=='add' else (i,):
                if j in values:values.remove(j)
                else:values.add(j)
        if not values:return self.node('zero')
        if len(values)==1:return next(iter(values))
        return self.node('add',*sorted(values))
    def times(self,a,b):
        if self.nodes[a][0]=='zero' or self.nodes[b][0]=='zero':return self.node('zero')
        if self.nodes[a][0]=='one':return b
        if self.nodes[b][0]=='one':return a
        if a==b:return self.node('square',a)
        return self.node('mul',*sorted((a,b)))
    def varpower(self,var,exponent):
        if exponent==0:return self.node('one')
        out=self.node('var',var);power=1
        while power*2<=exponent:out=self.node('square',out);power*=2
        return out if power==exponent else self.times(out,self.varpower(var,exponent-power))
    @lru_cache(None)
    def compile(self,poly):
        if not poly:return self.node('zero')
        if all(not any(v) for v in poly):return self.node('one')
        if all(not any(e%2 for e in v) for v in poly):return self.node('square',self.compile(frozenset(tuple(e//2 for e in v) for v in poly)))
        common=tuple(min(v[j] for v in poly) for j in range(self.n))
        if any(common):
            rest=self.compile(frozenset(tuple(e-c for e,c in zip(v,common)) for v in poly))
            factors=[self.varpower(j,c) for j,c in enumerate(common) if c]
            # Larger shared expressions are left intact; monomial factors
            # are individually Frobenius-linear when their exponent is a power of two.
            for factor in factors:rest=self.times(rest,factor)
            return rest
        best=None
        for var in range(self.n):
            if not any(v[var] for v in poly):continue
            groups={}
            for v in poly:
                e=v[var];w=(*v[:var],0,*v[var+1:]);groups.setdefault(e,set()).add(w)
            terms=[self.times(self.varpower(var,e),self.compile(frozenset(group))) for e,group in groups.items()]
            result=self.plus(terms);score=len(self.operations[result])
            if best is None or score<best[0]:best=(score,result)
        return best[1]


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--pairs',type=int,nargs='+',default=[2,3,4]);args=p.parse_args();reports=[]
    for pairs in args.pairs:
        polynomials=derive(pairs);compiler=Compiler(2*pairs)
        roots={name:[compiler.compile(v) for v in (value if isinstance(value,list) else [value])] for name,value in polynomials.items()}
        products=set().union(*(compiler.operations[i] for root in roots.values() for i in root))
        row=dict(pairs=pairs,terms={k:[len(v) for v in (val if isinstance(val,list) else [val])] for k,val in polynomials.items()},products=len(products),roots=roots,nodes=len(compiler.nodes))
        reports.append(row);print(json.dumps(row),flush=True)
    Path('results/continuant_polynomial_profiles.json').write_text(json.dumps(reports,indent=2)+'\n')
