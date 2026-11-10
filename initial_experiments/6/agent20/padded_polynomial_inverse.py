"""Common-monomial generalized inverse via a quantum-selected quotient prefix.

Promise: F=x^m Ftrue (degree n), G=x^m Gtrue (degree n-1), true degree
d=n-m. Unary flags[j]=[j>=d]. Generic active quotient pivots are nonzero.
Extract all n quotient pairs unconditionally. First fake pair is (0,c),
all later pairs and terminal constant are zero. Reverse only the real prefix
and perform the continuant normalization, obtaining x^m/Gtrue mod Ftrue.
An explicit inverse is provided for every input, including unsupported ones.
"""
from contextlib import contextmanager
import polynomial_inverse_bilinear as base
from compact_monomial_unpad import degree_count,rotate


def layout(F,G):
    heads=[]
    while G:
        heads.append((F[-2],F[-1]));F,G=G,F[:-2]
    return heads,F[0]


@contextmanager
def odd_degree(m,E,n,counter):
    with m.opened(E,'inverse_degree_odd',1) as temp:
        bit=temp[0]
        if n%2:E.b.x(bit)
        E.cx(counter[0],bit)
        yield bit
        E.cx(counter[0],bit)
        if n%2:E.b.x(bit)


def middle(m,E,heads,c,flags,counter,inverse=False):
    n=len(heads)
    def orient(back=False):
        banks=[[h[j] for h in heads] for j in (0,1)]
        class WordBuilder:
            def cswap(self,control,a,b):E.swaps(control,list(zip(a,b)))
        def rev():
            for bank in banks:
                for i in range(n//2):E.swaps(None,list(zip(bank[i],bank[-1-i])))
        if back:rev()
        for k in (reversed(range(len(counter))) if back else range(len(counter))):
            for bank in banks:rotate(WordBuilder(),bank,counter[k],1<<k,not back)
        if not back:rev()
    def select_c():E.swaps(flags[-1],list(zip(heads[-1][1],c)))
    def scale(back=False):
        with m.opened(E,'padded_inverse_c_squared',12) as c2:
            rt=[];rows=tuple(sum(((m.gmul(1<<j,1<<j)>>i)&1)<<j for j in range(12)) for i in range(12))
            m.lin_rows_terms(rt,rows,c,c2);E.mac(rt)
            for parity in ((1,0) if back else (0,1)):
                base.cf.scale(m,E,[w for pair in heads[parity::2] for w in pair],c2,bool(parity)^back)
            E.free(rt)
    def invert_c():
        with odd_degree(m,E,n,counter) as odd:
            with m.inv_word(E,c) as ci:E.swaps(odd,list(zip(c,ci)))
    def route_c(back=False):
        with m.opened(E,'padded_inverse_degree_match',1) as work:
            bit=work[0]
            for j in (reversed(range(n)) if back else range(n)):
                flag=flags[j]
                E.cx(flag,bit)
                if j:E.cx(flags[j-1],bit)
                E.swaps(bit,list(zip(c,heads[j][1])))
                if j:E.cx(flags[j-1],bit)
                E.cx(flag,bit)
    if inverse:
        route_c(True);invert_c();scale(True);select_c();orient(True)
    else:
        orient();select_c();scale();invert_c();route_c()


def apply(m,E,F,G,flags,inverse=False):
    heads,c=layout(F,G)
    with m.opened(E,'padded_inverse_completion_count',len(flags).bit_length()) as counter:
        degree_count(E.b,flags,counter)
        base.extract(m,E,F,G)
        middle(m,E,heads,c,flags,counter,inverse)
        base.extract(m,E,F,G,True)
        degree_count(E.b,flags,counter)
