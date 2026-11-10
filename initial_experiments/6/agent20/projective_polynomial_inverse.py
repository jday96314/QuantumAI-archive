"""Projective quotient reversal, preserving roots instead of a monic locator.

For F=c K(q), G=c K(tail q), reverse q and invert c. This maps
(F,G) to (F/c^2, G^-1 mod F), on regular Euclid charts. Its middle
operation is an involution; an explicit inverse covers armed/singular inputs.
The scale of F is harmless for a root-conditioned phase. Common x padding
uses the real quotient prefix and retains an explicit leading locator word.
Prior arithmetic is imported read-only. This is not a complete submission.
"""
from contextlib import contextmanager
import polynomial_inverse_bilinear as base
from low_workspace_inverse import bind as low_bind


def bind(m,limit=1004):
    low=low_bind(m)
    @contextmanager
    def retained(E,src):
        # Only used by source-preserving quotient division. Reciprocal swaps
        # below call low directly; retained intermediates cannot survive them.
        if E.b.active_qubits+66<=limit:
            with m.opened(E,'projective_retained_inverse',12) as out:
                with m.chain_live(E,src,out):yield out
        else:
            with low(E,src) as out:yield out
    return retained


def invert_scalar(m,E,c):
    from reciprocal_inplace import apply
    apply(m,E,c)


def reverse_heads(E,heads):
    for i in range(len(heads)//2):
        for a,b in zip(heads[i],heads[-1-i]):E.swaps(None,list(zip(a,b)))


def apply(m,E,F,G,inverse=False):
    heads,c=base.extract(m,E,F,G)
    reverse_heads(E,heads);invert_scalar(m,E,c)
    base.extract(m,E,F,G,True)


def unary_count(E,n,flags,counter):
    # Valid flags[j]=[j>=d]. Binary n-d is a linear function of this unary
    # encoding. On arbitrary flags this still defines a reversible XOR map.
    for k,q in enumerate(counter):
        for j,flag in enumerate(flags):
            if (((n-j)^(n-j-1))>>k)&1:E.cx(flag,q)


def middle(m,E,heads,c,flags,counter,inverse=False):
    n=len(heads)
    def route(back=False):
        with m.opened(E,'projective_degree_match',1) as work:
            bit=work[0]
            for j in (reversed(range(n)) if back else range(n)):
                E.cx(flags[j],bit)
                if j:E.cx(flags[j-1],bit)
                E.swaps(bit,list(zip(c,heads[j][1])))
                if j:E.cx(flags[j-1],bit)
                E.cx(flags[j],bit)
    def orient(back=False):
        # A partial shift puts reversed real heads in their prefix, and leaves
        # the sole fake scalar in a small set of suffix positions.
        banks=[[h[k] for h in heads] for k in (0,1)]
        if back:reverse_heads(E,heads)
        for k in (reversed(range(len(counter))) if back else range(len(counter))):
            distance=1<<k
            for bank in banks:
                indices=range(max(0,n-distance)) if back else reversed(range(max(0,n-distance)))
                for i in indices:E.swaps(counter[k],list(zip(bank[i],bank[i+distance])))
        if not back:reverse_heads(E,heads)
    def park(back=False):
        # Highest set bit l of m gives fake position n-1-(n mod l).
        # Merge equal positions; each control is XOR-linear in unary flags.
        groups={}
        for k in range(n.bit_length()):
            l=1<<k
            if l>n:continue
            j=n-1-(n%l);toggles=groups.setdefault(j,set())
            for index in [n-l,*([n-2*l] if 2*l<=n else [])]:
                if index in toggles:toggles.remove(index)
                else:toggles.add(index)
        with m.opened(E,'projective_fake_scalar_selector',1) as temp:
            bit=temp[0];items=list(groups.items())
            for j,indices in (reversed(items) if back else items):
                for i in sorted(indices):E.cx(flags[i],bit)
                E.swaps(bit,list(zip(c,heads[j][1])))
                for i in reversed(sorted(indices)):E.cx(flags[i],bit)
    if inverse:
        route(True);invert_scalar(m,E,c);park(True);orient(True)
    else:
        orient();park();invert_scalar(m,E,c);route()


def padded_apply(m,E,F,G,flags,inverse=False):
    import padded_polynomial_inverse as pad
    heads,c=pad.layout(F,G);n=len(heads)
    with m.opened(E,'projective_completion_count',n.bit_length()) as counter:
        unary_count(E,n,flags,counter)
        base.extract(m,E,F,G)
        middle(m,E,heads,c,flags,counter,inverse)
        base.extract(m,E,F,G,True)
        unary_count(E,n,flags,counter)


def gcd_constant(m,f,g):
    # Independent ordinary polynomial EEA on an unpadded supported chart.
    f=list(f);g=list(g)
    while g:
        if not g[-1]:raise ValueError('interior degree gap')
        r=f[:];gi=m.ginv(g[-1])
        for j in reversed(range(len(f)-len(g)+1)):
            q=m.gmul(r[j+len(g)-1],gi)
            for i,v in enumerate(g):r[j+i]^=m.gmul(q,v)
        r=r[:len(g)-1]
        if r and not r[-1]:raise ValueError('interior degree gap')
        f,g=g,r
    return f[0]
