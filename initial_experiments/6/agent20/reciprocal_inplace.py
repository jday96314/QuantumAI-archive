"""66-AND in-place GF4096 reciprocal by transporting the retained norm chain.

The ordinary retained inverse cannot survive swapping its source/reciprocal.
Here the chain is changed explicitly: N <-> N^-1, M -> M^-1=(M^3)^2,
M^3 -> M^-3=M^4 in GF8. Every transport is Clifford. The existing
measurement eraser then sees exactly the chain of the new reciprocal source.
Apache-2.0 field formulas and emitter are imported read-only.
"""
from tower_inverse import map_ops


def linear(m,E,word,rows,control=None):
    cols=[sum(((rows[i]>>j)&1)<<i for i in range(len(rows))) for j in range(len(rows))]
    for op,i,j in map_ops(cols):
        if op=='cx':
            if control is None:E.cx(word[i],word[j])
            else:E.b.ccx(control,word[i],word[j])
        else:E.swaps(control,[(word[i],word[j])])


def apply(m,E,src,control=None,before_swap=None):
    with m.opened(E,'inplace_reciprocal_output',12) as out:
        wins,rts=m.chain_fwd(E,src,out)
        n,mm,s,y=wins
        if before_swap is not None:before_swap(src,out)
        E.swaps(control,list(zip(src,out)))
        E.swaps(control,list(zip(n,y)))
        sq=m.SQ3R
        fourth=tuple(sum(((sum(((sq[i]>>k)&1)*((sq[k]>>j)&1) for k in range(3))%2))<<j for j in range(3)) for i in range(3))
        linear(m,E,s,sq,control);linear(m,E,mm,fourth,control)
        E.swaps(control,list(zip(mm,s)))
        m.chain_undo(E,wins,rts)
