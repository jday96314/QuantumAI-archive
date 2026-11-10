"""Undo core monomial completion padding, then pad by squarefree factors.

Distinct irreducible quadratics and at most one x keep the padded locator
squarefree. The cofactor pair may share those known factors: the later
Sylvester map requires squarefreeness of the locator, not coprimeness of this
pair. H is multiplied by the same known padding polynomial and can be divided
afterward using controlled linear polynomial maps.
"""
from compact_euclid import product
from quantum_coprime_padding import multiply, permute, apply as emit_ops,xor_constant,add_constant
from coprime_padding import QUADRATICS


def frames(words):
    bp, b, c = [], [words[-2]], words[:-2]
    out = []
    for k in range(1, len(words)//2):
        beta, alpha = list(reversed(c))[:2]
        out.append((k, bp, b, beta, alpha))
        bp, b, c = b, bp+[beta, alpha], c[:-2]
    return out, (bp, b)


def fake_undo_ops(bp, b, beta, alpha):
    ops = []
    for i, word in enumerate(bp):
        if i:
            ops.extend(('cx', x, y) for x, y in zip(b[i-1], word))
    if len(b) > 1:
        ops.extend(('cx', x, y) for x, y in zip(b[-2], beta))
    ops.extend(('cx', x, y) for x, y in zip(b[-1], alpha))
    ops.append(('x', beta[0]))  # the recovered dummy beta=1 becomes zero
    return ops


def common_padding_ops(bp, b, q, quadratic, constant=None):
    ops = []
    if quadratic:
        P = [constant, 1, 1]
        cp, cb = [q[0], q[1], *bp], [q[2], q[3], *b]
        multiply(ops, cp, P)
        multiply(ops, cb, P)
        permute(ops, cp+cb, bp+q[:2]+b+q[2:])
    else:
        cp, cb = [q[1], *bp], [q[0], *b]
        # Multiplication by x is just the embedding of a zero remainder.
        permute(ops, cp+cb, b+bp+q)
    return ops


def apply(m, E, words, flags, inverse=False):
    assert (len(words)//2) % 2 == 1, 'this schedule pads to an odd locator degree'
    fs, final = frames(words)

    def unroll(inv):
        for k, bp, b, beta, alpha in (fs if inv else reversed(fs)):
            ops=common_padding_ops(bp,b,[beta,alpha],False)
            emit_ops(E.b,ops,flags[k],not inv)

    def pad(inv):
        for index in (reversed(range(len(fs))) if inv else range(len(fs))):
            k, bp, b, beta, alpha = fs[index]
            if k % 2:
                # A two-step common factor acts on the two untouched zero
                # remainder words for each cofactor.
                q = [beta, alpha, fs[index+1][3], fs[index+1][4]]
                ops = common_padding_ops(bp, b, q, True, QUADRATICS[(k-1)//2])
                emit_ops(E.b, ops, flags[k], inv)
            else:
                # Only the first dummy step on an even completion degree
                # requires the one linear x. Later flags are already set.
                with m.opened(E, 'deferred_linear_control', 1) as (control,):
                    rt = m.one_rt((flags[k],), 0, (flags[k-1],), 1, control)
                    E.mac(rt)
                    ops = common_padding_ops(bp, b, [beta, alpha], False)
                    emit_ops(E.b, ops, control, inv)
                    E.free(rt)
    if inverse:
        pad(True)
        unroll(True)
    else:
        unroll(False)
        pad(False)
    return final


def padding_polynomial(weight, maximum):
    P = [1]
    degree = weight
    if (maximum-weight) % 2:
        P = [0, 1]
        degree += 1
    while degree < maximum:
        P = product(P, [QUADRATICS[(degree-1)//2], 1, 1])
        degree += 2
    return P


def direct_repad(m,E,words,flags,inverse=False):
    """Replace common x pairs directly by irreducible quadratics.

    If y=x^d f, multiplication by 1+x^-1+b*x^-2 replaces one x²
    with x²+x+b. On the fixed-degree monic coefficient banks this is an
    upper triangular affine map. Under flag[k], the still-unreplaced common
    x-adic factor has degree at least n-k, allowing a shortened target range.
    No cofactor unrolling, degree routing, or separate linear x is required.
    """
    n=len(words)//2
    assert n%2==1
    _,(bp,b)=frames(words)
    stages=list(range(1,n,2))
    for k in (reversed(stages) if inverse else stages):
        ops=[];constant=QUADRATICS[(k-1)//2]
        for bank in [bp,b]:
            for i in range(max(0,n-k-2),len(bank)):
                for shift,c in [(1,1),(2,constant)]:
                    src=i+shift
                    if src<len(bank):xor_constant(ops,c,bank[src],bank[i])
                    elif src==len(bank):add_constant(ops,c,bank[i])
        emit_ops(E.b,ops,flags[k],inverse)
    return bp,b
