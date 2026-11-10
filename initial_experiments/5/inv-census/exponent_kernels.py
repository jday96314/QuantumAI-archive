"""New borrowed-headroom and Walsh exponent phase kernels; no external writes."""
from fractions import Fraction
import itertools
import numpy as np
import kickmix as km


def borrowed_phase(b, y, x, flags, *, head_bits=10, denominator=37):
    """Require low 2*s input bits zero, s=flags[1:] <= head_bits.

    Only the first s high positions are output, the rest contain shifted dirty
    input bits and must not be phased. All input/word bits are restored.
    """
    core=len(y)
    bank=km.array([*y,*reversed(list(x[:head_bits]))])
    for i,q in enumerate(flags[1:]): b.left_rotate(bank,1<<i,control=q)
    for i,q in enumerate(bank[:core]):b.z_pow(q,Fraction(2)**(i-denominator+1))
    ctrl,tmp=b.alloc_qubits(2)
    for j,q in enumerate(bank[core:]):
        b.flip_if_greater_than(flags[1:],j,target=ctrl)
        b.init_and(ctrl,q,tmp)
        b.z_pow(tmp,Fraction(2)**(core+j-denominator+1))
        b.del_and(ctrl,q,tmp)
        b.flip_if_greater_than(flags[1:],j,target=ctrl)
    b.free(km.array([ctrl,tmp]))
    for i in reversed(range(len(flags)-1)):b.right_rotate(bank,1<<i,control=flags[i+1])


def walsh_coefficients():
    """Minimum support real Walsh extension on s=11..15.

    Any sparsest solution of 11 independent real constraints can be extended
    to a nonsingular eleven-column basis. Enumerate all 4368 such bases;
    floating point screens them and exact rational arithmetic certifies best.
    This does not optimize equivalent angles modulo whole turns.
    """
    h=np.array([[(-1)**((s&m).bit_count()%2) for m in range(16)] for s in range(11)],dtype=np.int64)
    rhs=2.**np.arange(11)
    best=None;valid=0;hist={}
    for support in itertools.combinations(range(16),11):
        a=h[:,support].astype(float)
        try:v=np.linalg.solve(a,rhs)
        except np.linalg.LinAlgError:continue
        if np.max(abs(a@v-rhs))>1e-7:continue
        if np.linalg.cond(a)>1e10:continue
        valid+=1;n=int(np.count_nonzero(abs(v)>1e-7));hist[n]=hist.get(n,0)+1
        # Coefficients are rational with small denominators; no near-zero
        # result is accepted without exact reconstruction of all constraints.
        vv=[Fraction(float(t)).limit_denominator(1<<16) for t in v]
        if any(sum(int(h[s,m])*vv[i] for i,m in enumerate(support))!=1<<s for s in range(11)):
            raise ArithmeticError("rational reconstruction failed")
        score=(sum(t!=0 for t in vv),sum(abs(t) for t in vv))
        if best is None or score<best[0]:best=(score,{m:t for m,t in zip(support,vv) if t})
    return best[1],dict(bases_considered=4368,nonsingular_bases=valid,support_histogram=hist,minimum_support=best[0][0])


def walsh_phase(b,y,scale,coefficients,*,denominator=37):
    """CNOT-only synthesis, including all exponent-only phase corrections.

    p*chi(t) = (p XOR t)-t. Combine the -t corrections across all y bits.
    """
    parity=b.alloc_qubits(1)[0]
    for mask,c in coefficients.items():
        if mask==0:
            for j,q in enumerate(y):b.z_pow(q,c*Fraction(2)**(j-denominator+1))
            continue
        controls=km.array([q for i,q in enumerate(scale) if mask>>i&1])
        b.cx(controls,parity)
        for j,q in enumerate(y):
            b.cx(parity,q);b.z_pow(q,c*Fraction(2)**(j-denominator+1));b.cx(parity,q)
        b.z_pow(parity,-c*((1<<len(y))-1)*Fraction(2)**(1-denominator))
        b.cx(controls,parity)
    b.free(km.array([parity]))


def circular_phase(b,y,flags,*,denominator=37):
    """No headroom: rotate within core, repair the wrapped bit weights.

    After rotl(core,s), bit j<s has weight 2^(core+j), rather than 2^j.
    Add that difference conditioned on s>j. All register bits are restored.
    """
    core=len(y)
    for i,q in enumerate(flags[1:]):b.left_rotate(y,1<<i,control=q)
    for j,q in enumerate(y):b.z_pow(q,Fraction(2)**(j-denominator+1))
    ctrl,tmp=b.alloc_qubits(2)
    for j,q in enumerate(y[:min(10,core)]):
        b.flip_if_greater_than(flags[1:],j,target=ctrl)
        b.init_and(ctrl,q,tmp)
        b.z_pow(tmp,((1<<core)-1)*Fraction(2)**(j-denominator+1))
        b.del_and(ctrl,q,tmp)
        b.flip_if_greater_than(flags[1:],j,target=ctrl)
    b.free(km.array([ctrl,tmp]))
    for i in reversed(range(len(flags)-1)):b.right_rotate(y,1<<i,control=flags[i+1])
