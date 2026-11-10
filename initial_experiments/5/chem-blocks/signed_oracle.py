"""Generalize the round-4 signed-gauge XOR oracle to 7/8/9-bit words.

Reuses local prior-round emit_qrom and ANF erasure without changing them.
"""
import numpy as np
import kickmix as km
from build_selectswap_chem import emit_qrom
from monomial_unlookup import erase

def build(delta, bits, lowbits):
    A,J=delta.shape
    n=(A-1).bit_length()
    b=km.CircuitBuilder()
    address=b.create_quantum_register(n,'selection')
    words=[b.create_quantum_register(bits,f'word_{j}') for j in range(J)]
    assert np.all((delta>>(bits-1)&1)==(delta>>(bits-2)&1))
    for w in words:b.cx(w[bits-2],w[bits-1])
    public=km.array([q for w in words for q in w[:bits-1]])
    width=len(public)
    packed=[sum((int(row[j])&((1<<(bits-1))-1))<<((bits-1)*j) for j in range(J)) for row in delta]
    if A%2:packed.append(packed[-1])
    first=packed[::2]
    diff=[x^y for x,y in zip(first,packed[1::2])]
    active=[i for i in range(width) if any(v>>i&1 for v in diff)]
    # Padding guarantees enough clean bank wires for the shared ANF monomials.
    m=n-1
    needed=sum((1<<t)-t-1 for t in [m//2,m-m//2])
    bank=b.alloc_qubits(max(len(active),needed))
    compressed=[sum(((v>>i)&1)<<j for j,i in enumerate(active)) for v in diff]
    table=[x|(v<<width) for x,v in zip(first,compressed)]
    emit_qrom(b,address[1:],km.array([*public,*bank]),table,lowbits)
    for q,i in zip(bank,active):b.ccx(address[0],q,public[i])
    erase(b,address[1:],bank,compressed+[0]*((1<<m)-len(compressed)))
    b.free(bank)
    for w in words:b.cx(w[bits-2],w[bits-1])
    return b.finish_circuit()

def bound(A,bits,J):
    """Worst-case T/Q for this explicit builder, including padded ANF bank.

    Every nontrivial binary branch costs one Toffoli; a chunk of L rows has
    at most L-1. Prefix AND setup+measurement repair costs 0 or 2*p-3.
    Erasure costs exactly the number of nonlinear half-address monomials.
    """
    n=(A-1).bit_length();m=n-1;L=(A+1)//2
    needed=sum((1<<t)-t-1 for t in [m//2,m-m//2])
    active=J*(bits-1);bank=max(active,needed)
    options=[]
    for low in range(1,m):
        p=m-low;chunks=(L+(1<<low)-1)//(1<<low)
        prefix=0 if p==1 else 2*p-3
        T=L-chunks+chunks*prefix+active+needed
        workspace=bank+max(low,max(1,p-1))
        Q=n+J*bits+workspace
        if Q<=128:options.append(dict(T_per_oracle=T,Q_oracle=Q,angle_temp=workspace,lowbits=low,ANF_monomials=needed))
    return min(options,key=lambda x:(x['T_per_oracle'],x['Q_oracle'])) if options else None
