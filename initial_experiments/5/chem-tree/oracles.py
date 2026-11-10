"""Signed ANGLE QROM, optional shared constant low grid-phase bits.

Adapted from completed chem-qfit/oracles.py and read-only round-4 Codex
monomial_unlookup/build_selectswap_chem. Public words retain the exact host
width/semantics; constant bits are emitted with reversible X gates.
"""
import numpy as np
import kickmix as km
from monomial_unlookup import erase
from build_selectswap_chem import emit_qrom

def build(delta,bits=8,low_constant_bits=0,lowbits=8):
    A,J=delta.shape;n=(A-1).bit_length();b=km.CircuitBuilder()
    address=b.create_quantum_register(n,'selection')
    words=[b.create_quantum_register(bits,f'word_{j}') for j in range(J)]
    assert np.all((delta>>(bits-1)&1)==(delta>>(bits-2)&1))
    mask=(1<<low_constant_bits)-1
    assert np.all((delta&mask)==(delta[0]&mask)[None,:])
    for j,word in enumerate(words):
        for k in range(low_constant_bits):
            if int(delta[0,j])>>k&1:b.x(word[k])
        b.cx(word[bits-2],word[bits-1])
    public=km.array([q for word in words for q in word[low_constant_bits:bits-1]])
    payload_bits=bits-1-low_constant_bits;width=len(public)
    packed=[sum(((int(row[j])>>low_constant_bits)&((1<<payload_bits)-1))<<(payload_bits*j) for j in range(J)) for row in delta]
    if A%2:packed.append(packed[-1])
    first=packed[::2];diff=[x^y for x,y in zip(first,packed[1::2])]
    active=[i for i in range(width) if any(v>>i&1 for v in diff)]
    split=(n-1)//2;needed=sum((1<<nn)-nn-1 for nn in [split,n-1-split])
    extra=b.alloc_qubits(max(len(active),needed))
    compressed=[sum(((v>>i)&1)<<j for j,i in enumerate(active)) for v in diff]
    table=[x|(v<<width) for x,v in zip(first,compressed)]
    emit_qrom(b,address[1:],km.array([*public,*extra]),table,min(lowbits,n-1))
    for q,i in zip(extra,active):b.ccx(address[0],q,public[i])
    erase(b,address[1:],extra,compressed+[0]*((1<<(n-1))-len(compressed)));b.free(extra)
    for word in words:b.cx(word[bits-2],word[bits-1])
    return b.finish_circuit()
