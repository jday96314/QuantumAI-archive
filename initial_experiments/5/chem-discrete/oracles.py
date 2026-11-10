"""Selected signed oracle with enough clean bank space for 11-bit addresses.

Adapted from round-4 Codex build_signed_gauge_chem.py. Padding is necessary
for monomial phase repair's two 5-bit banks (52 clean qubits, versus 49
payload qubits). Padding bits hold zero and do not alter the XOR interface.
"""
import numpy as np
import kickmix as km
from monomial_unlookup import erase
from build_selectswap_chem import emit_qrom

def build(delta,lowbits=8):
    A,J=delta.shape;n=(A-1).bit_length();b=km.CircuitBuilder()
    address=b.create_quantum_register(n,'selection');words=[b.create_quantum_register(8,f'word_{j}') for j in range(J)]
    assert np.all((delta>>7&1)==(delta>>6&1))
    for word in words:b.cx(word[6],word[7])
    public=km.array([q for word in words for q in word[:7]]);width=len(public)
    packed=[sum((int(row[j])&127)<<(7*j) for j in range(J)) for row in delta]
    if A%2:packed.append(packed[-1])
    first=packed[::2];diff=[x^y for x,y in zip(first,packed[1::2])];active=[i for i in range(width) if any(v>>i&1 for v in diff)]
    split=(n-1)//2;needed=sum((1<<nn)-nn-1 for nn in [split,n-1-split])
    extra=b.alloc_qubits(max(len(active),needed));compressed=[sum(((v>>i)&1)<<j for j,i in enumerate(active)) for v in diff]
    table=[x|(v<<width) for x,v in zip(first,compressed)]
    emit_qrom(b,address[1:],km.array([*public,*extra]),table,min(lowbits,n-1))
    for q,i in zip(extra,active):b.ccx(address[0],q,public[i])
    erase(b,address[1:],extra,compressed+[0]*((1<<(n-1))-len(compressed)));b.free(extra)
    for word in words:b.cx(word[6],word[7])
    return b.finish_circuit()
