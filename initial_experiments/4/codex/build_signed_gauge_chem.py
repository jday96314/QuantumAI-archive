"""Exact sign-gauge angle compression for arbitrary-target ANGLE XOR oracles.

All words are sign extended from seven bits. Conjugating public bit7 by bit6
lets the lookup use seven-bit payloads while retaining the full public words.
"""
import argparse,json
from pathlib import Path
import numpy as np
import kickmix as km
from angle_gauge_probe import encode
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas
from monomial_unlookup import erase
from build_selectswap_chem import emit_qrom


def build(delta,lowbits=None):
    A,J=delta.shape;n=(A-1).bit_length();builder=km.CircuitBuilder();address=builder.create_quantum_register(n,'selection');words=[builder.create_quantum_register(8,f'word_{j}') for j in range(J)]
    assert np.all((delta>>7&1)==(delta>>6&1))
    for word in words:builder.cx(word[6],word[7])
    public=km.array([q for word in words for q in word[:7]]);width=len(public)
    packed=[sum((int(row[j])&127)<<(7*j) for j in range(J)) for row in delta]
    if A%2:packed.append(packed[-1])
    first=packed[::2];diff=[a^b for a,b in zip(first,packed[1::2])];active=[i for i in range(width) if any(v>>i&1 for v in diff)]
    extra=builder.alloc_qubits(len(active));compressed=[sum(((v>>i)&1)<<j for j,i in enumerate(active)) for v in diff];table=[a|(v<<width) for a,v in zip(first,compressed)];targets=km.array([*public,*extra])
    if lowbits is None:
        with builder.unary_iteration(address[1:]) as cursor:
            for a,row in enumerate(table):
                cursor.move_to(a);bits=km.array([q for i,q in enumerate(targets) if row>>i&1])
                if len(bits):builder.cx(cursor.match_qubit,bits)
    else:emit_qrom(builder,address[1:],targets,table,lowbits)
    for q,i in zip(extra,active):builder.ccx(address[0],q,public[i])
    deletion=compressed+[0]*((1<<len(address[1:]))-len(compressed));erase(builder,address[1:],extra,deletion);builder.free(extra)
    for word in words:builder.cx(word[6],word[7])
    return builder.finish_circuit()


def main():
    p=argparse.ArgumentParser();p.add_argument('--source',default='candidates/jacobi_l1/submission.npz');p.add_argument('--out',default='candidates/signed_gauge7');p.add_argument('--lowbits',type=int);a=p.parse_args()
    z=dict(np.load(a.source));tree=z['tree_edges'];old=z['angle_words'];u=reconstruct_u_vectors(old/128,tree);new=encode(u,tree,8,True);v=reconstruct_u_vectors(new/128,tree);sign=np.where((u*v).sum(1)<0,-1.,1.)
    assert abs(v-sign[:,None]*u).max()<1e-12
    assert np.all((new>>7&1)==(new>>6&1))
    z['angle_words']=new;z['U']=v[76:].reshape(z['U'].shape);out=Path(a.out);out.mkdir(parents=True,exist_ok=True);np.savez(out/'submission.npz',**z)
    rows=[]
    for k,delta in enumerate(batch_deltas(new,7)):
        c=build(delta,a.lowbits);assert c.num_qubits<=128;c.to_file(str(out/f'angle_k_{k}.kmb'),format='kmb')
        rows.append(dict(k=k,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like))
    Path(f'results/{out.name}_build.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(rows),flush=True)


if __name__=='__main__':main()
