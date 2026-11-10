"""Exact Givens sign gauges and linear rank of select-copy differences.

The projector u*u^T is unchanged by the backward signed-hypot convention.
This diagnostic does not change the selected submission.
"""
import json
from pathlib import Path
import numpy as np
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas


def encode(vectors,tree,bits,signed):
    v=vectors.copy();words=np.zeros((len(v),len(tree)),dtype=np.int64);scale=1<<(bits-1)
    for t in reversed(range(len(tree))):
        p,c=tree[t];a,b=v[:,p].copy(),v[:,c].copy()
        sign=np.where(a<0,-1.,1.) if signed else np.ones(len(v))
        theta=np.arctan2(b*sign,a*sign);w=np.rint(theta/np.pi*scale).astype(np.int64)
        if signed:
            flip=w>=scale//2;w[flip]-=scale;sign[flip]*=-1
            flip=w<-(scale//2);w[flip]+=scale;sign[flip]*=-1
        words[:,t]=w%(2*scale);v[:,p]=sign*np.hypot(a,b);v[:,c]=0
    return words


def rank(rows):
    basis={}
    for x in rows:
        while x:
            bit=x.bit_length()-1
            if bit in basis:x^=basis[bit]
            else:basis[bit]=x;break
    return len(basis)


def probe():
    z=dict(np.load('candidates/jacobi_l1/submission.npz'));old=z['angle_words'];tree=z['tree_edges'];u=reconstruct_u_vectors(old/128,tree);records=[]
    for name,new in [('original',old),('positive',encode(u,tree,8,False)),('signed',encode(u,tree,8,True))]:
        v=reconstruct_u_vectors(new/128,tree);sign=np.where((u*v).sum(1)<0,-1.,1.);err=float(abs(v-sign[:,None]*u).max());assert err<1e-12,(name,err)
        signext=int(np.count_nonzero((new>>7&1)!=(new>>6&1)))
        zz=dict(z);zz['angle_words']=new;zz['U']=v[76:].reshape(z['U'].shape);np.savez(f'parameters/angle_gauge_{name}.npz',**zz)
        entry=dict(name=name,vector_max_error_up_to_sign=err,sign_extension_exceptions=signext,batches={})
        for batch in [7,8,9,10]:
            data=[]
            for k,delta in enumerate(batch_deltas(new,batch)):
                packed=[sum(int(row[j])<<(8*j) for j in range(batch)) for row in delta]
                if len(packed)%2:packed.append(packed[-1])
                diff=[x^y for x,y in zip(packed[::2],packed[1::2])];active=int(np.bitwise_or.reduce(np.array(diff,dtype=object))).bit_count();r=rank(diff)
                data.append(dict(k=k,active_bits=active,difference_rank=r,Q_floor=8*batch+r+10))
            entry['batches'][str(batch)]=data
        records.append(entry)
    Path('results/angle_gauge_probe.json').write_text(json.dumps(records,indent=2)+'\n')
    print(json.dumps(records),flush=True)


if __name__=='__main__':probe()
