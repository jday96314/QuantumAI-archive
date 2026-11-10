"""Search 128 exact affine Givens gauges and reversible Gray decoding.

At reverse edge (p,c), choosing the signed norm of p shifts the angle by
pi. Choose that sign so bit 7 equals any linear parity of bits 0..6.
Subsequent reverse edges absorb the signed norm. The final root sign is a
physically irrelevant global vector sign. No continuous refitting is used.
"""
import argparse
from common import *
from linear import basis
from build import build,family
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas

def encode(vectors,tree,mask):
    v=vectors.copy();words=np.zeros((len(v),len(tree)),dtype=np.int64)
    for t in reversed(range(len(tree))):
        p,c=tree[t];a,b=v[:,p].copy(),v[:,c].copy()
        w=np.rint(np.arctan2(b,a)*128/np.pi).astype(np.int64)%256
        parity=np.zeros(len(v),dtype=np.int64)
        for i in range(7):
            if mask>>i&1:parity^=w>>i&1
        flip=(w>>7)!=parity;w[flip]^=128
        words[:,t]=w;v[:,p]=np.where(flip,-1.,1.)*np.hypot(a,b);v[:,c]=0
    return words

def run():
    z=dict(np.load(INC/'submission.npz'));old=z['angle_words'];tree=z['tree_edges'];u=reconstruct_u_vectors(old/128,tree);rows=[]
    for mask in range(128):
        w=encode(u,tree,mask);v=reconstruct_u_vectors(w/128,tree)
        sign=np.where((u*v).sum(1)<0,-1.,1.);error=float(abs(v-sign[:,None]*u).max());assert error<1e-12
        deltas=batch_deltas(w,7);ranks=[]
        for delta in deltas:
            packed=[sum((int(x)&127)<<(7*j) for j,x in enumerate(row)) for row in delta]
            packed+=[packed[-1]]*(len(packed)%2)
            ranks.append(len(basis([x^y for x,y in zip(packed[::2],packed[1::2])])))
        rows.append(dict(mask=mask,projector_vector_error=error,ranks=ranks,changed_words=int(np.count_nonzero(w!=old))))
        if mask in [0,1,63,64,127]:
            out=HERE/'candidates'/f'gauge{mask}';out.mkdir(exist_ok=True);zz=dict(z);zz['angle_words']=w;zz['U']=v[76:].reshape(z['U'].shape);np.savez(out/'submission.npz',**zz)
            family(out/'submission.npz',out,lowbits=8,relation=mask)
    save('gauge_search.json',rows)
    family(INC/'submission.npz',HERE/'candidates/gray7',lowbits=8,gray=True)

if __name__=='__main__':run()
