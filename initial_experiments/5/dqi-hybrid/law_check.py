"""One million host-law syndrome cases for the overwrite representation.

This checks a local representation, not the full BM decoder. Locator operands
are independently randomized (including forced exceptional degrees), so the
claim is stronger than assuming generic successful BM states for this local
map, but it says nothing about the full decoder's phase-error rate.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import kickmix as km
from dqi.verify_submission_dqi import optimal_dqi_weights,FieldTables
from published_dqi import isomorphism
from import_published import load
from consume_component import field_table,reference

HERE=Path(__file__).resolve().parent


def check(shots,salt,L):
    started=time.time();m=load();mul=field_table(m)
    tables=FieldTables(km.GF2Field(12))
    forward,_,transform=isomorphism(m)
    conversion=np.array([transform(v,forward) for v in range(4096)],dtype=np.uint16)
    inverse=np.array([m.ginv(v) for v in range(4096)],dtype=np.uint16)
    probs=optimal_dqi_weights(4095,35,12,2016)**2
    seed=int.from_bytes(hashlib.sha256(salt.encode()).digest()[:16],'little')
    rng=np.random.default_rng(seed)
    weights_hist=np.zeros(36,dtype=np.int64);zeros=0;cases=0;failures=0;forced=0
    fingerprint=hashlib.sha256();degree_hist=np.zeros(L+1,dtype=np.int64)
    for start in range(0,shots,4096):
        count=min(4096,shots-start)
        weights=rng.choice(36,size=count,p=probs);weights_hist+=np.bincount(weights,minlength=36)
        pos=np.zeros((count,35),dtype=np.int64);val=np.zeros_like(pos)
        for i,w in enumerate(weights):
            pos[i,:w]=rng.choice(4095,size=int(w),replace=False)
            val[i,:w]=rng.integers(1,4096,size=int(w))
        used=val!=0
        S=np.zeros((count,L+1),dtype=np.uint16)
        # r=2L-1: ascending i<L uses S[r-i], target uses S[r-L].
        rows=list(range(2*L,L,-1))+[L]
        for i,row in enumerate(rows):
            exponents=(tables.log[val]+row*pos)%4095
            terms=np.where(used,tables.antilog[exponents],0)
            S[:,i]=conversion[np.bitwise_xor.reduce(terms,axis=1)]
        A=rng.integers(0,4096,(count,L+1),dtype=np.uint16)
        # Exercise every zero/nonzero locator-degree boundary independently
        # of its rarity in the syndrome law. S retains the exact host law.
        degrees=(np.arange(count)+start)%(L+1)
        A[np.arange(L+1)[None,:]>degrees[:,None]]=0
        degree_hist+=np.bincount(degrees,minlength=L+1)
        forced+=int(np.count_nonzero(degrees<L))
        out,park=reference(mul,A,S)
        other=np.bitwise_xor.reduce(mul[A[:,:-1],S[:,:-1]],axis=1)
        diagonal=A[:,-1];zero=diagonal==0;armed=np.where(zero,1,diagonal)
        recovered=mul[out[:,-1]^other,inverse[armed]]
        recovered=np.where(zero,park,recovered)
        failure=(recovered!=S[:,-1])|((~zero)&(park!=0))
        failures+=int(failure.sum());zeros+=int(zero.sum());cases+=count
        fingerprint.update(S.tobytes());fingerprint.update(A.tobytes());fingerprint.update(out.tobytes())
        if start//4096%32==0:print(json.dumps(dict(checked=cases,failures=failures,seconds=time.time()-started)),flush=True)
    report=dict(status='PASS' if failures==0 else 'FAIL',scope='local discrepancy overwrite and reconstruction, NOT full BM decode or host DQI acceptance',
        shots=cases,salt=salt,L=L,value_roundtrip_errors=failures,zero_pivot_cases=zeros,
        forced_exceptional_degree_cases=forced,degree_histogram=degree_hist.tolist(),
        law_weight_histogram=weights_hist.tolist(),exact_weight_probabilities=probs.tolist(),
        syndrome_law='host optimal_dqi_weights; uniform distinct positions; uniform nonzero values; native-to-upstream field isomorphism',
        operand_stream_sha256=fingerprint.hexdigest(),seconds=time.time()-started,
        naive_unparked_counterexample=dict(pivot=0,other_terms=7,inputs=[1,2],both_discrepancies=[7,7],
            implication='12 bits of retained information are necessary on the zero-pivot fiber when other operands remain fixed'))
    (HERE/f'results/law_{salt}.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)
    assert failures==0


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--shots',type=int,default=1000000);p.add_argument('--salt',default='law-independent-a');p.add_argument('--l',type=int,default=34);a=p.parse_args()
    check(a.shots,a.salt,a.l)
