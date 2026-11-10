"""Hamiltonian-preserving row orders and exact parity-cluster pairing.

Extend bucket labels by GF(2) features whose counts have even parity in
each existing bucket. Rows can then be paired inside every bucket. Such
features annihilate every pairwise difference and lower bank rank.
"""
import argparse,random
from functools import reduce
from operator import xor
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors,reconstruct_x,lambda_h2
from common import *
from linear import basis,nullspace
from build import family

def pack(words):return [sum((int(v)&127)<<(7*j) for j,v in enumerate(row)) for row in words]

def parity_order(words,oracle=0,trials=64):
    table=pack(batch_deltas(words,7)[oracle]);rng=random.Random(2026102402)
    groups=[list(range(76)),list(range(76,len(words)-1))];features=[];history=[]
    while True:
        constraints=[reduce(xor,(table[i] for i in g),0) for g in groups]
        ns=nullspace(constraints,49);oldrank=len(basis(features));choices=[]
        for t in range(trials):
            f=reduce(xor,(v for v in ns if rng.getrandbits(1)),0)
            if len(basis(features+[f]))==oldrank:continue
            gg=[[i for i in g if (table[i]&f).bit_count()%2==v] for g in groups for v in [0,1]]
            assert all(len(g)%2==0 for g in gg)
            choices.append((sum(bool(g) for g in gg),sum(len(g)==2 for g in gg),f,gg))
        if not choices:break
        _,_,f,groups=min(choices,key=lambda x:x[:3]);features.append(f)
        history.append(dict(features=len(features),nonempty_groups=sum(bool(g) for g in groups),nullity=len(ns)))
    # Keep the two physical sectors separate; no one-body row crosses 76.
    order=[i for g in groups for i in g if i<76]+[i for g in groups for i in g if i>=76]+[len(words)-1]
    assert sorted(order)==list(range(len(words)))
    for f in features:
        assert all(((table[order[i]]^table[order[i+1]])&f).bit_count()%2==0 for i in range(0,len(words)-1,2))
    return np.array(order),dict(features=features,history=history,oracle=oracle)

def greedy(words):
    # Greedy Hamming similarity over the complete 525-bit signed table.
    data=pack(words);order=[]
    for ids in [list(range(76)),list(range(76,len(words)))]:
        current=ids.pop(0);order.append(current)
        while ids:
            loc=min(range(len(ids)),key=lambda j:((data[current]^data[ids[j]]).bit_count(),ids[j]))
            current=ids.pop(loc);order.append(current)
    return np.array(order)

def transform(order,name,details):
    z=dict(np.load(INC/'submission.npz'));old=z['angle_words'];u=reconstruct_u_vectors(old/128,z['tree_edges']);x=reconstruct_x(u[76:].reshape(z['U'].shape),z['W'])
    order=np.asarray(order);assert np.all(order[:76]<76) and np.all(order[76:]>=76)
    z['angle_words']=old[order];two=order[76:]-76
    z['U']=z['U'][:,two,:];z['W']=z['W'][:,two,:]
    v=reconstruct_u_vectors(z['angle_words']/128,z['tree_edges']);xx=reconstruct_x(v[76:].reshape(z['U'].shape),z['W'])
    inverse=np.argsort(order)
    assert np.array_equal(z['angle_words'][inverse],old)
    assert np.array_equal(v[ inverse],u)
    out=HERE/'candidates'/name;out.mkdir(exist_ok=True);np.savez(out/'submission.npz',**z)
    data=dict(order=order.tolist(),inverse=inverse.tolist(),max_X_difference=float(abs(x-xx).max()),
              lambda_h2=lambda_h2(z['W'],z['H2_I']-z['W'].sum(1)),details=details)
    assert data['max_X_difference']<1e-13
    save(name+'_transform.json',data)
    family(out/'submission.npz',out,batch=7,select=1,lowbits=8,compress=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('method',choices=['parity','greedy','lex']);p.add_argument('--oracle',type=int,default=0);a=p.parse_args();w=np.load(INC/'submission.npz')['angle_words']
    if a.method=='parity':order,details=parity_order(w,a.oracle);name=f'parity{a.oracle}'
    elif a.method=='greedy':order=greedy(w);details={};name='greedy7'
    else:order=np.r_[np.lexsort(w[:76].T[::-1]),76+np.lexsort(w[76:].T[::-1])];details={};name='lex7'
    transform(order,name,details)
