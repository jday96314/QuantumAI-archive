"""Group rows in fours by linear output features; combine compilers.

Selected oracle gets an exact invertible output basis that makes three
difference columns vanish per feature. Other oracles use quadratic bank
compression. Counts divisible by four in both physical sectors guarantee
the grouping preserves the one-/two-body boundary.
"""
import argparse,random,time
from common import *
from linear import basis
from build import build as linear_build
from quadratic import build as nonlinear_build
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors,reconstruct_x

def search(words,k=0,max_features=3):
    d=batch_deltas(words,4)[k]
    packed=np.array([sum((int(x)&127)<<(7*j) for j,x in enumerate(row)) for row in d],dtype=np.uint64)
    features=[];labels=np.zeros(len(words),dtype=np.int64);rng=random.Random(2026102404+k);history=[]
    for count in range(max_features):
        for trial in range(250000):
            f=rng.getrandbits(28)
            if len(basis(features+[f]))==len(features):continue
            new=2*labels+(np.bitwise_count(packed&np.uint64(f))%2).astype(np.int64)
            c1=np.bincount(new[:76],minlength=1<<(count+1))
            c2=np.bincount(np.r_[new[76:],new[-1]],minlength=1<<(count+1))
            if np.any(c1%4) or np.any(c2%4):continue
            labels=new;features.append(f);history.append(dict(feature=f,trial=trial,onebody=c1.tolist(),twobody_with_padding=c2.tolist()));break
        else:break
    classes=range(1<<len(features));order=[]
    for label in classes:order.extend(np.flatnonzero(labels[:76]==label).tolist())
    last=int(labels[-1]);order2=[]
    for label in [x for x in classes if x!=last]+[last]:
        order2.extend((76+np.flatnonzero(labels[76:-1]==label)).tolist())
    order+=order2+[len(words)-1]
    matrix=list(features)
    for i in range(28):
        if len(basis(matrix+[1<<i]))>len(matrix):matrix.append(1<<i)
    assert len(matrix)==28
    return np.array(order),matrix,history

def family(oracle=0):
    source=INC/'submission.npz';z=dict(np.load(source));old=z['angle_words'];order,matrix,history=search(old,oracle)
    out=HERE/'candidates'/f'hybrid4_o{oracle}';out.mkdir(exist_ok=True);z['angle_words']=old[order]
    two=order[76:]-76;z['U']=z['U'][:,two,:];z['W']=z['W'][:,two,:];np.savez(out/'submission.npz',**z)
    u=reconstruct_u_vectors(old/128,z['tree_edges']);v=reconstruct_u_vectors(z['angle_words']/128,z['tree_edges']);inverse=np.argsort(order)
    assert np.array_equal(u,v[inverse]);xx=reconstruct_x(u[76:].reshape(z['U'].shape),np.load(source)['W']);yy=reconstruct_x(v[76:].reshape(z['U'].shape),z['W'])
    transform=dict(order=order.tolist(),inverse=inverse.tolist(),matrix=matrix,features=history,max_X_difference=float(abs(xx-yy).max()))
    save(out.name+'_transform.json',transform);rows=[]
    for k,d in enumerate(batch_deltas(z['angle_words'],4)):
        if k==oracle:c,info=linear_build(d,select=2,lowbits=7,mix=matrix)
        else:c,info=nonlinear_build(d,k)
        assert c.num_qubits<=128;c.to_file(str(out/f'angle_k_{k}.kmb'),format='kmb')
        row=dict(k=k,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,compiler='linear' if k==oracle else 'quadratic',info=info);rows.append(row)
        print(json.dumps(dict(k=k,T=row['T'],Q=row['Q'])),flush=True)
    result=dict(rows=rows,predicted=predicted(4,rows));save(out.name+'_build.json',result)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--oracle',type=int,default=0);a=p.parse_args();family(a.oracle)
