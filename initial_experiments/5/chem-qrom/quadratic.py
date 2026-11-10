"""Content-aware nonlinear compression of a QROM difference bank.

For 199 cared-for rows, retain 79 of 84 difference bits and interpolate
five missing bits as shared products of affine parities of the retained
bank. Fixed product left factors and output vectors make the fit exactly
linear over GF(2). No numerical approximation or physical refit occurs.
"""
import argparse,random,time
import kickmix as km
from common import *
from erase import erase
from build_selectswap_chem import emit_qrom
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas

def solve(rows,values):
    pivots={}
    for x,y in zip(rows,values):
        while x:
            p=x.bit_length()-1
            if p in pivots:
                xx,yy=pivots[p];x^=xx;y^=yy
            else:pivots[p]=(x,y);break
        if not x and y:return None,len(pivots)
    solution=0
    for p,(row,rhs) in sorted(pivots.items()):
        if ((row&solution).bit_count()^rhs)&1:solution^=1<<p
    assert all((row&solution).bit_count()%2==value for row,value in zip(rows,values))
    return solution,len(pivots)

def fit(bank,missing,d=5,seed=2026102403,max_products=12,trials=32):
    width=max(v.bit_length() for v in bank);n=width+1;aug=[v|(1<<width) for v in bank];rng=random.Random(seed);attempts=[]
    for count in range(7,max_products+1):
        for trial in range(trials):
            left=[rng.getrandbits(n) for _ in range(count)]
            output=[rng.randrange(1,1<<d) for _ in range(count)]
            rows=[];values=[]
            for x,y in zip(aug,missing):
                active=[k for k,l in enumerate(left) if (l&x).bit_count()%2]
                for j in range(d):
                    row=x<<(j*n)
                    for k in active:
                        if output[k]>>j&1:row^=x<<((d+k)*n)
                    rows.append(row);values.append(y>>j&1)
            solution,rank=solve(rows,values);attempts.append(dict(products=count,trial=trial,rank=rank,consistent=solution is not None))
            if solution is None:continue
            masks=[(solution>>(j*n))&((1<<n)-1) for j in range(d+count)]
            model=dict(width=width,linear=masks[:d],left=left,right=masks[d:],output=output,products=count,trial=trial,seed=seed,attempts=attempts)
            for x,y in zip(aug,missing):
                value=sum(((x&l).bit_count()%2)<<j for j,l in enumerate(model['linear']))
                for l,r,o in zip(left,model['right'],output):
                    if (x&l).bit_count()%2 and (x&r).bit_count()%2:value^=o
                assert value==y
            return model
    raise ValueError('No exact quadratic interpolation found')

def build(delta,k=0,drop=5,lowbits=7,model=None):
    A,J=delta.shape;assert J==4;n=(A-1).bit_length();width=7*J;b=km.CircuitBuilder()
    a=b.create_quantum_register(n,'selection');words=[b.create_quantum_register(8,f'word_{j}') for j in range(J)]
    assert np.all((delta>>7&1)==(delta>>6&1))
    for w in words:b.cx(w[6],w[7])
    public=km.array([q for w in words for q in w[:7]])
    packed=[sum((int(v)&127)<<(7*j) for j,v in enumerate(row)) for row in delta];packed+=[packed[-1]]*((-A)%4);coeff=[]
    for off in range(0,len(packed),4):
        x=packed[off:off+4];coeff.append([x[0],x[0]^x[1],x[0]^x[2],x[0]^x[1]^x[2]^x[3]])
    active=[(m,i) for m in range(1,4) for i in range(width) if any(row[m]>>i&1 for row in coeff)]
    missing=[(1,i) for i in range(drop)];kept=[t for t in active if t not in missing]
    bank=[sum(((row[m]>>i)&1)<<j for j,(m,i) in enumerate(kept)) for row in coeff]
    omitted=[sum(((row[m]>>i)&1)<<j for j,(m,i) in enumerate(missing)) for row in coeff]
    if model is None:model=fit(bank,omitted,drop,seed=2026102403+k)
    extra=b.alloc_qubits(len(kept));table=[row[0]|(v<<width) for row,v in zip(coeff,bank)]
    emit_qrom(b,a[2:],km.array([*public,*extra]),table,lowbits)
    both=b.alloc_qubits(1)[0];b.init_and(a[0],a[1],both);controls={1:a[0],2:a[1],3:both}
    for q,(m,i) in zip(extra,kept):b.ccx(controls[m],q,public[i])
    b.del_and(a[0],a[1],both);b.free(km.array([both]))
    decoded=b.alloc_qubits(drop);left,right=b.alloc_qubits(2)
    def parity(mask,target):
        b.cx(km.array([q for i,q in enumerate(extra) if mask>>i&1]),target)
        if mask>>len(extra)&1:b.x(target)
    for q,l in zip(decoded,model['linear']):parity(l,q)
    for l,r,o in zip(model['left'],model['right'],model['output']):
        parity(l,left);parity(r,right);targets=km.array([q for j,q in enumerate(decoded) if o>>j&1])
        if len(targets)>1:b.cx(targets[0],targets[1:])
        b.ccx(left,right,targets[0])
        if len(targets)>1:b.cx(targets[0],targets[1:])
        parity(r,right);parity(l,left)
    for j,q in enumerate(decoded):b.ccx(a[0],q,public[missing[j][1]])
    outcomes=b.alloc_bits(drop);b.hmr(decoded,outcomes)
    # Affine phase corrections and quadratic CZs need no quantum ANDs.
    for outcome,l in zip(outcomes,model['linear']):
        parity(l,left);b.push_condition(outcome);b.z(left);b.pop_condition();parity(l,left)
    flag=b.alloc_bits(1)[0]
    for l,r,o in zip(model['left'],model['right'],model['output']):
        selected=km.array([outcome for j,outcome in enumerate(outcomes) if o>>j&1])
        b.cx(selected,flag);parity(l,left);parity(r,right);b.push_condition(flag);b.cz(left,right);b.pop_condition();parity(r,right);parity(l,left);b.cx(selected,flag)
    b.bit_store(outcomes,False);b.free(outcomes);b.free(km.array([flag]));b.free(km.array([left,right]));b.free(decoded)
    info=erase(b,a[2:],extra,bank);b.free(extra)
    for w in words:b.cx(w[6],w[7])
    info.update(bank_columns=len(kept),omitted=missing,model=model)
    return b.finish_circuit(),info

def family(drop=5,lowbits=7):
    out=HERE/'candidates'/f'quadratic4_d{drop}_l{lowbits}';out.mkdir(exist_ok=True);source=INC/'submission.npz';(out/'submission.npz').write_bytes(source.read_bytes());rows=[]
    for k,d in enumerate(batch_deltas(np.load(source)['angle_words'],4)):
        t=time.time();c,info=build(d,k,drop,lowbits);assert c.num_qubits<=128
        c.to_file(str(out/f'angle_k_{k}.kmb'),format='kmb');row=dict(k=k,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,compression=info,seconds=time.time()-t);rows.append(row);print(json.dumps(dict(k=k,T=row['T'],Q=row['Q'],products=info['model']['products'],seconds=row['seconds'])),flush=True)
        save(out.name+'_build.json',dict(rows=rows,complete=len(rows)==len(batch_deltas(np.load(source)['angle_words'],4)),predicted=predicted(4,rows)))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--drop',type=int,default=5);p.add_argument('--lowbits',type=int,default=7);a=p.parse_args();family(a.drop,a.lowbits)
