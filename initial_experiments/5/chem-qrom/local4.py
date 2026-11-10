"""Four-row lookup with locally interpolated omitted difference columns.

Drop d columns belonging to the first low-bit coefficient. On each of
three high-address regions (64,64,71 rows), express them as GF(2) linear
combinations of the remaining 84-d columns. Reconstruct one parity at a
time and gate it into the public output. This exchanges a few Toffolis for
enough clean workspace to fit the four-row lookup under 128 qubits.
"""
import argparse
import kickmix as km
from common import *
from erase import erase
from linear import solve_many
from build_selectswap_chem import emit_qrom
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas

def build(delta,drop=5,lowbits=7):
    A,J=delta.shape;assert J==4
    n=(A-1).bit_length();width=7*J;b=km.CircuitBuilder()
    a=b.create_quantum_register(n,'selection');words=[b.create_quantum_register(8,f'word_{j}') for j in range(J)]
    assert np.all((delta>>7&1)==(delta>>6&1))
    for w in words:b.cx(w[6],w[7])
    public=km.array([q for w in words for q in w[:7]])
    packed=[sum((int(v)&127)<<(7*j) for j,v in enumerate(row)) for row in delta]
    packed+=[packed[-1]]*((-A)%4);coeff=[]
    for off in range(0,len(packed),4):
        x=packed[off:off+4];coeff.append([x[0],x[0]^x[1],x[0]^x[2],x[0]^x[1]^x[2]^x[3]])
    active=[(m,i) for m in range(1,4) for i in range(width) if any(row[m]>>i&1 for row in coeff)]
    missing=[(1,i) for i in range(drop)];assert all(t in active for t in missing)
    kept=[t for t in active if t not in missing];extra=b.alloc_qubits(len(kept))
    bank=[sum(((row[m]>>i)&1)<<j for j,(m,i) in enumerate(kept)) for row in coeff]
    omitted=[sum(((row[m]>>i)&1)<<j for j,(m,i) in enumerate(missing)) for row in coeff]
    regions=[(0,64),(64,128),(128,len(bank))]
    decoders=[solve_many(bank[lo:hi],omitted[lo:hi],len(kept),drop) for lo,hi in regions]
    table=[row[0]|(v<<width) for row,v in zip(coeff,bank)]
    emit_qrom(b,a[2:],km.array([*public,*extra]),table,lowbits)
    both=b.alloc_qubits(1)[0];b.init_and(a[0],a[1],both)
    controls={1:a[0],2:a[1],3:both}
    for q,(m,i) in zip(extra,kept):b.ccx(controls[m],q,public[i])
    b.del_and(a[0],a[1],both);b.free(km.array([both]))
    parity=b.alloc_qubits(1)[0]
    def apply(control,decoder):
        for j,mask in enumerate(decoder):
            inputs=km.array([q for i,q in enumerate(extra) if mask>>i&1])
            b.cx(inputs,parity);b.ccx(control,parity,public[missing[j][1]]);b.cx(inputs,parity)
    # Base decoder applies everywhere; two regional XOR corrections suffice.
    apply(a[0],decoders[2])
    flag,gate=b.alloc_qubits(2)
    b.x(a[9]);b.x(a[8]);b.init_and(a[9],a[8],flag);b.x(a[8]);b.x(a[9])
    for region in [0,1]:
        b.init_and(a[0],flag,gate)
        apply(gate,[x^y for x,y in zip(decoders[region],decoders[2])])
        b.del_and(a[0],flag,gate)
        if region==0:b.x(a[9]);b.cx(a[9],flag);b.x(a[9])
    b.x(a[9]);b.del_and(a[9],a[8],flag);b.x(a[9]);b.free(km.array([flag,gate,parity]))
    info=erase(b,a[2:],extra,bank);b.free(extra)
    for w in words:b.cx(w[6],w[7])
    info.update(bank_columns=len(kept),omitted=missing,decoders=decoders,regions=regions)
    return b.finish_circuit(),info

def family(drop=5,lowbits=7):
    out=HERE/'candidates'/f'local4_d{drop}_l{lowbits}';out.mkdir(exist_ok=True);source=INC/'submission.npz';(out/'submission.npz').write_bytes(source.read_bytes());rows=[]
    for k,d in enumerate(batch_deltas(np.load(source)['angle_words'],4)):
        c,info=build(d,drop,lowbits);assert c.num_qubits<=128
        c.to_file(str(out/f'angle_k_{k}.kmb'),format='kmb')
        rows.append(dict(k=k,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,interpolation=info))
    result=dict(rows=rows,predicted=predicted(4,rows));save(out.name+'_build.json',result);print(json.dumps(result),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--drop',type=int,default=5);p.add_argument('--lowbits',type=int,default=7);a=p.parse_args();family(a.drop,a.lowbits)
