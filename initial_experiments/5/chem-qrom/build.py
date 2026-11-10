"""Generalized select-by-ANF lookup with streamed measured erasure."""
import argparse
import kickmix as km
import numpy as np
from build_selectswap_chem import emit_qrom
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas
from common import *
from erase import erase
from linear import basis,coordinates


def build(delta,select=1,lowbits=4,care=True,compress=False,gray=False,relation=64,mix=None):
    A,J=delta.shape;n=(A-1).bit_length();b=km.CircuitBuilder()
    address=b.create_quantum_register(n,'selection')
    words=[b.create_quantum_register(8,f'word_{j}') for j in range(J)]
    parity=np.zeros_like(delta)
    for i in range(7):
        if relation>>i&1:parity^=delta>>i&1
    assert np.all((delta>>7&1)==parity)
    for word in words:
        for i in range(7):
            if relation>>i&1:b.cx(word[i],word[7])
    if gray:
        for word in words:
            for i in range(6):b.cx(word[i+1],word[i])
        delta=(delta&127)^((delta&127)>>1)
    public=km.array([q for word in words for q in word[:7]]);width=len(public)
    packed=[sum((int(v)&127)<<(7*j) for j,v in enumerate(row)) for row in delta]
    if mix:
        from linear import row_reduction
        undo=row_reduction(mix)
        for op,i,j in reversed(undo):
            if op=='cx':b.cx(public[j],public[i])
            else:b.swap(public[i],public[j])
        packed=[sum(((x&row).bit_count()%2)<<i for i,row in enumerate(mix)) for x in packed]
    K=1<<select
    packed+=[packed[-1]]*((-A)%K)
    coeff=[]
    for off in range(0,len(packed),K):
        row=packed[off:off+K]
        for i in range(select):
            for m in range(K):
                if m>>i&1:row[m]^=row[m^(1<<i)]
        coeff.append(row)
    active=[(m,i) for m in range(1,K) for i in range(width) if any(row[m]>>i&1 for row in coeff)]
    bank=[sum(((row[m]>>i)&1)<<j for j,(m,i) in enumerate(active)) for row in coeff]
    bs=basis(bank) if compress else [1<<i for i in range(len(active))]
    if compress:bank=coordinates(bank,bs)
    extra=b.alloc_qubits(len(bs))
    table=[row[0]|(v<<width) for row,v in zip(coeff,bank)]
    emit_qrom(b,address[select:],km.array([*public,*extra]),table,lowbits)
    controls={}
    for m in range(1,K):
        if m.bit_count()==1:controls[m]=address[m.bit_length()-1]
        else:
            bit=m&-m;rest=m^bit;q=b.alloc_qubits(1)[0]
            b.init_and(controls[rest],address[bit.bit_length()-1],q);controls[m]=q
    for j,q in enumerate(extra):
        for m in range(1,K):
            targets=km.array([public[i] for a,(mm,i) in enumerate(active) if mm==m and bs[j]>>a&1])
            if len(targets):
                if len(targets)>1:b.cx(targets[0],targets[1:])
                b.ccx(controls[m],q,targets[0])
                if len(targets)>1:b.cx(targets[0],targets[1:])
    for m in reversed(range(1,K)):
        if m.bit_count()>1:
            bit=m&-m;rest=m^bit;q=controls[m]
            b.del_and(controls[rest],address[bit.bit_length()-1],q);b.free(km.array([q]))
    info=erase(b,address[select:],extra,bank,care=care);b.free(extra)
    info['bank_rank']=len(bs);info['uncompressed_bank']=len(active)
    if mix:
        for op,i,j in undo:
            if op=='cx':b.cx(public[j],public[i])
            else:b.swap(public[i],public[j])
    if gray:
        for word in words:
            for i in reversed(range(6)):b.cx(word[i+1],word[i])
    for word in words:
        for i in range(7):
            if relation>>i&1:b.cx(word[i],word[7])
    return b.finish_circuit(),info


def family(source,out,batch=7,select=1,lowbits=4,care=True,compress=False,gray=False,relation=64):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    (out/'submission.npz').write_bytes(Path(source).read_bytes())
    words=np.load(source)['angle_words'];rows=[]
    for k,d in enumerate(batch_deltas(words,batch)):
        c,info=build(d,select,lowbits,care,compress,gray,relation)
        c.to_file(str(out/f'angle_k_{k}.kmb'),format='kmb')
        rows.append(dict(k=k,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,erasure=info))
    result=dict(batch=batch,select=select,lowbits=lowbits,care=care,compress=compress,gray=gray,relation=relation,rows=rows,predicted=predicted(batch,rows))
    save(out.name+'_build.json',result);print(json.dumps(result),flush=True)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',default=str(INC/'submission.npz'));p.add_argument('--out',default=str(HERE/'candidates/pruned7'));p.add_argument('--batch',type=int,default=7);p.add_argument('--select',type=int,default=1);p.add_argument('--lowbits',type=int,default=8);p.add_argument('--zero-padding',action='store_true');p.add_argument('--compress',action='store_true');p.add_argument('--gray',action='store_true');p.add_argument('--relation',type=int,default=64);a=p.parse_args()
    family(a.source,a.out,a.batch,a.select,a.lowbits,not a.zero_padding,a.compress,a.gray,a.relation)
