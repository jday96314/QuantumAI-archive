"""Two-entry select-swap ANGLE XOR oracles with measured difference banks.

For address (h,l), load D0(h) into the public output and D0(h)^D1(h)
into clean workspace. The controlled difference selects D_l. Measuring and
unlooking the difference cleans the bank without touching arbitrary output z.
"""
import argparse, io, json, zipfile
from pathlib import Path
import numpy as np
import kickmix as km
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas


def emit_qrom(b,address,target,table,lowbits):
    n=len(address);A=len(table);lowbits=min(lowbits,n);prefix=address[lowbits:]
    def targets(row):return km.array([q for i,q in enumerate(target) if int(row)>>i&1])
    def load(row,control):
        qs=targets(row)
        if len(qs):b.cx(control,qs)
    def equal(control,value,inverse=False):
        nn=len(prefix);zeros=km.array([q for i,q in enumerate(prefix) if not(value>>i&1)])
        b.x(zeros);work=b.alloc_qubits(max(0,nn-2))
        if inverse:
            m=b.alloc_bits(1);b.hmr(control,m[0]);b.push_condition(m[0])
            if nn==1:b.z(prefix[0])
            elif nn==2:b.cz(prefix[0],prefix[1])
            else:
                for i in range(nn-2):b.init_and(prefix[0] if i==0 else work[i-1],prefix[i+1],work[i])
                b.cz(work[-1],prefix[-1])
                for i in reversed(range(nn-2)):b.del_and(prefix[0] if i==0 else work[i-1],prefix[i+1],work[i])
            b.pop_condition();b.bit_store(m,False);b.free(m)
        elif nn==1:b.cx(prefix[0],control)
        else:
            for i in range(nn-2):b.init_and(prefix[0] if i==0 else work[i-1],prefix[i+1],work[i])
            b.init_and(prefix[0] if nn==2 else work[-1],prefix[-1],control)
            for i in reversed(range(nn-2)):b.del_and(prefix[0] if i==0 else work[i-1],prefix[i+1],work[i])
        b.free(work);b.x(zeros)
    for high in range(1<<len(prefix)):
        offset=high<<lowbits
        if offset>=A:break
        control=b.alloc_qubits(1)[0];equal(control,high);anc=b.alloc_qubits(lowbits-1)
        def node(c,level,lo,hi):
            if lo+offset>=A:return
            end=min(hi+offset,A)
            if hi-lo==1 or all(table[j]==table[lo+offset] for j in range(lo+offset,end)):
                load(table[lo+offset],c);return
            bit=address[lowbits-level-1]
            if hi-lo==2:
                row=table[lo+offset];load(row,c)
                if lo+offset+1<A:
                    diff=targets(row^table[lo+offset+1])
                    if len(diff):
                        if len(diff)>1:b.cx(diff[0],diff[1:])
                        b.ccx(c,bit,diff[0])
                        if len(diff)>1:b.cx(diff[0],diff[1:])
                return
            mid=(lo+hi)//2
            if mid+offset>=A:node(c,level+1,lo,mid);return
            t=anc[level];b.x(bit);b.init_and(c,bit,t);b.x(bit)
            node(t,level+1,lo,mid);b.cx(c,t);node(t,level+1,mid,hi)
            b.del_and(c,bit,t)
        node(control,0,0,1<<lowbits)
        b.free(anc);equal(control,high,True);b.free(km.array([control]))


def build(delta,brot,lowbits=6,monomial=False):
    A,J=delta.shape;n=(A-1).bit_length();b=km.CircuitBuilder()
    address=b.create_quantum_register(n,'selection')
    words=[b.create_quantum_register(brot,f'word_{j}') for j in range(J)]
    public=km.array([q for w in words for q in w]);width=len(public)
    rows=[sum(int(row[j])<<(j*brot) for j in range(J)) for row in delta]
    if A%2:rows.append(rows[-1])
    first=rows[::2];diff=[a^c for a,c in zip(first,rows[1::2])]
    active=[i for i in range(width) if any(v>>i&1 for v in diff)]
    extra=b.alloc_qubits(len(active));compressed=[sum(((v>>i)&1)<<j for j,i in enumerate(active)) for v in diff]
    table=[a|(v<<width) for a,v in zip(first,compressed)]
    emit_qrom(b,address[1:],km.array([*public,*extra]),table,lowbits)
    for q,i in zip(extra,active):b.ccx(address[0],q,public[i])
    deletion=compressed+[0]*((1<<len(address[1:]))-len(compressed))
    if monomial:
        from monomial_unlookup import erase
        erase(b,address[1:],extra,deletion)
    else:b.del_lookup(table=deletion,address=address[1:],target=extra)
    b.free(extra)
    return b.finish_circuit()


def main():
    p=argparse.ArgumentParser();p.add_argument('--archive',default='baseline/submission.zip')
    p.add_argument('--out',default='candidates/selectswap7');p.add_argument('--batch',type=int,default=7)
    p.add_argument('--lowbits',type=int,default=6);p.add_argument('--monomial',action='store_true');a=p.parse_args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(a.archive) as arc:npz=arc.read('dfthcsos/submission.npz')
    (out/'submission.npz').write_bytes(npz)
    with np.load(io.BytesIO(npz)) as z:words=z['angle_words']
    for k,delta in enumerate(batch_deltas(words,a.batch)):
        c=build(delta,8,a.lowbits,a.monomial);c.to_file(str(out/f'angle_k_{k}.kmb'),format='kmb')
        print(json.dumps(dict(k=k,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like)),flush=True)
        assert c.num_qubits<=128

if __name__=='__main__':main()
