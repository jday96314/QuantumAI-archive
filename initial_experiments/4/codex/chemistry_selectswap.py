"""Nine-word QROM batches within the 128-qubit ANGLE cap.

Pebble four high address bits separately. At the bottom, apply a shared
CCX to a whole output vector by conjugating one target with CNOTs; this
saves the last unary-iteration ancilla without increasing its Toffoli cost.
"""
import argparse,json,shutil
from pathlib import Path
import numpy as np
import kickmix as km
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas


def build(delta,brot,lowbits=6):
    delta=np.asarray(delta);A,J=delta.shape;n=(A-1).bit_length()
    b=km.CircuitBuilder()
    address=b.create_quantum_register(n,'selection')
    words=[b.create_quantum_register(brot,f'word_{j}') for j in range(J)]
    lowbits=min(lowbits,n);prefix=address[lowbits:]
    def targets(row):
        return km.array([words[j][i] for j in range(J) for i in range(brot) if int(row[j])>>i&1])
    def load(row,ctrl):
        ts=targets(row)
        if len(ts):b.cx(ctrl,ts)
    def equal(ctrl,value,inverse=False):
        n=len(prefix)
        zeros=km.array([q for i,q in enumerate(prefix) if not (value>>i&1)])
        b.x(zeros)
        work=b.alloc_qubits(max(0,n-2))
        if inverse:
            m=b.alloc_bits(1);b.hmr(ctrl,m[0]);b.push_condition(m[0])
            if n==1:b.z(prefix[0])
            elif n==2:b.cz(prefix[0],prefix[1])
            else:
                for i in range(n-2):b.init_and(prefix[0] if i==0 else work[i-1],prefix[i+1],work[i])
                b.cz(work[-1],prefix[-1])
                for i in reversed(range(n-2)):b.del_and(prefix[0] if i==0 else work[i-1],prefix[i+1],work[i])
            b.pop_condition();b.bit_store(m,False);b.free(m)
        else:
            if n==1:b.cx(prefix[0],ctrl)
            else:
                for i in range(n-2):b.init_and(prefix[0] if i==0 else work[i-1],prefix[i+1],work[i])
                b.init_and(prefix[0] if n==2 else work[-1],prefix[-1],ctrl)
                for i in reversed(range(n-2)):b.del_and(prefix[0] if i==0 else work[i-1],prefix[i+1],work[i])
        b.free(work);b.x(zeros)
    for high in range(1<<len(prefix)):
        offset=high<<lowbits
        if offset>=A:break
        ctrl=b.alloc_qubits(1)[0]
        equal(ctrl,high)
        anc=b.alloc_qubits(lowbits-1)
        def node(c,level,lo,hi):
            if lo+offset>=A:return
            end=min(hi+offset,A)
            if hi-lo==1 or np.all(delta[lo+offset:end]==delta[lo+offset]):
                load(delta[lo+offset],c);return
            bit=address[lowbits-level-1]
            if hi-lo==2:
                row0=delta[lo+offset]
                load(row0,c)
                if lo+offset+1<A:
                    diff=targets(row0^delta[lo+offset+1])
                    if len(diff):
                        if len(diff)>1:b.cx(diff[0],diff[1:])
                        b.ccx(c,bit,diff[0])
                        if len(diff)>1:b.cx(diff[0],diff[1:])
                return
            mid=(lo+hi)//2
            if mid+offset>=A:
                node(c,level+1,lo,mid);return
            t=anc[level]
            b.x(bit);b.init_and(c,bit,t);b.x(bit)
            node(t,level+1,lo,mid)
            b.cx(c,t)
            node(t,level+1,mid,hi)
            b.del_and(c,bit,t)
        node(ctrl,0,0,1<<lowbits)
        b.free(anc)
        equal(ctrl,high,True)
        b.free(ctrl)
    return b.finish_circuit()


def main():
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,default=Path('from_host/sample_submission/dfthcsos/submission.npz'))
    p.add_argument('--output',type=Path,required=True);p.add_argument('--brot',type=int,default=12);p.add_argument('--batch',type=int,default=9)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    z=dict(np.load(a.input));np.savez(a.output/'submission.npz',**z)
    for i,d in enumerate(batch_deltas(z['angle_words'],a.batch)):
        c=build(d,a.brot);c.to_file(a.output/f'angle_k_{i}.kmb',format='kmb')
        print(i,c.num_qubits,c.operation_statistics().num_toffoli_like,flush=True)
if __name__=='__main__':main()
