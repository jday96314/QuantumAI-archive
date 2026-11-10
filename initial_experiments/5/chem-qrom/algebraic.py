"""Direct multi-output ANF: shared monomials, no difference bank.

This is a separate lookup architecture from select/QROM. Monomials are
streamed and each controls the full associated XOR output vector.
"""
import argparse
import kickmix as km
from common import *
from erase import coefficients,closure
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas

def build(delta):
    A,J=delta.shape;n=(A-1).bit_length();b=km.CircuitBuilder()
    address=b.create_quantum_register(n,'selection')
    words=[b.create_quantum_register(8,f'word_{j}') for j in range(J)]
    target=km.array([q for w in words for q in w])
    table=[sum(int(x)<<(8*j) for j,x in enumerate(row)) for row in delta]
    coeff=coefficients(table,n);needed=closure(m for m,c in enumerate(coeff) if c)
    work=b.alloc_qubits(max(m.bit_count()-1 for m in needed))
    def apply(m,q):
        value=coeff[m]
        if value:b.cx(q,km.array([t for i,t in enumerate(target) if value>>i&1]))
    apply(0,True)
    def visit(mask,start,parent,depth):
        for j in range(start,n):
            mm=mask|(1<<j)
            if mm not in needed:continue
            if mask:
                q=work[depth];b.init_and(parent,address[j],q)
                apply(mm,q);visit(mm,j+1,q,depth+1);b.del_and(parent,address[j],q)
            else:apply(mm,address[j]);visit(mm,j+1,address[j],depth)
    visit(0,0,True,0);b.free(work)
    return b.finish_circuit()

def family(batch):
    source=INC/'submission.npz';out=HERE/'candidates'/f'anf{batch}';out.mkdir(exist_ok=True)
    (out/'submission.npz').write_bytes(source.read_bytes());rows=[]
    for k,d in enumerate(batch_deltas(np.load(source)['angle_words'],batch)):
        c=build(d);c.to_file(str(out/f'angle_k_{k}.kmb'),format='kmb')
        rows.append(dict(k=k,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like))
    data=dict(strategy='shared streamed ANF',rows=rows,predicted=predicted(batch,rows))
    save(out.name+'_build.json',data);print(json.dumps(data),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--batch',type=int,default=13);a=p.parse_args();family(a.batch)
