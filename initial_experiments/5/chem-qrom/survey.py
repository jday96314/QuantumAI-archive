"""Compile, natively simulate, and price lookup/batch alternatives."""
from common import *
from build import build,family
from algebraic import family as anf_family
from check import component,check
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas

def run():
    words=np.load(INC/'submission.npz')['angle_words'];grid=[]
    for j in range(3,9):
        for s in [1,2]:
            d=batch_deltas(words,j)[0]
            for l in range(3,10-s):
                c,e=build(d,s,l)
                row=dict(J=j,select=s,lowbits=l,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,erasure=e)
                # Over-cap circuits are transient research measurements,
                # never saved as valid submission candidates.
                row['native_component']=component(c,d,2026102405+len(grid));grid.append(row)
    save('batch_component_grid.json',grid)
    for batch,select in [(3,2),(4,1),(5,1),(6,1)]:
        out=HERE/'candidates'/f'select{1<<select}_j{batch}'
        family(INC/'submission.npz',out,batch,select,9-select)
        check(out)
    anf_family(7);check(HERE/'candidates/anf7')
    check(HERE/'candidates/baseline')

if __name__=='__main__':run()
