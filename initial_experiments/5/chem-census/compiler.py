"""Compose banked physical state with a fresh structured ANGLE compilation.

Read-only chem-qrom hybrid feature grouping/quadratic interpolation, not its
physical NPZ or precomputed fitted models. Source paths are in NOTICE.md.
"""
from census_common import *
import sys
sys.path.insert(0,str(ROOT/'initial_experiments/5/chem-qrom'))
from hybrid import search
from build import build as linear_build
from quadratic import build as nonlinear_build
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas
from dfthcsos.scoreinfra.circuits import load_circuits,submission_sha256

def main():
    z=dict(np.load(BANK/'submission.npz'));old=z['angle_words']
    order,matrix,history=search(old)
    z['angle_words']=old[order];two=order[76:]-76
    z['U']=z['U'][:,two,:];z['W']=z['W'][:,two,:]
    out=HERE/'candidates/structured';out.mkdir(exist_ok=True)
    np.savez(out/'submission.npz',**z)
    write('structured_transform.json',dict(order=order.tolist(),matrix=matrix,features=history))
    inp=load_inputs(out/'submission.npz')
    metrics=evaluate_words(inp,z['angle_words'],z['tree_edges'],8)
    rows=[];cs=[]
    for k,d in enumerate(batch_deltas(z['angle_words'],4)):
        began=time.time()
        if k==0:c,info=linear_build(d,select=2,lowbits=7,mix=matrix)
        else:c,info=nonlinear_build(d,k)
        assert c.num_qubits<=128
        c.to_file(str(out/f'angle_k_{k}.kmb'),format='kmb');cs.append(c)
        row=dict(k=k,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,info=info,seconds=time.time()-began)
        rows.append(row);print(json.dumps(dict(k=k,Q=row['Q'],T=row['T'],seconds=row['seconds'])),flush=True)
        write('structured_build.json',dict(rows=rows,metrics=metrics,complete=len(rows)==19))
    rc=resource_costs(inp,cs,8,2*sum(r['T'] for r in rows));cost=metrics['lambda_eff']*rc['toffoli']*rc['qubits']
    digest=submission_sha256(load_circuits(out),out/'submission.npz',z['tree_edges'])
    record=dict(rows=rows,metrics=metrics,resources=rc,conditional_score=cost,improvement=BASE_COST/cost,
        complete=True,digest=digest,files={p.name:sha(p) for p in sorted(out.iterdir())})
    write('structured_build.json',record);print(json.dumps({k:v for k,v in record.items() if k!='rows'}),flush=True)
if __name__=='__main__':main()
