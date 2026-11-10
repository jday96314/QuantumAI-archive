"""Host state-prep costs and selected complete signed oracle resource screen."""
import json
from core import *
from dfthcsos.verify import heuristic_costs
from dfthcsos.scoreinfra.costs_prepare import SpinFreeData_V1
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas
from oracles import build
rows=[]
shapes=[(1,b,350) for b in [511,543,575,600,639,671,700,719]]+[(2,447,c) for c in [160,175,190]]+[(2,479,c) for c in [160,175,190]]+[(2,351,175)]
for R,B,C in shapes:
    if R==1:
        src=INC/'submission.npz' if B==719 else HERE/f'candidates/pruned_B{B}.npz'
        if B==575:src=HERE/'candidates/pruned_B575.npz'
    else:src=ROOT/f'initial_experiments/4/claude/chem/snap_B{B}.npz' if B!=351 else ROOT/'initial_experiments/4/claude/chem/L_r2b351c175.npz'
    z,u,w,tree,words=load(src)
    if words is None:words=encode(u.reshape(-1,76),tree,8,True)
    if len(words)==R*B:words=np.concatenate([np.load(INC/'submission.npz')['angle_words'][:76],words])
    records=[]
    for delta in batch_deltas(words,7):
        circ=build(delta,8);records.append(dict(Q=circ.num_qubits,T=circ.operation_statistics().num_toffoli_like))
    h=heuristic_costs(76,R,B,C,8);model=SpinFreeData_V1(R,B,C,76,10,8)
    T=2*sum(x['T'] for x in records)+h['toffoli'];Q=h['persistent_qubits']+56+max(h['temporary_qubits'],max(x['Q'] for x in records)-((76+R*B-1).bit_length())-56)
    row=dict(R=R,B_per_block=B,total_vectors=R*B,C=C,k2=model.find_qroam_size_to_minimize_global_temp_ancilla(int(np.ceil(np.log2(B)))+11),host=h,oracles=records,T_static=T,Q=Q,lambda_required_567M=567e6/T/Q,estimated_at_incumbent_lambda=T*Q*48.59417497035283)
    rows.append(row);print(json.dumps(row),flush=True)
(HERE/'results/shape_costs.json').write_text(json.dumps(rows,indent=2)+'\n')
