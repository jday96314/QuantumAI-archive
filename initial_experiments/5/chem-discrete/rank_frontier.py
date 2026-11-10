"""Exact-word rank and host resource frontier; no continuous dictionary fits."""
import json,time
from core import *
from compact import Compact
from oracles import build
from dfthcsos.verify import heuristic_costs
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas
prior=ROOT/'initial_experiments/5/chem-qfit/candidates'
rows=[]
for tag in ['legacy_exact_split_B479','hardfreeze_B447','anneal_B671']:
    z,u,w,tree,words=load(prior/f'{tag}.npz');R,B,oldC=w.shape;m=Metric(tree)
    u=m.decode(m.tensor(words.reshape(R,B,75)));c=Compact(m,u)
    one=np.load(INC/'submission.npz')['angle_words'][:76]
    sw=encode(u.cpu().numpy().reshape(-1,76),tree,8,True)
    circuits=[build(delta) for delta in batch_deltas(np.concatenate([one,sw]),7)]
    at=2*sum(circ.operation_statistics().num_toffoli_like for circ in circuits);aq=max(ci.num_qubits for ci in circuits)
    for C in ([175,190,200,210,225,250,280,350] if R==2 else [350,360,375,400]):
        start=time.time();ww=torch.zeros((R,B,C),**m.d);ww[:,:,:min(oldC,C)]=m.tensor(w[:,:,:min(oldC,C)])
        f=c.factors(ww);f,err=c.project(f,30);ww=c.weights(f)
        exact=float(m.loss(u,ww).sqrt());assert abs(exact-err)<1e-7
        h=heuristic_costs(76,R,B,C,8);T=at+h['toffoli'];Q=h['persistent_qubits']+56+max(h['temporary_qubits'],aq-56-(76+R*B-1).bit_length())
        save(HERE/f'candidates/frontier_{tag}_C{C}.npz',u,ww,tree,words,err2=exact)
        row=dict(tag=tag,R=R,B=B,C=C,brot=8,error_2=exact,T=T,Q=Q,oracle_Q=aq,lambda_break_even=708547388.9553/T/Q,lambda_110=708547388.9553/1.1/T/Q,seconds=time.time()-start)
        rows.append(row);print(json.dumps(row),flush=True)
        (HERE/'results/rank_frontier.json').write_text(json.dumps(rows,indent=2))
