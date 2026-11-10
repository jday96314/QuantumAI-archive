"""Native calibration, exact threshold map, rank bounds and integer deletion.

Greedy deletion criterion follows chem-qfit/greedy_prune.py. Every saved
dictionary contains unchanged exact incumbent words. Fresh native signed
pair-space target is used, including all negative eigenvalues.
"""
from census_common import *
from dataclasses import replace
from dfthcsos.scoreinfra.costs_prepare import SpinFreeData_V1
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas
from oracles import build

def spectral(m,u,C):
    P=m.feats(u)[0];Q,R=torch.linalg.qr(P.T);H=Q.T@m.G@Q
    ev,V=torch.linalg.eigh((H+H.T)/2)
    outside=float((m.G-Q@H@Q.T).square().sum())
    tail=ev.clone();tail[-C:]=tail[-C:].clamp(max=0)
    bound=(outside+float(tail.square().sum()))**.5
    w=torch.linalg.solve(R,V[:,-C:]*ev[-C:].clamp(min=0).sqrt())[None]
    return bound,w,ev

def main():
    began=time.time();z,m,u,w,off,words=initial()
    inp=load_inputs(BANK/'submission.npz')
    native=evaluate_words(inp,z['angle_words'],m.tree,8)
    pair=quantities(m,u,w,off)
    assert abs(pair['error2']-native['error_2'])<2e-9
    write('calibration.json',dict(native=native,pair=pair,source_npz_sha256=sha(BANK/'submission.npz')))
    print(json.dumps(dict(kind='calibration',native=native,pair=pair)),flush=True)
    threshold=[]
    for B in range(479,720):
        h=heuristic_costs(76,1,B,350,8);model=SpinFreeData_V1(1,B,350,76,10,8)
        threshold.append(dict(B=B,prep=h['inner_state_prep_toffoli'],Q_J7=h['qubits']+56,
            k2=model.find_qroam_size_to_minimize_global_temp_ancilla(int(np.ceil(np.log2(B)))+11)))
    jumps=[dict(before=p,after=q,prep_difference=q['prep']-p['prep'],Q_difference=q['Q_J7']-p['Q_J7'])
        for p,q in zip(threshold,threshold[1:]) if p['prep']!=q['prep'] or p['Q_J7']!=q['Q_J7']]
    write('thresholds.json',dict(rows=threshold,jumps=jumps))
    ranks=[]
    with torch.no_grad():
        for C in [350,349,348,344,340,336,320,300,280,256,224]:
            bound,_,_=spectral(m,u,C);h=heuristic_costs(76,1,719,C,8)
            ranks.append(dict(C=C,rank_optimal_error2=bound,T=h['toffoli']+10604,Q=h['qubits']+56,
                conditional_cost=(h['toffoli']+10604)*(h['qubits']+56)*BASE_LAMBDA))
            print(json.dumps(dict(kind='rank',**ranks[-1])),flush=True)
    write('rank_capacity.json',dict(rows=ranks,note='Exact fixed-dictionary PSD rank-C optimum; not a shape-wide bound'))
    precision=[]
    for mode in ['onebody_even','all_even']:
        aw=z['angle_words'].copy();part=aw[:76] if mode=='onebody_even' else aw
        part[:]=np.rint(part/2).astype(int)*2%256
        mt=evaluate_words(inp,aw,m.tree,8)
        precision.append(dict(mode=mode,metrics=mt))
    write('precision.json',precision)
    targets={711,703,687,671,655,639,607,575,543,511};rows=[]
    fullP=m.feats(u)[0];fullX=w[0].T@fullP;indices=np.arange(719)
    # Single-row pruning, recomputing exact covariance and leverage each time.
    with torch.no_grad():
        while u.shape[1]>511:
            P=m.feats(u)[0];S=P@P.T;inv=torch.linalg.inv(S);K=w[0]@w[0].T;d=inv.diag()
            score=2*(K@S@K).diag()/d-(K.diag()/d).square()
            ix=torch.arange(u.shape[1],device=u.device);ix=ix[ix!=score.argmin()]
            u=u[:,ix];w=w[:,ix];indices=indices[ix.cpu().numpy()];words=words[ix.cpu().numpy()]
            w,err=m.project(u,w,1);B=u.shape[1]
            if B in targets:
                store(f'prune_B{B}',m,u,w,off,words,indices=indices)
                P=m.feats(u)[0]
                transfer=torch.linalg.lstsq(P.T,fullX.T).solution[None]
                store(f'transfer_B{B}',m,u,transfer,off,words,indices=indices)
                bound,_,_=spectral(m,u,350)
                # Full shared-span relaxation also permits indefinite/rank-free coefficients.
                Q,_=torch.linalg.qr(P.T);outside=float((m.G-Q@(Q.T@m.G@Q)@Q.T).norm())
                o=off.clone();m.project_gap(u,w,o)
                row=dict(B=B,rank_bound=bound,span_bound=outside,unconstrained=quantities(m,u,w,o),
                    transfer=quantities(m,u,transfer,off),seconds=time.time()-began)
                rows.append(row);print(json.dumps(dict(kind='prune',**row)),flush=True)
                store(f'projected_B{B}',m,u,w,o,words,indices=indices)
                write('continuation.json',dict(rows=rows,completed=B==511))
    # Compile complete existing oracle family for measured host resource points.
    costs=[]
    for B in [719,703,671,639,575,543,511]:
        aw=z['angle_words'] if B==719 else np.concatenate([z['angle_words'][:76],np.load(HERE/'parameters'/f'prune_B{B}.npz')['angle_words']])
        cs=[build(d) for d in batch_deltas(aw,7)]
        ww=np.zeros((1,B,350));this=replace(inp,w=ww)
        r=resource_costs(this,cs,8,2*sum(c.operation_statistics().num_toffoli_like for c in cs))
        row=dict(B=B,resources=r,conditional_cost=BASE_LAMBDA*r['toffoli']*r['qubits'],
            maximum_lambda_for_goal=TARGET/r['toffoli']/r['qubits'])
        costs.append(row);print(json.dumps(dict(kind='resources',**row)),flush=True)
    write('continuation_costs.json',costs)
if __name__=='__main__':main()
