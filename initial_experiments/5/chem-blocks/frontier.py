"""Exact host PREPARE/Givens + conservative, realizable signed-ANGLE bounds."""
import csv
import math
import numpy as np
from common import HERE, SOURCE, SPECTRUM, report
from signed_oracle import bound,build
from dfthcsos.scoreinfra.costs_prepare import SpinFreeData_V1
from dfthcsos.verify import heuristic_costs

def price(R,B,C,bits,J):
    A=76+R*B
    signed=bound(A,bits,J)
    n=(A-1).bit_length()
    # Ordinary unary QROM is a constructive fallback for address spaces whose
    # ANF difference-erasure bank is too wide, or for larger word batches.
    plain=dict(T_per_oracle=A-2,Q_oracle=2*n-1+J*bits,angle_temp=n-1,lowbits=0,ANF_monomials=0)
    options=[]
    if signed is not None:options.append(('signed_selectswap',signed))
    if plain['Q_oracle']<=128:options.append(('ordinary_QROM',plain))
    if not options:return
    model=SpinFreeData_V1(R,B,C,76,10,bits,optimize_qubits=True)
    h=heuristic_costs(76,R,B,C,bits)
    k2=model.find_qroam_size_to_minimize_global_temp_ancilla(math.ceil(math.log2(B))+11)
    K=math.ceil(75/J)
    def total(ob):
        T=h['toffoli']+2*K*ob['T_per_oracle']
        Q=h['persistent_qubits']+J*bits+max(h['temporary_qubits'],ob['angle_temp'])
        return T,Q
    name,oracle=min(options,key=lambda x:np.prod(total(x[1])))
    T,Q=total(oracle)
    return dict(R=R,B=B,C=C,brot=bits,J=J,k2=k2,oracle_type=name,**oracle,**h,T_total_bound=T,Q_total_bound=Q,
                conditional_cost_at_incumbent_lambda=48.59417497035283*T*Q,
                conditional_improvement=708547388.9553173/(48.59417497035283*T*Q),
                physical_feasibility='UNTESTED')

def main():
    sp=np.load(SPECTRUM);ev=sp['eigenvalues']
    # Necessary best unrestricted PSD rank bound; includes the negative tail.
    def rank_error(k):return float(np.linalg.norm(np.r_[np.minimum(ev[:k],0),ev[k:]]))
    rows=[]
    for R in [1,2,3,4,5,6,8,10,12,16,24,28,32,48,64]:
      for B in [12,16,24,32,40,48,60,64,76,88,96,112,128,144,160,175,176,180,192,224,240,255,256,257,272,320,351,352,384,448,512,640,719]:
       for C in sorted(set([16,24,32,44,60,64,76,88,96,112,128,144,175,176,180,224,256,300,350,B//2,B])):
        if C>B or R*C<250 or R*C>1800:continue
        err=rank_error(R*C)
        if err>.0200001:continue
        for bits in [7,8,9]:
         for J in range(4,17):
          x=price(R,B,C,bits,J)
          if x is not None:x['unrestricted_rank_error_bound']=err;rows.append(x)
    path=HERE/'results/frontier.csv'
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    chosen=[]
    for R,B,C in [(1,719,350),(2,351,175),(2,351,160),(4,180,88),(3,240,117),(4,192,96),(4,256,88),(4,257,88),(2,256,175),(3,256,117)]:
      for bits in [7,8,9]:
        xs=[price(R,B,C,bits,J) for J in range(4,17)]
        chosen.append(min((x for x in xs if x),key=lambda x:x['conditional_cost_at_incumbent_lambda']))
    best=[];pareto=[]
    for bits in [7,8,9]:
        xs=[x for x in rows if x['brot']==bits]
        best+=sorted(xs,key=lambda x:x['conditional_cost_at_incumbent_lambda'])[:10]
        bestQ=float('inf')
        for x in sorted(xs,key=lambda x:(x['T_total_bound'],x['Q_total_bound'])):
            if x['Q_total_bound']<bestQ:
                pareto.append(x);bestQ=x['Q_total_bound']
    # Check bound against actual compiled operation counts, with fully active
    # random tables; these are oracle controls, not physical candidates.
    controls=[];rng=np.random.default_rng(847291)
    for A,bits,J in [(795,8,7),(796,7,8),(778,8,7),(1100,9,6),(844,7,8)]:
        ob=bound(A,bits,J)
        d=rng.integers(0,1<<(bits-1),(A,J),dtype=np.int64)
        d|=((d>>(bits-2))&1)<<(bits-1)
        c=build(d,bits,ob['lowbits'])
        T=int(c.operation_statistics().num_toffoli_like)
        assert T<=ob['T_per_oracle'] and c.num_qubits<=ob['Q_oracle']
        controls.append(dict(A=A,brot=bits,J=J,measured_T=T,measured_Q=c.num_qubits,bound=ob))
    report(HERE/'results/frontier.json',dict(status='RESOURCE_FRONTIER_ONLY',rows=len(rows),
        rank_minimum=next(k for k in range(1,len(ev)) if rank_error(k)<=.02),
        assumptions='Exact unmodified host with optimize_qubits=True, chi10; conservative explicit oracle upper bounds. Costs use incumbent lambda conditionally; no physical feasibility or win claimed.',
        selected_shapes=chosen,resource_best=best,pareto_T_Q=pareto,compiled_oracle_controls=controls))

if __name__=='__main__':main()
