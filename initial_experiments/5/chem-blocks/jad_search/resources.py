"""Host-exact frontier plus constructive conservative signed ANGLE bounds."""
import argparse,csv,math
import numpy as np
from core import HERE,SPECTRUM,report,emit
from dfthcsos.verify import heuristic_costs
from dfthcsos.scoreinfra.costs_prepare import SpinFreeData_V1

_ev=np.load(SPECTRUM)['eigenvalues']
_positive=_ev[_ev>0]
_negative_squared=float(np.square(_ev[_ev<0]).sum())

def rank_bound(capacity):
    # Exact symmetric-feature Frobenius lower bound on any PSD Gram of this
    # rank. Negative reference eigenvalues are retained in the bound.
    return float(np.sqrt(_negative_squared+np.square(_positive[capacity:]).sum()))

def tree_cost(length,capacity):
    if length<=1:return 0
    if capacity==2:return 1
    h=capacity//2
    return tree_cost(length,h) if length<=h else 1+tree_cost(h,h)+tree_cost(length-h,h)

def bound(A,bits,J):
    n=(A-1).bit_length();m=n-1;active=J*(bits-1)
    needed=sum((1<<nn)-nn-1 for nn in [m//2,m-m//2])
    if needed>active:return None
    options=[]
    for low in range(1,m):
        high=m-low;L=(A+1)//2;pref=math.ceil(L/(1<<low))
        load=sum(tree_cost(min(1<<low,L-i*(1<<low)),1<<low) for i in range(pref))
        eq=pref*(max(0,high-1)+max(0,high-2))
        temp=active+max(low,high-1);q=n+bits*J+temp
        if q<=128:options.append(dict(angles_T_upper=2*math.ceil(75/J)*(load+eq+active+needed),angles_temp_upper=temp,oracle_Q_upper=q,lowbits=low))
    return min(options,key=lambda x:x['angles_T_upper']) if options else None

def row(R,B,C,bits,J):
    A=76+R*B;n=(A-1).bit_length()
    angles=bound(A,bits,J)
    if angles is not None:angles=dict(oracle_kind='signed_twoentry',**angles)
    # The host's ordinary QROM is also a valid construction. It may be
    # preferable when a larger J saves batches despite extra persistent width.
    direct=dict(oracle_kind='ordinary_qrom',angles_T_upper=2*math.ceil(75/J)*(A-2),angles_temp_upper=n-1,oracle_Q_upper=2*n-1+bits*J,lowbits=-1)
    if direct['oracle_Q_upper']<=128 and (angles is None or direct['angles_T_upper']<angles['angles_T_upper']):angles=direct
    if angles is None:return None
    h=heuristic_costs(76,R,B,C,bits);mod=SpinFreeData_V1(R,B,C,76,10,bits,optimize_qubits=True)
    k2=mod.find_qroam_size_to_minimize_global_temp_ancilla(math.ceil(math.log2(B))+11)
    T=h['toffoli']+angles['angles_T_upper'];Q=h['persistent_qubits']+J*bits+max(h['temporary_qubits'],angles['angles_temp_upper'])
    capacity=R*min(B,C);lower=rank_bound(capacity)
    return dict(R=R,B=B,C=C,brot=bits,J=J,k2=k2,rank_capacity=capacity,rank_error_lower_bound=lower,passes_necessary_rank_screen=lower<=.02,T_upper=T,Q_upper=Q,cost_conditional_upper=48.59417497035283*T*Q,**angles,**h)

def build(delta,bits,lowbits):
    import kickmix as km
    from build_selectswap_chem import emit_qrom
    from monomial_unlookup import erase
    A,J=delta.shape;n=(A-1).bit_length();b=km.CircuitBuilder()
    address=b.create_quantum_register(n,'selection');words=[b.create_quantum_register(bits,f'word_{j}') for j in range(J)]
    for w in words:b.cx(w[bits-2],w[bits-1])
    assert np.all((delta>>(bits-1)&1)==(delta>>(bits-2)&1))
    pub=km.array([q for w in words for q in w[:bits-1]]);width=len(pub)
    packed=[sum((int(row[j])&((1<<(bits-1))-1))<<((bits-1)*j) for j in range(J)) for row in delta]
    if A%2:packed.append(packed[-1])
    first=packed[::2];diff=[x^y for x,y in zip(first,packed[1::2])]
    active=[i for i in range(width) if any(v>>i&1 for v in diff)]
    m=n-1;needed=sum((1<<nn)-nn-1 for nn in [m//2,m-m//2])
    # Pad with zero lookup bits when a sparse actual table uses fewer bits
    # than the monomial erasure workspace. The resource bound already reserves
    # at least this many qubits.
    extra=b.alloc_qubits(max(len(active),needed));compressed=[sum(((v>>i)&1)<<j for j,i in enumerate(active)) for v in diff]
    table=[x|(v<<width) for x,v in zip(first,compressed)]
    emit_qrom(b,address[1:],km.array([*pub,*extra]),table,lowbits)
    for q,i in zip(extra,active):b.ccx(address[0],q,pub[i])
    erase(b,address[1:],extra,compressed+[0]*((1<<(n-1))-len(compressed)));b.free(extra)
    for w in words:b.cx(w[bits-2],w[bits-1])
    return b.finish_circuit()

def main():
    rows=[]
    for R in [1,2,3,4,5,6,8,10]:
     for B in [64,76,88,96,112,127,128,129,144,160,175,176,180,192,208,224,240,248,255,256,257,264,288,320,351,352,360,384,400,448,512,600,719]:
      if R*B>1000:continue
      for C in sorted({math.ceil(k/R) for k in [256,320,350,352,384,400,448]}):
       for bits in [7,8,9]:
        for J in range(4,16):
         r=row(R,B,C,bits,J)
         if r:rows.append(r)
    with (HERE/'results'/'frontier.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    pareto=[]
    for bits in [7,8,9]:
        best=math.inf
        for r in sorted([x for x in rows if x['brot']==bits and x['R']*x['C']>=350 and x['passes_necessary_rank_screen']],key=lambda x:(x['T_upper'],x['Q_upper'])):
            if r['Q_upper']<best:pareto.append(r);best=r['Q_upper']
    names=[];controls=[]
    for R,B,C in [(1,719,350),(2,351,175),(4,180,88),(4,176,88),(2,256,175),(2,257,175),(4,128,88),(4,129,88),(3,240,117)]:
      for bits in [7,8,9]:
        rr=[row(R,B,C,bits,J) for J in range(4,16)];rr=[r for r in rr if r]
        if rr:names.append(min(rr,key=lambda x:x['cost_conditional_upper']))
    # Compile dense independent signed tables at key cliffs, including unused
    # last words. Static counts upper-bound any sampled conditioned execution.
    rng=np.random.default_rng(1743)
    for A,bits,J in [(795,8,7),(796,7,8),(796,8,7),(588,8,7),(590,8,7),(796,9,6)]:
        bnd=bound(A,bits,J);d=rng.integers(0,1<<(bits-1),(A,J));d|=(d>>(bits-2)&1)<<(bits-1)
        c=build(d,bits,bnd['lowbits']);stat=c.operation_statistics();t=2*math.ceil(75/J)*stat.num_toffoli_like
        assert t<=bnd['angles_T_upper'] and c.num_qubits<=bnd['oracle_Q_upper'],(A,t,bnd)
        controls.append(dict(A=A,brot=bits,J=J,compiled_T_static_upper=t,compiled_Q=c.num_qubits,**bnd))
    data=dict(status='RESOURCE_SCREEN_ONLY',shape_count=len(rows),conditional_lambda_eff=48.59417497035283,physical_feasibility='UNTESTED; necessary PSD rank screen is not sufficient',necessary_psd_rank=next(k for k in range(len(_positive)) if rank_bound(k)<=.02),named=names,pareto=pareto,compiled_controls=controls,
        accounting='Exact host preparation and Givens costs, persistent plus maximum temporary; signed-gauge ANGLE full-activity upper bounds, all conditional Toffolis charged. No selected k2 override.')
    report('frontier.json',data)
    for r in names:emit({k:r[k] for k in ['R','B','C','brot','J','k2','T_upper','Q_upper','cost_conditional_upper']})
    emit(dict(compiled_controls=controls))

if __name__=='__main__':main()
