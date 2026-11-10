"""ANF measurement repair with one stored bank and one streamed bank.

The streamed bank visits each monomial once in a depth-first traversal. Its
parents stay live for free measurement cleanup. It uses the same AND count as
materializing both banks, but only logarithmic storage for the second bank.
"""
import argparse,json
from pathlib import Path
import kickmix as km


def erase(b,address,target,table,extra_outcomes=(),split=None):
    n=len(address);split=n//2 if split is None else split;other=n-split
    assert len(table)==1<<n and 0<=split<=n
    needed=(1<<split)-split-1+max(0,other-1)
    assert needed<=len(target),(needed,len(target))
    coeff=list(map(int,table))
    for i in range(n):
        for mask in range(1<<n):
            if mask>>i&1:coeff[mask]^=coeff[mask^(1<<i)]
    outcomes=b.alloc_bits(len(target));b.hmr(target,outcomes)
    outcomes=km.array([*outcomes,*extra_outcomes]);flag=b.alloc_bits(1)[0]
    low={0:True};records=[];at=0
    for mask in range(1,1<<split):
        bit=(mask&-mask).bit_length()-1;rest=mask^(1<<bit)
        if not rest:low[mask]=address[bit]
        else:
            q=target[at];at+=1;b.init_and(low[rest],address[bit],q)
            low[mask]=q;records.append((low[rest],address[bit],q))
    def phase(high,signal):
        for lm,left in low.items():
            c=coeff[lm|(high<<split)]
            if not c:continue
            selected=km.array([outcomes[j] for j in range(len(outcomes)) if c>>j&1])
            b.cx(selected,flag);b.push_condition(flag);b.cz(left,signal)
            b.pop_condition();b.cx(selected,flag)
    def visit(mask,start,signal,depth):
        phase(mask,signal)
        for i in range(start,other):
            if not mask:visit(1<<i,i+1,address[split+i],1)
            else:
                q=target[at+depth-1];b.init_and(signal,address[split+i],q)
                visit(mask|(1<<i),i+1,q,depth+1)
                b.del_and(signal,address[split+i],q)
    visit(0,0,True,0)
    for a,c,q in reversed(records):b.del_and(a,c,q)
    b.bit_store(outcomes,False);b.free(outcomes);b.free(km.array([flag]))


def main():
    import numpy as np
    p=argparse.ArgumentParser();p.add_argument('--maximum',type=int,default=11);a=p.parse_args()
    rng=np.random.default_rng(51739);rows=[]
    for n in range(2,a.maximum+1):
        low=n//2;needed=(1<<low)-low-1+max(0,n-low-1);width=max(8,needed)
        table=rng.integers(0,1<<min(width,61),1<<n,dtype=np.int64).tolist()
        b=km.CircuitBuilder();address=b.create_quantum_register(n,'a');target=b.alloc_qubits(width)
        erase(b,address,target,table);c=b.finish_circuit();shots=1<<n
        sim=km.Simulator(batch_size=shots,seed=n);sim.use_same_registers_as(c)
        bits=lambda v,k:((np.asarray(v,dtype=np.int64)[None,:]>>np.arange(k)[:,None])&1).astype(bool)
        sim.write_across_shots(address,bits(np.arange(shots),n));sim.write_across_shots(target,bits(table,width));sim.do(c)
        assert np.array_equal(sim.read_across_shots(address),bits(np.arange(shots),n))
        assert not sim.read_across_shots(target).any();assert not sim.read_phase_flipped_across_shots().any()
        r=dict(n=n,width=width,workspace=needed,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,shots=shots);rows.append(r);print(json.dumps(r),flush=True)
    Path('results/streamed_unlookup_check.json').write_text(json.dumps(dict(status='PASS',scope='all addresses, random tables, source restoration, zero output and phase',rows=rows),indent=2)+'\n')
if __name__=='__main__':main()
