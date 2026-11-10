"""Erase a table via bipartite algebraic-normal-form phase repair.

Measurement turns table erasure into a Boolean phase. All monomials within
both address halves are shared; classically controlled CZs join the halves.
"""
import kickmix as km


def erase(b,address,target,table,extra_outcomes=()):
    assert len(table)==1<<len(address)
    coeff=list(map(int,table));n=len(address);split=n//2
    for i in range(n):
        for mask in range(1<<n):
            if mask>>i&1:coeff[mask]^=coeff[mask^(1<<i)]
    outcomes=b.alloc_bits(len(target));b.hmr(target,outcomes)
    outcomes=km.array([*outcomes,*extra_outcomes])
    needed=sum((1<<nn)-nn-1 for nn in [split,n-split])
    assert needed<=len(target),(needed,len(target))
    workspace=iter(target);records=[]
    def bank(bits):
        mon={0:True}
        for mask in range(1,1<<len(bits)):
            low=mask&-mask;index=low.bit_length()-1;rest=mask^low
            if not rest:mon[mask]=bits[index]
            else:
                q=next(workspace);b.init_and(mon[rest],bits[index],q)
                mon[mask]=q;records.append((mon[rest],bits[index],q))
        return mon
    low=bank(address[:split]);high=bank(address[split:]);flag=b.alloc_bits(1)[0]
    for mask,c in enumerate(coeff):
        if not c:continue
        selected=km.array([outcomes[j] for j in range(len(outcomes)) if c>>j&1])
        b.cx(selected,flag);b.push_condition(flag)
        b.cz(low[mask&((1<<split)-1)],high[mask>>split])
        b.pop_condition();b.cx(selected,flag)
    for a,c,q in reversed(records):b.del_and(a,c,q)
    b.bit_store(outcomes,False);b.free(outcomes);b.free(km.array([flag]))


def test():
    import numpy as np
    rng=np.random.default_rng(71092)
    for n in range(2,11):
        width=max(8,(1<<(n//2))+(1<<(n-n//2))-n-2)
        table=rng.integers(0,1<<min(width,62),1<<n).tolist()
        b=km.CircuitBuilder();address=b.create_quantum_register(n,'a');target=b.alloc_qubits(width)
        erase(b,address,target,table);c=b.finish_circuit();N=1<<n
        sim=km.Simulator(batch_size=N,seed=n);sim.use_same_registers_as(c)
        bits=lambda xs,k:((np.array(xs,dtype=np.int64)[None,:]>>np.arange(k)[:,None])&1).astype(bool)
        sim.write_across_shots(address,bits(np.arange(N),n));sim.write_across_shots(target,bits(table,width));sim.do(c)
        assert not sim.read_phase_flipped_across_shots().any()
        assert not sim.read_across_shots(target).any()
        assert np.array_equal(sim.read_across_shots(address),bits(np.arange(N),n))
        print('PASS',n,width,c.operation_statistics().num_toffoli_like,flush=True)

if __name__=='__main__':test()
