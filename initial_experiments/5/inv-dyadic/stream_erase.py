"""ANF measurement-phase repair with one streamed bank of monomials.

Retain all monomials of the smaller address half. Visit monomials of the
larger half depth first and apply their full row of phase corrections before
erasing them. Gate count equals two resident banks; workspace is much smaller.
"""
import kickmix as km

def workspace_size(n):
    low=n//2;high=n-low
    return (1<<low)-low-1+max(0,high-1)

def erase(b,address,target,table,extra_outcomes=()):
    n=len(address);split=n//2;assert len(table)==1<<n
    coeff=list(map(int,table))
    for i in range(n):
        for mask in range(1<<n):
            if mask>>i&1:coeff[mask]^=coeff[mask^(1<<i)]
    outcomes=b.alloc_bits(len(target));b.hmr(target,outcomes);outcomes=km.array([*outcomes,*extra_outcomes])
    assert workspace_size(n)<=len(target)
    low={0:True};records=[];cursor=0
    for mask in range(1,1<<split):
        bit=mask&-mask;rest=mask^bit;q=address[bit.bit_length()-1]
        if not rest:low[mask]=q
        else:
            out=target[cursor];cursor+=1;b.init_and(low[rest],q,out)
            low[mask]=out;records.append((low[rest],q,out))
    flag=b.alloc_bits(1)[0]
    def apply(mask,high):
        for lm,lq in low.items():
            c=coeff[lm|(mask<<split)]
            if not c:continue
            selected=km.array([outcomes[j] for j in range(len(outcomes)) if c>>j&1])
            b.cx(selected,flag);b.push_condition(flag);b.cz(lq,high);b.pop_condition();b.cx(selected,flag)
    apply(0,True)
    highbits=address[split:]
    def visit(mask,start,parent,depth):
        for j in range(start,len(highbits)):
            mm=mask|(1<<j)
            if not mask:
                apply(mm,highbits[j]);visit(mm,j+1,highbits[j],depth)
            else:
                out=target[cursor+depth];b.init_and(parent,highbits[j],out)
                apply(mm,out);visit(mm,j+1,out,depth+1);b.del_and(parent,highbits[j],out)
    visit(0,0,True,0)
    for a,c,q in reversed(records):b.del_and(a,c,q)
    b.bit_store(outcomes,False);b.free(outcomes);b.free(km.array([flag]))

def check():
    import json,random,time
    from pathlib import Path
    import numpy as np
    from monomial_unlookup import erase as old_erase
    rng=random.Random(20261024);rows=[];t0=time.time()
    for n in range(3,12):
        width=max(8,workspace_size(n));table=[rng.getrandbits(width) for _ in range(1<<n)]
        b=km.CircuitBuilder();a=b.create_quantum_register(n,'a');t=b.alloc_qubits(width)
        erase(b,a,t,table);c=b.finish_circuit();count=1<<n
        sim=km.Simulator(batch_size=count,seed=n);sim.use_same_registers_as(c)
        sim.write_across_shots(a,np.array([[(v>>i)&1 for v in range(count)] for i in range(n)],dtype=bool))
        sim.write_across_shots(t,np.array([[(v>>i)&1 for v in table] for i in range(width)],dtype=bool));sim.do(c)
        assert not sim.read_across_shots(t).any()
        assert not sim.read_phase_flipped_across_shots().any()
        assert np.array_equal(sim.read_across_shots(a),np.array([[(v>>i)&1 for v in range(count)] for i in range(n)],dtype=bool))
        old_workspace=sum((1<<nn)-nn-1 for nn in [n//2,n-n//2])
        assert c.operation_statistics().num_toffoli_like==old_workspace
        rows.append(dict(address_bits=n,stream_workspace=width,bank_workspace=old_workspace,
            T=c.operation_statistics().num_toffoli_like,Q=c.num_qubits,shots=count))
    report=dict(status='PASS',verified_shots=sum(r['shots'] for r in rows),seconds=time.time()-t0,configurations=rows)
    (Path(__file__).resolve().parent/'results/stream_erase_component.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))

if __name__=='__main__':check()
