"""Variable-depth signed coordinates in input wires, with measured AND cleanup.

saved[j] is an unchanged copy of x[20-base-j]. The threshold flags encode
depth_flags[j] = (leaf_depth > base+j). Only K prefix bits are copied; no full
local-coordinate register is allocated. Each level uses one temporary AND.
"""
import kickmix as km

def transform(b,x,saved,depth_flags,base,inverse=False):
    K=len(saved);assert len(depth_flags)==K
    def level(j):
        pos=19-base-j;prev=saved[j];current=saved[j+1] if j+1<K else x[pos]
        if inverse:b.cx(depth_flags[j],x[pos])
        b.cx(current,prev);temp=b.alloc_qubits(1)[0]
        b.init_and(depth_flags[j],prev,temp)
        for bit in range(pos+1,21-base):b.cx(temp,x[bit])
        b.del_and(depth_flags[j],prev,temp);b.free(km.array([temp]));b.cx(current,prev)
        if not inverse:b.cx(depth_flags[j],x[pos])
    if not inverse:b.x(x[20-base])
    for j in (reversed(range(K)) if inverse else range(K)):level(j)
    if inverse:b.x(x[20-base])

def check():
    import json,time
    import numpy as np
    from pathlib import Path
    rng=np.random.default_rng(20261024);results=[];t0=time.time()
    bits=lambda xs,n:((xs[None,:]>>np.arange(n)[:,None])&1).astype(bool)
    for base in [5,6,7,8,9]:
      for K in range(5):
        b=km.CircuitBuilder();x=b.create_quantum_register(22,'x');f=b.create_quantum_register(K,'f');saved=b.alloc_qubits(K)
        for j,q in enumerate(saved):b.cx(x[20-base-j],q)
        transform(b,x,saved,f,base);forward=b.finish_circuit()
        # Fresh builder restores transformed coordinates, independently of forward.
        b=km.CircuitBuilder();xx=b.create_quantum_register(22,'x');ff=b.create_quantum_register(K,'f');ss=b.alloc_qubits(K)
        transform(b,xx,ss,ff,base,True)
        for j,q in enumerate(ss):b.cx(xx[20-base-j],q)
        inverse=b.finish_circuit()
        count=8192;values=rng.integers(2**21,2**22,count);depth=rng.integers(base,base+K+1,count)
        expected=values.copy();L=21-base;mask=(1<<L)-1
        for d in range(base,base+K+1):
            active=depth==d;low=21-d
            z=(values[active]&((1<<low)-1))-(1<<(low-1))
            expected[active]=(values[active]&~mask)|(z&mask)
        flagvalues=(1<<(depth-base))-1
        sim=km.Simulator(batch_size=count,seed=771+base+K);sim.use_same_registers_as(forward)
        sim.write_across_shots(x,bits(values,22))
        if K:sim.write_across_shots(f,bits(flagvalues,K))
        sim.do(forward)
        assert np.array_equal(sim.read_across_shots(x),bits(expected,22))
        assert np.array_equal(sim.read_across_shots(f),bits(flagvalues,K))
        assert not sim.read_phase_flipped_across_shots().any()
        sim.do(inverse)
        assert np.array_equal(sim.read_across_shots(x),bits(values,22))
        assert np.array_equal(sim.read_across_shots(f),bits(flagvalues,K))
        assert not sim.read_across_shots(km.array([km.q(i) for i in range(22+K,forward.num_qubits)])).any()
        assert not sim.read_phase_flipped_across_shots().any()
        results.append(dict(base=base,K=K,Q=forward.num_qubits,T_forward=forward.operation_statistics().num_toffoli_like,
            T_inverse=inverse.operation_statistics().num_toffoli_like,shots=count))
    report=dict(status='PASS',verified_shots=sum(r['shots'] for r in results),seconds=time.time()-t0,configurations=results)
    (Path(__file__).resolve().parent/'results/recenter_component.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))

if __name__=='__main__':check()
