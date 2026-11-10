"""Exact XOR-AND lookup networks with measured AND deletion.

Positive Davio decomposition f=f0 XOR x*(f0 XOR f1), sharing all affine
functions and repeated nonlinear nodes across output bits. Search concerns
resource counts only; no simulation seeds or salts enter synthesis.
"""
import itertools
import json
import random
from functools import lru_cache

import kickmix as km


class Network:
    def __init__(self, n, order):
        self.n=n; self.order=order; self.full=(1<<(1<<n))-1
        self.planes=[sum(((j>>i)&1)<<j for j in range(1<<n)) for i in range(n)]
        self.nodes=[]; self.products={}; self.memo={}
        for subset in range(1<<n):
            value=0; signal=0
            for i in range(n):
                if subset>>i&1: value^=self.planes[i]; signal^=1<<(i+1)
            self.memo[value]=signal; self.memo[value^self.full]=signal^1

    def function(self, value):
        if value in self.memo: return self.memo[value]
        for variable in self.order:
            step=1<<variable; plane=self.planes[variable]
            lo=value&(self.full^plane); hi=value&plane
            f0=lo|(lo<<step); f1=hi|(hi>>step)
            if f0!=f1: break
        else: raise AssertionError('nonconstant without a variable')
        base=self.function(f0); difference=self.function(f0^f1)
        if difference&1:
            base^=1<<(variable+1); difference^=1
        if not difference:
            self.memo[value]=base; self.memo[value^self.full]=base^1
            return base
        key=(variable,difference)
        if key not in self.products:
            self.products[key]=1<<(self.n+1+len(self.nodes))
            self.nodes.append(key)
        result=base^self.products[key]
        self.memo[value]=result; self.memo[value^self.full]=result^1
        return result

    def compile(self, table, width):
        self.outputs=[self.function(sum(((int(v)>>i)&1)<<j for j,v in enumerate(table))) for i in range(width)]
        return self

    def description(self):
        return dict(n=self.n,order=list(self.order),nodes=self.nodes,outputs=self.outputs)


@lru_cache(None)
def synthesize(table, width, trials=256):
    n=(len(table)-1).bit_length(); assert len(table)==1<<n
    rng=random.Random(20261025)
    orders=[tuple(range(n)),tuple(reversed(range(n)))]
    for _ in range(trials):
        order=list(range(n));rng.shuffle(order);orders.append(tuple(order))
    best=None
    for order in dict.fromkeys(orders):
        network=Network(n,order).compile(table,width)
        if best is None or len(network.nodes)<len(best.nodes): best=network
    return best


def prepare(b,address,target,table,workspace=(),limit=None):
    network=synthesize(tuple(map(int,table)),len(target))
    needed=len(network.nodes)+1
    borrowed=list(workspace[:needed])
    extra=needed-len(borrowed)
    if limit is not None and b.num_allocated_qubits+extra>limit:
        raise ValueError(f'Davio lookup needs {extra} extra wires at live {b.num_allocated_qubits}')
    owned=b.alloc_qubits(extra); work=km.array([*borrowed,*owned])
    temporary=work[-1]; values=km.array([*address,*work[:-1]])
    def linear(signal, target):
        if signal&1: b.x(target)
        signal>>=1
        for i in range(len(values)):
            if signal>>i&1: b.cx(values[i],target)
    for i,(variable,signal) in enumerate(network.nodes):
        linear(signal,temporary)
        b.init_and(address[variable],temporary,work[i])
        linear(signal,temporary)
    for signal,out in zip(network.outputs,target): linear(signal,out)
    for i,(variable,signal) in reversed(list(enumerate(network.nodes))):
        linear(signal,temporary)
        b.del_and(address[variable],temporary,work[i])
        linear(signal,temporary)
    b.free(owned)
    return dict(ANDs=len(network.nodes),workspace=needed)


if __name__=='__main__':
    import numpy as np
    from affine_invsqrt import table
    rows=table();reports={}
    for name,values,width in [('delayed2',(rows[:,0]>>21)&3,2),('curve_high3',rows[:,2]>>4,3)]:
        network=synthesize(tuple(map(int,values)),width,2048)
        b=km.CircuitBuilder();a=b.create_quantum_register(9,'a');t=b.create_quantum_register(width,'t')
        prepare(b,a,t,values)
        c=b.finish_circuit();xs=np.arange(512,dtype=np.uint64)
        bits=lambda v,w:((v[None,:]>>np.arange(w,dtype=np.uint64)[:,None])&1).astype(bool)
        for seed in (20261025,20261026):
            sim=km.Simulator(batch_size=512,seed=seed);sim.use_same_registers_as(c)
            incoming=(xs*3)%(1<<width)
            sim.write_across_shots(a,bits(xs,9));sim.write_across_shots(t,bits(incoming,width));sim.do(c)
            assert np.array_equal(sim.read_across_shots(t),bits(incoming^values.astype(np.uint64),width))
            assert np.array_equal(sim.read_across_shots(a),bits(xs,9))
            scratch=km.array([km.q(i) for i in range(9+width,c.num_qubits)])
            assert not sim.read_across_shots(scratch).any()
            assert not sim.read_phase_flipped_across_shots().any()
        reports[name]=dict(status='NATIVE_PASS',shots=1024,ANDs=len(network.nodes),
                           emitted_Toffolis=c.operation_statistics().num_toffoli_like,network=network.description())
        print(json.dumps({k:v for k,v in reports[name].items() if k!='network'}),flush=True)
    from pathlib import Path
    Path('results/davio_lookups.json').write_text(json.dumps(reports,indent=2)+'\n')
