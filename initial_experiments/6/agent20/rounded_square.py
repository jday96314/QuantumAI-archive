"""Carry-save XOR/AND network for round(u^2/64), 0<=u<=256.

The nine-bit promise coordinate uses the lower eight-bit square and a linear
endpoint bit. Outside the promise the reversible circuit still computes the
explicit total function round((u mod256)^2/64)+1024*bit8(u).
"""
from contextlib import contextmanager,ExitStack
from collections import Counter
from pathlib import Path
import argparse,json
import kickmix as km
import numpy as np


ZERO=(frozenset(),0)
ONE=(frozenset(),1)
def xor(*signals):
    terms=set();constant=0
    for ids,c in signals:terms.symmetric_difference_update(ids);constant^=c
    return frozenset(terms),constant


class Network:
    def __init__(self):self.nodes=[];self.intern={};self.ancestors={}
    def var(self,i):return frozenset({i}),0
    def product(self,a,b):
        if a==ZERO or b==ZERO:return ZERO
        if a==ONE:return b
        if b==ONE:return a
        if a==b:return a
        if xor(a,b)==ONE:return ZERO
        key=tuple(sorted((a,b),key=lambda v:(tuple(sorted(v[0])),v[1])))
        if key not in self.intern:
            index=8+len(self.nodes);self.intern[key]=index;self.nodes.append(key)
            parents=set(a[0])|set(b[0]);self.ancestors[index]=frozenset({index}).union(*(self.ancestors.get(i,{i}) for i in parents))
        return self.var(self.intern[key])
    def square(self):
        columns=[[] for _ in range(18)]
        for i in range(8):
            columns[2*i].append(self.var(i))
            for j in range(i+1,8):columns[i+j+1].append(self.product(self.var(i),self.var(j)))
        columns[5].append(ONE)
        for i in range(17):
            while len(columns[i])>=3:
                a,b,c=[columns[i].pop(0) for _ in range(3)]
                columns[i].append(xor(a,b,c));columns[i+1].append(xor(a,self.product(xor(a,b),xor(a,c))))
        output=[];carry=ZERO
        for i in range(18):
            signals=[v for v in [*columns[i],carry] if v!=ZERO]
            if len(signals)==3:
                a,b,c=signals;output.append(xor(a,b,c));carry=xor(a,self.product(xor(a,b),xor(a,c)))
            elif len(signals)==2:
                a,b=signals;output.append(xor(a,b));carry=self.product(a,b)
            elif signals:output.append(signals[0]);carry=ZERO
            else:output.append(ZERO);carry=ZERO
        return output[6:16]


class Emitter:
    def __init__(self,b,source,network,capacity):
        self.b,self.source,self.net,self.capacity=b,source,network,capacity
        self.work=b.alloc_qubits(2);self.cache={};self.pins=Counter();self.age={};self.tick=0;self.building=set();self.products=0;self.peak=0
    def touch(self,node):self.tick+=1;self.age[node]=self.tick
    def linear(self,s,target):
        terms,constant=s
        for i in terms:self.b.cx(self.source[i] if i<8 else self.cache[i][0],target)
        if constant:self.b.x(target)
    @contextmanager
    def signal(self,s):
        nodes=sorted(i for i in s[0] if i>=8)
        for i in nodes:self.ensure(i);self.pins[i]+=1;self.touch(i)
        try:yield s
        finally:
            for i in nodes:self.pins[i]-=1
    def remove(self,node):
        word=self.cache.pop(node);measurement=self.b.alloc_bits(1);self.b.hmr(word,measurement);self.b.free(word)
        left,right=self.net.nodes[node-8]
        with self.signal(left),self.signal(right):
            self.linear(left,self.work[0]);self.linear(right,self.work[1]);self.b.push_condition(measurement[0]);self.b.cz(*self.work);self.b.pop_condition()
            self.linear(right,self.work[1]);self.linear(left,self.work[0])
        self.b.bit_store(measurement,False);self.b.free(measurement)
    def room(self):
        while len(self.cache)>=self.capacity:
            choices=[i for i in self.cache if not self.pins[i] and not(self.net.ancestors[i]&self.building)]
            if not choices:return
            self.remove(min(choices,key=lambda i:(len(self.net.ancestors[i]),self.age[i])))
    def ensure(self,node):
        if node in self.cache:return
        self.building.add(node);left,right=self.net.nodes[node-8]
        with self.signal(left),self.signal(right):
            self.room();word=self.b.alloc_qubits(1);self.cache[node]=word;self.touch(node)
            self.linear(left,self.work[0]);self.linear(right,self.work[1]);self.b.init_and(*self.work,word[0]);self.products+=1
            self.linear(right,self.work[1]);self.linear(left,self.work[0]);self.peak=max(self.peak,len(self.cache))
            if self.products>100000:raise RuntimeError('bounded experiment: cache recomputation exceeded100000 ANDs')
        self.building.remove(node)
    def finish(self):
        while self.cache:self.remove(max(self.cache))
        self.b.free(self.work)


def apply(b,u,h,capacity=16):
    network=Network();outputs=network.square();emitter=Emitter(b,u,network,capacity)
    for signal,target in zip(outputs,h):
        with emitter.signal(signal):emitter.linear(signal,target)
    b.cx(u[8],h[10]);emitter.finish()
    return dict(and_nodes=len(network.nodes),executed_ands=emitter.products,peak_cache=emitter.peak)


def test(capacity):
    b=km.CircuitBuilder();u=b.create_quantum_register(9,'u');h=b.create_quantum_register(11,'h');extra=apply(b,u,h,capacity);c=b.finish_circuit()
    values=np.arange(512);rng=np.random.default_rng(2026103100+capacity);initial=rng.integers(0,2048,len(values));expected=initial^((((values&255)**2+32)>>6)+((values>>8)<<10))
    bits=lambda a,n:((a[None,:]>>np.arange(n)[:,None])&1).astype(bool)
    for seed in (20261031,20261101):
        sim=km.Simulator(batch_size=512,seed=seed);sim.use_same_registers_as(c);sim.write_across_shots(u,bits(values,9));sim.write_across_shots(h,bits(initial,11));sim.do(c)
        assert np.array_equal(sim.read_across_shots(u),bits(values,9))
        assert np.array_equal(sim.read_across_shots(h),bits(expected,11))
        assert not sim.read_phase_flipped_across_shots().any()
        assert not sim.read_across_shots(km.array([km.q(i) for i in range(20,c.num_qubits)])).any()
    row=dict(status='PASS_COMPONENT',capacity=capacity,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,shots=1024,**extra)
    c.to_file(str(Path('candidates')/f'rounded_square_cache{capacity}.kmb'),format='kmb');print(json.dumps(row),flush=True);return row


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--capacity',type=int,nargs='+',default=[12,16,24,32]);args=p.parse_args();rows=[]
    for capacity in args.capacity:
        rows.append(test(capacity));Path('results/rounded_square_components.json').write_text(json.dumps(rows,indent=2)+'\n')
