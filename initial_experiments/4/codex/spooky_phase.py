"""Bounded-cache AND graph phase emission with delayed HMR phase repair.

Eviction measures a computed bit and records its phase as a classical ghost.
Recomputation pays that ghost with a conditional Z. Final cleanup processes
nodes in reverse dependency order, so all ghosts are eventually discharged.
This compiler does not assume erased AND parents remain live.
"""
from collections import Counter
from fractions import Fraction
import kickmix as km


class CapacityError(ValueError):pass


def plan(net,budget,limit=100000):
    parents={net.n+i:sorted({j for s in pair for j in s.terms if j>=net.n}) for i,pair in enumerate(net.nodes)}
    live=set();ghost=set();pins=Counter();age={};tick=0;actions=[];computes=0;peak=0
    def touch(j):
        nonlocal tick
        tick+=1;age[j]=tick
    def room():
        if len(live)<budget:return
        choices=[j for j in live if not pins[j]]
        if not choices:raise CapacityError(f'all {budget} live nodes pinned')
        j=min(choices,key=lambda j:age[j]);live.remove(j);ghost.add(j);actions.append(('ghost',j))
    def ensure(j):
        nonlocal computes,peak
        if j<net.n:return
        if j in live:touch(j);return
        for p in parents[j]:ensure(p);pins[p]+=1
        room();actions.append(('compute',j));live.add(j);ghost.discard(j);touch(j);computes+=1;peak=max(peak,len(live))
        for p in parents[j]:pins[p]-=1
        if computes>limit:raise CapacityError(f'computation limit {limit}')
    for index,out in enumerate(net.outputs):
        actions.append(('phase_start',index))
        for j in sorted(out.terms):ensure(j);actions.append(('phase_xor',j))
        actions.append(('phase_finish',index))
        ghost.update(j for j in out.terms if j>=net.n and j not in live)
    while live or ghost:
        j=max(live|ghost);ensure(j);pins[j]+=1
        for p in parents[j]:ensure(p);pins[p]+=1
        actions.append(('clear',j));live.remove(j);ghost.discard(j)
        for p in parents[j]:pins[p]-=1
        pins[j]-=1
    return actions,dict(computes=computes,peak=peak,ghost_events=sum(a[0]=='ghost' for a in actions),actions=len(actions))


def emit(net,b,address,budget,actions=None):
    if actions is None:actions,stats=plan(net,budget)
    else:stats=None
    assert len(address)==net.n
    A,B,parity=b.alloc_qubits(3);live={};pending=set();ghost=b.alloc_bits(len(net.nodes));measurement=b.alloc_bits(1)[0]
    def linear(signal,target):
        if signal.const:b.x(target)
        for j in sorted(signal.terms):b.cx(address[j] if j<net.n else live[j],target)
    def phase_bit(measurement,j):
        if j<net.n:
            b.push_condition(measurement);b.z(address[j]);b.pop_condition()
        elif j in live:
            b.push_condition(measurement);b.z(live[j]);b.pop_condition()
        else:b.cx(measurement,ghost[j-net.n]);pending.add(j)
    for kind,j in actions:
        if kind=='compute':
            q=b.alloc_qubits(1)[0];a,c=net.nodes[j-net.n]
            linear(a,A);linear(c,B);b.init_and(A,B,q);linear(c,B);linear(a,A)
            if j in pending:
                b.push_condition(ghost[j-net.n]);b.z(q);b.pop_condition();b.bit_store(ghost[j-net.n],False);pending.remove(j)
            live[j]=q
        elif kind=='ghost':
            q=live.pop(j);b.hmr(q,measurement);b.cx(measurement,ghost[j-net.n]);b.bit_store(measurement,False);pending.add(j);b.free(km.array([q]))
        elif kind=='clear':
            q=live.pop(j);a,c=net.nodes[j-net.n]
            linear(a,A);linear(c,B);b.del_and(A,B,q);linear(c,B);linear(a,A);b.free(km.array([q]))
        elif kind=='phase_start':
            if net.outputs[j].const:b.x(parity)
        elif kind=='phase_xor':b.cx(address[j] if j<net.n else live[j],parity)
        elif kind=='phase_finish':
            b.z_pow(parity,Fraction(2)**(j+1-len(net.outputs)));b.hmr(parity,measurement)
            signal=net.outputs[j]
            if signal.const:
                b.push_condition(measurement);b.z(True);b.pop_condition()
            for v in sorted(signal.terms):phase_bit(measurement,v)
            b.bit_store(measurement,False)
        else:raise ValueError(kind)
    assert not live and not pending
    b.bit_store(ghost,False);b.free(ghost);b.free(km.array([measurement]));b.free(km.array([A,B,parity]));return stats


def main():
    import json,numpy as np
    from pathlib import Path
    from boolean_polynomial_phase import build,evaluate
    from invsqrt.verify_submission_invsqrt import half_turns_to_units
    from boolean_lookup import Network,Signal
    rng=np.random.default_rng(62059);reports=[]
    # A deep random network deliberately exceeds the cache, exercising ghosts.
    net=Network(8);signals=[Signal(frozenset([i])) for i in range(8)]
    for i in range(80):
        a=signals[int(rng.integers(len(signals)))];c=signals[int(rng.integers(len(signals)))];v=net.mul(a,c)
        signals.append(v^signals[int(rng.integers(len(signals)))])
    net.outputs=signals[-8:]
    for budget in [8,16,32,64]:
        try:actions,r=plan(net,budget)
        except CapacityError as e:reports.append(dict(budget=budget,status='CAPACITY',reason=str(e)));continue
        b=km.CircuitBuilder();inputs=b.create_quantum_register(net.n,'x');emit(net,b,inputs,budget,actions);c=b.finish_circuit();shots=256
        values=rng.integers(0,2,(net.n,shots),dtype=bool);want=evaluate(net,values)
        sim=km.Simulator(batch_size=shots,seed=budget);sim.use_same_registers_as(c);sim.write_across_shots(inputs,values);sim.do(c)
        assert np.array_equal(sim.read_across_shots(inputs),values)
        scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in set(inputs)]);assert not sim.read_across_shots(scratch).any()
        got=np.array([half_turns_to_units(sim.read_shot_phase(i)) for i in range(shots)],dtype=np.uint64)
        expected=want.astype(np.uint64)<<np.uint64(64-len(net.outputs));assert np.array_equal(got,expected),(budget,np.where(got!=expected)[0][:5])
        reports.append(dict(status='PASS',budget=budget,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,**r))
    Path('results/spooky_phase_check.json').write_text(json.dumps(reports,indent=2)+'\n');print(json.dumps(reports),flush=True)
    net,_,_=build();profiles=[]
    for budget in [16,24,32,40,48,64,80,96,128,192,256,512,800]:
        try:_,r=plan(net,budget,limit=20000);profiles.append(dict(status='PLANNED',budget=budget,Q_model=net.n+3+r['peak'],**r))
        except CapacityError as e:profiles.append(dict(status='CAPACITY',budget=budget,reason=str(e)))
        print(json.dumps(profiles[-1]),flush=True)
    Path('results/spooky_phase_polynomial_profiles.json').write_text(json.dumps(profiles,indent=2)+'\n')
if __name__=='__main__':main()
