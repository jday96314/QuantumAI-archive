"""Spooky AND scheduling with streamed linear forms.

Linear operands are accumulated one term at a time, then measured after use.
Their phase repairs distribute to live source bits or classical ghost masks.
This avoids pinning every AND in a wide XOR expression simultaneously.
"""
from collections import Counter
from functools import lru_cache
from fractions import Fraction
import kickmix as km
from boolean_lookup import Signal
from spooky_phase import CapacityError


def plan(net,budget,limit=20000,keep_actions=True):
    live=set();ghost=set();pins=Counter();age={};tick=0;computes=0;peak=0;temporaries=0;temp_id=net.n+len(net.nodes);actions=[];events=Counter()
    def action(*args):
        events[args[0]]+=1
        if keep_actions:actions.append(args)
    def resource():
        nonlocal peak
        peak=max(peak,len(live)+temporaries)
    @lru_cache(None)
    def height(j):
        if j<net.n:return 0
        return 1+max((height(p) for s in net.nodes[j-net.n] for p in s.terms),default=0)
    def touch(j):
        nonlocal tick
        tick+=1;age[j]=tick
    def repay_signal(s):
        ghost.update(j for j in s.terms if j>=net.n and j not in live)
    def room():
        if len(live)<budget:return
        choices=[j for j in live if not pins[j]]
        if not choices:raise CapacityError(f'all {budget} AND nodes pinned')
        j=min(choices,key=lambda j:age[j]);live.remove(j);ghost.add(j);action('ghost',j)
    def form(s):
        nonlocal temporaries,temp_id
        if len(s.terms)==1 and not s.const:
            j=next(iter(s.terms));ensure(j)
            if j>=net.n:pins[j]+=1
            return j,False,s
        terms=sorted(s.terms,key=lambda j:(-height(j),-j))
        if terms:ensure(terms[0])
        t=temp_id;temp_id+=1;temporaries+=1;resource();action('lin_start',t,s.const)
        for j in terms:ensure(j);action('lin_xor',t,j)
        return t,True,s
    def drop(ref):
        nonlocal temporaries
        j,temporary,s=ref
        if temporary:
            action('lin_drop',j,s);repay_signal(s);temporaries-=1
        elif j>=net.n:pins[j]-=1
    def operands(j):
        a,c=net.nodes[j-net.n]
        if max((height(p) for p in a.terms),default=0)<max((height(p) for p in c.terms),default=0):a,c=c,a
        return form(a),form(c)
    def ensure(j):
        nonlocal computes
        if j<net.n:return
        if j in live:touch(j);return
        a,c=operands(j);room();action('compute',j,a[0],c[0]);live.add(j);ghost.discard(j);touch(j);computes+=1;resource()
        drop(c);drop(a)
        if computes>limit:raise CapacityError(f'computation limit {limit}')
    for i,s in enumerate(net.outputs):
        ref=form(s);action('phase',ref[0],i);drop(ref)
    while live or ghost:
        j=max(live|ghost);ensure(j);pins[j]+=1;a,c=operands(j)
        action('clear',j,a[0],c[0]);live.remove(j);ghost.discard(j);pins[j]-=1
        drop(c);drop(a)
    assert not temporaries
    return actions,dict(computes=computes,workspace_peak=peak,ghost_events=events['ghost'],linear_measurements=events['lin_drop'],actions=sum(events.values()))


def emit(net,b,address,budget,actions=None):
    if actions is None:actions,stats=plan(net,budget)
    else:stats=None
    live={};temp={};pending=set();ghost=b.alloc_bits(len(net.nodes));measurement=b.alloc_bits(1)[0]
    def qubit(j):return address[j] if j<net.n else (live[j] if j<net.n+len(net.nodes) else temp[j])
    def repay_signal(s):
        if s.const:
            b.push_condition(measurement);b.z(True);b.pop_condition()
        for j in s.terms:
            if j<net.n or j in live:
                b.push_condition(measurement);b.z(qubit(j));b.pop_condition()
            else:b.cx(measurement,ghost[j-net.n]);pending.add(j)
        b.bit_store(measurement,False)
    for op in actions:
        kind,j,*args=op
        if kind=='compute':
            a,c=args;q=b.alloc_qubits(1)[0];b.init_and(qubit(a),qubit(c),q)
            if j in pending:
                b.push_condition(ghost[j-net.n]);b.z(q);b.pop_condition();b.bit_store(ghost[j-net.n],False);pending.remove(j)
            live[j]=q
        elif kind=='ghost':
            q=live.pop(j);b.hmr(q,measurement);b.cx(measurement,ghost[j-net.n]);b.bit_store(measurement,False);pending.add(j);b.free(km.array([q]))
        elif kind=='clear':
            a,c=args;q=live.pop(j);b.del_and(qubit(a),qubit(c),q);b.free(km.array([q]))
        elif kind=='lin_start':
            temp[j]=b.alloc_qubits(1)[0]
            if args[0]:b.x(temp[j])
        elif kind=='lin_xor':b.cx(qubit(args[0]),temp[j])
        elif kind=='lin_drop':
            q=temp.pop(j);b.hmr(q,measurement);repay_signal(args[0]);b.free(km.array([q]))
        elif kind=='phase':b.z_pow(qubit(j),Fraction(2)**(args[0]+1-len(net.outputs)))
        else:raise ValueError(kind)
    assert not pending and not live and not temp
    b.bit_store(ghost,False);b.free(ghost);b.free(km.array([measurement]));return stats


def main():
    import json,numpy as np
    from pathlib import Path
    from boolean_lookup import Network
    from boolean_polynomial_phase import build,evaluate
    from invsqrt.verify_submission_invsqrt import half_turns_to_units
    rng=np.random.default_rng(72051);net=Network(8);sig=[Signal(frozenset([i])) for i in range(8)]
    for _ in range(80):
        a=sig[int(rng.integers(len(sig)))];c=sig[int(rng.integers(len(sig)))];sig.append(net.mul(a,c)^sig[int(rng.integers(len(sig)))])
    net.outputs=sig[-8:];checks=[]
    for budget in [2,4,8,16,32]:
        try:actions,r=plan(net,budget,limit=20000)
        except CapacityError as e:checks.append(dict(status='CAPACITY',budget=budget,reason=str(e)));continue
        b=km.CircuitBuilder();address=b.create_quantum_register(net.n,'x');emit(net,b,address,budget,actions);c=b.finish_circuit()
        val=rng.integers(0,2,(net.n,256),dtype=bool);expected=evaluate(net,val);s=km.Simulator(batch_size=256,seed=budget);s.use_same_registers_as(c);s.write_across_shots(address,val);s.do(c)
        assert np.array_equal(s.read_across_shots(address),val)
        scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in set(address)]);assert not s.read_across_shots(scratch).any()
        phases=np.array([half_turns_to_units(s.read_shot_phase(i)) for i in range(256)],dtype=np.uint64)
        assert np.array_equal(phases,expected.astype(np.uint64)<<np.uint64(64-len(net.outputs)))
        checks.append(dict(status='PASS',budget=budget,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,**r));print(json.dumps(checks[-1]),flush=True)
    Path('results/spooky_linear_phase_check.json').write_text(json.dumps(checks,indent=2)+'\n')
    net,_,_=build();profiles=[]
    for budget in [4,8,12,16,24,32,48,64,96,128,256,512,800]:
        try:_,r=plan(net,budget,limit=20000,keep_actions=False);profiles.append(dict(status='PLANNED',budget=budget,Q_model=net.n+r['workspace_peak'],**r))
        except CapacityError as e:profiles.append(dict(status='CAPACITY',budget=budget,reason=str(e)))
        print(json.dumps(profiles[-1]),flush=True)
    Path('results/spooky_linear_phase_polynomial_profiles.json').write_text(json.dumps(profiles,indent=2)+'\n')
if __name__=='__main__':main()
