"""Replay exact selected scan functions, tracking native allocations.

Two field points suffice for scan-storage liveness: every point has the same
allocation/free sequence; only constant field maps and classical masks vary.
Full 4095-point native gate totals are measured separately in native_breakdown.
"""
from collections import Counter
from contextlib import contextmanager
import inspect
import json
from pathlib import Path
from types import SimpleNamespace
import time

import kickmix as km
from import_published import load
from published_dqi import isomorphism
from fast_dqi import FastBuilder
from profile_selected import HERE,R3,R4,TAG


class TrackingBuilder(km.CircuitBuilder):
    def __init__(self):
        super().__init__();self.live=set();self.stage='setup';self.peak=0;self.peaks={};self.events=[];self.roles={}
    def take(self,qs):
        ids={q.id for q in qs}
        assert not ids & self.live,ids & self.live
        self.live.update(ids);self.snapshot()
    def snapshot(self):
        self.peak=max(self.peak,len(self.live))
        old=self.peaks.get(self.stage)
        if old is None or len(self.live)>=old['live']:
            self.peaks[self.stage]=dict(live=len(self.live),by_role=dict(Counter(self.roles.get(q,'scan_scratch') for q in self.live)))
    def create_quantum_register(self,*args,**kwargs):
        qs=super().create_quantum_register(*args,**kwargs);self.take(qs);return qs
    def alloc_qubits(self,*args,**kwargs):
        qs=super().alloc_qubits(*args,**kwargs);self.take(qs);return qs
    def free(self,qs):
        items=[qs] if isinstance(qs,(km.q,km.b)) else qs
        ids={q.id for q in items if isinstance(q,km.q)}
        assert ids<=self.live,ids-self.live
        self.live.difference_update(ids)
        return super().free(qs)


def round3():
    import stream_forney as scan
    hook=json.loads((R3/f'parameters/{TAG}_hook.json').read_text())
    b=TrackingBuilder();b.roles={q:'syndrome_input_or_borrowed_template' for q in range(840)}
    for q in set(range(840,hook['maxq']))-set(hook['free']):b.roles[q]='persistent_template_history_controls'
    for word in hook['sigma']:
        for q in word:b.roles[q]='locator'
    words=[b.create_quantum_register(12,f'S{i}') for i in range(70)]
    masks=[b.create_classical_register(12,f'E{i}') for i in range(2)]
    q=[v for w in words for v in w]+list(b.alloc_qubits(hook['maxq']-840))
    m=load();forward,back,_=isomorphism(m)
    b.free(km.array([q[i] for i in hook['free']]))
    b.stage='scan_normalization';b.snapshot()
    originals={key:getattr(scan,key) for key in ['omega','word_mul','chien','prefix_zero']}
    depth=[0]
    def wrap(fn,key):
        def call(*args,**kwargs):
            old=b.stage
            if depth[0]==0:b.stage='scan/'+key;b.snapshot()
            depth[0]+=1
            try:return fn(*args,**kwargs)
            finally:depth[0]-=1;b.stage=old
        return call
    for key,fn in originals.items():setattr(scan,key,wrap(fn,key))
    try:scan.phase_kernel(b,hook,masks,q,forward,back)
    finally:
        for key,fn in originals.items():setattr(scan,key,fn)
    circuit=b.finish_circuit()
    return dict(family='round3',points=2,peak_live=b.peak,peak_scan_live=max(v['live'] for k,v in b.peaks.items() if k!='setup'),peak_native_wire_frontier=circuit.num_qubits,
                stage_peaks=b.peaks,end_live=len(b.live),scope='exact selected scan allocation sequence with two points; raw ids before lifetime compaction')


def round4():
    import compact_product_decoder as decoder
    from blocked_euclid_sequence import frames
    from cubic_inverse_product_kernel import apply as cubic
    flags=[]
    with (R4/'candidates/cubic_trace_product_n70_m4095.scan.kmx').open() as f:
        for line in f:
            if line.startswith('CSWAP'):break
            assert line.startswith('CX '),line
            source=int(line.split()[1][1:])
            if source not in flags:flags.append(source)
    assert len(flags)==35,flags
    _,(bp,bb,_,p)=frames([[12*i+j for j in range(12)] for i in range(70)],3)
    instances=[]
    class Tracker(TrackingBuilder):
        def __init__(self):
            super().__init__();instances.append(self)
            self.roles={q:'cofactor_previous' for word in bp for q in word}
            self.roles.update({q:'cofactor_current' for word in bb for q in word})
            self.roles.update({q:'scale' for q in p})
            self.roles.update({q:'completion_flag' for q in flags})
        def info_tag(self,name):
            self.stage=name;self.snapshot();return super().info_tag(name)
        def free(self,qs):
            result=super().free(qs)
            if self.stage=='COMPACT_SCAN_START':
                self.stage='scan_work';self.snapshot()
            return result
    proxy=SimpleNamespace(**{name:getattr(km,name) for name in dir(km) if not name.startswith('_')})
    proxy.CircuitBuilder=Tracker
    # Execute the exact selected source fragment without rebuilding the core.
    # No algebra, gate, allocation, or instruction is changed in the fragment.
    source=inspect.getsource(decoder.build)
    fragment=source[source.index('    forward,back,transform=isomorphism(m)'):source.index('    wrapper=output.with_suffix')]
    fragment='\n'.join(line[4:] if line.startswith('    ') else line for line in fragment.splitlines())
    env=dict(decoder.__dict__)
    env.update(m=load(),km=proxy,b=SimpleNamespace(next_qubit=961),N=70,M=2,flags=flags,bp=bp,bb=bb,p=p,
               begin=time.time(),virtual_lead=True,phase_kernel=cubic)
    exec(compile(fragment,str(R4/'compact_product_decoder.py')+' [unchanged scan fragment]','exec'),env)
    b=instances[0];circuit=env['kb'].finish_circuit()
    return dict(family='round4',points=2,peak_live=b.peak,peak_scan_live=max(v['live'] for k,v in b.peaks.items() if k not in ['setup','COMPACT_SCAN_START']),peak_native_wire_frontier=circuit.num_qubits,
                stage_peaks=b.peaks,end_live=len(b.live),flags=flags,
                scope='unchanged selected scan source fragment with two points; setup reserves core frontier then releases clean core scratch')


if __name__=='__main__':
    result=[round3(),round4()]
    (HERE/'results/scan_liveness.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
