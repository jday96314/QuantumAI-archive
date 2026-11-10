"""Replay the EXACT selected builders, changing bookkeeping only.

Round-3 capture and round-4 core tapes are compared byte for byte to the
selected tapes. No inherited file is written. Scan statistics are read from
the exact selected split native KMX, with outcome conditions distinguished
from classical input conditions by the native parser.
"""
from collections import Counter
from contextlib import contextmanager
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import kickmix as km
from import_published import load

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]
R3 = PROJECT / 'initial_experiments/3'
R4 = PROJECT / 'initial_experiments/4/codex'
TAG = 'stream10_offset3_select_noscale_compact'


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        while data := f.read(8 << 20):
            h.update(data)
    return h.hexdigest()


def summarize(b, family, regions, state=None):
    """Actual allocator events, with no assumed cleanliness of live data.

    Ownership reports physical storage. Public syndrome wires remain pinned
    even after their logical role changes. Claimed scratch and bank sets
    partition private storage; region end states add the logical layout.
    """
    live = set()
    owner = {}
    records = []
    at = 0
    for name, begin, end, layout, event_begin, event_end in regions:
        while at < event_begin:
            event, data, _, _ = b.wire_log[at]
            if event == 'alloc': live.add(data)
            elif event == 'free': live.remove(data); owner.pop(data, None)
            elif event == 'claim':
                label, ids = data
                if not label.startswith('~'): owner.update({q: label for q in ids})
            elif event == 'release':
                label, ids = data
                for q in ids:
                    if owner.get(q) == label: owner.pop(q)
            at += 1
        peak = None
        def snap():
            counts = Counter()
            for q in live:
                if q < 840: key = 'input_array'
                else: key = owner.get(q, 'fixed_or_anonymous')
                if key not in ['sig', 'bb', 'stk', 'input_array', 'fixed_or_anonymous']:
                    key = 'scratch/' + key
                counts[key] += 1
            return dict(live=len(live), physical_by_role=dict(counts))
        peak = snap()
        while at < event_end:
            event, data, _, _ = b.wire_log[at]
            if event == 'alloc': live.add(data)
            elif event == 'free': live.remove(data); owner.pop(data, None)
            elif event == 'claim':
                label, ids = data
                if not label.startswith('~'): owner.update({q: label for q in ids})
            elif event == 'release':
                label, ids = data
                for q in ids:
                    if owner.get(q) == label: owner.pop(q)
            if len(live) >= peak['live']: peak = snap()
            at += 1
        records.append(dict(region=name, begin=begin, end=end, peak=peak,
                            end_state=snap(), logical_layout=layout))
    return dict(family=family, peak_allocated=b.peak_qubits,
                maximum_wire_frontier=b.next_qubit,
                phase_statistics=b.phase_summary(), timeline=records)


def round3():
    import stream_capture as capture
    m = load()
    builders = []
    old_builder = m.Builder
    regions = []
    class Builder(old_builder):
        def __init__(self):
            super().__init__(); builders.append(self)
        @contextmanager
        def phase(self, name):
            start = len(self._items)
            event_start = len(self.wire_log)
            with super().phase(name): yield
            regions.append((name, start, len(self._items), {},event_start,len(self.wire_log)))
    m.Builder = Builder
    capture.load = lambda: m
    capture.HERE = HERE
    save = np.save
    comparison = {}
    def compare(path, rows):
        source = R3 / f'parameters/{TAG}_tape.npy'
        old = np.load(source, mmap_mode='r')
        assert rows.shape == old.shape
        for i in range(0, len(rows), 100000):
            assert np.array_equal(rows[i:i+100000], old[i:i+100000]), i
        comparison.update(equal=True, rows=len(rows), source=str(source),
                          source_sha256=digest(source))
    np.save = compare
    try: capture.capture(TAG)
    finally: np.save = save
    b = builders[0]
    report = summarize(b, 'round3_selected_BM_core', regions)
    report['tape_comparison'] = comparison
    hook = json.loads((HERE / f'parameters/{TAG}_hook.json').read_text())
    # Logical consumed-input/template aliasing is visible directly in the
    # selected Bank.mat implementation. Input storage stays physically pinned.
    for rec in report['timeline']:
        if '/r' in rec['region']:
            r = int(rec['region'].split('/r')[1])
            size = min(37, max(m.KMSV[r], m.KMMV[r]) + 2)
            rec['logical_layout'] = dict(
                locator_max_coefficients=m.IMAXV[r]+1,
                template_active_window=size,
                template_words_aliased_to_input=max(0, size-3),
                interpretation='window bounds, not a claim that all coefficients are nonzero')
    report['core_forward_T'] = sum(p['toffoli'] for p in report['phase_statistics']['phases'] if p['name'].startswith('bm/'))
    report['core_reverse_T'] = sum(p['toffoli'] for p in report['phase_statistics']['phases'] if p['name'].startswith('reverse/'))
    report['hook'] = {k: hook[k] for k in ['index', 'peak', 'maxq', 'toffoli', 'ops']}
    (HERE / 'results/round3_liveness.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['timeline','phase_statistics']}), flush=True)


class CoreDone(Exception): pass


def round4():
    import compact_product_decoder as decoder
    m = load()
    builders = []
    old_builder = m.Builder
    class Builder(old_builder):
        def __init__(self): super().__init__(); builders.append(self)
    m.Builder = Builder
    decoder.load = lambda: m
    regions = []
    sequence = Counter()
    originals = {key:getattr(decoder,key) for key in ['first_step','block','serial']}
    def wrap(fn, kind):
        def call(*args, **kwargs):
            b = args[1].b
            inverse = (args[4] if kind=='first' and len(args)>4 else
                       args[9] if kind=='block' and len(args)>9 else
                       args[8] if kind=='serial' and len(args)>8 else
                       kwargs.get('inverse',False))
            # block positional signature: m,E,math,Bp,B,C,pinv,block,inverse
            if kind=='block': inverse = args[8] if len(args)>8 else kwargs.get('inverse',False)
            direction='reverse' if inverse else 'forward'
            idx=sequence[(kind,direction)];sequence[(kind,direction)]+=1
            name=f'{direction}/{kind}{idx:02d}'
            start=len(b._items)
            event_start=len(b.wire_log)
            if kind=='first': layout=dict(cofactor_previous_words=0,cofactor_current_words=1,remainder_words=68,scale_words=1)
            else: layout=dict(cofactor_previous_words=len(args[3]),cofactor_current_words=len(args[4]),remainder_words=len(args[5]),scale_words=1)
            with b.phase(name): result=fn(*args,**kwargs)
            regions.append((name,start,len(b._items),layout,event_start,len(b.wire_log)))
            return result
        return call
    decoder.first_step=wrap(originals['first_step'],'first')
    decoder.block=wrap(originals['block'],'block')
    decoder.serial=wrap(originals['serial'],'serial')
    save=np.save
    comparison={}
    def compare(path, rows):
        import inspect
        local=inspect.currentframe().f_back.f_locals
        hook=local['hook'];b=local['b']
        source=R4/'candidates/cubic_trace_product_n70_m4095.npy'
        old=np.load(source,mmap_mode='r')
        assert rows.shape==old.shape, (rows.shape,old.shape)
        for i in range(0,len(rows),100000):
            assert np.array_equal(rows[i:i+100000],old[i:i+100000]),i
        comparison.update(equal=True,rows=len(rows),hook=hook,source=str(source),source_sha256=digest(source))
        regions.sort(key=lambda x:x[1])
        # Fill exact gaps (arming, final normalization, and flag cleanup).
        full=[];pos=0;event_pos=0
        for region in regions:
            if region[1]>pos or region[4]>event_pos:full.append(('control_or_scale',pos,region[1],{},event_pos,region[4]))
            full.append(region);pos=region[2];event_pos=region[5]
        if pos<len(rows) or event_pos<len(b.wire_log):full.append(('control_or_scale',pos,len(rows),{},event_pos,len(b.wire_log)))
        report=summarize(b,'round4_selected_compact_core',full)
        report['tape_comparison']=comparison
        def count(a):return int(np.count_nonzero(np.isin(a[:,0]&255,[6,7])))
        report['core_forward_T']=count(old[:hook])
        report['core_reverse_T']=count(old[hook:])
        report['hook_layout']={key:len(local[key]) for key in ['bp','bb','flags']}
        (HERE/'results/round4_liveness.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps({k:v for k,v in report.items() if k not in ['timeline','phase_statistics']}),flush=True)
        raise CoreDone()
    np.save=compare
    try:
        decoder.build(70,3,22,4095,HERE/'candidates/profile_unused.kmb',linear_padding=True,
                      virtual_zero=True,virtual_lead=True,low_inverse=True,lazy_flags=True)
    except CoreDone: pass
    finally:
        np.save=save
        for key,fn in originals.items():setattr(decoder,key,fn)


def native_scan():
    results=[]
    for family,path,result in [
        ('round3',R3/f'candidates/{TAG}_full.sweep.kmx',R3/f'results/{TAG}_verify.json'),
        ('round4',R4/'candidates/cubic_trace_product_n70_m4095.scan.kmx',R4/'results/cubic_trace_product_n70_m4095_verify.json')]:
        # Native loading accounts for the dynamic classical masks. These split
        # files use original physical ids, not the final compacted width.
        circuit=km.Circuit.from_file(str(path),format='kmx')
        static=dict(circuit.operation_statistics())
        full=json.loads(result.read_text())
        core=json.loads((HERE/f'results/{family}_liveness.json').read_text())
        core_t=core['core_forward_T']+core['core_reverse_T']
        results.append(dict(family=family,scan_source=str(path),scan_sha256=digest(path),
            scan_native_statistics=static,scan_static_T=circuit.operation_statistics().num_toffoli_like,
            forward_core_T=core['core_forward_T'],reverse_core_T=core['core_reverse_T'],
            sum_static_T=core_t+circuit.operation_statistics().num_toffoli_like,
            selected_max_T=full['max_toffoli_count'],selected_mean_T=full['toffoli_count'],
            selected_Q=full['num_qubits'],selected_R=full['repetition_factor'],selected_cost=full['score'],
            core_plus_scan_matches_selected=(core_t+circuit.operation_statistics().num_toffoli_like==full['max_toffoli_count'])))
        print(json.dumps(results[-1]),flush=True)
    (HERE/'results/native_breakdown.json').write_text(json.dumps(results,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['round3','round4','scan']);a=p.parse_args()
    {'round3':round3,'round4':round4,'scan':native_scan}[a.mode]()
