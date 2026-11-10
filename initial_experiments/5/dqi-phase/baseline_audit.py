"""Attribute selected round-4 bytes without rebuilding its locator recurrence."""
import hashlib
import json
from pathlib import Path
import numpy as np

OLD=Path(__file__).resolve().parent/'../../4/codex'


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        while data:=f.read(16<<20):h.update(data)
    return h.hexdigest()


def main():
    name='cubic_trace_product_n70_m4095'
    report=json.loads((OLD/'results'/f'{name}_verify.json').read_text())
    digest=sha(OLD/'candidates'/f'{name}.kmb');assert digest==report['circuit_sha256']
    rows=np.load(OLD/'candidates'/f'{name}.npy',mmap_mode='r')
    coreT=0;counts={}
    for i in range(0,len(rows),1<<20):
        kinds=rows[i:i+(1<<20),0]&255
        for k,n in zip(*np.unique(kinds,return_counts=True)):counts[int(k)]=counts.get(int(k),0)+int(n)
        coreT+=int(np.count_nonzero((kinds==6)|(kinds==7)))
    scan=88*4095
    # Toffoli-like CSWAP counts once. Read only the wrapper's scan section,
    # which includes both cofactor unpadding rotations and the point bodies.
    scan_counts={}
    with (OLD/'candidates'/f'{name}.scan.kmx').open('rb') as f:
        for line in f:
            kind=line.partition(b' ')[0].strip()
            if kind in [b'CCX',b'CCZ',b'CSWAP']:
                scan_counts[kind.decode()]=scan_counts.get(kind.decode(),0)+1
    wrapperT=sum(scan_counts.values());static=coreT+wrapperT
    assert static==report['max_toffoli_count'],(static,report['max_toffoli_count'])
    out=dict(status='PASS',control=str(OLD/'candidates'/f'{name}.kmb'),circuit_sha256=digest,
        actual_core_tape_operations=len(rows),core_kind_counts=counts,core_forward_and_reverse_T=coreT,
        scan_wrapper_kind_counts=scan_counts,scan_wrapper_T=wrapperT,point_body_T=88,positions=4095,
        all_point_bodies_static_T=scan,cofactor_unpadding_and_repadding_T=wrapperT-scan,
        static_total_T=static,native_Q=report['num_qubits'],native_mean_T=report['toffoli_count'],
        native_R=report['repetition_factor'],native_epsilon_plus=report['effective_error_rate'],native_cost=report['score'],
        native_shots=report['verified_shots'],native_value_errors=report['value_errors'],native_phase_errors=report['phase_errors'],
        scope='read-only attribution of selected round-4 candidate; existing native evidence, no new whole-decoder validation')
    Path('results/baseline_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)


if __name__=='__main__':main()
