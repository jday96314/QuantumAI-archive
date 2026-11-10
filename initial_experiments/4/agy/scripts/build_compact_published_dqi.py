"""Capture and build the 442k Toffoli Psi-route DQI circuit with 1429 qubits."""

import json
import subprocess
import time
from pathlib import Path
import numpy as np
import sys

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
AGY = ROOT / "initial_experiments/4/agy"
sys.path[:0] = [
    str(ROOT / "initial_experiments/3"),
    str(ROOT / "initial_experiments/2"),
    str(ROOT / "initial_experiments/3/vendor/uncompute"),
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
]

import kickmix as km
import published_dqi as p
from tower_inverse import map_ops, apply_map
import importlib.util

def capture():
    print("Loading discovered circuit...")
    SOURCE = ROOT / "initial_experiments/3/vendor/uncompute/tasks/dqi-reed-solomon/circuits/discovered/circuit.py"
    spec = importlib.util.spec_from_file_location("discovered_circuit", SOURCE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    OldBank = mod.BankReg
    bank_offset = 3

    class SyndBank(OldBank):
        __slots__ = ('slot',)
        def __init__(self, A, name, width=mod.B):
            super().__init__(A, name, width)
            self.slot = None
            if name == 'bb':
                self.slot = A.bb_count
                A.bb_count += 1
        def mat(self):
            if self.wires is None and self.slot is not None and self.slot >= bank_offset:
                A = self.A
                self.wires = list(A.syn[self.slot - bank_offset])
                A.borrowed.update(self.wires)
                self.round = A.round
                A.bank_log.setdefault(A.round, []).append(self)
                return self.wires
            return super().mat()

    mod.BankReg = SyndBank

    OldAlloc = mod.Alloc
    class SyndAlloc(OldAlloc):
        def __init__(self, b):
            super().__init__(b)
            self.bb_count = 0
            self.borrowed = set()
            self.syn = None

    mod.Alloc = SyndAlloc

    OldBM = mod.BM
    class SyndBM(OldBM):
        def __init__(self, E, A, syn, Dw, dW):
            A.syn = syn
            super().__init__(E, A, syn, Dw, dW)

    mod.BM = SyndBM

    mod.ZMASKI = [0] * 70
    hook = {}

    def sweep(b, E, A, omw, ebw, abw):
        hook.update(
            index=len(b._items),
            omega=[list(x) for x in omw],
            eo=[list(x) for x in ebw],
            ab=[list(x) for x in abw],
            free=list(b.free_qubits),
        )

    mod.emit_sweep = sweep
    mod.assert_metrics = lambda b, p: b._stream_totals()[0]

    print("Building circuit in Builder...")
    t0 = time.time()
    b, t = mod.build()
    print(f"Built tape in {time.time()-t0:.2f}s: {len(b._items)} ops, {t:,} Toffolis, {b.peak_qubits} peak qubits, {b.next_qubit} next qubit.")

    rows = np.asarray(b._items._rows, dtype=np.uint32).reshape(-1, 9)
    tape_path = AGY / "runs/published_compact_tape.npy"
    np.save(tape_path, rows)
    hook.update(
        peak=b.peak_qubits,
        toffoli=t,
        ops=len(rows),
        maxq=b.next_qubit,
        outcomes=b.next_out,
        bits=b.next_bit,
    )
    hook_path = AGY / "runs/published_compact_hook.json"
    hook_path.write_text(json.dumps(hook, indent=2))
    print(f"Saved tape ({tape_path.stat().st_size:,} bytes) and hook to runs/")
    return hook, tape_path

def build_kmb(m=4095, output_kmb=AGY / "runs/dqi_compact_psi.kmb"):
    pub = p.load()
    forward, back, transform = p.isomorphism(pub)
    p.ISOM_R = pub.LOG[forward[1]]
    hook = json.loads((AGY / "runs/published_compact_hook.json").read_text())
    
    print("Building wrapper circuit for runtime sweep...")
    b = km.CircuitBuilder()
    syn = [b.create_quantum_register(12, f'S{i}') for i in range(70)]
    masks = [b.create_classical_register(12, f'E{i}') for i in range(m)]
    q = [v for w in syn for v in w] + list(b.alloc_qubits(hook['maxq'] - 840))
    ops = map_ops(forward)
    for w in syn:
        apply_map(b, w, ops)
    b.info_tag('SPLIT_SWEEP_START')
    p.sweep(b, pub, hook, masks, q, back, transform, range(m))
    b.info_tag('SPLIT_SWEEP_END')
    for w in syn:
        apply_map(b, w, ops, True)
    c = b.finish_circuit()
    
    wrapper = output_kmb.with_suffix('.wrapper.kmx')
    c.to_file(str(wrapper), format='kmx')
    paths = [output_kmb.with_suffix('.' + x + '.kmx') for x in ['pre', 'sweep', 'post']]
    handles = [x.open('w') for x in paths]
    part = 0
    with wrapper.open() as f:
        for line in f:
            if 'SPLIT_SWEEP_START' in line:
                part = 1
            elif 'SPLIT_SWEEP_END' in line:
                part = 2
            else:
                handles[part].write(line)
    for f in handles:
        f.close()
        
    full = output_kmb.with_suffix('.kmx')
    print("Running translate_published...")
    tape_path = AGY / "runs/published_compact_tape.npy"
    subprocess.run([
        str(ROOT / "initial_experiments/3/translate_published"),
        str(tape_path),
        str(hook['index']),
        str(c.num_bits),
        *map(str, paths),
        str(full),
        str(hook['outcomes'])
    ], check=True)
    
    print(f"Loading full circuit from {full}...")
    c_full = km.Circuit.from_file(str(full))
    t_full = c_full.operation_statistics().num_toffoli_like
    print(f"Assembled full circuit: {c_full.num_qubits} qubits, {c_full.num_bits} bits, {t_full:,} Toffolis.")
    c_full.to_file(str(output_kmb), format='kmb')
    print(f"Wrote {output_kmb} ({output_kmb.stat().st_size:,} bytes)!")
    return output_kmb

def main():
    AGY.joinpath("runs").mkdir(parents=True, exist_ok=True)
    hook, tape_path = capture()
    out_kmb = build_kmb()
    print("Done building compact published DQI circuit!")

if __name__ == "__main__":
    main()
