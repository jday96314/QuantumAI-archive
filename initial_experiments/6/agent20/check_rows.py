"""Native checks against the unchanged independent integer phase model."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import kickmix as km
import numpy as np
from invsqrt.verify_submission_invsqrt import (
    ProblemSpecInvSqrt, validate_registers, to_bits, from_bits, half_turns_to_units,
)

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('reference_model', HERE.parent.parent/'5/inv-phase/model.py')
reference = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reference)


def check(path):
    circuit = km.Circuit.from_file(str(path))
    x, scratch = validate_registers(circuit, ProblemSpecInvSqrt.for_n_p(10))
    values = [1,2,3,4,958102,2965859,3133452]
    for k in range(22):
        for i in range(257):
            v = ((256+i)*(1<<k))//256
            values.extend(w for w in (v-1,v,v+1) if 1<=w<=3133452)
    rng = np.random.default_rng(20261025)
    values = np.unique(np.r_[values,rng.integers(1,3133453,32768)]).astype(np.uint64)
    xs = np.resize(values, ((len(values)+255)//256)*256)
    total_tof = 0
    for salt in ('agent20-karatsuba-rows-component-a', 'agent20-karatsuba-rows-component-b'):
        sim = km.Simulator(batch_size=256,counting_operations=True,
                           seed=int(hashlib.sha256(salt.encode()).hexdigest(),16))
        sim.use_same_registers_as(circuit)
        for offset in range(0,len(xs),256):
            xx = xs[offset:offset+256]
            sim.clear_for_shot()
            sim.write_across_shots(x,to_bits(xx,22))
            sim.do(circuit)
            assert np.array_equal(from_bits(sim.read_across_shots(x)),xx), 'input restoration'
            assert not sim.read_across_shots(scratch).any(), 'dirty scratch'
            actual = np.array([half_turns_to_units(sim.read_shot_phase(j)) for j in range(256)],dtype=np.uint64)
            expected = reference.evaluate(xx.astype(np.int64))
            assert np.array_equal(actual,expected), ('phase',salt,int(xx[np.flatnonzero(actual != expected)[0]]))
        total_tof += sim.operation_statistics().num_toffoli_like
    report = dict(status='EXACT_NATIVE_MODEL_MATCH', scope='development diagnostic; full verifier still required',
                  path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                  Q=circuit.num_qubits, static_Toffolis=circuit.operation_statistics().num_toffoli_like,
                  distinct_inputs=len(values),shots=2*len(xs),seeds=2,
                  mean_Toffolis=total_tof/(2*len(xs)))
    print(json.dumps(report),flush=True)
    (HERE/'results/rows_exact_native.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('path',type=Path,nargs='?',default=HERE/'candidates/rows_best.kmb')
    check(p.parse_args().path)
