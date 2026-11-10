"""Test offset=2, 3, 4, 5 in stream_model across 10,000 samples."""
import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3")
sys.path[:0] = [
    str(ROOT),
    str(ROOT / "vendor/uncompute"),
    str(ROOT.parent / "1"),
    str(ROOT.parent / "2"),
]

import kickmix as km
import dqi.verify_submission_dqi as v
from stream_model import decode, pub
from published_dqi import isomorphism
import numpy as np

spec = v.ProblemSpecDQI.benchmark_instance(num_shots=10000)
tables = v.FieldTables(km.GF2Field(12))
probs = v.optimal_dqi_weights(4095, 35, 12, 2016)**2
rng = np.random.default_rng(42)
s_samples, m_samples, ph_samples = v.sample_shots(spec, tables, rng, 10000, probs)

f, back, tr = isomorphism(pub)

for off in [2, 3, 4, 5]:
    failures = {}
    for i in range(10000):
        ss = [tr(int(x), f) for x in s_samples[i]]
        try:
            sig, aborted = decode(ss, offset=off)
            if aborted:
                reason = aborted[0]
                failures[reason] = failures.get(reason, 0) + 1
        except Exception as e:
            failures['assert_error'] = failures.get('assert_error', 0) + 1
    print(f"Offset {off}: Total failures out of 10000: {sum(failures.values())} | Reasons: {failures}")
