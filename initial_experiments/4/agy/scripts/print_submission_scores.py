#!/usr/bin/env python3
import sys
from pathlib import Path

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
]

from dfthcsos.verify import (
    load_tree_edges,
    load_inputs,
    load_circuits,
    evaluate_circuits,
    SUBMISSION_FILE
)

candidates = [
    ROOT / "initial_experiments/3/candidates/residual719_q8_initial_submission",
    ROOT / "initial_experiments/3/candidates/residual719_continue_submission",
    ROOT / "initial_experiments/3/candidates/residual719_highlr_submission",
    ROOT / "initial_experiments/3/candidates/residual719_lbfgs2_submission",
    ROOT / "initial_experiments/3/candidates/residual719_wb_submission",
    ROOT / "initial_experiments/3/candidates/lbfgs500_q10c375_submission",
    ROOT / "initial_experiments/3/candidates/global500_coord11_c375_exact_submission",
    ROOT / "initial_experiments/3/candidates/global500_coord11_c375_qr2_submission",
    ROOT / "initial_experiments/3/candidates/global500_coord11_c400_submission",
    ROOT / "initial_experiments/3/candidates/global500_q13_c400_qr2_submission",
]

base = 3355599588.3155375

for cdir in candidates:
    if not cdir.exists():
        continue
    try:
        edges = load_tree_edges(None, cdir)
        inputs = load_inputs(cdir / SUBMISSION_FILE)
        family = load_circuits(cdir, format="kmb")
        res = evaluate_circuits(inputs, family, edges, shots_per_address=100)
        score = float(res["score"])
        resrc = res["resources"]
        m = res["metrics"]
        t = int(resrc['toffoli'])
        q = int(resrc['qubits'])
        print(f"{cdir.name:45s} | Score: {score:14,.2f} ({base/score:6.3f}x) | T: {t:6,d} | Q: {q:4d} | err2: {m['error_2']*1000:6.3f} mHa")
    except Exception as e:
        print(f"{cdir.name:45s} | Failed: {e}")
