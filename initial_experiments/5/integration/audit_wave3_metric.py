"""Apply the unmodified host archive readers and score combiner to wave three.

Run with the existing torch2.13 Python (which includes pandas), after sourcing
integration/env.sh. Native verification is recorded separately and bound to
the extracted archive by package_submission.py; this script does not rerun it.
"""

import hashlib
import importlib.util
import json
import math
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / "wave3"


def read(path):
    return json.loads(path.read_text())


def sha(path):
    with path.open("rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


def main():
    metric_path = ROOT / "from_host/metric.py"
    before = sha(metric_path)
    manifest = read(OUT / "results/package_manifest.json")
    if before != manifest["host_sha256"]["from_host/metric.py"]:
        raise ValueError("Metric changed after packaging")
    spec = importlib.util.spec_from_file_location("wave3_unmodified_host_metric", metric_path)
    metric = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(metric)
    extracted = OUT / "submission"
    summary = read(OUT / "submission.summary.json")
    if sha(OUT / "submission.zip") != summary["sha256"]:
        raise ValueError("Archive changed")
    for name, registers in (("invsqrt", 1), ("dqi", 70 + 4095)):
        path = metric.find_one(extracted, f"{name}.km[xb]", f"{name} circuit")
        circuit = metric.read_circuit(path, registers)
        print(f"Official metric circuit reader: {name} PASS", flush=True)
        del circuit
    chem = metric.find_one(extracted, "dfthcsos", "dfthcsos directory", is_dir=True)
    paths = list(chem.glob("angle_k_*.km[xb]"))
    actual = {path.stem: path for path in paths}
    stems = [f"angle_k_{k}" for k in range(len(paths))]
    if not stems or len(actual) != len(paths) or sorted(actual) != sorted(stems):
        raise ValueError("Invalid chemistry circuit numbering")
    for stem in stems:
        metric.read_circuit(actual[stem], 128)
    arrays = metric.read_npz(metric.find_one(chem, "submission.npz", "chemistry NPZ"), {
        "U": ("f", ("R", "B", 76)),
        "W": ("f", ("R", "B", "C")),
        "H2_I": ("f", ("R", "C")),
        "H1_I": ("f", ()),
        "Hsym_lower": ("f", (76 * 77 // 2,)),
        "tree_edges": ("iu", (75, 2)),
    })
    print(f"Official metric chemistry readers: {len(stems)} circuits and {len(arrays)} arrays PASS", flush=True)
    tasks = list(summary["task_scores"])
    solution = pd.DataFrame({"row_id": list(range(len(tasks))), "benchmark": tasks,
                             "baseline_cost": [summary["baseline_costs"][task] for task in tasks]})
    native = pd.DataFrame({"row_id": list(range(len(tasks))),
                          "cost": [summary["task_scores"][task] for task in tasks],
                          "seconds": [summary["verification_runs"][task]["seconds"] for task in tasks]})
    score, details = metric.score(solution, native, "row_id", verification_shot_divisor=1)
    if not math.isclose(score, summary["overall_score"], rel_tol=1e-14):
        raise ValueError("Official score combiner disagrees with package summary")
    if sha(metric_path) != before:
        raise ValueError("Metric changed during audit")
    report = dict(status="PASS", archive_sha256=summary["sha256"], metric_sha256=before,
                  overall_score=score, task_details=details, angle_circuits=len(stems),
                  scope="unmodified metric find_one/read_circuit/read_npz and score; costs from separately hash-bound full native reports; no shortened-shot metric evaluation")
    (OUT / "results/metric_preflight_and_score.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
