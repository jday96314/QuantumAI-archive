"""Verify and score the current submission in initial_experiments/4/agy."""

import argparse
import json
import math
from pathlib import Path
import sys
import time

AGY_DIR = Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy")
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")

sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
]

BASELINES = {
    "invsqrt": 2076774.233715282,
    "dqi": 24439463583672.93,
    "dfthcsos": 3355599588.3155375,
}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--submission-dir", type=Path, default=AGY_DIR)
    parser.add_argument("--output-dir", type=Path, default=AGY_DIR / "logs")
    parser.add_argument("--tasks", nargs="+", default=["dfthcsos", "invsqrt", "dqi"])
    parser.add_argument("--salt", default="")
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    import kickmix as km
    from invsqrt.verify_submission_invsqrt import verify as inv_verify, ProblemSpecInvSqrt
    from dqi.verify_submission_dqi import verify as dqi_verify, ProblemSpecDQI
    from dfthcsos.verify import verify as chem_verify

    results = {}

    for task in args.tasks:
        print(f"\n==========================================")
        print(f"VERIFYING TASK: {task}")
        print(f"==========================================")
        start = time.monotonic()
        if task == "dfthcsos":
            chem_path = args.submission_dir / "dfthcsos"
            print(f"Running DFTHC-SOS verification on {chem_path}...")
            result = chem_verify(chem_path, salt=args.salt)
        elif task == "invsqrt":
            c_path = args.submission_dir / "invsqrt.kmb"
            print(f"Loading {c_path}...")
            c = km.Circuit.from_file(c_path, format="kmb")
            print(f"Running InvSqrt verification (200,540,928 shots, {args.workers} workers)...")
            result = inv_verify(c, ProblemSpecInvSqrt.for_n_p(10), salt=args.salt, max_workers=args.workers)
        elif task == "dqi":
            c_path = args.submission_dir / "dqi.kmb"
            print(f"Loading {c_path}...")
            c = km.Circuit.from_file(c_path, format="kmb")
            print(f"Running DQI verification (102,400 shots, {args.workers} workers)...")
            result = dqi_verify(c, ProblemSpecDQI.benchmark_instance(), salt=args.salt, max_workers=args.workers)
        else:
            raise ValueError(f"Unknown task: {task}")

        wall_sec = time.monotonic() - start
        result["wall_seconds"] = wall_sec
        results[task] = result

        out_file = args.output_dir / f"verify_{task}.json"
        out_file.write_text(json.dumps(result, indent=2) + "\n")

        status = result.get("status", "FAIL")
        score = result.get("score", float("nan"))
        base = BASELINES[task]
        ratio = base / score if score else 0.0

        print(f"\nTASK: {task}")
        print(f"STATUS: {status}")
        print(f"COST: {score:,.4f}")
        print(f"BASELINE: {base:,.4f}")
        print(f"IMPROVEMENT: {ratio:.4f}x")
        print(f"TIME: {wall_sec:.2f}s")

        if status != "PASS":
            print(f"FAILED on task {task}: {result.get('error')}")
            sys.exit(1)

    print("\n==========================================")
    print("ALL REQUESTED TASKS PASSED!")
    print("==========================================")

    # If all 3 tasks were verified or logs exist for all 3:
    all_scores = {}
    for task in ["invsqrt", "dqi", "dfthcsos"]:
        if task in results:
            all_scores[task] = results[task]["score"]
        else:
            lf = args.output_dir / f"verify_{task}.json"
            if lf.exists():
                ld = json.loads(lf.read_text())
                if ld.get("status") == "PASS":
                    all_scores[task] = ld["score"]

    if len(all_scores) == 3:
        ratios = {k: BASELINES[k] / all_scores[k] for k in all_scores}
        geom_mean = (ratios["invsqrt"] * ratios["dqi"] * ratios["dfthcsos"]) ** (1.0 / 3.0)
        print("\nFINAL OVERALL SUBMISSION SCORE:")
        print(f"  invsqrt : {all_scores['invsqrt']:16.2f} ({ratios['invsqrt']:.4f}x)")
        print(f"  dqi     : {all_scores['dqi']:16.2f} ({ratios['dqi']:.4f}x)")
        print(f"  dfthcsos: {all_scores['dfthcsos']:16.2f} ({ratios['dfthcsos']:.4f}x)")
        print(f"  --> OVERALL GEOMETRIC MEAN: {geom_mean:.6f}x <--")

        summary = {
            "status": "PASS",
            "overall_score": geom_mean,
            "improvements": ratios,
            "costs": all_scores,
            "baselines": BASELINES,
        }
        (args.output_dir / "submission_score_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

if __name__ == "__main__":
    main()
