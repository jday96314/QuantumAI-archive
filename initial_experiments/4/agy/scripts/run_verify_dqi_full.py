"""Run full verification of dqi.kmb on 102,400 shots and save to logs/verify_dqi_full.json."""

import json
import subprocess
import time

def main():
    cmd = [
        "/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy/scripts/run.sh",
        "/media/shared/ExpansionDrive2/QuantumAI/from_host/verifiers/dqi/verify_submission_dqi.py",
        "/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy/dqi.kmb",
        "--max-workers", "4",
    ]
    print("Running DQI full verification on 102,400 shots (4 workers)...")
    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    duration = time.time() - t0
    print(f"Finished in {duration:.1f}s, exit code: {proc.returncode}")
    
    if proc.returncode != 0:
        print("STDOUT:\n", proc.stdout)
        print("STDERR:\n", proc.stderr)
        return

    data = json.loads(proc.stdout)
    out_log = "/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy/logs/verify_dqi_full.json"
    with open(out_log, "w") as f:
        json.dump(data, f, indent=2)
    print(f"Saved {out_log}")
    print(f"Score: {data.get('score'):,.2f}")
    base = 24439463583672.93
    ratio = base / data.get("score")
    print(f"DQI Improvement Ratio: {ratio:.5f}x")

if __name__ == "__main__":
    main()
