"""Replay the new default stream on the original accepted gap decoder.

The new default run observed one phase error. This controlled comparison checks
that trimming preserves that counter and the exact gate saving, without using
a more favorable validation stream as the reported submission score.
"""

import fcntl
import json
from pathlib import Path
import subprocess
import sys
import time

from dqi.verify_submission_dqi import _initial_seed

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / "wave3"


def main():
    candidate = json.loads((OUT / "results/gap_trimmed_verify_default.json").read_text())
    previous = json.loads((HERE / "results/dqi_verify.json").read_text())
    seed = _initial_seed(candidate["circuit_sha256"], "")
    circuit = ROOT / "initial_experiments/5/dqi-gap/candidates/gap_best_n70_m4095.kmb"
    command = [sys.executable, str(ROOT / "from_host/verifiers/dqi/verify_submission_dqi.py"),
               str(circuit), "--max-workers", "4", "--seed", str(seed)]
    print("Replaying original gap decoder on the new candidate's fixed default stream.", flush=True)
    started = time.monotonic()
    with open("/tmp/r6-fullverify.lock", "a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        with (OUT / "results/gap_original_new_default.json").open("w") as out, (OUT / "results/gap_original_new_default.err").open("w") as err:
            subprocess.run(command, stdout=out, stderr=err, check=True)
    original = json.loads((OUT / "results/gap_original_new_default.json").read_text())
    assert original["status"] == "PASS" and original["verified_shots"] == 102400
    assert original["circuit_sha256"] == previous["circuit_sha256"]
    for key in ("value_errors", "phase_errors", "repetition_factor", "num_qubits"):
        assert original[key] == candidate[key], key
    difference = original["toffoli_count"] - candidate["toffoli_count"]
    assert abs(difference - 23352) < 1e-6, difference
    report = dict(status="PASS", scope="original decoder and trimmed decoder on the new default seed; no stream selection",
                  original_circuit_sha256=original["circuit_sha256"],
                  new_circuit_sha256=candidate["circuit_sha256"],
                  value_errors=candidate["value_errors"], phase_errors=candidate["phase_errors"],
                  repetition_factor=candidate["repetition_factor"], saved_mean_Toffolis=difference,
                  command=command, seconds=time.monotonic() - started)
    (OUT / "results/new_default_stream_comparison.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
