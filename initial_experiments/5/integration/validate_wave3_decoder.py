"""Accept the integrated decoder only after full native and restoration checks."""

import argparse
import fcntl
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROUND5 = HERE.parent
GAP = ROUND5 / "dqi-gap"
OUT = HERE / "wave3"
HOST = ROUND5.parents[1] / "from_host/verifiers/dqi/verify_submission_dqi.py"


def read(path):
    return json.loads(path.read_text())


def write(name, data):
    (OUT / "results" / name).write_text(json.dumps(data, indent=2) + "\n")


def sha(path):
    with path.open("rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


def residual_check(circuit_path):
    import kickmix as km
    import numpy as np
    from dqi.verify_submission_dqi import ProblemSpecDQI, to_bits, validate_registers

    cases = [json.loads(line) for line in (GAP / "parameters/residual_inputs.jsonl").read_text().splitlines()]
    reference = read(GAP / "results/best_residual_preservation.json")
    per = 128
    count = per * len(cases)
    rng = np.random.default_rng(2026102507)
    masks = rng.integers(0, 4096, (count, 4095), dtype=np.int64)
    syndromes = np.asarray([case["syndromes"] for case in cases for _ in range(per)], dtype=np.int64)
    truth = np.asarray([sum((v & int(masks[i, j])).bit_count()
                           for j, v in zip(cases[i // per]["positions"], cases[i // per]["values"])) % 2
                        for i in range(count)], dtype=bool)
    circuit = km.Circuit.from_file(str(circuit_path))
    syn, mask, anc = validate_registers(circuit, ProblemSpecDQI.benchmark_instance())
    simulator = km.Simulator(batch_size=count, seed=710521)
    simulator.use_same_registers_as(circuit)
    bits = to_bits(syndromes, 12)
    simulator.write_across_shots(km.array(syn), bits)
    simulator.write_across_shots(km.array(mask), to_bits(masks, 12))
    simulator.do(circuit)
    phase = simulator.read_phase_flipped_across_shots() ^ truth
    dirty = (simulator.read_across_shots(km.array(syn)) != bits).any(0) | simulator.read_across_shots(anc).any(0)
    rows = [dict(seed=case["seed"], shot=case["shot"], mask_cases=per,
                 truth_phase_errors=int(phase[i * per:(i + 1) * per].sum()),
                 value_errors=int(dirty[i * per:(i + 1) * per].sum()))
            for i, case in enumerate(cases)]
    simulator.do(circuit)
    twice = dict(phase_errors=int(simulator.read_phase_flipped_across_shots().sum()),
                 value_errors=int(((simulator.read_across_shots(km.array(syn)) != bits).any(0)
                                  | simulator.read_across_shots(anc).any(0)).sum()))
    matches = [row["truth_phase_errors"] for row in rows] == [row["truth_phase_errors"] for row in reference["cases"]]
    passed = not dirty.any() and not any(twice.values()) and matches
    report = dict(status="PASS" if passed else "FAIL", shots=count,
                  scope="same five known residual cases and masks; retained phase failures, clean inputs/scratch, twice-executed restoration",
                  circuit_sha256=sha(circuit_path), cases=rows, twice=twice,
                  reference_phase_error_counts_match=matches)
    write("gap_trimmed_residual.json", report)
    if not passed:
        raise ValueError("Residual preservation differs from the accepted gap decoder")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--residual-only", action="store_true")
    args = parser.parse_args()
    manifest = read(OUT / "parameters/gap_trimmed_manifest.json")
    circuit = OUT / manifest["circuit"]
    if sha(circuit) != manifest["circuit_sha256"]:
        raise ValueError("Candidate bytes changed")
    if args.residual_only:
        residual_check(circuit)
        return
    for source, digest in manifest["sources_sha256_relative_to_round5"].items():
        if sha(ROUND5 / source) != digest:
            raise ValueError(f"Frozen source changed: {source}")
    trimmed = ROUND5 / "dqi-core"
    accepted = read(trimmed / "parameters/trimmed_decoder_source_manifest.json")
    if sha(trimmed / "components/trimmed_middle.py") != accepted["sources"]["components/trimmed_middle.py"]:
        raise ValueError("Trimmed emitter differs from its original accepted version")
    if sha(HOST) != manifest["host_verifier_sha256"]:
        raise ValueError("Host verifier changed")
    from dqi.verify_submission_dqi import _initial_seed
    previous = read(HERE / "results/dqi_verify.json")
    paired_seed = _initial_seed(previous["circuit_sha256"], "")
    runs = [
        ("default", []),
        ("independent", ["--salt", "round5-wave3-trimmed-independent-20261025-a"]),
        ("paired", ["--seed", str(paired_seed)]),
    ]
    write("gap_trimmed_validation_plan.json", {
        "candidate_sha256": manifest["circuit_sha256"], "runs": runs,
        "default_score_reported": True,
        "paired_reference": "initial_experiments/5/integration/results/dqi_verify.json",
        "frozen_before_native_sampling": True,
    })
    reports = {}
    with open("/tmp/r6-fullverify.lock", "a") as lock:
        print("Waiting for shared full-verification lock.", flush=True)
        fcntl.flock(lock, fcntl.LOCK_EX)
        for name, options in runs:
            command = [sys.executable, str(HOST), str(circuit), "--max-workers", "4", *options]
            print(f"Starting full 102400-shot DQI run: {name}", flush=True)
            start = time.monotonic()
            with (OUT / f"results/gap_trimmed_verify_{name}.json").open("w") as out, (OUT / f"results/gap_trimmed_verify_{name}.err").open("w") as err:
                result = subprocess.run(command, stdout=out, stderr=err, cwd=OUT)
            if result.returncode:
                raise RuntimeError(f"Native {name} verification failed")
            report = read(OUT / f"results/gap_trimmed_verify_{name}.json")
            if report["status"] != "PASS" or report["verified_shots"] != 102400:
                raise ValueError(f"Invalid full-shot report: {name}")
            if report["circuit_sha256"] != manifest["circuit_sha256"] or report["num_qubits"] != previous["num_qubits"]:
                raise ValueError(f"Circuit identity or width mismatch: {name}")
            reports[name] = report
            write(f"gap_trimmed_verify_{name}_command.json", dict(command=command, seconds=time.monotonic() - start))
            print(f"{name}: PASS; value/phase errors {report['value_errors']}/{report['phase_errors']}; cost {report['score']}", flush=True)
        paired = reports["paired"]
        if any(paired[key] != previous[key] for key in ("value_errors", "phase_errors", "repetition_factor")):
            raise ValueError("Decoding differs on the original gap decoder's fixed stream")
        saved = read(OUT / "results/gap_trimmed_profile.json")["saved_core_T"]
        if abs((previous["toffoli_count"] - paired["toffoli_count"]) - saved) > 1e-6:
            raise ValueError("Paired native gate saving differs from exact emitted saving")
        for name in ("default", "independent"):
            if reports[name]["score"] >= previous["score"]:
                raise ValueError(f"Integrated decoder does not improve the prior cost: {name}")
        print("Checking 1024 independent adversarial cases and twice-execution restoration.", flush=True)
        command = [sys.executable, str(GAP / "verify_adversarial.py"), str(circuit),
                   "--shots", "1024", "--output", str(OUT / "results/gap_trimmed_adversarial.json")]
        with (OUT / "results/gap_trimmed_adversarial.log").open("w") as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, cwd=GAP, check=True)
        print("Checking known residual input/scratch and twice-execution restoration.", flush=True)
        with (OUT / "results/gap_trimmed_residual.log").open("w") as log:
            subprocess.run([sys.executable, str(Path(__file__).resolve()), "--residual-only"],
                           stdout=log, stderr=subprocess.STDOUT, cwd=OUT, check=True)
    if sha(HOST) != manifest["host_verifier_sha256"]:
        raise ValueError("Host verifier changed during validation")
    write("gap_trimmed_acceptance.json", dict(
        status="VERIFIED_IMPROVEMENT", circuit=str(circuit.relative_to(OUT)),
        circuit_sha256=manifest["circuit_sha256"], Q=reports["default"]["num_qubits"],
        default_score=reports["default"]["score"], previous_score=previous["score"],
        improvement=previous["score"] / reports["default"]["score"],
        saved_static_Toffolis=saved, paired_decoding_identical=True,
        native_runs={name: f"results/gap_trimmed_verify_{name}.json" for name, _ in runs},
        adversarial="results/gap_trimmed_adversarial.json", residual="results/gap_trimmed_residual.json"))
    print("Integrated trimmed gap decoder accepted after all native checks.", flush=True)


if __name__ == "__main__":
    main()
