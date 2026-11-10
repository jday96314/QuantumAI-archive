"""Package the verified round-five components and score the extracted archive.

Run with the Python environment selected by env.sh. All outputs stay beside
this script; inputs and the host verifiers are read-only.
"""

import argparse
import fcntl
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import time
import zipfile

import kickmix as km
import numpy as np
from dfthcsos.scoreinfra.circuits import submission_sha256


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ROUND4 = ROOT / "initial_experiments/4/codex"
ROUND5 = ROOT / "initial_experiments/5"
HOST = ROOT / "from_host"
RESULTS = HERE / "results"
EXTRACTED = HERE / "submission"
ARCHIVE = HERE / "submission.zip"
INVERSE = ROUND4 / "candidates/zero_head_affine.kmb"
DQI = ROUND5 / "dqi-gap/candidates/gap_best_n70_m4095.kmb"
CHEMISTRY = ROUND5 / "chem-audit/candidates/banked"
PREVIOUS_SUMMARY = ROUND4 / "submission.summary.json"
BASELINES = {
    "invsqrt": 2076774.233715,
    "dqi": 24439463583672.929688,
    "dfthcsos": 3355599588.315537,
}
EXPECTED_SHOTS = {"invsqrt": 200540928, "dqi": 102400, "dfthcsos": 50880}
EVIDENCE = {
    "invsqrt": ROUND4 / "results/zero_head_affine_verify.json",
    "dqi": ROUND5 / "dqi-gap/results/best_verify_default.json",
    "dfthcsos": ROUND5 / "chem-audit/results/banked_native.json",
}


def read_json(path):
    return json.loads(path.read_text())


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n")


def sha256(path):
    with path.open("rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


def chemistry_digest(directory, count):
    circuits = [km.Circuit.from_file(str(directory / f"angle_k_{k}.kmb"), format="kmb")
                for k in range(count)]
    with np.load(directory / "submission.npz", allow_pickle=False) as arrays:
        params = {key: arrays[key] for key in ("U", "W", "H2_I", "H1_I", "Hsym_lower")}
        return submission_sha256(circuits, params, arrays["tree_edges"])


def check_report(task, report):
    if report["status"] != "PASS" or report["verified_shots"] != EXPECTED_SHOTS[task]:
        raise ValueError(f"{task}: expected PASS at full default shots, got {report}")
    if not math.isfinite(report["score"]) or report["score"] <= 0:
        raise ValueError(f"{task}: invalid score")


def build():
    RESULTS.mkdir(parents=True, exist_ok=True)
    evidence = {task: read_json(path) for task, path in EVIDENCE.items()}
    for task, report in evidence.items():
        check_report(task, report)
    chem = CHEMISTRY
    files = {
        "invsqrt.kmb": INVERSE,
        "dqi.kmb": DQI,
        "dfthcsos/submission.npz": chem / "submission.npz",
    }
    for k in range(evidence["dfthcsos"]["num_circuits"]):
        files[f"dfthcsos/angle_k_{k}.kmb"] = chem / f"angle_k_{k}.kmb"
    if sha256(files["invsqrt.kmb"]) != evidence["invsqrt"]["circuit_sha256"]:
        raise ValueError("Inverse artifact does not match its native evidence")
    if sha256(files["dqi.kmb"]) != evidence["dqi"]["circuit_sha256"]:
        raise ValueError("DQI artifact does not match its native evidence")
    if chemistry_digest(chem, evidence["dfthcsos"]["num_circuits"]) != evidence["dfthcsos"]["submission_sha256"]:
        raise ValueError("Chemistry artifacts do not match their native evidence")
    for task, source in EVIDENCE.items():
        destination = RESULTS / "source_evidence" / f"{task}.json"
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
    gap_evidence = RESULTS / "source_evidence"
    if DQI != ROUND5 / "dqi-gap/candidates/gap_best_n70_m4095.kmb":
        gap_evidence = gap_evidence / "previous_gap_reference"
        gap_evidence.mkdir(parents=True, exist_ok=True)
    for name in ("best_validation_summary.json", "best_adversarial.json", "best_residual_preservation.json"):
        shutil.copyfile(ROUND5 / "dqi-gap/results" / name, gap_evidence / name)
    print("Source hashes and chemistry canonical digest match full native evidence.", flush=True)
    temporary = ARCHIVE.with_suffix(".zip.tmp")
    members = {}
    with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for name, source in sorted(files.items()):
            print(f"Packaging {name} ({source.stat().st_size:,} bytes)", flush=True)
            info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            with source.open("rb") as src, archive.open(info, "w", force_zip64=True) as dst:
                shutil.copyfileobj(src, dst, length=2**20)
            members[name] = {"source": str(source.relative_to(ROOT)), "sha256": sha256(source),
                             "bytes": source.stat().st_size}
    temporary.replace(ARCHIVE)
    # Extract exactly the archive that will be delivered, and bind every byte.
    with zipfile.ZipFile(ARCHIVE) as archive:
        if len(archive.namelist()) != len(files) or set(archive.namelist()) != set(files):
            raise ValueError("Unexpected archive layout")
        for name in sorted(files):
            destination = EXTRACTED / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            with archive.open(name) as src, destination.open("wb") as dst:
                shutil.copyfileobj(src, dst, length=2**20)
            if sha256(destination) != members[name]["sha256"]:
                raise ValueError(f"Archive byte mismatch: {name}")
    manifest = {
        "archive": "submission.zip", "archive_sha256": sha256(ARCHIVE),
        "archive_bytes": ARCHIVE.stat().st_size, "members": members,
        "source_evidence": {task: str(path.relative_to(ROOT)) for task, path in EVIDENCE.items()},
        "host_sha256": {str(path.relative_to(ROOT)): sha256(path) for path in (
            HOST / "metric.py", HOST / "verifiers/invsqrt/verify_submission_invsqrt.py",
            HOST / "verifiers/dqi/verify_submission_dqi.py", HOST / "verifiers/dfthcsos/verify.py")},
    }
    write_json(RESULTS / "package_manifest.json", manifest)
    print(f"Archive created and extracted: {ARCHIVE.stat().st_size:,} bytes.", flush=True)


def verify(workers, reuse_tasks=()):
    commands = {
        "invsqrt": [str(HOST / "verifiers/invsqrt/verify_submission_invsqrt.py"),
                    str(EXTRACTED / "invsqrt.kmb"), "--max-workers", str(workers)],
        "dfthcsos": [str(HOST / "verifiers/dfthcsos/verify.py"), str(EXTRACTED / "dfthcsos")],
        "dqi": [str(HOST / "verifiers/dqi/verify_submission_dqi.py"),
                str(EXTRACTED / "dqi.kmb"), "--max-workers", str(workers)],
    }
    timing = {}
    with open("/tmp/r6-fullverify.lock", "a") as lock:
        if set(commands) - set(reuse_tasks):
            print("Waiting for the shared full-verification lock.", flush=True)
            fcntl.flock(lock, fcntl.LOCK_EX)
        for task, command in commands.items():
            if task in reuse_tasks:
                # build() has already bound the input bytes/canonical chemistry
                # digest to this full-shot report. The final audit binds them
                # again after extraction from the delivered archive.
                report = read_json(EVIDENCE[task])
                check_report(task, report)
                shutil.copyfile(EVIDENCE[task], RESULTS / f"{task}_verify.json")
                timing[task] = {"mode": "reused_full_native_verification", "seconds": 0,
                                "source": str(EVIDENCE[task].relative_to(ROOT)),
                                "binding": "circuit/canonical chemistry digest checked before packaging and after ZIP extraction"}
                write_json(RESULTS / "verification_runs.json", timing)
                print(f"{task}: reused full-shot PASS bound to the exact packaged bytes.", flush=True)
                continue
            print(f"Starting full-default {task} verification from extracted ZIP.", flush=True)
            started = time.monotonic()
            with (RESULTS / f"{task}_verify.json").open("w") as out, (RESULTS / f"{task}_verify.err").open("w") as err:
                run = subprocess.run([sys.executable, *command], stdout=out, stderr=err, cwd=HERE)
            if run.returncode:
                raise RuntimeError(f"{task} verification failed; inspect results/{task}_verify.err and .json")
            report = read_json(RESULTS / f"{task}_verify.json")
            check_report(task, report)
            timing[task] = {"mode": "fresh_full_native_verification", "seconds": time.monotonic() - started,
                            "command": [sys.executable, *command], "seed": "default hash-derived; no override"}
            write_json(RESULTS / "verification_runs.json", timing)
            print(f"{task}: PASS, cost {report['score']:.10f}, {timing[task]['seconds']:.1f}s", flush=True)


def score_and_audit():
    manifest = read_json(RESULTS / "package_manifest.json")
    if sha256(ARCHIVE) != manifest["archive_sha256"]:
        raise ValueError("Archive changed after packaging")
    for name, expected in manifest["host_sha256"].items():
        if sha256(ROOT / name) != expected:
            raise ValueError(f"Host file changed: {name}")
    reports = {task: read_json(RESULTS / f"{task}_verify.json") for task in BASELINES}
    for task, report in reports.items():
        check_report(task, report)
    # Reopen the ZIP independently of extraction and native verification.
    hashes = {}
    with zipfile.ZipFile(ARCHIVE) as archive:
        if len(archive.namelist()) != len(manifest["members"]) or set(archive.namelist()) != set(manifest["members"]):
            raise ValueError("Archive members changed")
        for name, expected in manifest["members"].items():
            with archive.open(name) as member:
                hashes[name] = hashlib.file_digest(member, "sha256").hexdigest()
            if hashes[name] != expected["sha256"] or sha256(EXTRACTED / name) != hashes[name]:
                raise ValueError(f"Archive or extracted artifact changed: {name}")
    for task in ("invsqrt", "dqi"):
        if hashes[f"{task}.kmb"] != reports[task]["circuit_sha256"]:
            raise ValueError(f"{task}: native digest mismatch")
    digest = chemistry_digest(EXTRACTED / "dfthcsos", reports["dfthcsos"]["num_circuits"])
    if digest != reports["dfthcsos"]["submission_sha256"]:
        raise ValueError("Chemistry native digest mismatch")
    previous = read_json(PREVIOUS_SUMMARY)
    improvements = {task: BASELINES[task] / reports[task]["score"] for task in BASELINES}
    # Same geometric mean as from_host/metric.py:score, with the recorded local
    # reference costs. Kaggle's secret split salts and baselines are unavailable.
    overall = math.prod(improvements.values()) ** (1 / len(improvements))
    dqi_only = math.prod(BASELINES[task] / (reports[task]["score"] if task == "dqi"
                         else previous["task_scores"][task]) for task in BASELINES) ** (1 / 3)
    summary = {
        "status": "PASS", "overall_score": overall, "score_scope": "local default-seed full-shot score",
        "score_formula": "geometric_mean(baseline_cost / verified_task_cost); from_host/metric.py:score",
        "baseline_costs": BASELINES, "task_scores": {task: report["score"] for task, report in reports.items()},
        "improvements": improvements, "previous_overall_score": previous["overall_score"],
        "overall_improvement": overall / previous["overall_score"],
        "previous_summary": str(PREVIOUS_SUMMARY.relative_to(ROOT)),
        "previous_archive_sha256": previous["sha256"],
        "task_improvements": {task: previous["task_scores"][task] / reports[task]["score"]
                              for task in BASELINES},
        "dqi_improvement": previous["task_scores"]["dqi"] / reports["dqi"]["score"],
        "dqi_only_overall_score": dqi_only,
        "archive": "submission.zip", "sha256": manifest["archive_sha256"],
        "archive_bytes": manifest["archive_bytes"],
        "sources": {task: f"results/{task}_verify.json" for task in BASELINES},
        "verified_shots": EXPECTED_SHOTS,
        "verification_runs": read_json(RESULTS / "verification_runs.json"),
    }
    write_json(HERE / "submission.summary.json", summary)
    write_json(RESULTS / "archive_audit.json", {
        "status": "PASS", "archive_sha256": manifest["archive_sha256"], "members": hashes,
        "chemistry_digest": digest, "verified_shots": EXPECTED_SHOTS,
        "overall_score": overall,
        "scope": "ZIP reopened; every member and extracted file hashed; full native reports bound to circuit/canonical chemistry digests; host hashes unchanged",
    })
    print(json.dumps(summary, indent=2), flush=True)


def main():
    global HERE, RESULTS, EXTRACTED, ARCHIVE, DQI, CHEMISTRY, PREVIOUS_SUMMARY
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--audit-only", action="store_true", help="Recheck the existing archive and full native results")
    parser.add_argument("--output-dir", type=Path, default=HERE, help="Use a separate directory to preserve earlier archives")
    parser.add_argument("--chemistry", type=Path, default=CHEMISTRY)
    parser.add_argument("--chemistry-result", type=Path, default=EVIDENCE["dfthcsos"])
    parser.add_argument("--dqi", type=Path, default=DQI)
    parser.add_argument("--dqi-result", type=Path, default=EVIDENCE["dqi"])
    parser.add_argument("--previous-summary", type=Path, default=PREVIOUS_SUMMARY)
    parser.add_argument("--reuse-verified-task", action="append", choices=tuple(BASELINES), default=[],
                        help="Reuse full-shot input evidence after binding it to the exact archive bytes")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    HERE = args.output_dir.resolve()
    RESULTS, EXTRACTED, ARCHIVE = HERE / "results", HERE / "submission", HERE / "submission.zip"
    DQI, CHEMISTRY = args.dqi.resolve(), args.chemistry.resolve()
    PREVIOUS_SUMMARY = args.previous_summary.resolve()
    EVIDENCE["dfthcsos"] = args.chemistry_result.resolve()
    EVIDENCE["dqi"] = args.dqi_result.resolve()
    if not args.audit_only:
        build()
        verify(args.workers, args.reuse_verified_task)
    score_and_audit()


if __name__ == "__main__":
    main()
