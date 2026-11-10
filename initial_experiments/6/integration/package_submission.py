"""Combine the best verified round-six components without changing their bytes.

From the repository root:
  source initial_experiments/6/agent20/chem_env.sh
  "$EXP_PYTHON" initial_experiments/6/integration/package_submission.py
  "$EXP_CHEM_PYTHON" initial_experiments/6/agent20/audit_best_metric.py initial_experiments/6/integration

Full native reports are reused only after digest checks. No native sampling,
host verifier, source workspace, or saved milestone is changed.
"""

import importlib.util
from pathlib import Path
import shutil
import zipfile


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
AGENT14 = HERE.parent / "agent14"
AGENT20 = HERE.parent / "agent20"


def main():
    spec = importlib.util.spec_from_file_location(
        "round5_packager", ROOT / "initial_experiments/5/integration/package_submission.py"
    )
    package = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(package)
    package.HERE = HERE
    package.RESULTS = HERE / "results"
    package.EXTRACTED = HERE / "submission"
    package.ARCHIVE = HERE / "submission.zip"
    package.PREVIOUS_SUMMARY = AGENT14 / "submission.summary.json"
    package.EVIDENCE = {
        "invsqrt": AGENT20 / "results/public_phuc_invsqrt_default_verify.json",
        "dqi": AGENT14 / "results/raw_cf_destructive_n70_verify_default.json",
        "dfthcsos": AGENT20 / "results/paired_chemistry_default_verify.json",
    }
    independent = {
        "invsqrt": AGENT20 / "results/public_phuc_invsqrt_independent_verify.json",
        "dqi": AGENT14 / "results/raw_cf_destructive_n70_verify_independent.json",
        "dfthcsos": AGENT20 / "results/paired_chemistry_independent_verify.json",
    }
    reports = {task: package.read_json(path) for task, path in package.EVIDENCE.items()}
    files = {
        "invsqrt.kmb": AGENT20 / "public/phuc/invsqrt.kmb",
        "dqi.kmb": AGENT14 / "submission/dqi.kmb",
        "dfthcsos/submission.npz": AGENT20 / "candidates/paired_chemistry/submission.npz",
        "NOTICE.md": HERE / "NOTICE.md",
        "LICENSE": AGENT20 / "best/submission/LICENSE",
    }
    chemistry = AGENT20 / "candidates/paired_chemistry"
    for k in range(reports["dfthcsos"]["num_circuits"]):
        name = f"angle_k_{k}.kmb"
        files[f"dfthcsos/{name}"] = chemistry / name
    for task, report in reports.items():
        package.check_report(task, report)
    for task in ("invsqrt", "dqi"):
        assert package.sha256(files[f"{task}.kmb"]) == reports[task]["circuit_sha256"], task
    digest = package.chemistry_digest(chemistry, reports["dfthcsos"]["num_circuits"])
    assert digest == reports["dfthcsos"]["submission_sha256"]

    # Independently sampled reports must describe these same selected components.
    independent_reports = {}
    for task, source in independent.items():
        report = package.read_json(source)
        package.check_report(task, report)
        key = "submission_sha256" if task == "dfthcsos" else "circuit_sha256"
        assert report[key] == reports[task][key], task
        independent_reports[task] = report
        destination = package.RESULTS / "independent_evidence" / f"{task}.json"
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
    for task, source in package.EVIDENCE.items():
        destination = package.RESULTS / "source_evidence" / f"{task}.json"
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)

    # Keep the existing chemistry's exact physical parameters.
    assert package.sha256(files["dfthcsos/submission.npz"]) == package.sha256(
        AGENT14 / "submission/dfthcsos/submission.npz"
    )
    members = {}
    temporary = package.ARCHIVE.with_suffix(".zip.tmp")
    with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for name, source in sorted(files.items()):
            info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            with source.open("rb") as src, archive.open(info, "w", force_zip64=True) as dst:
                shutil.copyfileobj(src, dst, length=2**20)
            members[name] = {"source": str(source.relative_to(ROOT)),
                             "sha256": package.sha256(source), "bytes": source.stat().st_size}
    temporary.replace(package.ARCHIVE)
    with zipfile.ZipFile(package.ARCHIVE) as archive:
        assert len(archive.namelist()) == len(files) and set(archive.namelist()) == set(files)
        for name in sorted(files):
            destination = package.EXTRACTED / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            with archive.open(name) as src, destination.open("wb") as dst:
                shutil.copyfileobj(src, dst, length=2**20)
            assert package.sha256(destination) == members[name]["sha256"], name
    host_files = ("metric.py", "verifiers/invsqrt/verify_submission_invsqrt.py",
                  "verifiers/dqi/verify_submission_dqi.py", "verifiers/dfthcsos/verify.py")
    manifest = {
        "archive": "submission.zip", "archive_sha256": package.sha256(package.ARCHIVE),
        "archive_bytes": package.ARCHIVE.stat().st_size, "members": members,
        "source_evidence": {task: str(path.relative_to(ROOT)) for task, path in package.EVIDENCE.items()},
        "independent_evidence": {task: str(path.relative_to(ROOT)) for task, path in independent.items()},
        "host_sha256": {str((package.HOST / name).relative_to(ROOT)): package.sha256(package.HOST / name)
                        for name in host_files},
    }
    package.write_json(package.RESULTS / "package_manifest.json", manifest)
    package.verify(2, reuse_tasks=tuple(package.BASELINES))
    package.score_and_audit()
    independent_score = 1.0
    for task, report in independent_reports.items():
        independent_score *= package.BASELINES[task] / report["score"]
    package.write_json(package.RESULTS / "independent_audit.json", {
        "status": "PASS", "archive_sha256": manifest["archive_sha256"],
        "overall_score": independent_score ** (1 / 3),
        "task_scores": {task: report["score"] for task, report in independent_reports.items()},
        "verified_shots": package.EXPECTED_SHOTS,
        "source_evidence": manifest["independent_evidence"],
        "scope": "Existing full independently sampled native reports bound to the same packaged component digests; no new native run",
    })


if __name__ == "__main__":
    main()
