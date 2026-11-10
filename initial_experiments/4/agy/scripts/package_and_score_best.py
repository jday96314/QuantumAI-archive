"""Package and score a submission with the best-scoring components thus far.

Builds submission.zip in /media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy,
verifies all components from the archive, and computes the geometric mean score.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import shutil
import sys
import tempfile
import time
import zipfile

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

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()

def create_submission_zip(agy_dir: Path, zip_path: Path):
    print(f"Creating submission archive: {zip_path}")
    if zip_path.exists():
        zip_path.unlink()
    
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        # 1. invsqrt.kmb
        inv_src = agy_dir / "invsqrt.kmb"
        if not inv_src.exists():
            raise FileNotFoundError(f"Missing {inv_src}")
        print(f"  Adding {inv_src.name} ({inv_src.stat().st_size:,} bytes)...")
        z.write(inv_src, arcname="invsqrt.kmb")
        
        # 2. dqi.kmb
        dqi_src = agy_dir / "dqi.kmb"
        if not dqi_src.exists():
            raise FileNotFoundError(f"Missing {dqi_src}")
        print(f"  Adding {dqi_src.name} ({dqi_src.stat().st_size:,} bytes)...")
        z.write(dqi_src, arcname="dqi.kmb")
        
        # 3. dfthcsos/
        chem_dir = agy_dir / "dfthcsos"
        if not chem_dir.is_dir():
            raise FileNotFoundError(f"Missing {chem_dir}")
        for item in sorted(chem_dir.iterdir()):
            if item.is_file() and not item.name.startswith("."):
                arcname = f"dfthcsos/{item.name}"
                print(f"  Adding {arcname} ({item.stat().st_size:,} bytes)...")
                z.write(item, arcname=arcname)
                
    zip_size = zip_path.stat().st_size
    zip_sha = sha256_file(zip_path)
    print(f"Submission archive created: {zip_size:,} bytes, SHA-256: {zip_sha}")
    return zip_sha, zip_size

def verify_extracted_package(extract_dir: Path, workers: int = 8, salt: str = ""):
    import kickmix as km
    from invsqrt.verify_submission_invsqrt import verify as inv_verify, ProblemSpecInvSqrt
    from dqi.verify_submission_dqi import verify as dqi_verify, ProblemSpecDQI
    from dfthcsos.verify import verify as chem_verify

    results = {}

    # Task 1: DFTHC-SOS
    print("\n" + "=" * 60)
    print("VERIFYING TASK 1/3: DFTHC-SOS Chemistry")
    print("=" * 60)
    t0 = time.monotonic()
    chem_sub = extract_dir / "dfthcsos"
    res_chem = chem_verify(chem_sub, salt=salt)
    t_chem = time.monotonic() - t0
    res_chem["wall_seconds"] = t_chem
    results["dfthcsos"] = res_chem
    print(f"DFTHC-SOS status: {res_chem.get('status')} in {t_chem:.2f}s, cost: {res_chem.get('score', 0):,.2f}")
    if res_chem.get("status") != "PASS":
        raise RuntimeError(f"DFTHC-SOS verification failed: {res_chem.get('error')}")

    # Task 2: Inverse Square Root
    print("\n" + "=" * 60)
    print("VERIFYING TASK 2/3: Inverse Square Root (InvSqrt)")
    print("=" * 60)
    t0 = time.monotonic()
    inv_c = km.Circuit.from_file(extract_dir / "invsqrt.kmb", format="kmb")
    res_inv = inv_verify(inv_c, ProblemSpecInvSqrt.for_n_p(10), salt=salt, max_workers=workers)
    t_inv = time.monotonic() - t0
    res_inv["wall_seconds"] = t_inv
    results["invsqrt"] = res_inv
    print(f"InvSqrt status: {res_inv.get('status')} in {t_inv:.2f}s, cost: {res_inv.get('score', 0):,.2f}")
    if res_inv.get("status") != "PASS":
        raise RuntimeError(f"InvSqrt verification failed: {res_inv.get('error')}")

    # Task 3: DQI
    print("\n" + "=" * 60)
    print("VERIFYING TASK 3/3: Discrete Quantum Inference (DQI)")
    print("=" * 60)
    t0 = time.monotonic()
    dqi_c = km.Circuit.from_file(extract_dir / "dqi.kmb", format="kmb")
    res_dqi = dqi_verify(dqi_c, ProblemSpecDQI.benchmark_instance(), salt=salt, max_workers=workers)
    t_dqi = time.monotonic() - t0
    res_dqi["wall_seconds"] = t_dqi
    results["dqi"] = res_dqi
    print(f"DQI status: {res_dqi.get('status')} in {t_dqi:.2f}s, cost: {res_dqi.get('score', 0):,.2f}")
    if res_dqi.get("status") != "PASS":
        raise RuntimeError(f"DQI verification failed: {res_dqi.get('error')}")

    return results

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--salt", default="")
    args = parser.parse_args()

    zip_path = AGY_DIR / "submission.zip"
    zip_sha, zip_size = create_submission_zip(AGY_DIR, zip_path)

    # Extract to tempdir for strictly hermetic end-to-end verification
    extract_dir = Path("/tmp/verify_best_sub")
    if extract_dir.exists():
        shutil.rmtree(extract_dir)
    extract_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\nExtracting {zip_path} to {extract_dir}...")
    with zipfile.ZipFile(zip_path, "r") as z:
        z.extractall(extract_dir)

    results = verify_extracted_package(extract_dir, workers=args.workers, salt=args.salt)

    # Compute overall geometric mean score
    inv_cost = results["invsqrt"]["score"]
    dqi_cost = results["dqi"]["score"]
    chem_cost = results["dfthcsos"]["score"]

    inv_ratio = BASELINES["invsqrt"] / inv_cost
    dqi_ratio = BASELINES["dqi"] / dqi_cost
    chem_ratio = BASELINES["dfthcsos"] / chem_cost

    overall_score = (inv_ratio * dqi_ratio * chem_ratio) ** (1.0 / 3.0)

    summary = {
        "status": "PASS",
        "archive": str(zip_path),
        "archive_size_bytes": zip_size,
        "archive_sha256": zip_sha,
        "overall_score": overall_score,
        "improvements": {
            "invsqrt": inv_ratio,
            "dqi": dqi_ratio,
            "dfthcsos": chem_ratio,
        },
        "costs": {
            "invsqrt": inv_cost,
            "dqi": dqi_cost,
            "dfthcsos": chem_cost,
        },
        "baselines": BASELINES,
        "wall_seconds": {
            "invsqrt": results["invsqrt"]["wall_seconds"],
            "dqi": results["dqi"]["wall_seconds"],
            "dfthcsos": results["dfthcsos"]["wall_seconds"],
        },
        "component_details": {
            "invsqrt": {
                "qubits": results["invsqrt"]["num_qubits"],
                "toffoli_count": results["invsqrt"]["toffoli_count"],
                "toffoli_equivalent": results["invsqrt"]["toffoli_equivalent"],
                "rotation_count": results["invsqrt"]["rotation_count"],
                "rotation_eps": results["invsqrt"]["rotation_eps"],
                "verified_shots": results["invsqrt"]["verified_shots"],
                "worst_error_over_tolerance": results["invsqrt"]["worst_error_over_tolerance"],
            },
            "dqi": {
                "qubits": results["dqi"]["num_qubits"],
                "toffoli_count": results["dqi"]["toffoli_count"],
                "repetition_factor": results["dqi"]["repetition_factor"],
                "effective_error_rate": results["dqi"]["effective_error_rate"],
                "verified_shots": results["dqi"]["verified_shots"],
            },
            "dfthcsos": {
                "qubits": results["dfthcsos"]["num_qubits"],
                "toffoli_count": results["dfthcsos"]["toffoli_count"],
                "lambda_eff": results["dfthcsos"]["lambda_eff"],
                "verified_shots": results["dfthcsos"]["verified_shots"],
                "metrics": results["dfthcsos"].get("metrics", {}),
            }
        }
    }

    out_file = AGY_DIR / "logs/best_submission_score.json"
    out_file.write_text(json.dumps(summary, indent=2) + "\n")

    print("\n" + "=" * 60)
    print("SUBMISSION VERIFICATION & SCORING COMPLETE")
    print("=" * 60)
    print(f"Archive: {zip_path} ({zip_size:,} bytes)")
    print(f"SHA-256: {zip_sha}")
    print(f"\nTask Results:")
    print(f"  1. InvSqrt:  cost = {inv_cost:16,.2f} | speedup = {inv_ratio:8.4f}x (baseline: {BASELINES['invsqrt']:16,.2f})")
    print(f"  2. DQI:      cost = {dqi_cost:16,.2f} | speedup = {dqi_ratio:8.4f}x (baseline: {BASELINES['dqi']:16,.2f})")
    print(f"  3. Chemistry:cost = {chem_cost:16,.2f} | speedup = {chem_ratio:8.4f}x (baseline: {BASELINES['dfthcsos']:16,.2f})")
    print("-" * 60)
    print(f"OVERALL GEOMETRIC-MEAN SCORE: {overall_score:.6f}x")
    print("=" * 60)

if __name__ == "__main__":
    main()
