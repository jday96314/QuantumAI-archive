"""Copy /tmp/c24_submission.npz to dfthcsos/ and verify with verify.py."""

import json
import shutil
from pathlib import Path
import sys
sys.path.insert(0, "/media/shared/ExpansionDrive2/QuantumAI/from_host/verifiers")
from dfthcsos.verify import verify

def main():
    src_npz = Path("/tmp/c24_submission.npz")
    if not src_npz.exists():
        raise FileNotFoundError(f"{src_npz} not found.")

    out_dir = Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy/dfthcsos")
    out_dir.mkdir(parents=True, exist_ok=True)
    dst_npz = out_dir / "submission.npz"
    shutil.copy(src_npz, dst_npz)
    print(f"Copied {src_npz} -> {dst_npz}")

    # Ensure angle circuits exist in out_dir
    for k in range(10):
        src = f"/media/shared/ExpansionDrive2/QuantumAI/from_host/sample_submission/dfthcsos/angle_k_{k}.kmb"
        dst = out_dir / f"angle_k_{k}.kmb"
        if not dst.exists():
            shutil.copy(src, dst)

    print("Running official verify(dfthcsos/, format='kmb')...")
    res = verify(out_dir, format="kmb")
    print("Verification result:")
    for k, v in res.items():
        if k not in ("resources", "metrics"):
            print(f"  {k}: {v}")

    base_score = 3355599588.3155375
    new_score = res["score"]
    ratio = base_score / new_score
    print(f"\nDFTHC-SOS Improvement Ratio: {ratio:.4f}x (New Cost: {new_score:,.2f}, Base: {base_score:,.2f})")

    log_path = "/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy/logs/verify_dfthcsos.json"
    with open(log_path, "w") as f:
        json.dump(res, f, indent=2)
    print(f"Saved {log_path}")

if __name__ == "__main__":
    main()
