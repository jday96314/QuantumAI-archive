#!/usr/bin/env python3
import sys
from pathlib import Path
import json

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
]

from dfthcsos.verify import verify

candidates = [
    ROOT / "initial_experiments/3/candidates/grow511_q11_c350_submission",
    ROOT / "initial_experiments/3/candidates/grow725_q8_c350_submission",
    ROOT / "initial_experiments/3/candidates/residual719_q8_c350_submission",
    ROOT / "initial_experiments/3/candidates/global500_coord12_c350_qr2_submission",
    ROOT / "initial_experiments/3/candidates/global500_final_q13_c350_submission",
]

for cdir in candidates:
    if not cdir.exists():
        continue
    print(f"\n--- Checking {cdir.name} ---")
    try:
        res = verify(cdir)
        passed = res.get("passed")
        cost = res.get("cost")
        reason = res.get("reason")
        print(f"Passed: {passed} | Cost: {cost} | Reason: {reason}")
        if "metrics" in res:
            m = res["metrics"]
            print(f"  err2: {m.get('error_2', 0)*1000:.3f} mHa, gap: {m.get('E_gap')}, delta_e: {m.get('total_delta_e')}")
    except Exception as e:
        import traceback
        traceback.print_exc()
