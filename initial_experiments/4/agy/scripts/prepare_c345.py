"""Prepare C=345 candidate up to nudged.npz."""

import os, sys, time, json, subprocess
from pathlib import Path
import numpy as np

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
PY = "/mnt/data01/miniconda3/envs/torch2.13/bin/python"
ENV = dict(os.environ, LD_LIBRARY_PATH="/mnt/data01/miniconda3/envs/torch2.13/lib",
           PYTHONPATH=f"{ROOT}/from_host/kickmix:{ROOT}/from_host/verifiers:{ROOT}/initial_experiments/3:{ROOT}/initial_experiments/2")

def run(cmd):
    print(f"\n>>> Running: {cmd}")
    t0 = time.time()
    res = subprocess.run(cmd, shell=True, env=ENV, capture_output=True, text=True)
    dt = time.time() - t0
    if res.returncode != 0:
        print(f"FAILED (code {res.returncode}) in {dt:.1f}s:\n{res.stderr}\n{res.stdout}")
        raise RuntimeError(f"Command failed: {cmd}")
    print(f"Done in {dt:.1f}s:\n{res.stdout.strip()}")
    return res.stdout

C = 345
rdir = Path("runs/c345")
rdir.mkdir(parents=True, exist_ok=True)

# 1. Project Factor
f_proj = rdir / "proj.npz"
run(f"{PY} ../../3/project_factor.py --fit dfthcsos/submission.npz --out {f_proj} --columns {C}")

# 2. Fast Optimize Offsets
f_off = rdir / "off.npz"
run(f"{PY} scripts/fast_optimize_offsets.py --fit {f_proj} --out {f_off}")

# 3. Align Columns
f_aligned = rdir / "aligned.npz"
run(f"{PY} ../../3/align_columns.py --fit {f_off} --reference ../../3/candidates/mixing_mu1e7.npz --out {f_aligned} --physical")

# 4. Exact Jacobi Rotation
f_rot_in = rdir / "rot_in.bin"
f_rot_out = rdir / "rot_out.bin"
z = dict(np.load(f_aligned))
w = z["w"][0]
B = w.shape[0]
w_cm = np.ascontiguousarray(w.T, dtype=np.float64)
w_cm.tofile(f_rot_in)

print(f"\n>>> Running rotate_exact on w shape ({C}, {B})...")
run(f"../../3/rotate_exact {f_rot_in} {f_rot_out} {C} {B} 300000")

w_rot_cm = np.fromfile(f_rot_out, dtype=np.float64).reshape(C, B)
w_rot = np.ascontiguousarray(w_rot_cm.T)
z["w"] = w_rot[None, :, :]
f_rotated = rdir / "rotated.npz"
np.savez(f_rotated, **z)

# 5. Nudge Gap
f_nudged = rdir / "nudged.npz"
run(f"{PY} ../../3/nudge_gap.py --fit {f_rotated} --out {f_nudged} --target 3.9152 --cap 0.01999")
print("\n*** C=345 Preparation Complete! ***")
