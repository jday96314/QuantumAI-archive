"""Automated pipeline to produce a fully verified chemistry submission for target C."""

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

def run_pipeline(C=347, run_dir="runs/c347"):
    rdir = Path(run_dir)
    rdir.mkdir(parents=True, exist_ok=True)
    t_start = time.time()
    
    print(f"==================================================")
    print(f"Starting Chemistry Optimization Pipeline for C={C}")
    print(f"Working Directory: {rdir}")
    print(f"==================================================")
    
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
    w = z["w"][0] # (B, C)
    B = w.shape[0]
    # rotate_exact expects column-major order: A[c * B + k]
    w_cm = np.ascontiguousarray(w.T, dtype=np.float64) # (C, B)
    w_cm.tofile(f_rot_in)
    
    print(f"\n>>> Running rotate_exact on w shape ({C}, {B})...")
    run(f"../../3/rotate_exact {f_rot_in} {f_rot_out} {C} {B} 300000")
    
    w_rot_cm = np.fromfile(f_rot_out, dtype=np.float64).reshape(C, B)
    w_rot = np.ascontiguousarray(w_rot_cm.T) # (B, C)
    z["w"] = w_rot[None, :, :]
    f_rotated = rdir / "rotated.npz"
    np.savez(f_rotated, **z)
    
    # 5. Nudge Gap
    f_nudged = rdir / "nudged.npz"
    run(f"{PY} ../../3/nudge_gap.py --fit {f_rotated} --out {f_nudged} --target 3.9152 --cap 0.01999")
    
    # 6. Constrained Factor (fine-tuning)
    f_factor = rdir / "factor.npz"
    run(f"{PY} ../../3/constrained_factor.py --fit {f_nudged} --out {f_factor} --steps 4000 --lr 2e-6 --cap 0.01999 --whitened")
    
    # 7. Quantize One-Body Angles
    f_onebody = rdir / "onebody.npz"
    run(f"{PY} ../../3/quantize_onebody.py --fit {f_factor} --out {f_onebody} --target 0.012 --steps 1500")
    
    # 8. Build Circuits
    sub_dir = rdir / "submission"
    run(f"{PY} ../../3/build_global_chem.py --fit {f_onebody} --out {sub_dir} --batch 13")
    
    # 9. Verify on Host
    print(f"\n>>> Running Host Verification on {sub_dir}...")
    f_verify = rdir / "verify.json"
    run(f"{PY} {ROOT}/from_host/verifiers/dfthcsos/verify.py {sub_dir} --out {f_verify}")
    
    with open(f_verify) as fp:
        v_res = json.load(fp)
        
    print("\n" + "="*50)
    print("VERIFICATION RESULT:")
    print(f"Status: {v_res.get('status')}")
    print(f"Cost:   {v_res.get('score'):,.2f}")
    print(f"Toffoli:{v_res.get('toffoli_count')}")
    print(f"Qubits: {v_res.get('num_qubits')}")
    print(f"lambda: {v_res.get('lambda_eff'):.6f}")
    print(f"err_1:  {v_res.get('metrics', {}).get('error_1')*1000:.3f} mHa")
    print(f"err_2:  {v_res.get('metrics', {}).get('error_2')*1000:.3f} mHa")
    print(f"E_gap:  {v_res.get('metrics', {}).get('E_gap'):.6f} Ha")
    print(f"Total pipeline time: {time.time()-t_start:.1f}s")
    print("="*50)
    return v_res

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--C", type=int, default=347)
    ap.add_argument("--dir", default="runs/c347")
    a = ap.parse_args()
    run_pipeline(a.C, a.dir)
