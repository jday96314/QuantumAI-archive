"""Build and verify a candidate factor checkpoint end-to-end."""

import sys, time, json, os, subprocess
from pathlib import Path
import numpy as np

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
PY = "/mnt/data01/miniconda3/envs/torch2.13/bin/python"
ENV = dict(os.environ, LD_LIBRARY_PATH="/mnt/data01/miniconda3/envs/torch2.13/lib",
           PYTHONPATH=f"{ROOT}/from_host/kickmix:{ROOT}/from_host/verifiers:{ROOT}/initial_experiments/3:{ROOT}/initial_experiments/2")

def evaluate_checkpoint(factor_path, work_dir="runs/test_eval"):
    wdir = Path(work_dir)
    wdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    
    # 1. Quantize One-Body Angles
    f_onebody = wdir / "onebody.npz"
    print(f"\n>>> Running quantize_onebody on {factor_path}...")
    cmd1 = f"{PY} ../../3/quantize_onebody.py --fit {factor_path} --out {f_onebody} --target 0.012 --steps 1000"
    res1 = subprocess.run(cmd1, shell=True, env=ENV, capture_output=True, text=True)
    if res1.returncode != 0:
        print("quantize_onebody failed:\n", res1.stderr)
        return None
    print("quantize_onebody done in", time.time() - t0, "s")
    
    # Fix brot and angle_words shape for build_global_chem
    z = dict(np.load(f_onebody))
    words = z["angle_words"].reshape(-1, 75)
    if len(words) == 795:
        z["angle_words"] = words[76:]
    z["brot"] = 8
    np.savez(f_onebody, **z)
    
    # 2. Build Circuits
    sub_dir = wdir / "submission"
    print(f"\n>>> Building circuits in {sub_dir}...")
    cmd2 = f"{PY} ../../3/build_global_chem.py --fit {f_onebody} --out {sub_dir} --batch 13"
    res2 = subprocess.run(cmd2, shell=True, env=ENV, capture_output=True, text=True)
    if res2.returncode != 0:
        print("build_global_chem failed:\n", res2.stderr)
        return None
    print("build_global_chem done!")
    
    # 3. Host Verification
    print(f"\n>>> Running Host Verification...")
    f_verify = wdir / "verify.json"
    cmd3 = f"{PY} {ROOT}/from_host/verifiers/dfthcsos/verify.py {sub_dir} > {f_verify}"
    res3 = subprocess.run(cmd3, shell=True, env=ENV)
    
    with open(f_verify) as fp:
        res = json.load(fp)
        
    print("\n" + "="*50)
    print(f"VERIFICATION RESULT FOR {factor_path}:")
    print(f"Status:   {res.get('status')}")
    print(f"Score:    {res.get('score'):,.2f}")
    print(f"Toffoli:  {res.get('toffoli_count')}")
    print(f"Qubits:   {res.get('num_qubits')}")
    print(f"lambda:   {res.get('lambda_eff'):.6f}")
    print(f"error_1:  {res.get('metrics', {}).get('error_1')*1000:.3f} mHa")
    print(f"error_2:  {res.get('metrics', {}).get('error_2')*1000:.3f} mHa")
    print(f"E_gap:    {res.get('metrics', {}).get('E_gap'):.6f} Ha")
    print(f"Incumbent: 757,163,489.50")
    if res.get('score', 1e12) < 757163489.50:
        diff = 757163489.50 - res.get('score')
        print(f"*** BEAT INCUMBENT BY {diff:,.2f} POINTS! ***")
    print("="*50)
    return res

if __name__ == "__main__":
    factor = sys.argv[1] if len(sys.argv) > 1 else "runs/c347/factor_16k.npz"
    work = sys.argv[2] if len(sys.argv) > 2 else "runs/c347_eval"
    evaluate_checkpoint(factor, work)
