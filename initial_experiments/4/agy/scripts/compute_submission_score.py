"""Compute overall score from verification JSONs of the three benchmarks."""

import json
import math

BASELINES = {
    "invsqrt": 2076774.233715282,
    "dqi": 24439463583672.93,
    "dfthcsos": 3355599588.3155375,
}

def main():
    try:
        with open("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy/logs/verify_invsqrt.json") as f:
            v_invsqrt = json.load(f)
    except Exception as e:
        v_invsqrt = {"status": "FAIL", "error": str(e)}

    try:
        with open("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy/logs/verify_dqi_full.json") as f:
            v_dqi = json.load(f)
    except Exception as e:
        v_dqi = {"status": "FAIL", "error": str(e)}

    try:
        with open("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy/logs/verify_dfthcsos.json") as f:
            v_dfthc = json.load(f)
    except Exception as e:
        v_dfthc = {"status": "FAIL", "error": str(e)}

    print("=== Verification Summary ===")
    ratios = {}
    for name, res in [("invsqrt", v_invsqrt), ("dqi", v_dqi), ("dfthcsos", v_dfthc)]:
        status = res.get("status", "UNKNOWN")
        score = res.get("score")
        base = BASELINES[name]
        if status == "PASS" and score:
            ratio = base / score
            ratios[name] = ratio
            print(f"{name:10s}: PASS | Cost: {score:,.2f} | Baseline: {base:,.2f} | Improvement: {ratio:.4f}x")
        else:
            print(f"{name:10s}: {status} | Error: {res.get('error', 'N/A')}")

    if len(ratios) == 3:
        geom_mean = (ratios["invsqrt"] * ratios["dqi"] * ratios["dfthcsos"]) ** (1.0 / 3.0)
        print(f"\n>>> OVERALL GEOMETRIC MEAN SCORE: {geom_mean:.4f}x <<<")
    else:
        print(f"\nMissing some passing benchmarks ({len(ratios)}/3 available).")

if __name__ == "__main__":
    main()
