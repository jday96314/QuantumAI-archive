"""Find all verified components and results across initial_experiments to identify the best-scoring components thus far."""

from pathlib import Path
import json

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
BASELINES = {
    "invsqrt": 2076774.233715282,
    "dqi": 24439463583672.93,
    "dfthcsos": 3355599588.3155375,
}

def scan_results():
    all_results = {"invsqrt": [], "dqi": [], "dfthcsos": []}
    
    # Scan for json files containing 'score' or 'status'
    for exp_dir in [ROOT / "initial_experiments/4/agy", ROOT / "initial_experiments/3", ROOT / "initial_experiments/2", ROOT / "initial_experiments/1"]:
        if not exp_dir.exists():
            continue
        for p in exp_dir.rglob("*.json"):
            if "package" in p.name or "node_modules" in str(p) or ".agent-home" in str(p):
                continue
            try:
                data = json.loads(p.read_text())
                if isinstance(data, dict):
                    status = data.get("status")
                    score = data.get("score")
                    if status == "PASS" and isinstance(score, (int, float)):
                        # Detect task
                        if "n_p" in data or "invsqrt" in p.name.lower() or "qubits" in data and "rotations" in data:
                            all_results["invsqrt"].append((score, str(p), data))
                        elif "shots" in data and data.get("shots") == 102400 or "dqi" in p.name.lower() or "stream" in p.name.lower():
                            all_results["dqi"].append((score, str(p), data))
                        elif "shots" in data and data.get("shots") == 50880 or "chem" in p.name.lower() or "dfthcsos" in p.name.lower() or "lambda_eff" in data:
                            all_results["dfthcsos"].append((score, str(p), data))
                    # Also check submission summary format
                    if "task_scores" in data:
                        for k, v in data["task_scores"].items():
                            if k in all_results and isinstance(v, (int, float)):
                                all_results[k].append((v, str(p), data))
            except Exception:
                pass

    print("=== SUMMARY OF BEST VERIFIED COMPONENTS FOUND ===")
    for task, items in all_results.items():
        if not items:
            print(f"Task {task}: No verified results found")
            continue
        # Sort by score ascending (lower is better)
        items.sort(key=lambda x: x[0])
        print(f"\nTask {task} (Baseline: {BASELINES[task]:,.2f}):")
        seen = set()
        count = 0
        for score, path, data in items:
            key = round(score, 2)
            if key in seen:
                continue
            seen.add(key)
            ratio = BASELINES[task] / score
            print(f"  Score: {score:,.4f} ({ratio:.4f}x) from {path}")
            count += 1
            if count >= 5:
                break

if __name__ == "__main__":
    scan_results()
