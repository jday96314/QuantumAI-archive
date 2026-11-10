import os
import glob
import json

base_dirs = [
    "/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/1",
    "/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/2",
    "/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3",
    "/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy",
]

baselines = {
    "invsqrt": 2076774.2337,
    "dqi": 24439463583672.93,
    "dfthcsos": 3355599588.3155,
}

best_by_task = {
    "invsqrt": [],
    "dqi": [],
    "dfthcsos": [],
}

for bdir in base_dirs:
    for root, dirs, files in os.walk(bdir):
        for f in files:
            if f.endswith(".json") and ("verify" in f or "score" in f or "summary" in f):
                p = os.path.join(root, f)
                try:
                    with open(p, "r") as fp:
                        data = json.load(fp)
                except Exception:
                    continue
                
                # Check format
                # Format 1: host verifier direct output
                # {"verdict": "PASS", "task": ..., "cost": ...}
                # Format 2: summary output
                # {"task_scores": {"invsqrt": ..., ...}}
                if isinstance(data, dict):
                    if data.get("verdict") == "PASS":
                        task = data.get("task")
                        cost = data.get("cost")
                        if task in best_by_task and cost is not None:
                            best_by_task[task].append((cost, p, data))
                    elif "task_scores" in data:
                        for task, cost in data["task_scores"].items():
                            if task in best_by_task and cost is not None:
                                best_by_task[task].append((cost, p, data))

for task in ["invsqrt", "dqi", "dfthcsos"]:
    print(f"=== TASK: {task} (Baseline: {baselines[task]}) ===")
    entries = sorted(best_by_task[task], key=lambda x: x[0])
    seen = set()
    count = 0
    for cost, path, data in entries:
        ratio = baselines[task] / cost
        cost_key = round(cost, 2)
        if cost_key not in seen:
            seen.add(cost_key)
            print(f"Cost: {cost:18.4f} | Ratio: {ratio:8.4f}x | Path: {path}")
            count += 1
            if count >= 5:
                break
    print()
