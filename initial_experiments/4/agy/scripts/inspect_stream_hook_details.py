"""Inspect keys and register mappings in stream10_offset3_select_noscale_compact_hook.json."""

from pathlib import Path
import json

hook_path = Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/parameters/stream10_offset3_select_noscale_compact_hook.json")
data = json.loads(hook_path.read_text())

print("Hook keys:", list(data.keys()))
print("maxq:", data.get("maxq"))
print("index:", data.get("index"))
print("outcomes:", data.get("outcomes"))

if "sigma" in data:
    sig = data["sigma"]
    print(f"sigma: {len(sig)} words, first word: {sig[0]}, last word: {sig[-1]}")
    all_sig = [x for w in sig for x in w]
    print(f"sigma min wire: {min(all_sig)}, max wire: {max(all_sig)}")

if "free" in data:
    fr = data["free"]
    print(f"free: {len(fr)} wires, min wire: {min(fr) if fr else None}, max wire: {max(fr) if fr else None}")

if "phases" in data:
    ph = data["phases"]
    print(f"phases: {len(ph)} phases")
