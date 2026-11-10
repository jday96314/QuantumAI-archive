"""Count the existing gap schedule with wave-three trimmed middle products.

This is a resource probe, not native acceptance. Only integration/wave3/results
is written; the original decoder and accepted artifacts remain read-only.
"""

import importlib.util
import json
import os
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROUND5 = HERE.parent
GAP = ROUND5 / "dqi-gap"
sys.path.insert(0, str(GAP))
# The inherited classical field reference reads its tables during import.
os.chdir(GAP)
import hybrid_decoder as decoder
from blocked_monic_euclid import schedule
from profile_hybrid import module as count_module


def frozen_frames(words, width):
    counts = json.loads((GAP / "parameters/best_schedule.json").read_text())["counts"]
    previous, current, tail, pivot = [], [words[-2]], words[:-2], words[-1]
    k, stages = 1, []
    for count in counts:
        stages.append((k, count, previous, current, tail, pivot))
        _, (previous, current) = schedule(previous, current, list(reversed(tail))[:2 * count])
        tail = tail[:-2 * count]
        k += count
    assert k == 35 and not tail
    return stages, (previous, current, tail, pivot)


def main():
    destination = HERE / "wave3"
    (destination / "results").mkdir(parents=True, exist_ok=True)
    os.chdir(destination)
    path = ROUND5 / "dqi-core/components/trimmed_middle.py"
    spec = importlib.util.spec_from_file_location("wave3_trimmed_middle", path)
    trimmed = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(trimmed)
    decoder.BoundedArithmetic = trimmed.TrimmedArithmetic
    decoder.load = count_module
    decoder.frames = frozen_frames
    decoder.build.profile_only = True
    report = decoder.build(70, 3, 100, 4095, Path("candidates/gap_trimmed_profile.kmb"),
                           linear_padding=True, virtual_zero=True, virtual_lead=True,
                           low_inverse=True, lazy_flags=True)
    original = json.loads((GAP / "parameters/best_candidate_manifest.json").read_text())
    original_core = 824740
    report.update(scope="resource-only; no native acceptance", original_core_T=original_core,
                  original_Q=original["Q"], saved_core_T=original_core - report["T"])
    Path("results/gap_trimmed_profile.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: value for key, value in report.items() if key != "rows"}), flush=True)


if __name__ == "__main__":
    main()
