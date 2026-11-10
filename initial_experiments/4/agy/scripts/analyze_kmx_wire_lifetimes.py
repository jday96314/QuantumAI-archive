#!/usr/bin/env python3
import re
from pathlib import Path

kmx_path = Path("/tmp/pub_temp.kmx")
if not kmx_path.exists():
    print("Not found:", kmx_path)
    exit(1)

q_pattern = re.compile(r'\bq(\d+)\b')

born = {}
died = {}
total_lines = 0

with open(kmx_path) as f:
    for line_idx, line in enumerate(f):
        total_lines += 1
        if line.startswith("INFO_TAG") or line.startswith("REGISTER") or line.startswith("APPEND_TO_REGISTER"):
            continue
        for m in q_pattern.finditer(line):
            qid = int(m.group(1))
            if qid not in born:
                born[qid] = line_idx
            died[qid] = line_idx

print(f"Total lines: {total_lines}")
print(f"Total unique wires: {len(born)}")
pinned = 840
ancilla_wires = [q for q in born if q >= pinned]
print(f"Total ancilla wires (>= 840): {len(ancilla_wires)}")

# Check how many ancilla wires die near the end vs early
near_end_count = sum(1 for q in ancilla_wires if died[q] > total_lines - 100000)
print(f"Ancillas dying in last 100k lines: {near_end_count} / {len(ancilla_wires)}")

# Show distribution of died lines for ancillas
died_percentiles = [died[q] / total_lines for q in ancilla_wires]
import numpy as np
pcts = [0, 10, 25, 50, 75, 90, 95, 99, 100]
print("Percentile of wire death in circuit:")
for p, v in zip(pcts, np.percentile(died_percentiles, pcts)):
    print(f"  {p:3d}%: {v*100:5.1f}% of circuit")
