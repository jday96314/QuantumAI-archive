# Follow-up: ANGLE pruning probe, then package a submission

1. **Probe (a few hours):** for the banked candidate's 11 oracles, measure how much lookup work could be skipped if aligned address subtrees with identical XOR-difference words were collapsed. Measure under the current address order and under optimized address permutations (permuting dictionary rows together with their weights is physically free). Report the projected Toffoli saving per oracle and in total. Build and native-verify a pruned family only if it projects ≥3% of chemistry T.
2. **Package:** assemble a complete archive from the best accepted artifacts:
   - chemistry: `candidates/banked`, or a better verified result from step 1;
   - inverse: the round-4 Codex selection;
   - DQI: the round-4 Codex selection, unless `initial_experiments/5/dqi-core/RESULTS.md` says "ACCEPTED: trimmed decoder" with full native runs.
   Run the metric's packaging and validation checks (`from_host/metric.py`) and record the local combined score. Write the archive here.
GPU free; CPU cap 2.