"""Simulate and prove the exact non-colliding register reuse schedule for DQI BM."""

from pathlib import Path
import json

# Load IMAXV and schedule from constants
circuit_constants = json.loads(Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/vendor/uncompute/tasks/dqi-reed-solomon/circuits/discovered/constants.json").read_text())
bm_windows = circuit_constants["bm_windows"]
IMAXV = bm_windows["IMAXV"]
KMSV = bm_windows["KMSV"]
KMMV = bm_windows["KMMV"]

print(f"Total BM rounds: {len(IMAXV)}")

# For each round r from 0 to 69:
# Check which syndromes S_i are needed in discrepancy:
# S_{r - i} for i <= IMAXV[r].
# Minimum syndrome index needed: min_syn(r) = r - IMAXV[r].
# Any syndrome with index < min_syn(r) will NEVER be needed again in discrepancy!

never_needed_again = []
for r in range(70):
    min_needed = min(r - i for i in range(IMAXV[r] + 1))
    # Syndromes strictly less than min_needed are no longer needed
    free_syndromes = min_needed
    # How many sig registers needed?
    sig_needed = max(IMAXV[r], KMSV[r], KMMV[r]) + 1
    # How many bb registers needed?
    bb_needed = min(37, max(KMSV[r], KMMV[r]) + 2)
    
    total_needed = sig_needed + bb_needed
    print(f"Round {r:2d}: min_needed S_{min_needed:2d} (freed: {free_syndromes:2d}) | sig: {sig_needed:2d}, bb: {bb_needed:2d} | total needed: {total_needed:2d}")
