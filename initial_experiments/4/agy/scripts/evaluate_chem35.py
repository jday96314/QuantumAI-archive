"""Evaluate chem35c45_rot.npz using the host verifier evaluate_words."""

import sys
from pathlib import Path
import numpy as np

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
]
from dfthcsos.verify import load_inputs, evaluate_words, heuristic_costs

cand_path = ROOT / "initial_experiments/3/candidates/chem35c45_rot.npz"
z = dict(np.load(cand_path))
U = z["U"]
W = z["W"] # here W is w (R, B, C)
H2_I = z["H2_I"]
R, B, N = U.shape
C = W.shape[2]
tree = z["tree_edges"]
words = z["angle_words"]

fact = dict(U=U, W=W, H2_I=H2_I, H1_I=z["H1_I"], Hsym_lower=z["Hsym_lower"])
inputs = load_inputs(fact)

metrics = evaluate_words(inputs, words, tree, brot=12)
err1 = metrics["error_1"] * 1000
err2 = metrics["error_2"] * 1000
gap = metrics["E_gap"]
lam_eff = metrics["lambda_eff"]
print(f"chem35c45_rot (brot=12):")
print(f"  err1: {err1:.3f} mHa, err2: {err2:.3f} mHa, gap: {gap:.3f} Ha, lam_eff: {lam_eff:.3f}")
print(f"  gap_in_domain: {metrics['gap_in_domain']}")

base = 3355599588.3155375
h = heuristic_costs(N, R, B, C, 12)
J = 10
K = 8
A = N + R * B
sel = (A - 1).bit_length()
ang_toff = 2 * K * (A - 2)
T = ang_toff + h["toffoli"]
Q = h["persistent_qubits"] + (J * 12) + max(h["temporary_qubits"], sel - 1)
cost = lam_eff * T * Q
ratio = base / cost
print(f"  T={T:,}, Q={Q:,} -> Cost: {cost:,.2f} ({ratio:.4f}x)")
