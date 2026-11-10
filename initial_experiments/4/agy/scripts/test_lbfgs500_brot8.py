#!/usr/bin/env python3
import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
    str(ROOT / "initial_experiments/3"),
]
import numpy as np
from compress_chem import encode
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors, evaluate_dfthc_vectors
from dfthcsos.verify import load_inputs, evaluate_words

sub_npz = ROOT / "initial_experiments/3/candidates/lbfgs500_q10c375_submission/submission.npz"
z = dict(np.load(sub_npz))
U = z["U"]
W = z["W"]
H2_I = z["H2_I"]
R, B, N = U.shape
C = W.shape[2]
tree = z["tree_edges"]

print(f"lbfgs500: R={R}, B={B}, C={C}")

fact = dict(U=U, W=W, H2_I=H2_I, H1_I=z["H1_I"], Hsym_lower=z["Hsym_lower"])
inputs = load_inputs(fact)

# Original q10
m10 = evaluate_words(inputs, z["angle_words"], tree, brot=10)
print(f"q10 original -> err1: {m10['error_1']*1000:.3f} mHa, err2: {m10['error_2']*1000:.3f} mHa, gap: {m10['E_gap']:.3f}")

# Re-encode at brot=8
words_8 = encode(z["angle_words"][:, :], tree, 8) if False else None
# Encode directly from U
u_flat = U.reshape(-1, N)
two_8 = encode(u_flat, tree, 8)
# For one-body:
one_8 = encode(reconstruct_u_vectors(z["angle_words"][:76] / (1 << 9), tree), tree, 8)
words_8 = np.r_[one_8, two_8]

m8 = evaluate_words(inputs, words_8, tree, brot=8)
print(f"q8 direct   -> err1: {m8['error_1']*1000:.3f} mHa, err2: {m8['error_2']*1000:.3f} mHa, gap: {m8['E_gap']:.3f}")
