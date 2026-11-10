"""Step-by-step debug of C=24 evaluation."""

import sys
sys.path.insert(0, "/media/shared/ExpansionDrive2/QuantumAI/from_host/verifiers")
import numpy as np
from dfthcsos.verify import load_inputs, evaluate_words

def main():
    inputs = load_inputs("/tmp/c24_submission.npz")
    sub = np.load("/media/shared/ExpansionDrive2/QuantumAI/from_host/sample_submission/dfthcsos/submission.npz")
    tree_edges = np.asarray(sub["tree_edges"], dtype=np.int64)
    base_words = np.asarray(sub["angle_words"], dtype=np.int64)

    print(f"Inputs w shape: {inputs.w.shape}, W shape: {inputs.W.shape}")
    res = evaluate_words(inputs, base_words, tree_edges, 12)
    print("Results with base_words:")
    for k, v in res.items():
        print(f"  {k}: {v}")

if __name__ == "__main__":
    main()
