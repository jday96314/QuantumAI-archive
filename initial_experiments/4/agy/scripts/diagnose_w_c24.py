"""Diagnose w_c24 vs w_base."""

import numpy as np

def main():
    c24 = np.load("/tmp/c24_submission.npz")
    base = np.load("/media/shared/ExpansionDrive2/QuantumAI/from_host/sample_submission/dfthcsos/submission.npz")
    w_c24 = np.asarray(c24["W"])
    w_base = np.asarray(base["W"])

    print("w_c24 shape:", w_c24.shape, "norm:", np.linalg.norm(w_c24))
    print("w_base shape:", w_base.shape, "norm:", np.linalg.norm(w_base))
    
    # What was the initial rotated w before optimization?
    # Let's check how much w_c24 differs from w_base:
    print("Max w_c24:", np.max(np.abs(w_c24)))
    print("Max w_base:", np.max(np.abs(w_base)))

if __name__ == "__main__":
    main()
