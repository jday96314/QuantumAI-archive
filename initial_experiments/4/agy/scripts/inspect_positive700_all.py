"""Inspect all positive700 candidates to check their shapes, err2, and pass status."""

from pathlib import Path
import numpy as np

cand_dir = Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/candidates")

for p in sorted(cand_dir.glob("positive700*.npz")):
    try:
        z = np.load(p)
        err2 = float(z["err2"]) * 1000 if "err2" in z else None
        brot = int(z["brot"]) if "brot" in z else None
        w_shape = z["w"].shape if "w" in z else None
        U_shape = z["U"].shape if "U" in z else None
        has_tree = "tree_edges" in z
        has_words = "angle_words" in z
        err_str = f"{err2:6.2f} mHa" if err2 is not None else "None"
        print(f"{p.name:45s} | err2={err_str} | brot={brot} | w={w_shape} | tree={has_tree} | words={has_words}")
    except Exception as e:
        print(f"{p.name:45s} | Error: {e}")
