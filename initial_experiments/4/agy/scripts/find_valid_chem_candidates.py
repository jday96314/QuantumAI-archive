import os
import glob
import numpy as np

candidates = glob.glob("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/candidates/*.npz")

valid = []
for p in candidates:
    try:
        z = np.load(p)
        err2 = None
        if "err2" in z:
            err2 = float(z["err2"])
        elif "metrics" in z and "error_2" in z["metrics"]:
            err2 = float(z["metrics"]["error_2"])
        
        # Check shape of U and w
        U_shape = z["U"].shape if "U" in z else None
        w_shape = z["w"].shape if "w" in z else None
        
        if err2 is not None and err2 <= 0.0205: # Within tolerance or close
            valid.append((err2, p, U_shape, w_shape, z.get("brot")))
    except Exception:
        pass

valid.sort(key=lambda x: x[0])
print(f"Found {len(valid)} valid or nearly valid candidates (err2 <= 20.5 mHa):")
for err2, p, U_shape, w_shape, brot in valid[:25]:
    fname = os.path.basename(p)
    print(f"err2: {err2*1000:6.3f} mHa | U: {str(U_shape):15s} | w: {str(w_shape):15s} | brot: {brot} | {fname}")
