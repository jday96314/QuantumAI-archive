"""Fast offset optimizer using Scipy L-BFGS-B with analytical gradient."""

import sys, time, json, argparse
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize

p = argparse.ArgumentParser()
p.add_argument("--fit", required=True)
p.add_argument("--out", required=True)
p.add_argument("--reference", default="dfthcsos/submission.npz")
a = p.parse_args()

HERE = Path(__file__).resolve().parent
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
ref = np.load(ROOT / "initial_experiments/3/parameters/chem_reference.npz")

z = dict(np.load(a.fit))
u = z["U"][0] # (B, 76)
w = z["w"][0] # (B, C)
B, C = w.shape
n = 76

X = np.einsum("bp,bq,bc->cpq", u, u, w, optimize=True) # (C, 76, 76)
xf = X.reshape(C, -1)
tr = np.trace(X, axis1=1, axis2=2)

h0 = ref["h1_thc"] - 0.5 * int(ref["electrons"]) * ref["S"] + float(ref["beta_1"]) * np.eye(n)
constant = -22140.27335 - float(ref["ecore"]) + int(ref["electrons"]) * float(ref["beta_1"])

def fun(o, return_grad=True):
    hh = h0 - (o @ xf).reshape(n, n)
    ev, ve = eigh(0.5 * (hh + hh.T))
    r = np.sqrt(ev * ev + 1e-14)
    gap = constant + r.sum() - ev.sum() + 0.5 * np.sum(o * o)
    if return_grad:
        grad = o + tr - xf @ ((ve * (ev / r)) @ ve.T).reshape(-1)
        return gap, grad
    lam_h1 = r.sum()
    wb = o + w.sum(0)
    lam_h2 = 0.25 * np.sum((np.abs(w).sum(0) + np.abs(wb))**2)
    Lambda = lam_h1 + lam_h2
    return gap, Lambda

# Initialize from reference submission W if available, or zero
if Path(a.reference).exists():
    ref_sub = np.load(a.reference)
    W_init_full = ref_sub["H2_I"][0] - ref_sub["W"][0].sum(0)
    W_init = W_init_full[:C] if len(W_init_full) >= C else np.pad(W_init_full, (0, C - len(W_init_full)))
else:
    W_init = np.zeros(C)

t0 = time.time()
print(f"Optimizing offset W for B={B}, C={C} using Scipy L-BFGS-B...")
res = minimize(fun, W_init, jac=True, method="L-BFGS-B", options=dict(maxiter=200, ftol=1e-12, gtol=1e-7))
W_opt = res.x

gap, Lambda = fun(W_opt, return_grad=False)
rad = max(0.0, gap * (2.0 * Lambda - gap))
lam_eff = np.sqrt(rad)
print(f"Optimization finished in {time.time()-t0:.2f}s:")
print(f"  gap = {gap:.6f} Ha, Lambda = {Lambda:.2f}, lambda_eff = {lam_eff:.4f}")

z["W"] = W_opt[None, :]
np.savez(a.out, **z)
print(f"Saved optimized offset to {a.out}")
