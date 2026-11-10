"""Optimize 1-body offset vector W alone for residual719_q8_c340.npz to achieve gap >= 3.9151 Ha and minimize lambda_eff."""

import sys, time
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
    str(ROOT / "initial_experiments/3"),
]
import numpy as np
import torch
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
from dfthcsos.verify import load_inputs, evaluate_words

device = torch.device("cuda")

z = np.load("../../3/candidates/residual719_q8_c340.npz")
U_raw = z["U"]
w_raw = z["w"]
R, B, N = U_raw.shape
C = w_raw.shape[-1]
tree = z["tree_edges"]
brot = int(z["brot"])
words2 = z["angle_words"].reshape(-1, 75)

ref_dict = np.load(ROOT / "initial_experiments/3/parameters/chem_reference.npz")
S_ref = np.asarray(ref_dict["S"], dtype=np.float64)
beta_1 = float(ref_dict["beta_1"])
h1_thc = np.asarray(ref_dict["h1_thc"], dtype=np.float64)
ecore = float(ref_dict["ecore"])
electrons = int(ref_dict["electrons"])
energy_upper = -22140.27335

# Reconstruct u from words2
U_q = reconstruct_u_vectors(words2 / 128.0, tree).reshape(R, B, N)
u_norm = U_q / np.linalg.norm(U_q, axis=2, keepdims=True)

# Reconstruct X_c
u_t = torch.tensor(u_norm[0], dtype=torch.float64, device=device) # (B, N)
w_t = torch.tensor(w_raw[0], dtype=torch.float64, device=device) # (B, C)
# X_c: sum_b w_bc u_bp u_bq
X_c = torch.einsum("bp,bq,bc->cpq", u_t, u_t, w_t)
# Symmetrize X_c from tril
X_tril = torch.tril(X_c)
X_c_sym = X_tril + X_tril.transpose(-1, -2) - torch.diagonal(X_tril, dim1=-2, dim2=-1).diag_embed()

h1_fixed = torch.tensor(h1_thc - 0.5 * electrons * S_ref + beta_1 * np.eye(N), dtype=torch.float64, device=device)

# Initial W: baseline has 0.5 * sum(W^2) approx 540.35
# Initialize W_c proportionally to w_t.sum(0)
w_sum = w_t.sum(0)
W_var = torch.full((C,), np.sqrt(1080.0 / C), dtype=torch.float64, device=device, requires_grad=True)

opt = torch.optim.Adam([W_var], lr=0.01)

print(f"Optimizing W for C={C}...")
t0 = time.time()
best_lam = float("inf")
best_W = None

for step in range(1000):
    opt.zero_grad()
    h1_p = h1_fixed - torch.einsum("c,cpq->pq", W_var, X_c_sym)
    ev = torch.linalg.eigvalsh(0.5 * (h1_p + h1_p.T))
    lam_h1 = torch.sum(torch.abs(ev))
    tr_h1 = torch.sum(ev)
    
    e_sos = ecore - electrons * beta_1 - (lam_h1 - tr_h1) - 0.5 * torch.sum(W_var**2)
    e_gap = energy_upper - e_sos
    
    w_abs = torch.sum(torch.abs(w_t), dim=0)
    lam_h2 = 0.25 * torch.sum((w_abs + torch.abs(W_var + w_sum))**2)
    Lambda = lam_h1 + lam_h2
    
    gap_pen = torch.relu(3.9152 - e_gap)**2 * 1e4 + torch.relu(e_gap - 4.05)**2 * 1e4
    rad = torch.clamp(e_gap * (2.0 * Lambda - e_gap), min=1.0)
    lam_eff = torch.sqrt(rad)
    
    loss = lam_eff + gap_pen
    loss.backward()
    opt.step()
    
    if e_gap.item() >= 3.9151 and e_gap.item() <= 4.05:
        if lam_eff.item() < best_lam:
            best_lam = lam_eff.item()
            best_W = W_var.detach().clone()
            
    if (step + 1) % 100 == 0:
        print(f"Step {step+1:4d} | gap = {e_gap.item():6.4f} Ha | lam_eff = {lam_eff.item():6.2f} | best = {best_lam:6.2f}")

W_opt = best_W.cpu().numpy()[None, :]
print(f"\nOptimization complete in {time.time()-t0:.2f}s. Best lambda_eff = {best_lam:.3f}")

# Save candidate
out_path = "dfthcsos/residual719_c340_fitted.npz"
out = dict(z)
for k in ["onebody_angle_words", "onebody_brot", "error_1_optimized"]:
    out.pop(k, None)
out.update(
    U=U_raw,
    w=w_raw,
    W=W_opt,
    err2=float(z["err2"]),
)
np.savez(out_path, **out)
print(f"Saved {out_path} with C={C}, best lambda_eff = {best_lam:.3f}!")
