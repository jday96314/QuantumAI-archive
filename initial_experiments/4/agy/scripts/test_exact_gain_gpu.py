"""Test exact formula for Delta(u) (gain in captured Frobenius norm squared) on GPU."""

import time, torch, numpy as np

device = torch.device("cuda")
sp = np.load("../../3/parameters/reference_spectrum.npz")
ii = torch.tensor(sp["ii"], device=device)
jj = torch.tensor(sp["jj"], device=device)
fac = torch.tensor(sp["fac"], dtype=torch.float64, device=device)
ev = torch.tensor(sp["eigenvalues"], dtype=torch.float64, device=device)
V = torch.tensor(sp["eigenvectors"], dtype=torch.float64, device=device)
G = (V * torch.clamp(ev, min=0.0)) @ V.T
total_norm_sq = torch.sum(torch.clamp(ev, min=0.0)**2).item()

# Current dictionary: say 400 random unit vectors
B = 400
u = torch.randn(B, 76, dtype=torch.float64, device=device)
u = u / torch.norm(u, dim=-1, keepdim=True)
P = u[:, ii] * u[:, jj] * fac
S = P @ P.T
se, sv = torch.linalg.eigh(S)
pos = se > 1e-10
sv = sv[:, pos]
se = se[pos]
iroot = (sv / torch.sqrt(se)) @ sv.T
Q = iroot @ P
H = Q @ G @ Q.T
E0 = torch.sum(H**2).item()

A = Q @ G # shape (B, 2926)

# Test candidate batch of M atoms
M = 2048
u_cand = torch.randn(M, 76, dtype=torch.float64, device=device)
u_cand = u_cand / torch.norm(u_cand, dim=-1, keepdim=True)
p_cand = u_cand[:, ii] * u_cand[:, jj] * fac # (M, 2926)

torch.cuda.synchronize()
t0 = time.time()

# Compute gain for all M candidates in parallel
# q = Q @ p_cand.T -> (B, M)
q = Q @ p_cand.T
norm_p_sq = torch.sum(p_cand**2, dim=1) # (M,)
norm_q_sq = torch.sum(q**2, dim=0) # (M,)
d = torch.clamp(norm_p_sq - norm_q_sq, min=1e-10) # (M,)

# Ap = p_cand @ A.T -> (M, B)
Ap = p_cand @ A.T
# Hq = H @ q -> (B, M)
Hq = H @ q
# r = Ap.T - Hq -> (B, M)
r = Ap.T - Hq
norm_r_sq = torch.sum(r**2, dim=0) # (M,)

# pGp = diag(p_cand @ G @ p_cand.T)
pGp = torch.sum((p_cand @ G) * p_cand, dim=1) # (M,)
# qAp = sum_b q_bm * Ap_mb
qAp = torch.sum(q * Ap.T, dim=0) # (M,)
# qHq = sum_b q_bm * Hq_bm
qHq = torch.sum(q * Hq, dim=0) # (M,)
s = pGp - 2.0 * qAp + qHq # (M,)

gains = 2.0 * norm_r_sq / d + (s / d)**2 # (M,)

torch.cuda.synchronize()
t1 = time.time()

print(f"Computed exact gain Delta(u) for {M} candidate atoms in {(t1 - t0)*1000:.2f} ms!")
print(f"Max gain: {torch.max(gains).item():.6f}, Min gain: {torch.min(gains).item():.6f}")

# Verify exactness for the best candidate
best_idx = torch.argmax(gains).item()
best_u = u_cand[best_idx:best_idx+1]
P_new = torch.cat([P, best_u[:, ii] * best_u[:, jj] * fac], dim=0)
S_new = P_new @ P_new.T
se_n, sv_n = torch.linalg.eigh(S_new)
pos_n = se_n > 1e-10
iroot_n = (sv_n[:, pos_n] / torch.sqrt(se_n[pos_n])) @ sv_n[:, pos_n].T
Q_new = iroot_n @ P_new
H_new = Q_new @ G @ Q_new.T
E_new = torch.sum(H_new**2).item()
actual_gain = E_new - E0

print(f"Predicted gain: {gains[best_idx].item():.6f}")
print(f"Actual gain:    {actual_gain:.6f}")
print(f"Relative error: {abs(gains[best_idx].item() - actual_gain) / actual_gain:.2e}")
