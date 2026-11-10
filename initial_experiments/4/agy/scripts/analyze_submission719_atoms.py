"""Analyze the 719 atoms of submission719: overlaps, leverage, importance."""

import numpy as np
import torch

device = torch.device("cuda")
sp = np.load("../../3/parameters/reference_spectrum.npz")
ii = torch.tensor(sp["ii"], device=device)
jj = torch.tensor(sp["jj"], device=device)
fac = torch.tensor(sp["fac"], dtype=torch.float64, device=device)
ev = torch.tensor(sp["eigenvalues"], dtype=torch.float64, device=device)
V = torch.tensor(sp["eigenvectors"], dtype=torch.float64, device=device)
G = (V * torch.clamp(ev, min=0.0)) @ V.T
total_norm_sq = torch.sum(torch.clamp(ev, min=0.0)**2).item()

sub = np.load("dfthcsos/submission.npz")
U = torch.tensor(sub["U"][0], dtype=torch.float64, device=device) # (719, 76)
B = U.shape[0]

# 1. Mutual coherence
dots = torch.abs(U @ U.T)
dots.fill_diagonal_(0.0)
max_dots = torch.max(dots, dim=1).values
print(f"Mutual coherence |u_i^T u_j| across 719 atoms:")
print(f"  Max overlap: {torch.max(max_dots).item():.4f}")
print(f"  Overlaps > 0.99: {torch.sum(max_dots > 0.99).item()} atoms")
print(f"  Overlaps > 0.95: {torch.sum(max_dots > 0.95).item()} atoms")
print(f"  Overlaps > 0.90: {torch.sum(max_dots > 0.90).item()} atoms")
print(f"  Overlaps > 0.80: {torch.sum(max_dots > 0.80).item()} atoms")

# 2. Base projection floor
P = U[:, ii] * U[:, jj] * fac
S = P @ P.T
se, sv = torch.linalg.eigh(S)
pos = se > 1e-10
sv = sv[:, pos]
se = se[pos]
iroot = (sv / torch.sqrt(se)) @ sv.T
Q = iroot @ P
H = Q @ G @ Q.T
inside_sq = torch.sum(H**2).item()
floor_719 = np.sqrt(max(0.0, total_norm_sq - inside_sq))
print(f"Base floor (719 atoms): {floor_719*1000:.3f} mHa")

# 3. Leave-one-out floor increase for each of the 719 atoms (exact in batch)
# For each atom b, what is the floor if we drop b?
# We can do this efficiently:
print("Evaluating leave-one-out floor impact for all 719 atoms...")
loos = []
for b in range(min(50, B)): # sample first 50
    keep = [i for i in range(B) if i != b]
    Pk = P[keep]
    Sk = Pk @ Pk.T
    sek, svk = torch.linalg.eigh(Sk)
    posk = sek > 1e-10
    svk = svk[:, posk]
    sek = sek[posk]
    irootk = (svk / torch.sqrt(sek)) @ svk.T
    Qk = irootk @ Pk
    Hk = Qk @ G @ Qk.T
    in_k = torch.sum(Hk**2).item()
    fl_k = np.sqrt(max(0.0, total_norm_sq - in_k))
    loos.append(fl_k)

loos = np.array(loos)
print(f"Sample of 50 leave-one-out floors:")
print(f"  Min floor increase: {(np.min(loos) - floor_719)*1000:.4f} mHa (floor = {np.min(loos)*1000:.3f} mHa)")
print(f"  Max floor increase: {(np.max(loos) - floor_719)*1000:.4f} mHa (floor = {np.max(loos)*1000:.3f} mHa)")
print(f"  Median floor increase: {(np.median(loos) - floor_719)*1000:.4f} mHa")
