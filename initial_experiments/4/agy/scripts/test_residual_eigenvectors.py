"""Test computing leading eigenvectors of residual tensor R on GPU."""

import time, torch, numpy as np

device = torch.device("cuda")
sp = np.load("../../3/parameters/reference_spectrum.npz")
ii = torch.tensor(sp["ii"], device=device)
jj = torch.tensor(sp["jj"], device=device)
fac = torch.tensor(sp["fac"], dtype=torch.float64, device=device)
ev = torch.tensor(sp["eigenvalues"], dtype=torch.float64, device=device)
V = torch.tensor(sp["eigenvectors"], dtype=torch.float64, device=device)
G = (V * torch.clamp(ev, min=0.0)) @ V.T

# Random basis Q of 450 atoms
Q = torch.randn(450, 2926, dtype=torch.float64, device=device)
Q, _ = torch.linalg.qr(Q.T)
Q = Q.T

torch.cuda.synchronize()
t0 = time.time()

# Residual R = (I - Q^T Q) G (I - Q^T Q)
QG = Q @ G
R = G - Q.T @ QG - QG.T @ Q + Q.T @ (QG @ Q.T) @ Q
ev_r, v_r = torch.linalg.eigh(R)
top_ev = ev_r[-1].item()
top_vec = v_r[:, -1]

# Unpack top_vec to 76x76 matrix M
M = torch.zeros((76, 76), dtype=torch.float64, device=device)
M[ii, jj] = top_vec / fac
M[jj, ii] = top_vec / fac
w_m, v_m = torch.linalg.eigh(M) # 76 eigenvectors in R^76!

torch.cuda.synchronize()
t1 = time.time()

print(f"Computed top mode of R and its 76 orbital eigenvectors in {(t1 - t0)*1000:.1f} ms!")
print(f"Top eigenvalue of residual tensor R: {top_ev:.6f}")
print(f"Top 5 eigenvalues of unpacked orbital matrix M: {w_m[-5:].cpu().numpy().round(4)}")
