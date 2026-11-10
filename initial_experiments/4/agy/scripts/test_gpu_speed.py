"""Test GPU speed of global_fit on RTX 4090."""
import sys
import time
import torch
import numpy as np

print("CUDA available:", torch.cuda.is_available())
if torch.cuda.is_available():
    print("Device:", torch.cuda.get_device_name(0))

sp = np.load("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/parameters/reference_spectrum.npz")
ev, V = sp["eigenvalues"], sp["eigenvectors"]
ii, jj = sp["ii"], sp["jj"]
fac = sp["fac"]

B = 300
N = 76
C = 250

device = "cuda"
dtype = torch.float64

# Start from random orthonormal vectors or top eigenvectors
rng = np.random.default_rng(42)
u0 = rng.normal(size=(B, N))
u0 /= np.linalg.norm(u0, axis=1, keepdims=True)

u = torch.tensor(u0, device=device, dtype=dtype, requires_grad=True)
G = torch.tensor((V * np.maximum(ev, 0)) @ V.T, device=device, dtype=dtype)
ti = torch.tensor(ii, device=device)
tj = torch.tensor(jj, device=device)
ff = torch.tensor(fac, device=device, dtype=dtype)
eye = torch.eye(B, device=device, dtype=dtype)

opt = torch.optim.Adam([u], lr=0.002)

start = time.time()
for step in range(100):
    opt.zero_grad()
    un = u / u.norm(dim=1, keepdim=True)
    P = un[:, ti] * un[:, tj] * ff
    S = (un @ un.T).square()
    H = (P @ G) @ P.T
    Sr = S + 1e-6 * eye
    L = torch.linalg.cholesky(Sr)
    K = torch.cholesky_solve(torch.cholesky_solve(H, L).T, L).T
    K = 0.5 * (K + K.T)
    KS = K @ Sr
    loss = G.square().sum() - 2 * (K * H).sum() + (KS * KS.T).sum()
    loss.backward()
    opt.step()

elapsed = time.time() - start
print(f"100 steps took {elapsed:.2f} seconds ({elapsed/100*1000:.1f} ms/step)!")
with torch.no_grad():
    physical = K @ S
    err2 = float((G.square().sum() - 2 * (K * H).sum() + (physical * physical.T).sum()).clamp(min=0).sqrt())
    print(f"err2 after 100 steps: {err2*1000:.2f} mHa")
