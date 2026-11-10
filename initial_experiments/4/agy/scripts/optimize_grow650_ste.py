"""Run STE Adam on the 650 discrete atoms of grow650 to drive floor below 20 mHa."""

import time, torch, numpy as np
import sys
from pathlib import Path

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
    str(ROOT / "initial_experiments/3"),
]
from compress_chem import encode
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors

device = torch.device("cuda")
sp = np.load("../../3/parameters/reference_spectrum.npz")
ii = torch.tensor(sp["ii"], device=device)
jj = torch.tensor(sp["jj"], device=device)
fac = torch.tensor(sp["fac"], dtype=torch.float64, device=device)
ev = torch.tensor(sp["eigenvalues"], dtype=torch.float64, device=device)
V = torch.tensor(sp["eigenvectors"], dtype=torch.float64, device=device)
G = (V * torch.clamp(ev, min=0.0)) @ V.T
total_norm_sq = torch.sum(torch.clamp(ev, min=0.0)**2).item()
tree = np.load("../../3/candidates/chem35c45_rot.npz")["tree_edges"]

stages = []
for j, (pa, ch) in enumerate(tree):
    pa, ch = int(pa), int(ch)
    if not stages or any(pa in (x, y) or ch in (x, y) for _, x, y in stages[-1]):
        stages.append([])
    stages[-1].append((j, pa, ch))
stages = [tuple(torch.tensor(x, device=device) for x in zip(*s)) for s in stages]

bits = 8
alpha = np.pi / (1 << (bits - 1))

z = np.load("../../3/candidates/grow650_q8_c350.npz")
words = z["angle_words"].reshape(-1, 75)
B = len(words)
print(f"Loaded grow650 with B={B} atoms.")

theta = torch.tensor(words, dtype=torch.float64, device=device, requires_grad=True)

def decode(th):
    angles = (th + (th.round() - th).detach()) * alpha
    un = torch.zeros((len(th), 76), dtype=torch.float64, device=device)
    un[:, 0] = 1.0
    for ts, ps, cs in stages:
        ang = angles[:, ts]
        c, s = ang.cos(), ang.sin()
        up, uc = un[:, ps], un[:, cs]
        new = un.clone()
        new[:, ps] = c * up - s * uc
        new[:, cs] = s * up + c * uc
        un = new
    return un

# Initial discrete evaluation
with torch.no_grad():
    un0 = decode(theta.round())
    P0 = un0[:, ii] * un0[:, jj] * fac
    S0 = P0 @ P0.T + 1e-6 * torch.eye(B, dtype=torch.float64, device=device)
    L0 = torch.linalg.cholesky(S0)
    H0 = P0 @ G @ P0.T
    Y0 = torch.linalg.solve_triangular(L0, H0, upper=False)
    Y0 = torch.linalg.solve_triangular(L0, Y0.T, upper=False).T
    in0 = torch.sum(Y0**2).item()
    fl0 = np.sqrt(max(0.0, total_norm_sq - in0))
    print(f"Initial discrete floor: {fl0*1000:.4f} mHa")

opt = torch.optim.Adam([theta], lr=0.05)
sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=200, eta_min=0.002)

t0 = time.time()
best_floor = fl0
best_words = theta.round().detach().cpu().numpy().astype(np.int64)

for step in range(200):
    opt.zero_grad()
    un = decode(theta)
    P = un[:, ii] * un[:, jj] * fac
    S = P @ P.T + 1e-6 * torch.eye(B, dtype=torch.float64, device=device)
    L = torch.linalg.cholesky(S)
    H = P @ G @ P.T
    Y = torch.linalg.solve_triangular(L, H, upper=False)
    Y = torch.linalg.solve_triangular(L, Y.T, upper=False).T
    inside_sq = torch.sum(Y**2)
    loss = -inside_sq
    loss.backward()
    opt.step()
    sched.step()
    
    with torch.no_grad():
        un_e = decode(theta.round())
        Pe = un_e[:, ii] * un_e[:, jj] * fac
        Se = Pe @ Pe.T + 1e-6 * torch.eye(B, dtype=torch.float64, device=device)
        Le = torch.linalg.cholesky(Se)
        He = Pe @ G @ Pe.T
        Ye = torch.linalg.solve_triangular(Le, He, upper=False)
        Ye = torch.linalg.solve_triangular(Le, Ye.T, upper=False).T
        in_e = torch.sum(Ye**2).item()
        fl_e = np.sqrt(max(0.0, total_norm_sq - in_e))
        if fl_e < best_floor:
            best_floor = fl_e
            best_words = theta.round().detach().cpu().numpy().astype(np.int64)
            
        if (step + 1) % 10 == 0:
            print(f"Step {step+1:3d} | Discrete Floor: {fl_e*1000:7.4f} mHa (best: {best_floor*1000:7.4f} mHa) | Elapsed: {time.time()-t0:.1f}s")
            
        if fl_e * 1000 < 19.80:
            print(f"\n*** SUCCESS: Reached floor {fl_e*1000:.4f} mHa < 19.80 mHa at step {step+1}! ***")
            break

out_path = "/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy/dfthcsos/grow650_ste_best.npz"
U_best = reconstruct_u_vectors(best_words / 128.0, tree)
np.savez(
    out_path,
    U=U_best[None, :, :],
    angle_words=best_words,
    tree_edges=tree,
    brot=8,
    floor=best_floor,
)
print(f"Saved optimized dictionary to {out_path} with floor = {best_floor*1000:.4f} mHa")
