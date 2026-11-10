"""Run 300 steps of STE Adam optimization on B=450 to test convergence floor."""

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
tree = np.load("../../3/candidates/chem35c45_rot.npz")["tree_edges"]

stages = []
for j, (pa, ch) in enumerate(tree):
    pa, ch = int(pa), int(ch)
    if not stages or any(pa in (x, y) or ch in (x, y) for _, x, y in stages[-1]):
        stages.append([])
    stages[-1].append((j, pa, ch))
stages = [tuple(torch.tensor(x, device=device) for x in zip(*s)) for s in stages]

B = 450
bits = 8
alpha = np.pi / (1 << (bits - 1))

sub = np.load("dfthcsos/submission.npz")
words = sub["angle_words"][76:76+B]
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

opt = torch.optim.Adam([theta], lr=0.1)
sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=300, eta_min=0.005)

t0 = time.time()
best_floor = float("inf")
best_theta = None

for step in range(300):
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
        un_eval = decode(theta.round())
        Pe = un_eval[:, ii] * un_eval[:, jj] * fac
        Se = Pe @ Pe.T + 1e-6 * torch.eye(B, dtype=torch.float64, device=device)
        Le = torch.linalg.cholesky(Se)
        He = Pe @ G @ Pe.T
        Ye = torch.linalg.solve_triangular(Le, He, upper=False)
        Ye = torch.linalg.solve_triangular(Le, Ye.T, upper=False).T
        in_e = torch.sum(Ye**2).item()
        fl_e = np.sqrt(max(0.0, total_norm_sq - in_e))
        if fl_e < best_floor:
            best_floor = fl_e
            best_theta = theta.round().cpu().numpy().astype(np.int64)
        if (step + 1) % 25 == 0:
            print(f"Step {step+1:3d} | Discrete Floor: {fl_e*1000:7.3f} mHa (best: {best_floor*1000:7.3f} mHa) | Elapsed: {time.time()-t0:.1f}s")

print(f"\nFinal Best Discrete Floor: {best_floor*1000:.3f} mHa")
np.savez("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy/dfthcsos/ste450_best.npz", angle_words=best_theta, tree_edges=tree, brot=8, floor=best_floor)
