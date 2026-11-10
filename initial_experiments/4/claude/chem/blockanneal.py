"""Continuation from continuous Givens angles to the brot-bit lattice for block factorizations.

Block generalization of experiment 3's anneal_angles.py: Adam on continuous angle
words with a lattice penalty whose strength grows geometrically; the block PSD
weights (rank <= C) are re-solved each step by warm-started block Gauss-Seidel.
"""
import argparse, json, time, sys
from pathlib import Path
import numpy as np, torch

p = argparse.ArgumentParser()
p.add_argument('--init', required=True); p.add_argument('--out', required=True)
p.add_argument('--brot', type=int, default=8); p.add_argument('--steps', type=int, default=12000)
p.add_argument('--lr', type=float, default=.02); p.add_argument('--ridge', type=float, default=1e-5)
p.add_argument('--columns', type=int, required=True); p.add_argument('--strength', type=float, default=100.)
p.add_argument('--start-strength', type=float, default=1e-4); p.add_argument('--sweeps', type=int, default=2)
a = p.parse_args()
E3 = Path('/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3')
z = np.load(a.init); U0 = z['U']; R, B, n = U0.shape
tree = z['tree_edges'] if 'tree_edges' in z else np.load(E3 / 'candidates/chem35c45_rot.npz')['tree_edges']
v = U0.reshape(R * B, n).copy(); th = np.zeros((R * B, 75))
for j in reversed(range(75)):
    pa, ch = tree[j]; th[:, j] = np.arctan2(v[:, ch], v[:, pa]); v[:, pa] = np.hypot(v[:, pa], v[:, ch]); v[:, ch] = 0
alpha = np.pi / (1 << (a.brot - 1)); d = dict(device='cuda', dtype=torch.float64)
theta = torch.tensor(th / alpha, requires_grad=True, **d)
sp = np.load(E3 / 'parameters/reference_spectrum.npz'); ev, V = sp['eigenvalues'], sp['eigenvectors']
G = torch.tensor((V * np.maximum(ev, 0)) @ V.T, **d); const = G.square().sum()
ti = torch.tensor(sp['ii'], device='cuda'); tj = torch.tensor(sp['jj'], device='cuda'); fac = torch.tensor(sp['fac'], **d)
eye = torch.eye(B, **d); startvec = torch.zeros((R * B, n), **d); startvec[:, 0] = 1
stages = []
for j, (pa, ch) in enumerate(tree):
    pa, ch = int(pa), int(ch)
    if not stages or any(pa in (x, y) or ch in (x, y) for _, x, y in stages[-1]):
        stages.append([])
    stages[-1].append((j, pa, ch))
stages = [tuple(torch.tensor(x, device='cuda') for x in zip(*s)) for s in stages]
K = [torch.zeros(B, B, **d) for _ in range(R)]
if 'w' in z and np.asarray(z['w']).shape[-1] == a.columns:
    w0 = torch.tensor(z['w'], **d); K = [w0[r] @ w0[r].T for r in range(R)]


def decode(t):
    u = startvec
    for ts, ps, cs in stages:
        ang = t[:, ts] * alpha; c, s = ang.cos(), ang.sin(); up, uc = u[:, ps], u[:, cs]
        new = u.clone(); new[:, ps] = c * up - s * uc; new[:, cs] = s * up + c * uc; u = new
    return u.reshape(R, B, n)


def state(t, update=True):
    u = decode(t); P = u[:, :, ti] * u[:, :, tj] * fac
    with torch.no_grad():
        S = [[(u[r] @ u[s].T).square() for s in range(R)] for r in range(R)]
        H = [P[r] @ G @ P[r].T for r in range(R)]
        if update:
            for _ in range(a.sweeps):
                for r in range(R):
                    M = H[r].clone()
                    for s in range(R):
                        if s != r:
                            M -= S[r][s] @ K[s] @ S[s][r]
                    L = torch.linalg.cholesky(S[r][r] + a.ridge * eye)
                    hh = torch.linalg.solve_triangular(L, M, upper=False)
                    hh = torch.linalg.solve_triangular(L, hh.T, upper=False).T
                    e, vv = torch.linalg.eigh(.5 * (hh + hh.T))
                    f = torch.linalg.solve_triangular(L.T, vv[:, -a.columns:] * e[-a.columns:].clamp(min=0).sqrt(), upper=True)
                    K[r] = f @ f.T
    lin = sum((K[r] * (P[r] @ G @ P[r].T)).sum() for r in range(R))
    quad = 0
    for r in range(R):
        for s in range(R):
            Srs = (u[r] @ u[s].T).square(); quad = quad + ((K[r] @ Srs) * (K[s] @ Srs.T).T).sum()
    loss = const - 2 * lin + quad
    return loss, u


opt = torch.optim.Adam([theta], lr=a.lr); best = float('inf'); t0 = time.time()
for step in range(a.steps):
    opt.zero_grad(); loss, u = state(theta)
    strength = 10 ** (np.log10(a.start_strength) + (np.log10(a.strength) - np.log10(a.start_strength)) * step / max(1, a.steps - 1))
    penalty = ((1 - torch.cos(2 * np.pi * theta)) / (4 * np.pi * np.pi)).sum()
    (loss + strength * float(loss.detach()) / theta.numel() * penalty).backward(); opt.step()
    opt.param_groups[0]['lr'] = a.lr * (.05 + .95 * .5 * (1 + np.cos(np.pi * step / a.steps)))
    if step % 100 == 0 or step == a.steps - 1:
        with torch.no_grad():
            Ksave = [k.clone() for k in K]
            words = theta.round(); hl, uq = state(words); hard = float(hl.clamp(min=0).sqrt())
            if hard < best:
                best = hard; ws = []
                for r in range(R):
                    e, vv = torch.linalg.eigh(K[r]); ws.append((vv[:, -a.columns:] * e[-a.columns:].clamp(min=0).sqrt()).flip(1).cpu().numpy())
                np.savez(a.out, U=uq.cpu().numpy(), w=np.stack(ws), W=np.zeros((R, a.columns)),
                         angle_words=(words.cpu().numpy().astype(np.int64) % (1 << a.brot)), tree_edges=tree, brot=a.brot, err2=hard)
            K[:] = Ksave
        print(json.dumps(dict(step=step, soft_error=float(loss.detach().clamp(min=0).sqrt()), hard_error=hard, best=best,
                              strength=strength, penalty=float(penalty.detach()), s=round(time.time() - t0))), flush=True)
