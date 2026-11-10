"""Variable-projection L-BFGS fit of block DF-THC vectors.

For each evaluation the block PSD weights (rank <= C) are re-solved for the
current vectors (warm-started block coordinate descent); the gradient with
respect to the vectors uses those weights fixed (envelope theorem).
"""
import argparse, json, time
import numpy as np, torch

p = argparse.ArgumentParser()
p.add_argument('--init', required=True)
p.add_argument('--out', required=True)
p.add_argument('--C', type=int, required=True)
p.add_argument('--iters', type=int, default=200)
p.add_argument('--inner', type=int, default=20)
p.add_argument('--cycles', type=int, default=2)
p.add_argument('--ridge', type=float, default=1e-9)
a = p.parse_args()
sp = np.load('/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/parameters/reference_spectrum.npz')
ev, V = sp['eigenvalues'], sp['eigenvectors']
d = dict(device='cuda', dtype=torch.float64)
G = torch.tensor((V * np.maximum(ev, 0)) @ V.T, **d); const = G.square().sum()
ti = torch.tensor(sp['ii'], device='cuda'); tj = torch.tensor(sp['jj'], device='cuda'); ff = torch.tensor(sp['fac'], **d)
z = np.load(a.init)
U = torch.tensor(z['U'], **d)
if U.dim() == 2:
    U = U[None]
R, B, n = U.shape; C = a.C
eye = torch.eye(B, **d)
Kw = [None] * R
if 'w' in z:
    w0 = torch.tensor(z['w'], **d).reshape(R, B, -1)
    Kw = [w0[r] @ w0[r].T for r in range(R)]


def feats(u):
    un = u / u.norm(dim=-1, keepdim=True)
    return un, un[..., ti] * un[..., tj] * ff


def solve_k(un, P, warm):
    with torch.no_grad():
        S = [[(un[r] @ un[s].T).square() for s in range(R)] for r in range(R)]
        H = [P[r] @ G @ P[r].T for r in range(R)]
        K = [k.clone() if k is not None else torch.zeros(B, B, **d) for k in warm]
        Ls = [torch.linalg.cholesky(S[r][r] + a.ridge * eye) for r in range(R)]
        for _ in range(a.cycles if R > 1 else 1):
            for r in range(R):
                M = H[r].clone()
                for s in range(R):
                    if s != r:
                        M -= S[r][s] @ K[s] @ S[s][r]
                L = Ls[r]
                hh = torch.linalg.solve_triangular(L, M, upper=False)
                hh = torch.linalg.solve_triangular(L, hh.T, upper=False).T
                e, v = torch.linalg.eigh(.5 * (hh + hh.T))
                f = torch.linalg.solve_triangular(L.T, v[:, -C:] * e[-C:].clamp(min=0).sqrt(), upper=True)
                K[r] = f @ f.T
        return K


def loss_fixed(un, P, K):
    lin = sum((K[r] * (P[r] @ G @ P[r].T)).sum() for r in range(R))
    quad = 0
    for r in range(R):
        for s in range(R):
            S_rs = (un[r] @ un[s].T).square()
            quad = quad + ((K[r] @ S_rs) * (K[s] @ S_rs.T).T).sum()
    return const - 2 * lin + quad


u = U.clone().requires_grad_(True)
best = float('inf'); t0 = time.time()
opt = torch.optim.LBFGS([u], lr=1, max_iter=a.inner, history_size=30, line_search_fn='strong_wolfe')
for it in range(a.iters):
    def closure():
        global Kw
        opt.zero_grad()
        un, P = feats(u)
        Kw = solve_k(un.detach(), P.detach(), Kw)
        l = loss_fixed(un, P, Kw); l.backward(); return l
    opt.step(closure)
    with torch.no_grad():
        un, P = feats(u); Kw = solve_k(un, P, Kw); err = float(loss_fixed(un, P, Kw).clamp(min=0).sqrt())
    if err < best:
        best = err; ws = []
        for r in range(R):
            e, v = torch.linalg.eigh(Kw[r]); ws.append((v[:, -C:] * e[-C:].clamp(min=0).sqrt()).flip(1).cpu().numpy())
        np.savez(a.out, U=un.cpu().numpy(), w=np.stack(ws), err=err)
    print(json.dumps(dict(it=it, err=err, best=best, s=round(time.time() - t0))), flush=True)
