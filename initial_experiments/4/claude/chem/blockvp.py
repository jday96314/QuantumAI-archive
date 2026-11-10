"""Variable-projection Adam fit of block THC vectors (experiment 3's global_fit.py for R > 1).

Every step re-solves the optimal block-diagonal symmetric weight matrices K_rho
(full rank, small ridge) by warm-started block Gauss-Seidel on the normal
equations, then takes an Adam step on the vectors with K fixed (envelope
gradient). --columns optionally projects each K_rho to its top PSD rank.
"""
import argparse, json, time
import numpy as np, torch

p = argparse.ArgumentParser()
p.add_argument('--init', required=True)
p.add_argument('--out', required=True)
p.add_argument('--steps', type=int, default=5000)
p.add_argument('--lr', type=float, default=.002)
p.add_argument('--ridge', type=float, default=1e-7)
p.add_argument('--sweeps', type=int, default=2)
p.add_argument('--columns', type=int)
p.add_argument('--anneal', type=float, default=0., help='fraction of steps over which rank decays from B to columns')
p.add_argument('--psd', action='store_true', help='fit the positive spectral part (default)')
a = p.parse_args()
sp = np.load('/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/parameters/reference_spectrum.npz')
ev, V = sp['eigenvalues'], sp['eigenvectors']
d = dict(device='cuda', dtype=torch.float64)
G = torch.tensor((V * np.maximum(ev, 0)) @ V.T, **d); const = G.square().sum()
ti = torch.tensor(sp['ii'], device='cuda'); tj = torch.tensor(sp['jj'], device='cuda'); ff = torch.tensor(sp['fac'], **d)
z = np.load(a.init)
U0 = torch.tensor(z['U'], **d)
R, B, n = U0.shape
u = U0.clone().requires_grad_(True)
eye = torch.eye(B, **d)
K = [torch.zeros(B, B, **d) for _ in range(R)]
FAC = [None] * R
if 'w' in z:
    w0 = torch.tensor(z['w'], **d)
    K = [w0[r] @ w0[r].T for r in range(R)]
opt = torch.optim.Adam([u], lr=a.lr)


def solve(un, P):
    with torch.no_grad():
        S = [[(un[r] @ un[s].T).square() for s in range(R)] for r in range(R)]
        H = [P[r] @ G @ P[r].T for r in range(R)]
        Ls = [torch.linalg.cholesky(S[r][r] + a.ridge * eye) for r in range(R)]
        for _ in range(a.sweeps if R > 1 else 1):
            for r in range(R):
                M = H[r].clone()
                for s in range(R):
                    if s != r:
                        M -= S[r][s] @ K[s] @ S[s][r]
                if a.columns:
                    L = Ls[r]
                    hh = torch.linalg.solve_triangular(L, M, upper=False)
                    hh = torch.linalg.solve_triangular(L, hh.T, upper=False).T
                    e, v = torch.linalg.eigh(.5 * (hh + hh.T))
                    f = torch.linalg.solve_triangular(L.T, v[:, -CUR:] * e[-CUR:].clamp(min=0).sqrt(), upper=True)
                    K[r] = f @ f.T; FAC[r] = f
                else:
                    Kr = torch.cholesky_solve(torch.cholesky_solve(M, Ls[r]).T, Ls[r]).T
                    K[r] = .5 * (Kr + Kr.T)
        return S, H


def error_of(S, H, Ks):
    lin = sum((Ks[r] * H[r]).sum() for r in range(R))
    quad = sum(((Ks[r] @ S[r][s]) * (Ks[s] @ S[s][r]).T).sum() for r in range(R) for s in range(R))
    return float((const - 2 * lin + quad).clamp(min=0).sqrt())


def truncated(S, H):
    """Rank-C PSD projection of each block (whitened), refit block-coordinate."""
    Kc = []
    for r in range(R):
        L = torch.linalg.cholesky(S[r][r] + a.ridge * eye)
        M = H[r].clone()
        for s in range(R):
            if s != r:
                M -= S[r][s] @ K[s] @ S[s][r]
        hh = torch.linalg.solve_triangular(L, M, upper=False)
        hh = torch.linalg.solve_triangular(L, hh.T, upper=False).T
        e, v = torch.linalg.eigh(.5 * (hh + hh.T))
        f = torch.linalg.solve_triangular(L.T, v[:, -CUR:] * e[-CUR:].clamp(min=0).sqrt(), upper=True)
        Kc.append(f)
    return Kc


best = float('inf'); t0 = time.time()
CUR = a.columns
for step in range(a.steps):
    if a.columns and a.anneal > 0:
        frac = min(1., step / (a.anneal * a.steps)); CUR = int(round(B - (B - a.columns) * frac))
    opt.zero_grad()
    un = u / u.norm(dim=-1, keepdim=True)
    P = un[:, :, ti] * un[:, :, tj] * ff
    S, H = solve(un.detach(), P.detach())
    lin = sum((K[r] * (P[r] @ G @ P[r].T)).sum() for r in range(R))
    quad = 0
    for r in range(R):
        for s in range(R):
            Srs = (un[r] @ un[s].T).square()
            quad = quad + ((K[r] @ Srs) * (K[s] @ Srs.T).T).sum()
    loss = const - 2 * lin + quad
    loss.backward(); opt.step()
    if step % 100 == 0 or step == a.steps - 1:
        err = float(loss.detach().clamp(min=0).sqrt())
        rec = dict(step=step, error_2=err, s=round(time.time() - t0))
        if a.columns:
            fs = FAC
        score = rec.get('rank_error', err)
        rec['rank'] = CUR
        if score < best and (not a.columns or CUR == a.columns):
            best = score
            out = dict(U=un.detach().cpu().numpy(), err=score)
            if a.columns:
                out['w'] = np.stack([f.cpu().numpy() for f in fs])
            np.savez(a.out, **out)
        rec['best'] = best
        print(json.dumps(rec), flush=True)
