"""Continuation from a global (R=1) dictionary to R blocks by gradual per-block pruning.

Every block starts with the full dictionary and a co-clustered share of the columns,
which reproduces the global fit. Each round drops each block's least-used vectors and
re-optimizes vectors by variable-projection Adam (rank-C PSD block weights).
"""
import argparse, json, time
import numpy as np, torch
from sklearn.cluster import SpectralCoclustering

p = argparse.ArgumentParser()
p.add_argument('--init', required=True); p.add_argument('--out', required=True)
p.add_argument('--R', type=int, default=2); p.add_argument('--B', type=int, required=True)
p.add_argument('--drop', type=int, default=24); p.add_argument('--steps', type=int, default=300)
p.add_argument('--final-steps', type=int, default=3000)
p.add_argument('--lr', type=float, default=3e-4); p.add_argument('--ridge', type=float, default=1e-7)
p.add_argument('--sweeps', type=int, default=2)
p.add_argument('--criterion', choices=['diag', 'leverage'], default='diag')
a = p.parse_args()
sp = np.load('/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/parameters/reference_spectrum.npz')
ev, V = sp['eigenvalues'], sp['eigenvectors']
d = dict(device='cuda', dtype=torch.float64)
G = torch.tensor((V * np.maximum(ev, 0)) @ V.T, **d); const = G.square().sum()
ti = torch.tensor(sp['ii'], device='cuda'); tj = torch.tensor(sp['jj'], device='cuda'); ff = torch.tensor(sp['fac'], **d)
z = np.load(a.init); U0 = z['U'].reshape(-1, 76); W0 = z['w'].reshape(len(U0), -1)
R = a.R; Ctot = W0.shape[1]; C = Ctot // R
A = W0 ** 2
cl = SpectralCoclustering(n_clusters=R, random_state=0).fit(A + 1e-12).column_labels_
# balance column groups to exactly C each
aff = np.stack([A[:, cl == k].sum(0) for k in range(R)], 1) if False else None
order = np.argsort(cl, kind='stable'); groups = [[] for _ in range(R)]
for c in order:
    k = cl[c]
    if len(groups[k]) >= C:
        k = int(np.argmin([len(g) for g in groups]))
    groups[k].append(c)
groups = [np.array(g[:C]) for g in groups]
u = torch.tensor(np.stack([U0] * R), **d)
K = []
for r in range(R):
    wr = torch.tensor(W0[:, groups[r]], **d); K.append(wr @ wr.T)


def run(u, K, steps, lr):
    Bc = u.shape[1]; eye = torch.eye(Bc, **d)
    u = u.clone().requires_grad_(True); opt = torch.optim.Adam([u], lr=lr); err = None
    for step in range(steps):
        opt.zero_grad()
        un = u / u.norm(dim=-1, keepdim=True); P = un[:, :, ti] * un[:, :, tj] * ff
        with torch.no_grad():
            S = [[(un[r] @ un[s].T).square() for s in range(R)] for r in range(R)]
            H = [P[r] @ G @ P[r].T for r in range(R)]
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
                    f = torch.linalg.solve_triangular(L.T, vv[:, -C:] * e[-C:].clamp(min=0).sqrt(), upper=True)
                    K[r] = f @ f.T
        lin = sum((K[r] * (P[r] @ G @ P[r].T)).sum() for r in range(R))
        quad = 0
        for r in range(R):
            for s in range(R):
                Srs = (un[r] @ un[s].T).square(); quad = quad + ((K[r] @ Srs) * (K[s] @ Srs.T).T).sum()
        loss = const - 2 * lin + quad; loss.backward(); opt.step()
        err = float(loss.detach().clamp(min=0).sqrt())
    with torch.no_grad():
        un = u / u.norm(dim=-1, keepdim=True)
    return un.detach(), K, err


t0 = time.time()
u, K, err = run(u, K, 50, a.lr)
print(json.dumps(dict(B=u.shape[1], err=err, s=round(time.time() - t0))), flush=True)
while u.shape[1] > a.B:
    newB = max(a.B, u.shape[1] - a.drop)
    nu, nK = [], []
    for r in range(R):
        # importance: diagonal of K weighted by feature Gram diagonal (=1 for unit vectors)
        if a.criterion == 'leverage':
            Sr = (u[r] @ u[r].T).square(); Q = torch.linalg.inv(Sr + 1e-9 * torch.eye(len(Sr), **d))
            N = K[r]; qd = Q.diagonal().clamp(min=1e-30)
            imp = 2 * (N @ Sr @ N).diagonal() / qd - (N.diagonal() / qd).square()
        else:
            imp = K[r].diagonal().abs() + (K[r].square().sum(1)).sqrt()
        keep = torch.argsort(imp, descending=True)[:newB].sort().values
        nu.append(u[r][keep]); nK.append(K[r][keep][:, keep])
    u = torch.stack(nu); K = nK
    u, K, err = run(u, K, a.steps if u.shape[1] > a.B else a.final_steps, a.lr)
    print(json.dumps(dict(B=u.shape[1], err=err, s=round(time.time() - t0))), flush=True)
    ws = []
    for r in range(R):
        e, vv = torch.linalg.eigh(K[r]); ws.append((vv[:, -C:] * e[-C:].clamp(min=0).sqrt()).flip(1).cpu().numpy())
    np.savez(a.out, U=u.cpu().numpy(), w=np.stack(ws), err=err)
