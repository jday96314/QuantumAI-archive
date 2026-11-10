"""Float fit of a block DF-THC shape (R,B,C) to the PSD pair-space target.

error_2 = || sum_{rho,c} X X^T - G ||_F in the symmetric pair metric.
Adam on (U,w) with periodic exact block-coordinate refits of w for fixed U.
"""
import argparse, json, time
import numpy as np, torch

p = argparse.ArgumentParser()
p.add_argument('--R', type=int, required=True)
p.add_argument('--B', type=int, required=True)
p.add_argument('--C', type=int, required=True)
p.add_argument('--steps', type=int, default=20000)
p.add_argument('--lr', type=float, default=3e-3)
p.add_argument('--project-every', type=int, default=1000)
p.add_argument('--cycles', type=int, default=2)
p.add_argument('--seed', type=int, default=0)
p.add_argument('--init')
p.add_argument('--out', required=True)
a = p.parse_args()
torch.manual_seed(a.seed); np.random.seed(a.seed)
sp = np.load('/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/parameters/reference_spectrum.npz')
ev, V = sp['eigenvalues'], sp['eigenvectors']
d = dict(device='cuda', dtype=torch.float32)
G64 = torch.tensor((V * np.maximum(ev, 0)) @ V.T, device='cuda', dtype=torch.float64)
G = G64.float()
ii = torch.tensor(sp['ii'], device='cuda'); jj = torch.tensor(sp['jj'], device='cuda')
fac = torch.tensor(sp['fac'], **d)
R, B, C, n = a.R, a.B, a.C, 76


def feats(u):
    un = u / u.norm(dim=-1, keepdim=True)
    return un[..., ii] * un[..., jj] * fac.to(u.dtype)


def project(u, w, cycles):
    """Exact PSD rank-C least squares for each block's H given the others."""
    with torch.no_grad():
        F = feats(u.double())
        cont = []
        for r in range(R):
            Xr = w[r].double().T @ F[r]
            cont.append(Xr.T @ Xr)
        tot = sum(cont)
        for _ in range(cycles):
            for r in range(R):
                Gr = G64 - (tot - cont[r])
                Q, Rm = torch.linalg.qr(F[r].T)
                Hp = Q.T @ Gr @ Q; Hp = (Hp + Hp.T) / 2
                e, v = torch.linalg.eigh(Hp)
                idx = torch.argsort(e, descending=True)[:C]
                root = v[:, idx] * e[idx].clamp(min=0).sqrt()
                wr = torch.linalg.solve(Rm, root)
                w[r] = wr.float()
                Xr = wr.T @ F[r]; new = Xr.T @ Xr
                tot = tot - cont[r] + new; cont[r] = new
        return float((tot - G64).square().sum().sqrt())


if a.init:
    z = np.load(a.init)
    u = torch.tensor(z['U'], **d).reshape(R, B, n); w = torch.tensor(z['w'], **d).reshape(R, B, C)
else:
    u = torch.randn(R, B, n, **d); w = torch.zeros(R, B, C, **d)
print(json.dumps(dict(projected=project(u, w, a.cycles + 2))), flush=True)
u.requires_grad_(True); w.requires_grad_(True)
opt = torch.optim.Adam([u, w], lr=a.lr); best = 1e9; t0 = time.time()
for step in range(a.steps):
    lr = a.lr * (.02 + .98 * .5 * (1 + np.cos(np.pi * step / a.steps)))
    if step and step % a.project_every == 0:
        pe = project(u.detach(), w.data, a.cycles)
        opt = torch.optim.Adam([u, w], lr=lr)
    for g in opt.param_groups: g['lr'] = lr
    opt.zero_grad()
    X = torch.bmm(w.transpose(1, 2), feats(u)).reshape(-1, len(ii))
    e2 = (X.T @ X - G).square().sum(); e2.backward()
    err = float(e2.detach().sqrt())
    if err < best and step % 50 == 0:
        best = err
        np.savez(a.out, U=(u / u.norm(dim=-1, keepdim=True)).detach().cpu().numpy().astype(float),
                 w=w.detach().cpu().numpy().astype(float), err=err)
    if step % 1000 == 0:
        print(json.dumps(dict(step=step, err=err, best=best, s=round(time.time() - t0))), flush=True)
    opt.step()
print(json.dumps(dict(final=best, R=R, B=B, C=C)), flush=True)
