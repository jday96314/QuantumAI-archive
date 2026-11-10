"""Minimize lambda_eff of a fixed quantized block dictionary inside the error and gap gates.

Block generalization of experiment 3's constrained_factor.py (augmented Lagrangian on
error_2 <= cap and E_gap >= floor), optimizing weights w (R,B,C) and offsets W (R,C).
"""
import argparse, json, time
from pathlib import Path
import numpy as np, torch

p = argparse.ArgumentParser()
p.add_argument('--fit', required=True)
p.add_argument('--out', required=True)
p.add_argument('--steps', type=int, default=6000)
p.add_argument('--lr', type=float, default=1e-5)
p.add_argument('--cap', type=float, default=.019995)
p.add_argument('--gapfloor', type=float, default=3.9151)
p.add_argument('--offset-lr-factor', type=float, default=3.)
a = p.parse_args()
E3 = Path('/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3')
z = dict(np.load(a.fit)); u = z['U']; w0 = z['w']; W0 = z['W']
R, B, C = w0.shape; n = 76
sp = np.load(E3 / 'parameters/reference_spectrum.npz'); ref = np.load(E3 / 'parameters/chem_reference.npz')
ii, jj, fac = sp['ii'], sp['jj'], sp['fac']
d = dict(device='cuda', dtype=torch.float64)
Gfull = torch.tensor((sp['eigenvectors'] * sp['eigenvalues']) @ sp['eigenvectors'].T, **d)
const = Gfull.square().sum()
ut = torch.tensor(u, **d); ut = ut / ut.norm(dim=-1, keepdim=True)
P = ut[:, :, ii] * ut[:, :, jj] * torch.tensor(fac, **d)  # R,B,pairs
Hs = [P[r] @ Gfull @ P[r].T for r in range(R)]
S = [[(ut[r] @ ut[s].T).square() for s in range(R)] for r in range(R)]
uu = torch.einsum('rbp,rbq->rbpq', ut, ut)  # R,B,n,n
w = torch.tensor(w0, **d, requires_grad=True)
off = torch.tensor(W0, **d, requires_grad=True)
h0 = torch.tensor(ref['h1_thc'] - .5 * int(ref['electrons']) * ref['S'] + float(ref['beta_1']) * np.eye(n), **d)
constant = -22140.27335 - float(ref['ecore']) + int(ref['electrons']) * float(ref['beta_1'])


def vals():
    K = [w[r] @ w[r].T for r in range(R)]
    lin = sum((K[r] * Hs[r]).sum() for r in range(R))
    quad = 0
    for r in range(R):
        for s in range(R):
            quad = quad + ((K[r] @ S[r][s]) * (K[s] @ S[s][r]).T).sum()
    e2 = const - 2 * lin + quad
    X = torch.einsum('rbc,rbpq->rcpq', w, uu)
    h = h0 - torch.einsum('rc,rcpq->pq', off, X)
    e = torch.linalg.eigvalsh(.5 * (h + h.T)); rr = torch.sqrt(e * e + 1e-14)
    gap = constant + rr.sum() - e.sum() + .5 * off.square().sum()
    eps = 1e-7
    lam = rr.sum() + .25 * (torch.sqrt(w * w + eps * eps).sum(1) + torch.sqrt((off + w.sum(1)).square() + eps * eps)).square().sum()
    return e2, gap, lam, gap * (2 * lam - gap)


opt = torch.optim.Adam([{'params': [w], 'lr': a.lr}, {'params': [off], 'lr': a.lr * a.offset_lr_factor}])
rho = 2.; mue = 0.; mug = 0.; best = float('inf'); t0 = time.time()


def save(er):
    out = dict(z)
    for k in ['onebody_angle_words', 'onebody_brot', 'error_1_optimized']:
        out.pop(k, None)
    out.update(w=w.detach().cpu().numpy(), W=off.detach().cpu().numpy(), err2=er)
    tmp = Path(a.out + '.tmp.npz'); np.savez(tmp, **out); tmp.replace(a.out)


for step in range(a.steps):
    opt.zero_grad()
    e2, gap, lam, score = vals()
    ge = e2 / (a.cap * a.cap) - 1; gg = (a.gapfloor + .0004 - gap) / 3.9155
    pen = lambda g, mu: (torch.relu(mu + rho * g).square() - mu * mu) / (2 * rho)
    loss = score / 3000 + pen(ge, mue) + pen(gg, mug)
    loss.backward(); opt.step()
    factor = .1 + .9 * (1 - step / a.steps)
    opt.param_groups[0]['lr'] = a.lr * factor; opt.param_groups[1]['lr'] = a.offset_lr_factor * a.lr * factor
    if step % 20 == 0:
        with torch.no_grad():
            e2, gap, lam, score = vals()
        er = float(e2.clamp(min=0).sqrt()); gv = float(gap); sv = float(score.clamp(min=0).sqrt())
        if er <= a.cap + 4e-6 and gv >= a.gapfloor and sv < best:
            best = sv; save(er)
        mue = max(0, mue + rho * (float(e2) / (a.cap * a.cap) - 1)); mug = max(0, mug + rho * (a.gapfloor + .0004 - gv) / 3.9155)
        if step % 200 == 0:
            print(json.dumps(dict(step=step, error_2=er, gap=gv, Lambda=float(lam), lambda_eff=sv, best=best,
                                  mue=mue, mug=mug, s=round(time.time() - t0))), flush=True)
print(json.dumps(dict(final=best)), flush=True)
