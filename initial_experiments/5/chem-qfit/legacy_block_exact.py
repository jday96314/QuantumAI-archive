# Adapted from round-4 Claude block_coordinate.py; original control remains read-only.
# Changes: full signed target, three CPU threads, exclude null-angle proposals.
"""Integer Givens-angle descent for a block DF-THC factorization (R blocks).

Generalizes experiment 3's global_coordinate.py to R > 1. For fixed vectors the
block PSD weights K_rho (rank <= C each) are solved by block-coordinate whitened
least squares (variable projection). Moves change groups of angle words by +-1
along the negative gradient and are accepted only when the projected loss drops.
"""
import argparse, json, time, sys
from pathlib import Path
import numpy as np, torch

sys.path.insert(0, '/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3')
from compress_chem import encode  # noqa: E402

p = argparse.ArgumentParser()
p.add_argument('--init', required=True)
p.add_argument('--out', required=True)
p.add_argument('--brot', type=int, default=8)
p.add_argument('--steps', type=int, default=2000)
p.add_argument('--ridge', type=float, default=1e-5)
p.add_argument('--tree', help='numpy file with Givens tree edges')
p.add_argument('--columns', type=int, required=True, help='rank C per block')
p.add_argument('--moves', type=int, default=64)
p.add_argument('--cycles', type=int, default=2)
p.add_argument('--padding', type=float, default=.003)
p.add_argument('--pair-candidates', type=int, default=1024)
a = p.parse_args(); torch.set_num_threads(3)

z = np.load(a.init)
U0 = z['U']; R, B, n = U0.shape
tree = np.load(a.tree) if a.tree else (z['tree_edges'] if 'tree_edges' in z else
        np.load('/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/candidates/chem35c45_rot.npz')['tree_edges'])
T = len(tree)
sp = np.load('/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/parameters/reference_spectrum.npz')
ev, V, ii, jj, fac = sp['eigenvalues'], sp['eigenvectors'], sp['ii'], sp['jj'], sp['fac']
d = dict(device='cuda', dtype=torch.float64)
G = torch.tensor((V * ev) @ V.T, **d); const = G.square().sum()
ti = torch.tensor(ii, device='cuda'); tj = torch.tensor(jj, device='cuda'); ff = torch.tensor(fac, **d)
if 'angle_words' in z and int(z.get('brot', 0)) == a.brot:
    words = np.asarray(z['angle_words']).reshape(R * B, T)
else:
    words = encode(U0.reshape(R * B, n), tree, a.brot)
theta = torch.tensor(words, **d, requires_grad=True)
alpha = np.pi / (1 << (a.brot - 1))
stages = []
for t, (pa, ch) in enumerate(tree):
    pa, ch = int(pa), int(ch)
    if not stages or any(pa in (x, y) or ch in (x, y) for _, x, y in stages[-1]):
        stages.append([])
    stages[-1].append((t, pa, ch))
stages = [tuple(torch.tensor(x, device='cuda') for x in zip(*s)) for s in stages]
startvec = torch.zeros((R * B, n), **d); startvec[:, 0] = 1
eye = torch.eye(B, **d)
Kcur = [None] * R
if 'w' in z and np.asarray(z['w']).shape[-1] == a.columns:
    _w = torch.tensor(z['w'], **d); Kcur = [_w[r] @ _w[r].T for r in range(R)]


def decode(th):
    ang = th * alpha; u = startvec
    for ts, ps, cs in stages:
        an = ang[:, ts]; co, sn = an.cos(), an.sin(); up, uc = u[:, ps], u[:, cs]
        new = u.clone(); new[:, ps] = co * up - sn * uc; new[:, cs] = sn * up + co * uc; u = new
    return u.reshape(R, B, n)


def solve_k(u, warm):
    """Block-coordinate PSD rank-C weights for fixed vectors; returns K list and pieces."""
    with torch.no_grad():
        P = u[:, :, ti] * u[:, :, tj] * ff
        S = [[(u[r] @ u[s].T).square() for s in range(R)] for r in range(R)]
        H = [P[r] @ G @ P[r].T for r in range(R)]
        K = [k.clone() if k is not None else torch.zeros(B, B, **d) for k in warm]
        Ls = [torch.linalg.cholesky(S[r][r] + a.ridge * eye) for r in range(R)]
        for _ in range(a.cycles):
            for r in range(R):
                M = H[r].clone()
                for s in range(R):
                    if s != r:
                        M -= S[r][s] @ K[s] @ S[s][r]
                L = Ls[r]
                hh = torch.linalg.solve_triangular(L, M, upper=False)
                hh = torch.linalg.solve_triangular(L, hh.T, upper=False).T
                e, v = torch.linalg.eigh(.5 * (hh + hh.T))
                f = torch.linalg.solve_triangular(L.T, v[:, -a.columns:] * e[-a.columns:].clamp(min=0).sqrt(), upper=True)
                K[r] = f @ f.T
        return K


def loss_of(u, K):
    P = u[:, :, ti] * u[:, :, tj] * ff
    lin = sum((K[r] * (P[r] @ G @ P[r].T)).sum() for r in range(R))
    quad = 0
    for r in range(R):
        for s in range(R):
            S_rs = (u[r] @ u[s].T).square()
            quad = quad + ((K[r] @ S_rs) * (K[s] @ S_rs.T).T).sum()
    return (const - 2 * lin + quad)


def evaluate(th, warm):
    u = decode(th)
    K = solve_k(u.detach(), warm)
    return loss_of(u, K), K


def curvature(K, g=None, pad=None):
    pad = a.padding if pad is None else pad
    """Per-angle predicted curvature of a +-1 word move (experiment 3's estimate, per block)."""
    hs = []
    with torch.no_grad():
        for r in range(R):
            rows = slice(r * B, (r + 1) * B)
            th = theta[rows]
            v = startvec[rows].clone(); J = torch.zeros((B, T, n), **d)
            for t, (pp, cc) in enumerate(tree):
                pp, cc = int(pp), int(cc); ang = th[:, t] * alpha; co, sn = ang.cos(), ang.sin()
                jp, jc = J[:, :, pp].clone(), J[:, :, cc].clone()
                J[:, :, pp] = co[:, None] * jp - sn[:, None] * jc; J[:, :, cc] = sn[:, None] * jp + co[:, None] * jc
                vp, vc = v[:, pp].clone(), v[:, cc].clone(); v[:, pp] = co * vp - sn * vc; v[:, cc] = sn * vp + co * vc
                J[:, t, pp] = -alpha * v[:, cc]; J[:, t, cc] = alpha * v[:, pp]
            S = (v @ v.T).square(); L = torch.linalg.cholesky(S + a.ridge * eye)
            raw = 2 * J.square().sum(2)
            cross = 2 * (J @ v.T) * (v @ v.T)[:, None, :]
            flat = cross.reshape(B * T, B)
            projected = (flat * torch.cholesky_solve(flat.T, L).T).sum(1).reshape(B, T)
            eff = (raw - projected).clamp(min=0) + pad * raw
            yn = (K[r] @ S * K[r]).sum(1)
            h = (4 * eff * yn[:, None]).clamp(min=1e-18); hs.append(h)
            if g is not None and a.pair_candidates:
                gr = g[rows]
                barrier = .5 * h - gr.abs(); ratio = barrier / h
                cands = [barrier.flatten().argsort()[:a.pair_candidates], ratio.flatten().argsort()[:a.pair_candidates // 2],
                         (gr.square() / h).flatten().argsort(descending=True)[:a.pair_candidates // 2]]
                ids = torch.unique(torch.cat(cands)); rws = ids // T
                deriv = J.reshape(B * T, n)[ids]; uv = v[rws]; ju = deriv @ uv.T
                gram = 2 * (deriv @ deriv.T) * (uv @ uv.T) + 2 * ju * ju.T
                cr = flat[ids]; gram -= cr @ torch.cholesky_solve(cr.T, L)
                yy = K[r] @ S @ K[r]; hh = 4 * gram * yy[rws[:, None], rws[None, :]]
                signs = -gr.flatten()[ids].sign(); bb = .5 * h.flatten()[ids] - gr.flatten()[ids].abs()
                inter = (hh * signs[:, None] * signs[None, :]).cpu().numpy(); np.fill_diagonal(inter, 0)
                field = bb.cpu().numpy().copy(); sel = np.zeros(len(ids), bool)
                for _ in range(max(512, 2 * a.moves)):
                    gain = field * (1 - 2 * sel.astype(int))
                    if sel.sum() >= max(32, a.moves): gain[~sel] = np.inf
                    j = int(gain.argmin())
                    if gain[j] >= -1e-16: break
                    dirn = -1 if sel[j] else 1; sel[j] = not sel[j]; field += dirn * inter[:, j]
                chosen = (ids[torch.tensor(np.flatnonzero(sel), device='cuda')] + r * B * T).cpu().tolist() if sel.any() else []
                PAIRS.append(chosen)
    return torch.cat(hs)


PAIRS = []
best = float('inf'); start = time.time(); accepted = 0; stale = 0
loss, Kcur = evaluate(theta, Kcur)
for step in range(a.steps):
    theta.grad = None
    loss, Kcur = evaluate(theta, Kcur)
    loss.backward(); g = theta.grad.detach().clone()
    value = float(loss.detach()); err = float(np.sqrt(max(value, 0)))
    with torch.no_grad():
        if err < best:
            best = err
            u = decode(theta).cpu().numpy()
            np.savez(a.out, U=u, K=np.stack([k.cpu().numpy() for k in Kcur]), err2=err,
                     angle_words=theta.cpu().numpy().astype(np.int64) % (1 << a.brot), tree_edges=tree, brot=a.brot)
        PAIRS.clear(); h = curvature(Kcur, g); scores = -g.abs() + .5 * h
        scores[h < 1e-16] = torch.inf
        order = torch.argsort(scores.flatten())
        rows_seen = set(); choices = []
        for x in order[:max(512, 16 * a.moves)].cpu().tolist():
            if float(scores.view(-1)[x]) >= 0:
                break
            row = x // T
            if row not in rows_seen:
                choices.append(x); rows_seen.add(row)
            if len(choices) >= a.moves:
                break
        changed = False
        sizes = sorted({a.moves >> k for k in range(a.moves.bit_length())} | {1}, reverse=True)
        together = [x for pr in PAIRS for x in pr]
        groups = ([together] if together else []) + [pr for pr in PAIRS if pr and len(PAIRS) > 1] + [choices[:k] for k in sizes if 0 < k <= len(choices)] + [[x] for x in order[1:33].cpu().tolist()]
        for grp in groups:
            k = len(grp); ids = torch.tensor(grp, device='cuda')
            flat = theta.view(-1); direction = -g.flatten()[ids].sign()
            flat[ids] += direction
            nv = float(evaluate(theta, Kcur)[0])
            if nv < value - 1e-14:
                changed = True; accepted += k; break
            flat[ids] -= direction
        stale = 0 if changed else stale + 1
    if step % 5 == 0:
        print(json.dumps(dict(step=step, error_2=err, best=best, accepted=accepted, changed=changed,
                              seconds=round(time.time() - start))), flush=True)
    if stale >= 3:
        break
z2 = dict(np.load(a.out)); Ks = z2.pop('K'); ws = []
for r in range(R):
    e, v = np.linalg.eigh(Ks[r]); ix = np.argsort(e)[::-1]; e = np.maximum(e[ix], 0); v = v[:, ix]
    ws.append(v[:, :a.columns] * np.sqrt(e[:a.columns]))
z2.update(w=np.stack(ws), W=np.zeros((R, a.columns))); np.savez(a.out, **z2)
print(json.dumps(dict(final=best)), flush=True)
