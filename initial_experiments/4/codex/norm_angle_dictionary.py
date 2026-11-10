"""Integer-angle search with factor columns projected in their existing gauge.

Unlike pure two-body variable projection, this retains the physical target
factor matrices and their column orientation. Changing an atom re-expresses
those targets in the new dictionary before evaluating normalization, one-body
shifts and the actual reference residual. No new column rotation is introduced.
Outputs are research fits, requiring one-body rebuilding and native validation.
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
import torch


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--source', default='candidates/jacobi_l1/submission.npz')
    p.add_argument('--output', default='candidates/norm_angles719.npz')
    p.add_argument('--steps', type=int, default=200)
    p.add_argument('--trials', type=int, default=128)
    p.add_argument('--moves', type=int, default=16)
    p.add_argument('--mue', type=float, default=.62)
    p.add_argument('--mug', type=float, default=.6)
    p.add_argument('--lattice', action='store_true')
    p.add_argument('--rows', type=int, default=128)
    p.add_argument('--max-step', type=int, default=8)
    p.add_argument('--project-gap', action='store_true')
    a = p.parse_args()
    torch.set_num_threads(2)
    z = dict(np.load(a.source))
    sp = np.load('../../3/parameters/reference_spectrum.npz')
    ref = np.load('../../3/parameters/chem_reference.npz')
    tree = z['tree_edges']
    native = 'w' not in z
    w0 = z['W'][0] if native else z['w'][0]
    off0 = z['H2_I'][0] - w0.sum(0) if native else z['W'][0]
    d = dict(device='cuda', dtype=torch.float64)
    theta = torch.tensor(z['angle_words'][76:] if native else z['angle_words'], requires_grad=True, **d)
    B, T = theta.shape
    ii = torch.tensor(sp['ii'], device='cuda')
    jj = torch.tensor(sp['jj'], device='cuda')
    fac = torch.tensor(sp['fac'], **d)
    G = torch.tensor((sp['eigenvectors'] * sp['eigenvalues']) @ sp['eigenvectors'].T, **d)
    constant = G.square().sum()
    off = torch.tensor(off0, **d)
    h0 = torch.tensor(ref['h1_thc'] - .5 * int(ref['electrons']) * ref['S'] + float(ref['beta_1']) * np.eye(76), **d)
    gap_constant = -22140.27335 - float(ref['ecore']) + int(ref['electrons']) * float(ref['beta_1'])
    start = torch.zeros((B, 76), **d)
    start[:, 0] = 1
    alpha = np.pi / 128
    stages = []
    for t, (pa, ch) in enumerate(tree):
        pa, ch = int(pa), int(ch)
        if not stages or any(pa in (x, y) or ch in (x, y) for _, x, y in stages[-1]):
            stages.append([])
        stages[-1].append((t, pa, ch))
    stages = [tuple(torch.tensor(x, device='cuda') for x in zip(*stage)) for stage in stages]

    def decode():
        u = start
        for ts, ps, cs in stages:
            ang = theta[:, ts] * alpha
            co, sn = ang.cos(), ang.sin()
            up, uc = u[:, ps], u[:, cs]
            new = u.clone()
            new[:, ps] = co * up - sn * uc
            new[:, cs] = sn * up + co * uc
            u = new
        return u

    def jacobian():
        u = start.clone()
        J = torch.zeros((B, T, 76), **d)
        for t, (pa, ch) in enumerate(tree):
            pa, ch = int(pa), int(ch)
            ang = theta[:, t] * alpha
            co, sn = ang.cos(), ang.sin()
            up, uc = u[:, pa].clone(), u[:, ch].clone()
            jp, jc = J[:, :, pa].clone(), J[:, :, ch].clone()
            J[:, :, pa] = co[:, None] * jp - sn[:, None] * jc
            J[:, :, ch] = sn[:, None] * jp + co[:, None] * jc
            u[:, pa] = co * up - sn * uc
            u[:, ch] = sn * up + co * uc
            J[:, t, pa] = -alpha * u[:, ch]
            J[:, t, ch] = alpha * u[:, pa]
        return J

    with torch.no_grad():
        u0 = decode()
        p0 = u0[:, ii] * u0[:, jj] * fac
        target = torch.tensor(w0.T, **d) @ p0

    def project_gap(X):
        # The error_2 gate is independent of offsets. Project the offsets by
        # minimum-norm Newton corrections onto the SOS gap margin, using the
        # exact eigenvalue derivative. This remains classical fit optimization.
        with torch.no_grad():
            for _ in range(5):
                packed = (off @ X.detach()) / fac
                hc = torch.zeros((76, 76), **d)
                hc[ii, jj] = packed
                hc = hc + hc.T - torch.diag(hc.diag())
                eig, vec = torch.linalg.eigh(h0 - hc)
                gap = gap_constant + eig.abs().sum() - eig.sum() + .5 * off.square().sum()
                delta = 3.9155 - float(gap)
                if abs(delta) < 1e-10:
                    break
                negative = vec[:, eig < 0]
                projector = negative @ negative.T
                dual = projector[ii, jj] * fac
                grad = off + 2 * (X.detach() @ dual)
                off.add_(delta * grad / grad.square().sum().clamp(min=1e-20))

    def state():
        u = decode()
        P = u[:, ii] * u[:, jj] * fac
        S = P @ P.T
        w = torch.linalg.solve(S, P @ target.T)
        X = w.T @ P
        if a.project_gap:
            project_gap(X)
        # ||X^T X-G||^2 via smaller factors; keep the full reference spectrum.
        err2 = constant - 2 * (X * (X @ G)).sum() + (X @ X.T).square().sum()
        packed = (off @ X) / fac
        hc = torch.zeros((76, 76), **d)
        hc[ii, jj] = packed
        hc = hc + hc.T - torch.diag(hc.diag())
        eig = torch.linalg.eigvalsh(h0 - hc)
        abs_eig = torch.sqrt(eig.square() + 1e-14)
        gap = gap_constant + abs_eig.sum() - eig.sum() + .5 * off.square().sum()
        lam = abs_eig.sum() + .25 * (torch.sqrt(w.square() + 1e-14).sum(0) + torch.sqrt((off + w.sum(0)).square() + 1e-14)).square().sum()
        score2 = gap * (2 * lam - gap)
        ge = err2 / .019995**2 - 1
        gg = (3.9155 - gap) / 3.9155
        rho = 4.
        penalty = lambda g, mu: (torch.relu(mu + rho * g).square() - mu**2) / (2 * rho)
        loss = score2 / 3000 + penalty(ge, a.mue) + penalty(gg, a.mug)
        return loss, err2, gap, lam, score2, u, w

    best = float('inf')
    bestloss = float('inf')
    began = time.time()
    total = 0

    def save(path, st):
        out = dict(U=st[5].detach().cpu().numpy()[None], w=st[6].detach().cpu().numpy()[None],
                   W=off.detach().cpu().numpy()[None], angle_words=theta.detach().cpu().numpy().astype(np.int64) % 256,
                   tree_edges=tree, brot=np.array(8), err2=float(st[1].clamp(min=0).sqrt()))
        path = Path(path)
        tmp = path.with_name(path.stem + '.tmp.npz')
        np.savez(tmp, **out)
        tmp.replace(path)

    for step in range(a.steps):
        theta.grad = None
        st = state()
        st[0].backward()
        gradient = theta.grad.detach().clone()
        value = float(st[0].detach())
        with torch.no_grad():
            # Use gradient orders at several scales so a few highly weighted
            # coordinates cannot monopolize the nonlinear trial budget.
            order = gradient.abs().flatten().argsort(descending=True).tolist()
            choices = []
            if a.lattice:
                from lattice_babai import proposals
                J = jacobian()
                u, w = st[5], st[6]
                P = u[:, ii] * u[:, jj] * fac
                S = P @ P.T
                K = w @ w.T
                cross = 2 * (J @ u.T) * (u @ u.T)[:, None, :]
                flat = cross.reshape(B*T, B)
                solved = torch.linalg.solve(S, flat.T).T.reshape(B, T, B)
                ynorm = (K @ S * K).sum(1)
                metric_scale = 4 * a.mue / .019995**2
                diag = (metric_scale * ynorm[:, None] * (2*J.square().sum(2) - (cross*solved).sum(2))).clamp(min=1e-24)
                ids = (gradient.square()/diag).sum(1).argsort(descending=True)[:a.rows]
                for row in ids.tolist():
                    raw = 2 * J[row] @ J[row].T
                    gram = raw - cross[row] @ solved[row].T
                    H = metric_scale * ynorm[row] * (gram + .001 * raw)
                    for predicted, delta in proposals(H.cpu().numpy(), gradient[row].cpu().numpy(), max_step=a.max_step):
                        choices.append((predicted, row, delta))
                choices.sort(key=lambda x:x[0])
            else:
                for index in order:
                    row, col = divmod(index, T)
                    delta = np.zeros(T)
                    delta[col] = -1 if float(gradient[row, col]) > 0 else 1
                    choices.append((-float(gradient[row, col].abs()), row, delta))
            seen = {}
            trials = 0
            accepted = 0
            for predicted, row, delta in choices:
                if seen.get(row, 0) >= 2:
                    continue
                seen[row] = seen.get(row, 0) + 1
                delta = torch.tensor(delta, **d)
                oldoff = off.clone()
                theta[row] += delta
                nxt = state()
                trials += 1
                if float(nxt[0]) < value - 1e-12:
                    value = float(nxt[0])
                    st = nxt
                    accepted += 1
                    total += 1
                else:
                    theta[row] -= delta
                    off.copy_(oldoff)
                if trials >= a.trials or accepted >= a.moves:
                    break
            st = state()
            value = float(st[0])
            er = float(st[1].clamp(min=0).sqrt())
            gap = float(st[2])
            score = float(st[4].clamp(min=0).sqrt())
            if value < bestloss:
                bestloss = value
                save(str(Path(a.output).with_suffix('')) + '_bestloss.npz', st)
            if er <= .019999 and gap >= 3.9151 and score < best:
                best = score
                save(a.output, st)
            print(json.dumps(dict(step=step, error_2=er, gap=gap, Lambda=float(st[3]), lambda_eff=score,
                                  loss=value, best_valid=best, accepted=accepted, total_moves=total,
                                  trials=trials, seconds=time.time()-began)), flush=True)
            if not accepted:
                break
    print(json.dumps(dict(status='FINISHED', scope='projected-column-gauge integer-angle research',
                          best_valid=best, best_loss=bestloss, seconds=time.time()-began)), flush=True)


if __name__ == '__main__':
    main()
