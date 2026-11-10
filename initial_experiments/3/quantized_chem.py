"""Fit DF-THC weights and Givens angles with quantization in the forward pass.

Runs with the torch environment. Export chem_reference.npz using the host
reader first; verification remains separate and uses the unchanged verifier.
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
import torch


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--src', required=True)
    p.add_argument('--ref', default='parameters/chem_reference.npz')
    p.add_argument('--out', required=True)
    p.add_argument('--brot', type=int, default=11)
    p.add_argument('--src-brot', type=int, default=12)
    p.add_argument('--steps', type=int, default=4000)
    p.add_argument('--warmup', type=int, default=500)
    p.add_argument('--angle-lr', type=float, default=0.015)
    p.add_argument('--weight-lr', type=float, default=1e-4)
    a = p.parse_args()
    torch.set_num_threads(4)
    torch.backends.cuda.matmul.allow_tf32 = False
    z = dict(np.load(a.src))
    ref = np.load(a.ref)
    R, B, C = z['W'].shape
    n = 76
    edges = z['tree_edges']
    # Parallel disjoint rotations at each depth of the supplied tree.
    stages = []
    for t, (pa, ch) in enumerate(edges):
        pa, ch = int(pa), int(ch)
        if not stages or any(pa in (p0, c0) or ch in (p0, c0)
                             for _, p0, c0 in stages[-1]):
            stages.append([])
        stages[-1].append((t, pa, ch))
    stages = [tuple(torch.tensor(v, device='cuda') for v in zip(*st)) for st in stages]
    ii, jj = np.tril_indices(n)
    ix = torch.tensor(ii, device='cuda')
    jx = torch.tensor(jj, device='cuda')
    mul = torch.tensor(np.where(ii == jj, 1, np.sqrt(2)), device='cuda', dtype=torch.float32)
    flat = ii * n + jj
    gt = ref['G'][flat[:, None], flat[None, :]]
    G = torch.tensor(gt, device='cuda') * mul[:, None] * mul[None, :]
    words = z['angle_words'][n:] * (2. ** (a.brot - a.src_brot))
    theta = torch.tensor(words, device='cuda', dtype=torch.float32, requires_grad=True)
    w = torch.tensor(z['W'], device='cuda', dtype=torch.float32, requires_grad=True)
    eye0 = torch.zeros(R * B, n, device='cuda')
    eye0[:, 0] = 1

    def decode(t):
        rounded = t + (t.round() - t).detach()
        angles = rounded * (np.pi / (1 << (a.brot - 1)))
        u = eye0
        for ts, ps, cs in stages:
            ang = angles[:, ts]
            co, sn = ang.cos(), ang.sin()
            up, uc = u[:, ps], u[:, cs]
            nxt = u.clone()
            nxt[:, ps] = co * up - sn * uc
            nxt[:, cs] = sn * up + co * uc
            u = nxt
        return u

    def loss():
        u = decode(theta)
        P = ((u[:, ix] * u[:, jx]) * mul).reshape(R, B, -1)
        X = torch.bmm(w.transpose(1, 2), P).reshape(R * C, -1)
        D = X.T @ X - G
        return D.square().sum()

    opt = torch.optim.Adam([{'params': [theta], 'lr': a.angle_lr},
                            {'params': [w], 'lr': a.weight_lr}])
    best = float('inf')
    t0 = time.time()
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    for s in range(a.steps):
        opt.zero_grad()
        e2 = loss()
        val = float(e2.detach().sqrt())
        if val < best:
            best = val
            # Snapshot the exact weights and words that produced this loss.
            zz = dict(z)
            ww = w.detach().cpu().numpy().astype(np.float64)
            old_offset = z['H2_I'] - z['W'].sum(1)
            zz['W'] = ww
            zz['H2_I'] = old_offset + ww.sum(1)
            zz['U'] = decode(theta).detach().cpu().numpy().reshape(R, B, n).astype(np.float64)
            zz['angle_words'] = np.concatenate([
                np.rint(z['angle_words'][:n] * 2. ** (a.brot - a.src_brot)).astype(np.int64) % (1 << a.brot),
                theta.detach().round().cpu().numpy().astype(np.int64) % (1 << a.brot)])
            zz['brot'] = np.array(a.brot)
            np.savez(out, **zz)
        if s % 100 == 0 or s == a.steps - 1:
            print(json.dumps(dict(step=s, error_2=val, best=best,
                                  seconds=round(time.time()-t0, 1))), flush=True)
        e2.backward()
        opt.param_groups[0]['lr'] = 0 if s < a.warmup else a.angle_lr * (0.2 + 0.8 * (1 - s/a.steps))
        opt.param_groups[1]['lr'] = a.weight_lr * (0.1 + 0.9 * (1 - s/a.steps))
        opt.step()


if __name__ == '__main__':
    main()
