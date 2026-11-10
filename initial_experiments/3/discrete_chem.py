"""Greedy quantized-angle refinement of the two-body tensor residual."""
import argparse
import json
import time
from pathlib import Path

import numpy as np
import torch


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--src', required=True)
    p.add_argument('--ref', default='initial_experiments/3/parameters/chem_reference.npz')
    p.add_argument('--out', required=True)
    p.add_argument('--brot', type=int, default=12)
    p.add_argument('--steps', type=int, default=3000)
    p.add_argument('--block', type=int, default=64)
    a = p.parse_args()
    torch.set_num_threads(4)
    torch.backends.cuda.matmul.allow_tf32 = False
    z = dict(np.load(a.src))
    ref = np.load(a.ref)
    R, B, C = z['W'].shape
    N = 76
    words = z['angle_words'][N:].astype(np.int32)
    stages = []
    for t, (pa, ch) in enumerate(z['tree_edges']):
        pa, ch = int(pa), int(ch)
        if not stages or any(pa in (p0, c0) or ch in (p0, c0)
                             for _, p0, c0 in stages[-1]):
            stages.append([])
        stages[-1].append((t, pa, ch))
    stages = [tuple(torch.tensor(v, device='cuda') for v in zip(*st)) for st in stages]
    ii, jj = np.tril_indices(N)
    ix = torch.tensor(ii, device='cuda')
    jx = torch.tensor(jj, device='cuda')
    mul = torch.tensor(np.where(ii == jj, 1, np.sqrt(2)), device='cuda', dtype=torch.float32)
    flat = ii * N + jj
    gt = ref['G'][flat[:, None], flat[None, :]]
    G = torch.tensor(gt, device='cuda') * mul[:, None] * mul[None, :]
    w = torch.tensor(z['W'], device='cuda', dtype=torch.float32)
    eye0 = torch.zeros(R * B, N, device='cuda')
    eye0[:, 0] = 1

    def decode(t):
        u = eye0
        angles = t * (np.pi / (1 << (a.brot - 1)))
        for ts, ps, cs in stages:
            ang = angles[:, ts]
            co, sn = ang.cos(), ang.sin()
            up, uc = u[:, ps], u[:, cs]
            nxt = u.clone()
            nxt[:, ps] = co * up - sn * uc
            nxt[:, cs] = sn * up + co * uc
            u = nxt
        return u

    def residual(t):
        u = decode(t)
        P = ((u[:, ix] * u[:, jx]) * mul).reshape(R, B, -1)
        X = torch.bmm(w.transpose(1, 2), P).reshape(R * C, -1)
        D = X.T @ X - G
        return D.square().sum()

    current = torch.tensor(words, device='cuda', dtype=torch.float32)
    best = float(residual(current).sqrt())
    print(json.dumps(dict(step=0, best=best)), flush=True)
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    def save():
        zz = dict(z)
        zz['angle_words'] = np.concatenate([z['angle_words'][:N],
                                            current.cpu().numpy().astype(np.int64) % (1 << a.brot)])
        zz['U'] = decode(current).detach().cpu().numpy().reshape(R, B, N).astype(np.float64)
        np.savez(out, **zz)

    save()
    block = a.block
    start = time.time()
    stagnant = 0
    for step in range(1, a.steps + 1):
        t = current.detach().requires_grad_(True)
        val = residual(t)
        grad = torch.autograd.grad(val, t)[0]
        # Roundoff near the periodic boundary does not alter Givens rotations.
        flatg = grad.flatten()
        k = min(block, len(flatg))
        ids = flatg.abs().topk(k).indices
        trial = current.clone().flatten()
        trial[ids] -= flatg[ids].sign()
        trial = trial.reshape_as(current)
        with torch.no_grad():
            score = float(residual(trial).sqrt())
        if score < best - 1e-9:
            current = trial.detach()
            best = score
            stagnant = 0
            if block < a.block:
                block = min(a.block, block * 2)
            save()
        else:
            stagnant += 1
            block = max(1, block // 2)
        if step % 50 == 0 or step == a.steps:
            print(json.dumps(dict(step=step, best=best, block=block,
                                  stagnant=stagnant, seconds=round(time.time()-start, 1))), flush=True)
        if stagnant >= 100:
            break


if __name__ == '__main__':
    main()
