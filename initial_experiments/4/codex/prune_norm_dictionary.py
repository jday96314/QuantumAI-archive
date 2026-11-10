"""Prune a factorization using exact projected row-removal model errors.

Each deletion projects that row's physical rank-one feature onto all remaining
features, preserving factor rank and column gauge. Scores include the full
reference residual, not just the projection distance. This is research data.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import torch


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--source', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--basis', type=int, default=719)
    p.add_argument('--protect-from', type=int, help='retain original row IDs at or above this index')
    p.add_argument('--criterion', choices=['error', 'normalization'], default='error')
    a = p.parse_args()
    torch.set_num_threads(2)
    z = dict(np.load(a.source))
    sp = np.load('../../3/parameters/reference_spectrum.npz')
    u = z['U'][0]
    w = z['w'][0]
    words = z['angle_words'].reshape(-1, 75)
    d = dict(device='cuda', dtype=torch.float64)
    P = torch.tensor(u[:, sp['ii']] * u[:, sp['jj']] * sp['fac'], **d)
    G = torch.tensor((sp['eigenvectors'] * sp['eigenvalues']) @ sp['eigenvectors'].T, **d)
    constant = G.square().sum()
    S = P @ P.T
    H = P @ G @ P.T
    inverse = torch.linalg.inv(S)
    w = torch.tensor(w, **d)
    active = list(range(len(u)))
    records = []

    def error(weights, gram, cross):
        K = weights @ weights.T
        KS = K @ gram
        return constant - 2 * (K * cross).sum() + (KS * KS.T).sum()

    initial_error = float(error(w, S, H).clamp(min=0).sqrt())
    while len(active) > a.basis:
        K = w @ w.T
        Y = K @ S @ K
        Z = inverse @ H @ K
        T = inverse @ H @ inverse
        dd = inverse.diag().reciprocal()
        qq = K.diag()
        increment = -2 * dd * Y.diag() + 4 * dd * Z.diag() + dd.square() * (qq.square() - 2 * qq * T.diag())
        selection = increment.clone()
        if a.criterion == 'normalization':
            offset = torch.tensor(z['W'][0], **d)
            before = error(w, S, H)
            for lo in range(0, len(active), 16):
                inds = torch.arange(lo, min(lo+16, len(active)), device='cuda')
                # Keep a zero row at the removed index while scoring; this
                # gives exactly the same column L1 norm as deleting that row.
                changes = -inverse[:, inds].T / inverse.diag()[inds, None]
                ww = w[None] + changes[:, :, None] * w[inds, None, :]
                lam2 = .25 * (ww.abs().sum(1) + (offset[None] + ww.sum(1)).abs()).square().sum(1)
                ge = (before + increment[inds]) / .019995**2 - 1
                penalty = (torch.relu(.62 + 4 * ge).square() - .62**2) / 8
                selection[inds] = 2 * 3.9155 * lam2 / 3000 + penalty
        if a.protect_from is not None:
            protected = torch.tensor([j >= a.protect_from for j in active], device='cuda')
            selection[protected] = float('inf')
        j = int(selection.argmin())
        keep = torch.tensor([i for i in range(len(active)) if i != j], device='cuda')
        # P_j = -sum_{i != j} inverse[i,j]/inverse[j,j] P_i + r_j.
        coordinates = -inverse[keep, j] / inverse[j, j]
        nextw = w[keep] + coordinates[:, None] * w[j]
        nextS = S[keep][:, keep]
        nextH = H[keep][:, keep]
        before = error(w, S, H)
        after = error(nextw, nextS, nextH)
        # Independently evaluate the selected projected model; reject a sign
        # or transpose mistake before writing any candidate.
        discrepancy = abs(float(after - before - increment[j]))
        assert discrepancy < 2e-9, discrepancy
        records.append(dict(basis=len(active)-1, removed=active[j],
                            error_2=float(after.clamp(min=0).sqrt()),
                            predicted_squared_error_increment=float(increment[j]),
                            formula_discrepancy=discrepancy))
        inverse = inverse[keep][:, keep] - torch.outer(inverse[keep, j], inverse[j, keep]) / inverse[j, j]
        S, H, w = nextS, nextH, nextw
        active.pop(j)
        print(json.dumps(records[-1]), flush=True)
    out = dict(z)
    for key in ['onebody_angle_words', 'onebody_brot', 'error_1_optimized']:
        out.pop(key, None)
    out.update(U=u[active][None], w=w.cpu().numpy()[None], angle_words=words[active], err2=records[-1]['error_2'])
    np.savez(a.output, **out)
    record = dict(status='PASS', scope='classical model row pruning, not native verification',
                  source=a.source, output=a.output, criterion=a.criterion, initial_error_2=initial_error, rows=records)
    Path('results/' + Path(a.output).stem + '_prune.json').write_text(json.dumps(record, indent=2) + '\n')


if __name__ == '__main__':
    main()
