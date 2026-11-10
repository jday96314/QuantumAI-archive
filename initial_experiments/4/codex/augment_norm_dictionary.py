"""Column-generation probes for the chemistry dictionary's weighted L1 norm.

For each factor column, S^-1 sign(w) gives a matrix dual to the current
rank-one dictionary. Its extreme eigenvectors identify violated atomic-norm
constraints. Quantize those directions, score the actual decoded atoms, and
append them with zero weights, preserving the initial physical model exactly.
All output is research data; no submission or verifier is modified.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import torch

from angle_gauge_probe import encode
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--source', default='candidates/jacobi_l1/submission.npz')
    p.add_argument('--output', default='parameters/norm_augmented735.npz')
    p.add_argument('--extra', type=int, default=16)
    a = p.parse_args()
    torch.set_num_threads(2)
    z = dict(np.load(a.source))
    u, w = z['U'][0], z['W'][0]
    tree, words = z['tree_edges'], z['angle_words'][76:]
    W = z['H2_I'][0] - w.sum(0)
    ii, jj = np.tril_indices(76)
    fac = np.where(ii == jj, 1., np.sqrt(2.))
    P = u[:, ii] * u[:, jj] * fac
    d = dict(device='cuda', dtype=torch.float64)
    pt, wt = torch.tensor(P, **d), torch.tensor(w, **d)
    S = pt @ pt.T
    dual = torch.linalg.solve(S, wt.sign())
    packed = (dual.T @ pt).cpu().numpy() / fac
    matrices = np.zeros((w.shape[1], 76, 76))
    matrices[:, ii, jj] = packed
    matrices += matrices.transpose(0, 2, 1).copy()
    matrices[:, np.arange(76), np.arange(76)] *= .5
    ev, vec = torch.linalg.eigh(torch.tensor(matrices, **d))
    # Both signs of the extremal eigenvalue can violate the atomic constraint.
    pool = vec[:, :, [0, -1]].permute(0, 2, 1).reshape(-1, 76).cpu().numpy()
    pool_words = np.unique(encode(pool, tree, 8, True), axis=0)
    pool_u = reconstruct_u_vectors(pool_words / 128., tree)
    cp = torch.tensor(pool_u[:, ii] * pool_u[:, jj] * fac, **d)
    cross = pt @ cp.T
    coordinates = torch.linalg.solve(S, cross)
    residual = (cp.square().sum(1) - (cross * coordinates).sum(0)).clamp(min=1e-14)
    certificate = coordinates.T @ wt.sign()
    l1 = wt.abs().sum(0) + torch.tensor(abs(z['H2_I'][0]), **d)
    violations = (certificate.abs() - 1).clamp(min=0) * l1
    # Relative to a fixed physical perturbation budget, the tangent gain is
    # controlled by the component outside the existing dictionary span.
    score = violations.square().sum(1) / residual
    overlap = (torch.tensor(pool_u, **d) @ torch.tensor(pool_u, **d).T).square()
    chosen = []
    rows = []
    for _ in range(a.extra):
        j = int(score.argmax())
        chosen.append(j)
        rows.append(dict(index=j, score=float(score[j]),
                         outside_squared_norm=float(residual[j]),
                         largest_dual_violation=float(certificate[j].abs().max())))
        # Avoid spending most slots on near-identical column witnesses.
        score *= (1 - overlap[j]).clamp(min=0)
        score[chosen] = -1
    out = dict(U=np.r_[u, pool_u[chosen]][None],
               w=np.r_[w, np.zeros((a.extra, w.shape[1]))][None], W=W[None],
               angle_words=np.r_[words, pool_words[chosen]], tree_edges=tree,
               brot=np.array(8), onebody_angle_words=z['angle_words'][:76],
               onebody_brot=np.array(8))
    output = Path(a.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez(output, **out)
    # The new rows have exactly zero weight, so these algebraic checks do not
    # depend on a cancellation-prone Frobenius residual calculation.
    assert np.array_equal(out['w'][0, :len(w)], w)
    assert not out['w'][0, len(w):].any()
    assert np.array_equal(out['U'][0, :len(u)], u)
    assert np.array_equal(out['W'][0], W)
    record = dict(status='PASS', scope='zero-weight extension preserves initial model; research only',
                  source=a.source, output=str(output), basis=len(u) + a.extra,
                  candidate_count=len(pool_u), atoms=rows)
    Path('results/' + output.stem + '_build.json').write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record), flush=True)


if __name__ == '__main__':
    main()
