"""Truncate each ensemble to its leading vector-space principal components."""
import argparse
from pathlib import Path

import numpy as np

from dfthcsos.verify import load_inputs
from dfthcsos.scoreinfra.dfthc_math import (
    build_h1_prime, reconstruct_u_vectors, reconstruct_x,
)


def encode(vectors, tree, brot):
    v = vectors.copy()
    th = np.zeros((len(v), len(tree)))
    for t in reversed(range(len(tree))):
        p, c = tree[t]
        th[:, t] = np.arctan2(v[:, c], v[:, p])
        v[:, p] = np.hypot(v[:, p], v[:, c])
        v[:, c] = 0
    return np.mod(np.rint(th / np.pi * (1 << (brot - 1))).astype(np.int64), 1 << brot)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--columns', type=int, default=50)
    p.add_argument('--brot', type=int, default=12)
    a = p.parse_args()
    z = dict(np.load(a.source))
    old = load_inputs(a.source)
    R, B, C = old.w.shape
    n = old.num_orbitals
    if not 0 < a.columns <= C:
        raise ValueError('column count must lie between one and the source count')
    tree = z['tree_edges']
    words = z['angle_words'].copy()
    U = reconstruct_u_vectors(words[n:] / (1 << (a.brot - 1)), tree).reshape(R, B, n)
    Xold = reconstruct_x(U, old.w).reshape(R, C, n * n)
    rotations = np.empty((R, C, a.columns))
    for r in range(R):
        eigenvalues, eigenvectors = np.linalg.eigh(Xold[r] @ Xold[r].T)
        rotations[r] = eigenvectors[:, np.argsort(eigenvalues)[-a.columns:][::-1]]
    w = np.einsum('rbc,rcd->rbd', old.w, rotations)
    W = np.einsum('rc,rcd->rd', old.W, rotations)
    X = reconstruct_x(U, w)
    h = build_h1_prime(old.ref['h1_thc'], old.S, old.beta_1,
                       old.ref['electrons'], W, X)
    _, eigenvectors = np.linalg.eigh(0.5 * (h + h.T))
    words[:n] = encode(eigenvectors.T, tree, a.brot)
    z.update(W=w, H2_I=W + w.sum(1), angle_words=words)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez(a.output, **z)
    print(f'R={R} B={B} C={a.columns}')


if __name__ == '__main__':
    main()
