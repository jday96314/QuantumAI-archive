"""THC-structured (R=1) refit of the FeMoCo two-body tensor on GPU.

Target: sum_c X_c (x) X_c ~= tei + (S (x) I + I (x) S)/2, X_c = sum_b w_bc u_b u_b^T.
Pipeline: DF eigen-factors Y_c of the packed target -> variable-projection fit of
shared unit vectors u_b (w eliminated by ridge least squares) -> full quartic
L-BFGS fine-tune of (u, w). Reports error_2 (Frobenius, Hartree) per (B, C).

Usage: python thc_fit_gpu.py --B 256 --C 320 [--out file.npz]
"""

from __future__ import annotations

import argparse
import time

import numpy as np
import torch

N = 76
dev = torch.device("cuda")
torch.backends.cuda.matmul.allow_tf32 = False


def load_target():
    cache = np.load("/tmp/femoco_ref_cache.npz", allow_pickle=False)
    tei = np.asarray(cache["tei"], dtype=np.float64)
    S = np.asarray(cache["S"], dtype=np.float64)
    rows, cols = np.tril_indices(N)
    wt = np.where(rows == cols, 1.0, np.sqrt(2.0))
    h2 = tei.copy()
    d = np.arange(N)
    h2[:, :, d, d] += 0.5 * S[:, :, None]
    h2[d, d, :, :] += 0.5 * S[None, :, :]
    M = h2[rows[:, None], cols[:, None], rows[None, :], cols[None, :]] * (wt[:, None] * wt[None, :])
    return (torch.tensor(M, dtype=torch.float64, device=dev),
            torch.tensor(rows, device=dev), torch.tensor(cols, device=dev),
            torch.tensor(wt, dtype=torch.float64, device=dev))


def packed(u, rows, cols, wt):
    un = u / torch.norm(u, dim=-1, keepdim=True)
    return un[..., rows] * un[..., cols] * wt


def baseline_init(B: int, seed: int) -> np.ndarray:
    sub = np.load("/media/shared/ExpansionDrive2/QuantumAI/from_host/sample_submission/dfthcsos/submission.npz")
    U = np.asarray(sub["U"], dtype=np.float64).reshape(-1, N)
    U /= np.linalg.norm(U, axis=1, keepdims=True)
    w = np.asarray(sub["W"], dtype=np.float64).reshape(-1, sub["W"].shape[2])
    imp = np.sum(w * w, axis=1)
    order = np.argsort(-imp)
    chosen = []
    for i in order:
        v = U[i]
        if all(abs(v @ U[j]) < 0.95 for j in chosen):
            chosen.append(i)
        if len(chosen) == B:
            break
    rng = np.random.default_rng(seed)
    out = U[chosen]
    if len(chosen) < B:
        out = np.vstack([out, rng.normal(size=(B - len(chosen), N))])
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--B", type=int, default=256)
    ap.add_argument("--C", type=int, default=320)
    ap.add_argument("--vp-iters", type=int, default=300)
    ap.add_argument("--ft-iters", type=int, default=600)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=str, default="")
    a = ap.parse_args()
    t0 = time.time()
    M, rows, cols, wt = load_target()
    ev, evec = torch.linalg.eigh(M)
    idx = torch.argsort(ev, descending=True)
    ev, evec = ev[idx], evec[:, idx]
    print(f"neg eig mass: {torch.sqrt(torch.sum(ev[ev < 0] ** 2)).item():.5f}; "
          f"DF tail@C={a.C}: {torch.sqrt(torch.sum(ev[a.C:] ** 2)).item():.5f}", flush=True)
    Y = (evec[:, : a.C] * torch.sqrt(torch.clamp(ev[: a.C], min=0))).T  # (C, P)

    u = torch.tensor(baseline_init(a.B, a.seed), dtype=torch.float64, device=dev).requires_grad_(True)

    def solve_w(UB, ridge=1e-8):
        G = UB @ UB.T
        G = G + ridge * torch.trace(G) / G.shape[0] * torch.eye(G.shape[0], dtype=G.dtype, device=dev)
        return torch.linalg.solve(G, UB @ Y.T)  # (B, C)

    opt = torch.optim.LBFGS([u], lr=1.0, max_iter=a.vp_iters, history_size=50, line_search_fn="strong_wolfe")

    def vp_closure():
        opt.zero_grad()
        UB = packed(u, rows, cols, wt)
        w = solve_w(UB)
        loss = torch.sum((w.T @ UB - Y) ** 2)
        loss.backward()
        return loss

    opt.step(vp_closure)
    with torch.no_grad():
        UB = packed(u, rows, cols, wt)
        w0 = solve_w(UB)
        X = w0.T @ UB
        e2 = torch.norm(X.T @ X - M).item()
    print(f"[VP] B={a.B} C={a.C}: factor resid={vp_closure().sqrt().item():.5f}, error_2={e2:.5f} ({time.time()-t0:.0f}s)", flush=True)

    w = w0.clone().detach().requires_grad_(True)
    u2 = u.detach().clone().requires_grad_(True)
    opt2 = torch.optim.LBFGS([u2, w], lr=1.0, max_iter=a.ft_iters, history_size=50, line_search_fn="strong_wolfe")

    def ft_closure():
        opt2.zero_grad()
        X = w.T @ packed(u2, rows, cols, wt)
        loss = torch.sum((X.T @ X - M) ** 2)
        loss.backward()
        return loss

    opt2.step(ft_closure)
    e2 = ft_closure().sqrt().item()
    with torch.no_grad():
        l1 = torch.sum(torch.abs(w), dim=0)
        lam2_lower = 0.25 * torch.sum((l1 + torch.abs(torch.sum(w, dim=0))) ** 2).item()
    print(f"[FT] B={a.B} C={a.C}: error_2={e2:.5f}, lambda_h2(W=0)={lam2_lower:.2f} ({time.time()-t0:.0f}s)", flush=True)
    if a.out:
        un = (u2 / torch.norm(u2, dim=-1, keepdim=True)).detach().cpu().numpy()
        np.savez(a.out, u=un, w=w.detach().cpu().numpy(), error_2=e2)


if __name__ == "__main__":
    main()
