"""Exact Greedy Backward Elimination of incumbent 719 8-bit atoms on GPU."""

import time, torch, numpy as np
import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
    str(ROOT / "initial_experiments/3"),
]
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors

device = torch.device("cuda")
sp = np.load(ROOT / "initial_experiments/3/parameters/reference_spectrum.npz")
ii = torch.tensor(sp["ii"], device=device)
jj = torch.tensor(sp["jj"], device=device)
fac = torch.tensor(sp["fac"], dtype=torch.float64, device=device)
ev = torch.tensor(sp["eigenvalues"], dtype=torch.float64, device=device)
V = torch.tensor(sp["eigenvectors"], dtype=torch.float64, device=device)
G = (V * torch.clamp(ev, min=0.0)) @ V.T
total_norm_sq = torch.sum(torch.clamp(ev, min=0.0)**2).item()

sub = np.load("dfthcsos/submission.npz")
tree = sub["tree_edges"]
words = sub["angle_words"].reshape(-1, 75)
words_2body = words[76:]
U_q = reconstruct_u_vectors(words_2body / 128.0, tree)
all_words = list(words_2body)
all_U = list(torch.tensor(U_q / np.linalg.norm(U_q, axis=1, keepdims=True), dtype=torch.float64, device=device))

def eval_floor_exact(U_list):
    U_curr = torch.stack(U_list, dim=0)
    P = U_curr[:, ii] * U_curr[:, jj] * fac
    S = P @ P.T + 1e-6 * torch.eye(len(U_list), dtype=torch.float64, device=device)
    se, sv = torch.linalg.eigh(S)
    pos = se > 1e-10
    iroot = (sv[:, pos] / torch.sqrt(se[pos])) @ sv[:, pos].T
    Q = iroot @ P
    H = Q @ G @ Q.T
    E = torch.sum(H**2).item()
    floor = np.sqrt(max(0.0, total_norm_sq - E))
    return floor, P, S, Q, H

floor, P, S, Q, H = eval_floor_exact(all_U)
print(f"Starting B={len(all_words)}: Floor = {floor*1000:.4f} mHa")

t0 = time.time()
target_B = 680

# At each step, evaluate Sherman-Morrison LOO proxy for all current atoms,
# select top 8 candidates with smallest proxy loss, and pick the one that gives smallest exact floor.
while len(all_U) > target_B:
    B = len(all_U)
    Si = torch.linalg.inv(S)
    cap = Si @ H
    K = cap @ Si
    d = Si.diag()
    approx_loss = (2.0 * torch.sum(cap * K.T, dim=1) - (torch.sum(K * cap.T, dim=1))) / d
    cands = torch.argsort(approx_loss)[:8].cpu().numpy()
    
    best_cand = None
    best_cand_floor = 1e9
    for c in cands:
        keep = [i for i in range(B) if i != c]
        Uk = [all_U[i] for i in keep]
        fk, _, _, _, _ = eval_floor_exact(Uk)
        if fk < best_cand_floor:
            best_cand_floor = fk
            best_cand = c
            
    # Remove best_cand
    keep = [i for i in range(B) if i != best_cand]
    all_words = [all_words[i] for i in keep]
    all_U = [all_U[i] for i in keep]
    floor, P, S, Q, H = eval_floor_exact(all_U)
    
    if len(all_U) % 5 == 0 or len(all_U) <= target_B + 2:
        print(f"B={len(all_U):3d}: Floor = {floor*1000:7.4f} mHa | Elapsed: {time.time()-t0:.1f}s")
