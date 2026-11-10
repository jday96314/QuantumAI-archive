#!/usr/bin/env python3
import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
]
from dfthcsos.verify import heuristic_costs
import math

base = 3355599588.3155375
N = 76
lam = 48.90

candidates = []
for R in [1, 2, 3, 4]:
    for C in range(75, 400, 25):
        if R * C < 300: # Need total rank >= 300 for err2 <= 20 mHa
            continue
        for B in range(C, 450, 25):
            for brot in [8]:
                h = heuristic_costs(N, R, B, C, brot)
                J = 10
                K = brot
                A = N + R * B
                sel = math.ceil(math.log2(A))
                ang_toff = 2 * K * (A - 2)
                T = ang_toff + h["toffoli"]
                Q = h["persistent_qubits"] + (J * brot) + max(h["temporary_qubits"], sel - 1)
                cost = lam * T * Q
                ratio = base / cost
                candidates.append((ratio, cost, R, B, C, brot, T, Q))

candidates.sort(key=lambda x: -x[0])
print(f"Top 20 viable configurations with R*C >= 300:")
print(f"{'Ratio':7s} | {'Cost':13s} | {'R':2s} {'B':4s} {'C':4s} {'brot':4s} | {'T':6s} | {'Q':4s}")
print("-" * 55)
for r, cost, R, B, C, brot, T, Q in candidates[:20]:
    print(f"{r:6.2f}x | {cost:11,.0f} | {R:2d} {B:4d} {C:4d} {brot:4d} | {T:6d} | {Q:4d}")
