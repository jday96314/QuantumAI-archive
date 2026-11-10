#!/usr/bin/env python3
import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
]
from dfthcsos.verify import heuristic_costs

base = 3355599588.3155375
N = 76
lam = 48.90

results = []
for R in [1, 2, 4]:
    for B in range(64, 513, 16):
        for C in range(64, B + 1, 16):
            for brot in [8]:
                h = heuristic_costs(N, R, B, C, brot)
                J = 10
                K = brot
                A = N + R * B
                sel = (A - 1).bit_length()
                ang_toff = 2 * K * (A - 2)
                T = ang_toff + h["toffoli"]
                Q = h["persistent_qubits"] + (J * brot) + max(h["temporary_qubits"], sel - 1)
                cost = lam * T * Q
                ratio = base / cost
                if 10.0 <= ratio <= 30.0:
                    results.append((ratio, cost, R, B, C, brot, T, Q))

results.sort(key=lambda x: -x[0])
print(f"Found {len(results)} configurations with ratio between 10x and 30x:")
print(f"{'Ratio':7s} | {'Cost':13s} | {'R':2s} {'B':4s} {'C':4s} {'brot':4s} | {'T':6s} | {'Q':4s}")
print("-" * 55)
for r, cost, R, B, C, brot, T, Q in results[:25]:
    print(f"{r:6.2f}x | {cost:11,.0f} | {R:2d} {B:4d} {C:4d} {brot:4d} | {T:6d} | {Q:4d}")
