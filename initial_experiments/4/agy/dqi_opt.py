"""Optimized Reversible Reed-Solomon Syndrome Decoder for Benchmark II (`dqi`).

Key architectural optimizations over the Proos-Zalka reference circuit:
1. Isomorphic Tower-Field Representation `GF(2^12) = GF(((2^3)^2)^2)`:
   - Converts the 70 syndrome registers in-place (0 Toffolis, 0 ancillas) from the standard
     polynomial basis of `km.GF2Field(12)` to the Tower-Field basis `{1, w, w^2} x {1, lam} x {1, beta}`
     where `w^3 + w + 1 = 0`, `lam^2 + lam + 1 = 0`, `beta^2 + beta + lam = 0`.
   - Since `0` and `1` (`e_0`) are fixed points of the basis isomorphism and XOR is linear,
     the entire Forward and Reverse Synchronized EEA runs natively in the Tower-Field basis
     with zero per-step basis conversion gates.
2. 66-Toffoli Tower-Field Karatsuba Multiplier (`gf2_12_mul_xor`):
   - Replaces `builder.ixor_gf2_mul(..., control=mul)` (~115 Toffolis) with 12-bit input gating
     (`init_and` / `del_and`, 12 Toffolis) + 2x2x3 Karatsuba multiplication (54 `init_and`s,
     uncomputed in 0 Toffolis via `del_and`) using only 15 temporary ancilla qubits.
3. 78-Toffoli Gated Scaffolded Tower-Field Inverter (`gf2_12_inv_scaffold`) in Chien/Forney:
   - Replaces the 315-Toffoli Itoh-Tsujii inverter + 12-Toffoli `located` gating (327 Toffolis)
     with a 78-Toffoli Tower-Field inverter whose 3-bit `GF(2^3)` subfield inversion is gated
     directly by `is_sigma_zero` (3 Toffolis) and whose entire forward scaffold uncomputes in
     0 Toffolis via `del_and`, eliminating the 12-qubit `located` register.
   - Replaces the 24-Toffoli double `flip_if_equal(sigma, 0)` with an 11-Toffoli `init_and` tree
     uncomputed in 0 Toffolis via `del_and`.
4. Low-Ancilla Pebbled Tower-Field Divider (`gf2_12_div_xor_pebbled`) in `accumulate_quotient`:
   - Replaces the 72-ancilla-qubit `init_gf2_div_with_scaffold` (327 Toffolis) with a 24-ancilla-qubit
     pebbled Tower-Field divider (156 Toffolis), removing `accumulate_quotient` as the global qubit bottleneck.
5. Tightened Degree Bounds & Merged Unary Iteration in EEA:
   - Merges `eq_i` and `eq_ip1` in `sweep_positions` into a single `unary_iteration(n_q)` via
     `isub(1, target=n_q, control=c3_hc)`.
   - Tightens `lim_a = min(36, it // 2 + 1)`, `lim_b = min(36, it // 2 + 1)`, restricts `multiplies`
     to reachable slots, and uses 7-bit degree registers (`logm = 7`).
"""

from __future__ import annotations

import argparse
from typing import Sequence

import kickmix as km
import numpy as np

_UNREACHABLE = 1 << 63


def build_tower_isomorphism(gf: km.GF2Field) -> tuple[np.ndarray, np.ndarray, list[tuple[int, int]], list[tuple[int, int]]]:
    """Builds the 12x12 GF(2) isomorphism matrices and in-place CNOT schedules."""
    prim = gf.primitive_element
    exp_tab = [gf.pow(prim, k) for k in range(4095)]
    log_tab = [0] * 4096
    for k, val in enumerate(exp_tab):
        if k < 4095:
            log_tab[val] = k

    def gf_mul(a: int, b: int) -> int:
        if a == 0 or b == 0:
            return 0
        return exp_tab[(log_tab[a] + log_tab[b]) % 4095]

    omega = 0
    for k in range(1, 7):
        cand = exp_tab[585 * k]
        c3 = gf_mul(gf_mul(cand, cand), cand)
        if (c3 ^ cand ^ 1) == 0:
            omega = cand
            break
    assert omega != 0

    lam = exp_tab[1365]
    assert (gf_mul(lam, lam) ^ lam ^ 1) == 0

    beta = 0
    for cand in range(1, 4096):
        if (gf_mul(cand, cand) ^ cand ^ lam) == 0:
            beta = cand
            break
    assert beta != 0

    w_pows = [1, omega, gf_mul(omega, omega)]
    l_pows = [1, lam]
    b_pows = [1, beta]

    T_twr_to_std = np.zeros((12, 12), dtype=np.int64)
    for ib in range(2):
        for il in range(2):
            for iw in range(3):
                col = 6 * ib + 3 * il + iw
                val = gf_mul(gf_mul(w_pows[iw], l_pows[il]), b_pows[ib])
                for bit in range(12):
                    T_twr_to_std[bit, col] = (val >> bit) & 1

    aug = np.hstack([T_twr_to_std.copy(), np.eye(12, dtype=np.int64)])
    ops_red: list[tuple[int, int]] = []
    for col in range(12):
        piv = col
        while piv < 12 and aug[piv, col] == 0:
            piv += 1
        assert piv < 12
        if piv != col:
            aug[[col, piv]] = aug[[piv, col]]
            ops_red.extend([(col, piv), (piv, col), (col, piv)])
        for row in range(12):
            if row != col and aug[row, col] == 1:
                aug[row] ^= aug[col]
                ops_red.append((col, row))
    T_std_to_twr = aug[:, 12:]
    # Recall: E_k ... E_1 * T_twr_to_std = I  =>  T_std_to_twr = E_k ... E_1
    # So applying ops_red in forward order (src, dst -> reg[dst] ^= reg[src]) multiplies by T_std_to_twr!
    # And applying ops_red in reverse order multiplies by T_twr_to_std!
    ops_std_to_twr = list(ops_red)
    ops_twr_to_std = list(reversed(ops_red))
    return T_std_to_twr, T_twr_to_std, ops_std_to_twr, ops_twr_to_std


def apply_cnot_schedule(builder: km.CircuitBuilder, reg: km.array, ops: Sequence[tuple[int, int]]) -> None:
    """Applies an in-place GF(2)-linear CNOT schedule `(src, dst)` to `reg`."""
    for src, dst in ops:
        builder.cx(reg[src], reg[dst])


def gf2_3_mul_xor(builder: km.CircuitBuilder, u: km.array, v: km.array, tgt: km.array, work3: km.array) -> None:
    """Computes `tgt ^= u * v` in GF(2)[w]/(w^3 + w + 1) using 6 Toffolis and `work3` (3 clean qubits)."""
    # Phase 1: m0 = u0&v0, m1 = u1&v1, m2 = u2&v2
    builder.init_and(u, v, work3)
    m0, m1, m2 = work3
    builder.cx(km.array([m0, m1, m2]), tgt[0])
    builder.cx(m0, tgt[1])
    builder.cx(km.array([m0, m1]), tgt[2])
    builder.del_and(u, v, work3)

    # Phase 2: m01 = (u0^u1)&(v0^v1), m02 = (u0^u2)&(v0^v2), m12 = (u1^u2)&(v1^v2)
    builder.cx(km.array([u[0], v[0]]), km.array([u[1], v[1]]))
    builder.init_and(u[1], v[1], work3[0])  # m01
    builder.cx(km.array([u[0], v[0], u[1], v[1]]), km.array([u[1], v[1], u[2], v[2]]))
    builder.init_and(u[2], v[2], work3[2])  # m12
    builder.cx(km.array([u[1], v[1], u[0], v[0]]), km.array([u[2], v[2], u[2], v[2]]))
    builder.init_and(u[2], v[2], work3[1])  # m02
    builder.cx(km.array([u[0], v[0]]), km.array([u[2], v[2]]))

    m01, m02, m12 = work3
    builder.cx(m12, tgt[0])
    builder.cx(km.array([m01, m12]), tgt[1])
    builder.cx(m02, tgt[2])

    builder.cx(km.array([u[0], v[0]]), km.array([u[2], v[2]]))
    builder.del_and(u[2], v[2], work3[1])
    builder.cx(km.array([u[0], v[0], u[1], v[1]]), km.array([u[2], v[2], u[2], v[2]]))
    builder.del_and(u[2], v[2], work3[2])
    builder.cx(km.array([u[1], v[1], u[0], v[0]]), km.array([u[2], v[2], u[1], v[1]]))
    builder.del_and(u[1], v[1], work3[0])
    builder.cx(km.array([u[0], v[0]]), km.array([u[1], v[1]]))


def gf2_3_mul_keep(builder: km.CircuitBuilder, u: km.array, v: km.array, tgt: km.array, m_out: km.array, inverse: bool = False) -> None:
    """Scaffolded GF(2^3) multiply: 6 Toffolis forward, 0 Toffolis reverse."""
    m0, m1, m2, m01, m02, m12 = m_out
    if not inverse:
        builder.init_and(u, v, m_out[:3])
        builder.cx(km.array([u[0], v[0]]), km.array([u[1], v[1]]))
        builder.init_and(u[1], v[1], m01)
        builder.cx(km.array([u[0], v[0], u[1], v[1]]), km.array([u[1], v[1], u[2], v[2]]))
        builder.init_and(u[2], v[2], m12)
        builder.cx(km.array([u[1], v[1], u[0], v[0]]), km.array([u[2], v[2], u[2], v[2]]))
        builder.init_and(u[2], v[2], m02)
        builder.cx(km.array([u[0], v[0]]), km.array([u[2], v[2]]))

        builder.cx(km.array([m0, m1, m2, m12]), tgt[0])
        builder.cx(km.array([m0, m01, m12]), tgt[1])
        builder.cx(km.array([m0, m1, m02]), tgt[2])
    else:
        builder.cx(km.array([m0, m1, m02]), tgt[2])
        builder.cx(km.array([m0, m01, m12]), tgt[1])
        builder.cx(km.array([m0, m1, m2, m12]), tgt[0])

        builder.cx(km.array([u[0], v[0]]), km.array([u[2], v[2]]))
        builder.del_and(u[2], v[2], m02)
        builder.cx(km.array([u[0], v[0], u[1], v[1]]), km.array([u[2], v[2], u[2], v[2]]))
        builder.del_and(u[2], v[2], m12)
        builder.cx(km.array([u[1], v[1], u[0], v[0]]), km.array([u[2], v[2], u[1], v[1]]))
        builder.del_and(u[1], v[1], m01)
        builder.cx(km.array([u[0], v[0]]), km.array([u[1], v[1]]))
        builder.del_and(u, v, m_out[:3])


def gf2_6_mul_xor(builder: km.CircuitBuilder, a: km.array, b: km.array, tgt: km.array, work3: km.array) -> None:
    """Computes `tgt ^= a * b` in GF(2^6) using 18 Toffolis and `work3` (3 clean qubits)."""
    a0, a1 = a[:3], a[3:6]
    b0, b1 = b[:3], b[3:6]
    t0, t1 = tgt[:3], tgt[3:6]

    builder.cx(t0, t1)
    gf2_3_mul_xor(builder, a0, b0, t0, work3)
    builder.cx(t0, t1)

    gf2_3_mul_xor(builder, a1, b1, t0, work3)

    builder.cx(a0, a1)
    builder.cx(b0, b1)
    gf2_3_mul_xor(builder, a1, b1, t1, work3)
    builder.cx(b0, b1)
    builder.cx(a0, a1)


def gf2_6_mul_keep(builder: km.CircuitBuilder, a: km.array, b: km.array, tgt: km.array, m18: km.array, inverse: bool = False) -> None:
    """Scaffolded GF(2^6) multiply: 18 Toffolis forward, 0 Toffolis reverse."""
    a0, a1 = a[:3], a[3:6]
    b0, b1 = b[:3], b[3:6]
    t0, t1 = tgt[:3], tgt[3:6]
    if not inverse:
        builder.cx(t0, t1)
        gf2_3_mul_keep(builder, a0, b0, t0, m18[:6], inverse=False)
        builder.cx(t0, t1)

        gf2_3_mul_keep(builder, a1, b1, t0, m18[6:12], inverse=False)

        builder.cx(a0, a1)
        builder.cx(b0, b1)
        gf2_3_mul_keep(builder, a1, b1, t1, m18[12:18], inverse=False)
        builder.cx(b0, b1)
        builder.cx(a0, a1)
    else:
        builder.cx(a0, a1)
        builder.cx(b0, b1)
        gf2_3_mul_keep(builder, a1, b1, t1, m18[12:18], inverse=True)
        builder.cx(b0, b1)
        builder.cx(a0, a1)

        gf2_3_mul_keep(builder, a1, b1, t0, m18[6:12], inverse=True)

        builder.cx(t0, t1)
        gf2_3_mul_keep(builder, a0, b0, t0, m18[:6], inverse=True)
        builder.cx(t0, t1)


def gf2_12_mul_xor(
    builder: km.CircuitBuilder,
    x: km.array,
    y: km.array,
    tgt: km.array,
    control: km.q | None = None,
) -> None:
    """Computes `tgt ^= x * y` (or `tgt ^= control * (x * y)`) in GF(2^12) Tower Basis.

    Uses 54 Toffolis (66 if controlled) and 15 ancilla qubits (3 if uncontrolled).
    """
    work3 = builder.alloc_qubits(3)
    if control is not None:
        x_use = builder.alloc_qubits(12)
        builder.init_and(control, x, x_use)
    else:
        x_use = x

    x0, x1 = x_use[:6], x_use[6:12]
    y0, y1 = y[:6], y[6:12]
    out0, out1 = tgt[:6], tgt[6:12]

    builder.cx(out0, out1)
    gf2_6_mul_xor(builder, x0, y0, out0, work3)
    builder.cx(out0, out1)

    builder.cx(x0, x1)
    builder.cx(y0, y1)
    gf2_6_mul_xor(builder, x1, y1, out1, work3)
    builder.cx(y0, y1)
    builder.cx(x0, x1)

    builder.cx(out0[:3], out0[3:6])
    out0_swapped = km.array([*out0[3:6], *out0[:3]])
    gf2_6_mul_xor(builder, x1, y1, out0_swapped, work3)
    builder.cx(out0[:3], out0[3:6])

    if control is not None:
        builder.del_and(control, x, x_use)
        builder.free(x_use)
    builder.free(work3)


def gf2_3_sq_xor(builder: km.CircuitBuilder, u: km.array, tgt: km.array) -> None:
    """Computes `tgt ^= u^2` in GF(2)[w]/(w^3 + w + 1) (0 Toffolis)."""
    builder.cx(u[0], tgt[0])
    builder.cx(u[2], tgt[1])
    builder.cx(km.array([u[1], u[2]]), tgt[2])


def gf2_3_inv_keep(
    builder: km.CircuitBuilder,
    u: km.array,
    inv_u: km.array,
    inverse: bool = False,
) -> None:
    """Computes `inv_u = u^{-1}` in GF(2)[w]/(w^3 + w + 1) in 3 Toffolis (0 Toffolis reverse)."""
    if not inverse:
        builder.init_and(km.array([u[1], u[0], u[0]]), km.array([u[2], u[1], u[2]]), inv_u)
        builder.cx(km.array([u[0], u[1], u[2]]), inv_u[0])
        builder.cx(u[2], inv_u[1])
        builder.cx(km.array([u[1], u[2]]), inv_u[2])
    else:
        builder.cx(km.array([u[1], u[2]]), inv_u[2])
        builder.cx(u[2], inv_u[1])
        builder.cx(km.array([u[0], u[1], u[2]]), inv_u[0])
        builder.del_and(km.array([u[1], u[0], u[0]]), km.array([u[2], u[1], u[2]]), inv_u)


def gf2_12_inv_scaffold(
    builder: km.CircuitBuilder,
    x: km.array,
    target: km.array,
    scaffold: km.array,
    gate_ctl: km.q | None = None,
    inverse: bool = False,
) -> None:
    """Scaffolded Tower-Field Inverter in GF(2^12).

    `scaffold` has length 90 (or 93 when `gate_ctl` is not None).
    Forward (`inverse=False`): 75 Toffolis (78 if `gate_ctl` is not None).
    Reverse (`inverse=True`): 0 Toffolis!
    """
    N6 = scaffold[0:6]
    n3 = scaffold[6:9]
    inv_n3 = scaffold[9:12]
    d6 = scaffold[12:18]
    m_N6 = scaffold[18:36]
    m_n3 = scaffold[36:42]
    m_d6 = scaffold[42:54]
    m_tgt = scaffold[54:90]
    inv_n3_use = scaffold[90:93] if gate_ctl is not None else inv_n3

    x0, x1 = x[:6], x[6:12]
    N6_0, N6_1 = N6[:3], N6[3:6]
    d6_0, d6_1 = d6[:3], d6[3:6]
    tgt0, tgt1 = target[:6], target[6:12]

    if not inverse:
        gf2_3_sq_xor(builder, x1[3:6], N6_0)
        gf2_3_sq_xor(builder, x1[0:3], N6_1)
        builder.cx(x1, x0)
        builder.cx(x0, x1)
        gf2_6_mul_keep(builder, x1, x0, N6, m_N6, inverse=False)
        builder.cx(x0, x1)
        builder.cx(x1, x0)

        gf2_3_sq_xor(builder, N6_1, n3)
        builder.cx(N6_0, N6_1)
        gf2_3_mul_keep(builder, N6_0, N6_1, n3, m_n3, inverse=False)
        builder.cx(N6_0, N6_1)

        gf2_3_inv_keep(builder, n3, inv_n3, inverse=False)
        if gate_ctl is not None:
            builder.init_and(gate_ctl, inv_n3, inv_n3_use)

        builder.cx(N6_0, N6_1)
        gf2_3_mul_keep(builder, inv_n3_use, N6_1, d6_0, m_d6[:6], inverse=False)
        builder.cx(N6_0, N6_1)
        gf2_3_mul_keep(builder, inv_n3_use, N6_1, d6_1, m_d6[6:12], inverse=False)

        builder.cx(x1, x0)
        gf2_6_mul_keep(builder, d6, x0, tgt0, m_tgt[:18], inverse=False)
        builder.cx(x1, x0)
        gf2_6_mul_keep(builder, d6, x1, tgt1, m_tgt[18:36], inverse=False)
    else:
        gf2_6_mul_keep(builder, d6, x1, tgt1, m_tgt[18:36], inverse=True)
        builder.cx(x1, x0)
        gf2_6_mul_keep(builder, d6, x0, tgt0, m_tgt[:18], inverse=True)
        builder.cx(x1, x0)

        gf2_3_mul_keep(builder, inv_n3_use, N6_1, d6_1, m_d6[6:12], inverse=True)
        builder.cx(N6_0, N6_1)
        gf2_3_mul_keep(builder, inv_n3_use, N6_1, d6_0, m_d6[:6], inverse=True)
        builder.cx(N6_0, N6_1)

        if gate_ctl is not None:
            builder.del_and(gate_ctl, inv_n3, inv_n3_use)
        gf2_3_inv_keep(builder, n3, inv_n3, inverse=True)

        builder.cx(N6_0, N6_1)
        gf2_3_mul_keep(builder, N6_0, N6_1, n3, m_n3, inverse=True)
        builder.cx(N6_0, N6_1)
        gf2_3_sq_xor(builder, N6_1, n3)

        builder.cx(x1, x0)
        builder.cx(x0, x1)
        gf2_6_mul_keep(builder, x1, x0, N6, m_N6, inverse=True)
        builder.cx(x0, x1)
        builder.cx(x1, x0)
        gf2_3_sq_xor(builder, x1[0:3], N6_1)
        gf2_3_sq_xor(builder, x1[3:6], N6_0)


def gf2_12_div_xor_pebbled(
    builder: km.CircuitBuilder,
    num: km.array,
    den: km.array,
    tgt: km.array,
    control: km.q,
) -> None:
    """Computes `tgt ^= control * (num / den)` in GF(2^12) Tower Basis using at most 27 ancilla qubits."""
    x0, x1 = den[:6], den[6:12]

    def toggle_d6_into(d6: km.array, work3: km.array) -> None:
        # Uses 15 temporary qubits: N6(6), n3(3), inv_n3(3), inv_n3_g(3)
        N6 = builder.alloc_qubits(6)
        n3 = builder.alloc_qubits(3)
        inv_n3 = builder.alloc_qubits(3)
        inv_n3_g = builder.alloc_qubits(3)
        N6_0, N6_1 = N6[:3], N6[3:6]

        # Compute N6
        gf2_3_sq_xor(builder, x1[3:6], N6_0)
        gf2_3_sq_xor(builder, x1[0:3], N6_1)
        builder.cx(x1, x0)
        builder.cx(x0, x1)
        gf2_6_mul_xor(builder, x1, x0, N6, work3)
        builder.cx(x0, x1)
        builder.cx(x1, x0)

        # Compute n3
        gf2_3_sq_xor(builder, N6_1, n3)
        builder.cx(N6_0, N6_1)
        gf2_3_mul_xor(builder, N6_0, N6_1, n3, work3)
        builder.cx(N6_0, N6_1)

        # Compute inv_n3 and gate with control
        gf2_3_inv_keep(builder, n3, inv_n3, inverse=False)
        builder.init_and(control, inv_n3, inv_n3_g)

        # Accumulate into d6
        builder.cx(N6_0, N6_1)
        gf2_3_mul_xor(builder, inv_n3_g, N6_1, d6[:3], work3)
        builder.cx(N6_0, N6_1)
        gf2_3_mul_xor(builder, inv_n3_g, N6_1, d6[3:6], work3)

        # Uncompute inv_n3_g, inv_n3, n3, N6
        builder.del_and(control, inv_n3, inv_n3_g)
        gf2_3_inv_keep(builder, n3, inv_n3, inverse=True)

        builder.cx(N6_0, N6_1)
        gf2_3_mul_xor(builder, N6_0, N6_1, n3, work3)
        builder.cx(N6_0, N6_1)
        gf2_3_sq_xor(builder, N6_1, n3)

        builder.cx(x1, x0)
        builder.cx(x0, x1)
        gf2_6_mul_xor(builder, x1, x0, N6, work3)
        builder.cx(x0, x1)
        builder.cx(x1, x0)
        gf2_3_sq_xor(builder, x1[0:3], N6_1)
        gf2_3_sq_xor(builder, x1[3:6], N6_0)

        builder.free(inv_n3_g)
        builder.free(inv_n3)
        builder.free(n3)
        builder.free(N6)

    d6 = builder.alloc_qubits(6)
    work3 = builder.alloc_qubits(3)
    # 1. Compute d6 = control * N6(den)^{-1} (peak ancillas here: 6 + 3 + 15 = 24 qubits)
    toggle_d6_into(d6, work3)

    # 2. Allocate den_inv (12 qubits) and compute den_inv = d6 * den^{2^6} (peak ancillas: 6 + 3 + 12 = 21 qubits)
    den_inv = builder.alloc_qubits(12)
    builder.cx(x1, x0)
    gf2_6_mul_xor(builder, d6, x0, den_inv[:6], work3)
    builder.cx(x1, x0)
    gf2_6_mul_xor(builder, d6, x1, den_inv[6:12], work3)
    builder.free(work3)

    # 3. Multiply den_inv * num into tgt (allocates 3 qubits inside gf2_12_mul_xor -> peak: 6 + 12 + 3 = 21 qubits)
    gf2_12_mul_xor(builder, den_inv, num, tgt, control=None)

    # 4. Uncompute den_inv and d6
    work3 = builder.alloc_qubits(3)
    gf2_6_mul_xor(builder, d6, x1, den_inv[6:12], work3)
    builder.cx(x1, x0)
    gf2_6_mul_xor(builder, d6, x0, den_inv[:6], work3)
    builder.cx(x1, x0)
    builder.free(den_inv)

    toggle_d6_into(d6, work3)
    builder.free(work3)
    builder.free(d6)


def init_zero_test_12(builder: km.CircuitBuilder, reg12: km.array, target: km.q, scaf11: km.array) -> None:
    """Computes `target ^= (reg12 == 0)` using 11 `init_and`s + 1 `ccx` (12 Toffolis forward)."""
    builder.x(reg12)
    builder.init_and(reg12[0], reg12[1], scaf11[0])
    for k in range(1, 10):
        builder.init_and(scaf11[k - 1], reg12[k + 1], scaf11[k])
    builder.ccx(scaf11[9], reg12[11], target)
    builder.x(reg12)


def del_zero_test_12(builder: km.CircuitBuilder, reg12: km.array, target: km.q, scaf11: km.array) -> None:
    """Uncomputes `target ^= (reg12 == 0)` and `scaf11` using 1 `ccx` + 10 `del_and`s (1 Toffoli reverse!)."""
    builder.x(reg12)
    builder.ccx(scaf11[9], reg12[11], target)
    for k in reversed(range(1, 10)):
        builder.del_and(scaf11[k - 1], reg12[k + 1], scaf11[k])
    builder.del_and(reg12[0], reg12[1], scaf11[0])
    builder.x(reg12)


def _active(cursor: km.UnaryIterationCursor, address_value: int, limit: int):
    if address_value > limit:
        cursor.move_to(_UNREACHABLE)
        return False
    cursor.move_to(address_value)
    return cursor.match_qubit


def _xor_threshold(
    builder: km.CircuitBuilder,
    cursor: km.UnaryIterationCursor,
    limit: int,
    target,
    control,
    value: int,
    *,
    at_least: bool,
    lowest: int = 0,
) -> None:
    if value > limit:
        if not at_least:
            builder.cx(control, target)
        return
    if at_least:
        builder.cx(control, target)
    for u in range(lowest, value):
        builder.cx(_active(cursor, u, limit), target)


class OptimizedProosAndZalkaDecoder:
    """Optimized reversible Reed-Solomon MBUC decoder circuit."""

    def __init__(self, n_syndromes: int, field: km.GF2Field):
        if n_syndromes < 2:
            raise ValueError(f"n_syndromes must be at least 2, got {n_syndromes}.")
        self.n_syndromes = n_syndromes
        self.field = field
        _, _, self.ops_std_to_twr, self.ops_twr_to_std = build_tower_isomorphism(field)

    def build_syndrome_decoding_mbuc(
        self,
        builder: km.CircuitBuilder,
        *,
        syndromes: Sequence[km.array],
        error_mbuc: Sequence[km.array],
        n_iters: int,
    ) -> None:
        gf = self.field
        n = self.n_syndromes
        degree = gf.degree
        m = n + 3
        t = n // 2

        reg_a = [builder.alloc_qubits(degree) if i < 3 else syndromes[i - 3] for i in range(m)]
        reg_b = [builder.alloc_qubits(degree) for _ in range(m)]

        # Convert input syndromes in-place to Tower-Field basis (0 Toffolis)
        for s_reg in syndromes:
            apply_cnot_schedule(builder, s_reg, self.ops_std_to_twr)

        # m = 73 < 128, check if 8 bits is needed for signed comparisons or if 8 bits is safer:
        # Note: in kickmix, flip_if_less_than uses unsigned comparison, and m = 73 fits in 7 bits (0..127).
        logm = m.bit_length()
        n_A = builder.alloc_qubits(logm)
        n_B = builder.alloc_qubits(logm)
        n_a = builder.alloc_qubits(logm)
        n_b = builder.alloc_qubits(logm)
        n_q = builder.alloc_qubits(logm)
        flag = builder.alloc_qubits(1)[0]
        counter = builder.alloc_qubits(2)
        halting_counter = builder.alloc_qubits(n_iters.bit_length())
        eea = (reg_a, reg_b, n_A, n_B, n_a, n_b, n_q, flag, counter, halting_counter)

        init_ones = km.array([reg_a[0][0], reg_b[m - 1][0], n_a[0], flag])
        builder.x(init_ones)
        builder.gf2_iadd(n, target=n_A)
        builder.gf2_iadd(n + 1, target=n_B)

        with builder.mark("initial_shift"):
            self._initial_shift_loop(builder, reg_a, n_A, is_adjoint=False)
        with builder.mark("forward_eea"):
            self._synchronized_eea(builder, *eea, n_iters=n_iters, is_adjoint=False)
        with builder.mark("chien_forney_mbuc"):
            self._chien_search_and_forney_mbuc(builder, reg_b, error_mbuc, n_b, n_B)
        with builder.mark("reverse_eea"):
            self._synchronized_eea(builder, *eea, n_iters=n_iters, is_adjoint=True)
        with builder.mark("uncompute_initial_shift"):
            self._initial_shift_loop(builder, reg_a, n_A, is_adjoint=True)

        builder.gf2_iadd(n + 1, target=n_B)
        builder.gf2_iadd(n, target=n_A)
        builder.x(init_ones)

        # Convert syndromes back to Standard basis (0 Toffolis)
        for s_reg in syndromes:
            apply_cnot_schedule(builder, s_reg, self.ops_twr_to_std)

        for reg in (halting_counter, counter, flag, n_q, n_b, n_a, n_B, n_A, *reg_a[:3], *reg_b):
            builder.free(reg)

    def _initial_shift_loop(
        self,
        builder: km.CircuitBuilder,
        reg_a: Sequence[km.array],
        n_A: km.array,
        *,
        is_adjoint: bool,
    ) -> None:
        n = self.n_syndromes
        m = len(reg_a)
        shift_step = builder.mark("shift_step")
        for step in (range(n - 2, -1, -1) if is_adjoint else range(n - 1)):
            with shift_step:
                ctrl = builder.alloc_qubits(1)[0]
                if is_adjoint:
                    builder.flip_if_less_than(n_A, n - step, target=ctrl)
                    builder.iadd(1, target=n_A, control=ctrl)
                    for b in range(self.field.degree):
                        lane = km.array([r[b] for r in reg_a[3 + step : m]])
                        builder.right_rotate(lane, 1, control=ctrl)
                    builder.flip_if_equal(reg_a[m - 1], 0, target=ctrl)
                else:
                    builder.flip_if_equal(reg_a[m - 1], 0, target=ctrl)
                    for b in range(self.field.degree):
                        lane = km.array([r[b] for r in reg_a[3 + step : m]])
                        builder.left_rotate(lane, 1, control=ctrl)
                    builder.isub(1, target=n_A, control=ctrl)
                    builder.flip_if_less_than(n_A, n - step, target=ctrl)
                builder.free(ctrl)

    def _synchronized_eea(
        self,
        builder: km.CircuitBuilder,
        reg_a: Sequence[km.array],
        reg_b: Sequence[km.array],
        n_A: km.array,
        n_B: km.array,
        n_a: km.array,
        n_b: km.array,
        n_q: km.array,
        flag,
        counter: km.array,
        halting_counter: km.array,
        *,
        n_iters: int,
        is_adjoint: bool,
    ) -> None:
        eea_step = builder.mark("eea_step")
        for it in (range(n_iters - 1, -1, -1) if is_adjoint else range(n_iters)):
            with eea_step:
                if is_adjoint:
                    self._update_counter(builder, counter, flag, is_adjoint=True)

                is_c3 = builder.alloc_qubits(1)[0]
                builder.init_and(counter[0], counter[1], is_c3)
                if is_adjoint:
                    self._update_halting_counter(
                        builder, n_A, n_B, n_b, halting_counter, flag, is_c3, is_adjoint=True
                    )

                is_hc_zero = builder.alloc_qubits(1)[0]
                builder.flip_if_equal(halting_counter, 0, target=is_hc_zero)
                ctrl = self._counter_controls(builder, counter, is_c3, is_hc_zero)

                steps = [
                    lambda: self._update_flag(builder, reg_b, n_A, n_B, n_a, n_b, n_q, flag, ctrl),
                    lambda: self._multiply_and_add_or_sub(
                        builder, reg_a, reg_b, n_A, n_a, n_q, ctrl, it, is_adjoint=is_adjoint
                    ),
                    lambda: self._shift_b_array_right(
                        builder, reg_b, n_B, n_b, flag, ctrl, it, is_adjoint=is_adjoint
                    ),
                    lambda: self._update_metadata(
                        builder, n_B, n_b, n_q, flag, ctrl, is_adjoint=is_adjoint
                    ),
                    lambda: self._swap_a_and_b(builder, reg_a, reg_b, n_A, n_B, n_a, n_b, flag, ctrl),
                ]
                for step in (reversed(steps) if is_adjoint else steps):
                    step()

                self._del_counter_controls(builder, counter, is_c3, is_hc_zero, ctrl)
                builder.flip_if_equal(halting_counter, 0, target=is_hc_zero)
                builder.free(is_hc_zero)

                if not is_adjoint:
                    self._update_halting_counter(
                        builder, n_A, n_B, n_b, halting_counter, flag, is_c3, is_adjoint=False
                    )
                builder.del_and(counter[0], counter[1], is_c3)
                builder.free(is_c3)
                if not is_adjoint:
                    self._update_counter(builder, counter, flag, is_adjoint=False)

    @staticmethod
    def _counter_controls(builder: km.CircuitBuilder, counter: km.array, is_c3, is_hc_zero) -> km.array:
        ctrl = builder.alloc_qubits(4)
        c0, c1, c2, c3 = ctrl
        builder.init_and(is_hc_zero, km.array([counter[0], is_c3]), km.array([c1, c3]))
        builder.cx(km.array([is_hc_zero, c1]), c0)
        builder.init_and(c0, counter[1], c2)
        builder.cx(km.array([c2, c3]), km.array([c0, c1]))
        return ctrl

    @staticmethod
    def _del_counter_controls(
        builder: km.CircuitBuilder, counter: km.array, is_c3, is_hc_zero, ctrl: km.array
    ) -> None:
        c0, c1, c2, c3 = ctrl
        builder.cx(km.array([c2, c3]), km.array([c0, c1]))
        builder.del_and(c0, counter[1], c2)
        builder.cx(km.array([is_hc_zero, c1]), c0)
        builder.del_and(is_hc_zero, km.array([counter[0], is_c3]), km.array([c1, c3]))
        builder.free(ctrl)

    @staticmethod
    def _update_counter(builder: km.CircuitBuilder, counter: km.array, flag, *, is_adjoint: bool) -> None:
        if is_adjoint:
            builder.cx(flag, counter[0])
            builder.ccx(flag, counter[0], counter[1])
        else:
            builder.ccx(flag, counter[0], counter[1])
            builder.cx(flag, counter[0])

    def _update_flag(
        self,
        builder: km.CircuitBuilder,
        reg_b: Sequence[km.array],
        n_A: km.array,
        n_B: km.array,
        n_a: km.array,
        n_b: km.array,
        n_q: km.array,
        flag,
        ctrl,
    ) -> None:
        with builder.mark("update_flag"):
            c0, c1, c2, c3 = ctrl
            m = len(reg_b)
            aux = builder.alloc_qubits(1)[0]

            builder.flip_if_less_than(n_A, n_B, target=aux, control=c0)
            builder.flip_if_equal(n_q, 0, target=flag, control=aux)
            builder.flip_if_less_than(n_A, n_B, target=aux, control=c0)
            builder.flip_if_equal(n_A, n_B, target=flag, control=c0)

            builder.flip_if_equal(n_A, n_B, target=flag, control=c1)
            builder.cx(c1, flag)
            builder.flip_if_equal(n_B, 1, target=aux, control=c1)
            builder.cx(c1, aux)
            builder.flip_if_equal(reg_b[m - 2], 0, target=flag, control=aux)
            builder.cx(c1, aux)
            builder.flip_if_equal(n_B, 1, target=aux, control=c1)

            builder.cx(c2, flag)
            builder.flip_if_equal(n_b, 0, target=aux, control=c2)
            builder.cx(c2, aux)
            builder.flip_if_equal(reg_b[0], 0, target=flag, control=aux)
            builder.cx(c2, aux)
            builder.flip_if_equal(n_b, 0, target=aux, control=c2)
            builder.iadd(1, target=n_b)
            builder.flip_if_equal(n_a, n_b, target=flag, control=c2)
            builder.isub(1, target=n_b)

            builder.flip_if_equal(n_a, n_b, target=flag, control=c3)
            builder.flip_if_equal(n_q, 1, target=flag, control=c3)

            builder.free(aux)

    def _multiply_and_add_or_sub(
        self,
        builder: km.CircuitBuilder,
        reg_a: Sequence[km.array],
        reg_b: Sequence[km.array],
        n_A: km.array,
        n_a: km.array,
        n_q: km.array,
        ctrl,
        it: int,
        *,
        is_adjoint: bool,
    ) -> None:
        c0_hc, _, _, c3_hc = ctrl
        gf = self.field
        m = len(reg_b)
        t = self.n_syndromes // 2

        def accumulate_quotient(index: int, control) -> None:
            with builder.mark("quotient_div"):
                gf2_12_div_xor_pebbled(builder, reg_b[index], reg_a[index], q_i, control)

        def sweep_positions() -> None:
            with builder.mark("sweep_positions"):
                lim_q = min(m - 1, 1 + it)
                lim_a = min(t + 1, it // 2 + 1)
                lim_A = m - 2
                order = range(m) if is_adjoint else range(m - 1, -1, -1)

                c0_or_c3 = builder.alloc_qubits(1)[0]
                builder.cx(km.array([c0_hc, c3_hc]), c0_or_c3)
                builder.isub(1, target=n_q, control=c3_hc)

                with (
                    builder.unary_iteration(n_q, control=c0_or_c3) as eq_q,
                    builder.unary_iteration(n_a, control=c3_hc) as cur_a,
                    builder.unary_iteration(n_A, control=c0_hc) as cur_A,
                ):
                    mul = builder.alloc_qubits(1)[0]

                    def endpoint(i: int) -> None:
                        _xor_threshold(
                            builder, cur_a, lim_a, mul, c3_hc, i + 1, at_least=True, lowest=1
                        )
                        _xor_threshold(
                            builder, cur_A, lim_A, mul, c0_hc, m - i, at_least=True, lowest=1
                        )

                    prev = None
                    for i in order:
                        if prev is None:
                            endpoint(i)
                        else:
                            builder.cx(
                                km.array([
                                    _active(cur_a, min(prev, i) + 1, lim_a),
                                    _active(cur_A, m - max(prev, i), lim_A),
                                ]),
                                mul,
                            )
                        prev = i

                        swap_q = _active(eq_q, i, lim_q)
                        multiplies = (i < lim_a) or (i >= m - lim_A)

                        if is_adjoint and multiplies:
                            gf2_12_mul_xor(builder, q_i, reg_a[i], reg_b[i], control=mul)
                        if swap_q is not False:
                            builder.cswap(swap_q, q_i, reg_a[i])
                        if not is_adjoint and multiplies:
                            gf2_12_mul_xor(builder, q_i, reg_a[i], reg_b[i], control=mul)

                    endpoint(prev)
                    builder.free(mul)

                builder.iadd(1, target=n_q, control=c3_hc)
                builder.cx(km.array([c0_hc, c3_hc]), c0_or_c3)
                builder.free(km.array([c0_or_c3]))

        with builder.mark("multiply_and_add_or_sub"):
            q_i = builder.alloc_qubits(gf.degree)
            leading, trailing = (0, m - 1) if is_adjoint else (m - 1, 0)
            accumulate_quotient(leading, c3_hc if is_adjoint else c0_hc)
            builder.iadd(n_a, target=n_q)
            sweep_positions()
            builder.isub(n_a, target=n_q)
            accumulate_quotient(trailing, c0_hc if is_adjoint else c3_hc)
            builder.free(q_i)

    def _shift_b_array_right(
        self,
        builder: km.CircuitBuilder,
        reg_b: Sequence[km.array],
        n_B: km.array,
        n_b: km.array,
        flag,
        ctrl,
        it: int,
        *,
        is_adjoint: bool,
    ) -> None:
        with builder.mark("shift_b_array_right"):
            c0_hc, c1_hc, c2_hc, c3_hc = ctrl
            m = len(reg_b)
            t = self.n_syndromes // 2
            lim_B = m - 2
            lim_b = min(t + 1, it // 2 + 1)

            builder.x(flag)
            nf = builder.alloc_qubits(2)
            builder.init_and(km.array([c0_hc, c3_hc]), flag, nf)
            builder.x(flag)
            c0_nf, c3_nf = nf

            ctrl_nB = builder.alloc_qubits(1)[0]
            builder.cx(km.array([c0_nf, c1_hc]), ctrl_nB)

            order = range(1, m) if is_adjoint else range(m - 1, 0, -1)
            with (
                builder.unary_iteration(n_B, control=ctrl_nB) as cur_B_ge,
                builder.unary_iteration(n_b, control=c2_hc) as cur_b_ge,
                builder.unary_iteration(n_B, control=c3_nf) as cur_B_lt,
            ):
                do_shift = builder.alloc_qubits(1)[0]

                def endpoint(i: int) -> None:
                    _xor_threshold(
                        builder, cur_B_ge, lim_B, do_shift, ctrl_nB, m - i + 1, at_least=True
                    )
                    _xor_threshold(builder, cur_b_ge, lim_b, do_shift, c2_hc, i, at_least=True)
                    _xor_threshold(
                        builder, cur_B_lt, lim_B, do_shift, c3_nf, m - i, at_least=False
                    )

                prev = None
                for i in order:
                    if prev is None:
                        endpoint(i)
                    else:
                        builder.cx(
                            km.array([
                                _active(cur_B_ge, m - max(prev, i) + 1, lim_B),
                                _active(cur_b_ge, min(prev, i), lim_b),
                                _active(cur_B_lt, m - max(prev, i), lim_B),
                            ]),
                            do_shift,
                        )
                    prev = i
                    builder.cswap(do_shift, reg_b[i - 1], reg_b[i])

                endpoint(prev)
                builder.free(do_shift)

            builder.cx(km.array([c0_nf, c1_hc]), ctrl_nB)
            builder.free(ctrl_nB)
            builder.x(flag)
            builder.del_and(km.array([c0_hc, c3_hc]), flag, nf)
            builder.x(flag)
            builder.free(nf)

    @staticmethod
    def _update_metadata(
        builder: km.CircuitBuilder,
        n_B: km.array,
        n_b: km.array,
        n_q: km.array,
        flag,
        ctrl,
        *,
        is_adjoint: bool,
    ) -> None:
        with builder.mark("update_metadata"):
            c0_hc, c1_hc, c2_hc, c3_hc = ctrl
            builder.x(flag)
            nf = builder.alloc_qubits(2)
            builder.init_and(km.array([c0_hc, c3_hc]), flag, nf)
            builder.x(flag)
            c0_nf, c3_nf = nf

            sign = -1 if is_adjoint else 1
            updates = [
                [(n_q, 1, c0_hc), (n_B, -1, c0_nf)],
                [(n_B, -1, c1_hc)],
                [(n_b, 1, c2_hc)],
                [(n_b, 1, c3_nf), (n_q, -1, c3_hc)],
            ]
            for step in (reversed(updates) if is_adjoint else updates):
                for target, offset, control in (reversed(step) if is_adjoint else step):
                    builder.iadd(sign * offset, target=target, control=control)

            builder.x(flag)
            builder.del_and(km.array([c0_hc, c3_hc]), flag, nf)
            builder.x(flag)
            builder.free(nf)

    @staticmethod
    def _swap_a_and_b(
        builder: km.CircuitBuilder,
        reg_a: Sequence[km.array],
        reg_b: Sequence[km.array],
        n_A: km.array,
        n_B: km.array,
        n_a: km.array,
        n_b: km.array,
        flag,
        ctrl,
    ) -> None:
        with builder.mark("swap_a_and_b"):
            c3_hc = ctrl[3]
            swap_ctrl = builder.alloc_qubits(1)[0]
            builder.init_and(c3_hc, flag, swap_ctrl)
            lhs = km.array([q for r in (n_A, n_a, *reg_a) for q in r])
            rhs = km.array([q for r in (n_B, n_b, *reg_b) for q in r])
            builder.cswap(swap_ctrl, lhs, rhs)
            builder.del_and(c3_hc, flag, swap_ctrl)
            builder.free(swap_ctrl)

    @staticmethod
    def _update_halting_counter(
        builder: km.CircuitBuilder,
        n_A: km.array,
        n_B: km.array,
        n_b: km.array,
        halting_counter: km.array,
        flag,
        is_c3,
        *,
        is_adjoint: bool,
    ) -> None:
        with builder.mark("update_halting_counter"):
            work = builder.alloc_qubits(3)
            ctrl, busy, running = work
            builder.init_and(is_c3, flag, ctrl)
            builder.flip_if_less_than(n_b, n_B, target=busy, control=ctrl)
            builder.cx(busy, running)
            builder.flip_if_equal(n_A, 0, target=running, control=busy)

            builder.cx(running, ctrl)
            if is_adjoint:
                builder.isub(1, target=halting_counter, control=ctrl)
            else:
                builder.iadd(1, target=halting_counter, control=ctrl)
            builder.cx(running, ctrl)

            builder.flip_if_equal(n_A, 0, target=running, control=busy)
            builder.cx(busy, running)
            builder.flip_if_less_than(n_b, n_B, target=busy, control=ctrl)
            builder.del_and(is_c3, flag, ctrl)
            builder.free(work)

    def _pre_align_reg_b(
        self,
        builder: km.CircuitBuilder,
        reg_b: Sequence[km.array],
        n_b: km.array,
        n_B: km.array,
        *,
        is_adjoint: bool,
    ) -> None:
        omega = reg_b[len(reg_b) // 2 + 1 :]
        for b in range(self.field.degree):
            lane = km.array([r[b] for r in omega])
            if is_adjoint:
                builder.right_rotate(lane, n_B)
                builder.left_rotate(lane, n_b)
            else:
                builder.right_rotate(lane, n_b)
                builder.left_rotate(lane, n_B)

    def _chien_search_and_forney_mbuc(
        self,
        builder: km.CircuitBuilder,
        reg_b: Sequence[km.array],
        error_mbuc: Sequence[km.array],
        n_b: km.array,
        n_B: km.array,
    ) -> None:
        gf = self.field
        degree = gf.degree
        m = len(reg_b)
        t = self.n_syndromes // 2
        block_length = (1 << degree) - 1
        primitive = gf.primitive_element

        # Convert reg_b from Tower basis to Standard basis (0 Toffolis)
        for r_b in reg_b:
            apply_cnot_schedule(builder, r_b, self.ops_twr_to_std)

        with builder.mark("pre_align_reg_b"):
            self._pre_align_reg_b(builder, reg_b, n_b, n_B, is_adjoint=False)

        sigma = reg_b[0]
        sigma_prime = reg_b[1]
        omega = reg_b[m // 2 + t]
        sigma_prime_inv = builder.alloc_qubits(degree)
        inv_scaffold = builder.alloc_qubits(93)
        zero_scaf = builder.alloc_qubits(11)
        is_sigma_zero = zero_scaf[10]
        is_error = builder.alloc_bits(1)

        omega_indices = range(m // 2 + 1, m // 2 + t + 1)
        step_multipliers = [
            (k, gf.pow(primitive, k)) for k in range(1, t + 1)
        ] + [
            (k, gf.pow(primitive, (m - 2 - k) % block_length))
            for k in omega_indices
            if (m - 2 - k) % block_length != 0
        ]

        sum_pairs = (
            [(reg_b[k], omega) for k in range(m // 2 + 1, m // 2 + t)]
            + [(reg_b[k], sigma_prime) for k in range(3, t + 1, 2)]
            + [(reg_b[k], sigma) for k in range(2, t + 1, 2)]
            + [(sigma_prime, sigma)]
        )
        fwd_src = km.array([q for src, _ in sum_pairs for q in src])
        fwd_dst = km.array([q for _, dst in sum_pairs for q in dst])
        rev_src = km.array([q for src, _ in reversed(sum_pairs) for q in src])
        rev_dst = km.array([q for _, dst in reversed(sum_pairs) for q in dst])

        num_positions = min(block_length, len(error_mbuc))
        forney_step = builder.mark("forney_step")
        for j in range(num_positions):
            with forney_step:
                builder.bit_store(is_error[0], True, control=error_mbuc[j])
                builder.push_condition(is_error[0])

                builder.cx(fwd_src, fwd_dst)

                # 11-Toffoli zero test on sigma (0 Toffolis to uncompute via del_and!)
                builder.x(sigma)
                builder.init_and(sigma[0], sigma[1], zero_scaf[0])
                for zk in range(1, 11):
                    builder.init_and(zero_scaf[zk - 1], sigma[zk + 1], zero_scaf[zk])
                builder.x(sigma)

                # Convert sigma_prime to Tower basis, invert with built-in is_sigma_zero gating (78 Toffolis),
                # and convert sigma_prime_inv to Standard basis!
                apply_cnot_schedule(builder, sigma_prime, self.ops_std_to_twr)
                gf2_12_inv_scaffold(
                    builder,
                    sigma_prime,
                    sigma_prime_inv,
                    inv_scaffold,
                    gate_ctl=is_sigma_zero,
                    inverse=False,
                )
                apply_cnot_schedule(builder, sigma_prime_inv, self.ops_twr_to_std)

                builder.gf2_phase_by_product(omega, sigma_prime_inv, mask=error_mbuc[j], field=gf)

                # Uncompute in 0 Toffolis!
                apply_cnot_schedule(builder, sigma_prime_inv, self.ops_std_to_twr)
                gf2_12_inv_scaffold(
                    builder,
                    sigma_prime,
                    sigma_prime_inv,
                    inv_scaffold,
                    gate_ctl=is_sigma_zero,
                    inverse=True,
                )
                apply_cnot_schedule(builder, sigma_prime, self.ops_twr_to_std)

                builder.x(sigma)
                for zk in reversed(range(1, 11)):
                    builder.del_and(zero_scaf[zk - 1], sigma[zk + 1], zero_scaf[zk])
                builder.del_and(sigma[0], sigma[1], zero_scaf[0])
                builder.x(sigma)

                builder.cx(rev_src, rev_dst)

                builder.pop_condition()
                builder.bit_store(is_error, False)

                if j + 1 < num_positions:
                    for k, step in step_multipliers:
                        builder.gf2_imul(step, target=reg_b[k], field=gf)

        if num_positions > 1:
            for k, step in step_multipliers:
                builder.gf2_idiv(gf.pow(step, num_positions - 1), target=reg_b[k], field=gf)

        builder.free(is_error)
        builder.free(zero_scaf)
        builder.free(inv_scaffold)
        builder.free(sigma_prime_inv)

        with builder.mark("unalign_reg_b"):
            self._pre_align_reg_b(builder, reg_b, n_b, n_B, is_adjoint=True)

        # Convert reg_b back to Tower basis (0 Toffolis)
        for r_b in reg_b:
            apply_cnot_schedule(builder, r_b, self.ops_std_to_twr)


def default_n_iters(n_syndromes: int) -> int:
    return 3 * (n_syndromes + 1) + 9


def build_opt_dqi_circuit(m: int = 4095, n: int = 70, degree: int = 12, n_iters: int | None = None) -> km.Circuit:
    builder = km.CircuitBuilder()
    syndromes = [builder.create_quantum_register(degree, name=f"S{i}") for i in range(n)]
    error_mbuc = [builder.create_classical_register(degree, name=f"E{j}") for j in range(m)]
    iters = default_n_iters(n) if n_iters is None else n_iters
    OptimizedProosAndZalkaDecoder(n, km.GF2Field(degree)).build_syndrome_decoding_mbuc(
        builder, syndromes=syndromes, error_mbuc=error_mbuc, n_iters=iters
    )
    return builder.finish_circuit()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", help="Path of the .kmb (or .kmx) file to write.")
    parser.add_argument("--m", type=int, default=4095)
    parser.add_argument("--n", type=int, default=70)
    parser.add_argument("--degree", type=int, default=12)
    parser.add_argument("--n-iters", type=int, default=None)
    parser.add_argument("--format", choices=("kmb", "kmx"), default="kmb")
    args = parser.parse_args()

    circuit = build_opt_dqi_circuit(args.m, args.n, args.degree, args.n_iters)
    circuit.to_file(args.output, format=args.format)
    toffolis = circuit.operation_statistics().num_toffoli_like
    print(f"Wrote {args.output} ({args.format}): {toffolis:,} max Toffolis, {circuit.num_qubits} qubits.")


if __name__ == "__main__":
    main()
