"""Construct and verify the GF(2^12) = GF(((2^3)^2)^2) Tower-Field Multiplier and Inverter in kickmix."""

from __future__ import annotations

import kickmix as km
import numpy as np


def build_tower_isomorphism(gf: km.GF2Field) -> tuple[np.ndarray, np.ndarray, list[int], list[int]]:
    """Builds the 12x12 GF(2)-linear basis change matrices T_std_to_twr and T_twr_to_std."""
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
    for col in range(12):
        piv = col
        while piv < 12 and aug[piv, col] == 0:
            piv += 1
        assert piv < 12
        if piv != col:
            aug[[col, piv]] = aug[[piv, col]]
        for row in range(12):
            if row != col and aug[row, col] == 1:
                aug[row] ^= aug[col]
    T_std_to_twr = aug[:, 12:]
    return T_std_to_twr, T_twr_to_std, exp_tab, log_tab


def gf2_3_mul_xor(builder: km.CircuitBuilder, u: km.array, v: km.array, tgt: km.array, work6: km.array) -> None:
    """Computes `tgt ^= u * v` in GF(2)[omega]/(omega^3 + omega + 1) using 6 Toffolis and `work6` (6 clean qubits)."""
    builder.init_and(u, v, work6[:3])  # m0, m1, m2
    builder.cx(km.array([u[0], v[0]]), km.array([u[1], v[1]]))
    builder.init_and(u[1], v[1], work6[3])  # m01
    builder.cx(km.array([u[0], v[0], u[1], v[1]]), km.array([u[1], v[1], u[2], v[2]]))
    builder.init_and(u[2], v[2], work6[5])  # m12
    builder.cx(km.array([u[1], v[1], u[0], v[0]]), km.array([u[2], v[2], u[2], v[2]]))
    builder.init_and(u[2], v[2], work6[4])  # m02
    builder.cx(km.array([u[0], v[0]]), km.array([u[2], v[2]]))

    m0, m1, m2, m01, m02, m12 = work6
    builder.cx(km.array([m0, m1, m2, m12]), tgt[0])
    builder.cx(km.array([m0, m01, m12]), tgt[1])
    builder.cx(km.array([m0, m1, m02]), tgt[2])

    builder.cx(km.array([u[0], v[0]]), km.array([u[2], v[2]]))
    builder.del_and(u[2], v[2], work6[4])
    builder.cx(km.array([u[0], v[0], u[1], v[1]]), km.array([u[2], v[2], u[2], v[2]]))
    builder.del_and(u[2], v[2], work6[5])
    builder.cx(km.array([u[1], v[1], u[0], v[0]]), km.array([u[2], v[2], u[1], v[1]]))
    builder.del_and(u[1], v[1], work6[3])
    builder.cx(km.array([u[0], v[0]]), km.array([u[1], v[1]]))
    builder.del_and(u, v, work6[:3])


def gf2_3_mul_keep(builder: km.CircuitBuilder, u: km.array, v: km.array, tgt: km.array, m_out: km.array, inverse: bool = False) -> None:
    """Forward (`inverse=False`): computes `m_out` (6 qubits) via `init_and` and XORs `u * v` into `tgt` (6 Toffolis).
    Reverse (`inverse=True`): uncomputes `tgt ^= u * v` and `m_out` via `del_and` (0 Toffolis!).
    """
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


def gf2_6_mul_xor(builder: km.CircuitBuilder, a: km.array, b: km.array, tgt: km.array, work6: km.array) -> None:
    """Computes `tgt ^= a * b` in GF(2^6) using 18 Toffolis."""
    a0, a1 = a[:3], a[3:6]
    b0, b1 = b[:3], b[3:6]
    t0, t1 = tgt[:3], tgt[3:6]

    builder.cx(t0, t1)
    gf2_3_mul_xor(builder, a0, b0, t0, work6)
    builder.cx(t0, t1)

    gf2_3_mul_xor(builder, a1, b1, t0, work6)

    builder.cx(a0, a1)
    builder.cx(b0, b1)
    gf2_3_mul_xor(builder, a1, b1, t1, work6)
    builder.cx(b0, b1)
    builder.cx(a0, a1)


def gf2_6_mul_keep(builder: km.CircuitBuilder, a: km.array, b: km.array, tgt: km.array, m18: km.array, inverse: bool = False) -> None:
    """Forward (`inverse=False`): computes `tgt ^= a * b` in GF(2^6) keeping 18 AND-ancillas `m18` (18 Toffolis).
    Reverse (`inverse=True`): uncomputes `tgt ^= a * b` and `m18` via `del_and` (0 Toffolis!).
    """
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


def gf2_3_sq_xor(builder: km.CircuitBuilder, u: km.array, tgt: km.array) -> None:
    """Computes `tgt ^= u^2` in GF(2)[omega]/(omega^3 + omega + 1) (0 Toffolis): (u0, u2, u1 ^ u2)."""
    builder.cx(u[0], tgt[0])
    builder.cx(u[2], tgt[1])
    builder.cx(km.array([u[1], u[2]]), tgt[2])


def gf2_3_inv_keep(
    builder: km.CircuitBuilder,
    u: km.array,
    inv_u: km.array,
    p3: km.array,
    inverse: bool = False,
) -> None:
    """Computes `inv_u = u^{-1}` in GF(2)[omega]/(omega^3 + omega + 1) using 3 `init_and`s into `p3` (3 Toffolis).
    Uncomputes in 0 Toffolis when `inverse=True`.
    Exact ANF of u^{-1} = u^6 in GF(2)[omega]/(omega^3 + omega + 1):
      Let p01 = u0 & u1, p02 = u0 & u2, p12 = u1 & u2.
      Truth table of u^{-1}:
        0=(0,0,0) -> (0,0,0)
        1=(1,0,0) -> (1,0,0)
        2=(0,1,0)=w -> w^6 = 1+w^2 = (1,0,1)
        3=(1,1,0)=1+w=w^3 -> w^4 = w+w^2 = (0,1,1)
        4=(0,0,1)=w^2 -> w^5 = 1+w+w^2 = (1,1,1)
        5=(1,0,1)=1+w^2=w^6 -> w = (0,1,0)
        6=(0,1,1)=w+w^2=w^4 -> w^3 = 1+w = (1,1,0)
        7=(1,1,1)=1+w+w^2=w^5 -> w^2 = (0,0,1)
      Therefore:
        inv0 = u0 ^ u1 ^ u2 ^ p02
        inv1 = u2 ^ p01 ^ p12
        inv2 = u1 ^ u2 ^ p01 ^ p02
    """
    p01, p02, p12 = p3
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

    `scaffold` has length 96 (or 99 when `gate_ctl` is not None):
      - N6 (6), n3 (3), inv_n3 (3), d6 (6) = 18 intermediate value qubits
      - m_N6 (18), m_n3 (6), p_n3 (3), m_d6 (12), m_tgt (36) = 75 AND-ancilla qubits
      - if `gate_ctl` is not None: `inv_n3_gated` (3)
    Forward (`inverse=False`): 75 Toffolis (78 if `gate_ctl` is not None).
    Reverse (`inverse=True`): 0 Toffolis!
    """
    N6 = scaffold[0:6]
    n3 = scaffold[6:9]
    inv_n3 = scaffold[9:12]
    d6 = scaffold[12:18]
    m_N6 = scaffold[18:36]
    m_n3 = scaffold[36:42]
    p_n3 = scaffold[42:45]
    m_d6 = scaffold[45:57]
    m_tgt = scaffold[57:93]
    inv_n3_use = scaffold[93:96] if gate_ctl is not None else inv_n3

    x0, x1 = x[:6], x[6:12]
    N6_0, N6_1 = N6[:3], N6[3:6]
    d6_0, d6_1 = d6[:3], d6[3:6]
    tgt0, tgt1 = target[:6], target[6:12]

    if not inverse:
        # 1. N6 = x0 * (x0 ^ x1) ^ lambda * x1^2  (18 Toffolis)
        gf2_3_sq_xor(builder, x1[3:6], N6_0)
        gf2_3_sq_xor(builder, x1[0:3], N6_1)
        builder.cx(x1, x0)  # x0 now holds x0 ^ x1
        builder.cx(x0, x1)  # x1 now holds orig x0
        gf2_6_mul_keep(builder, x1, x0, N6, m_N6, inverse=False)
        # Note: leave x0 = orig_x0 ^ orig_x1 and x1 = orig_x0? Wait: restore x0, x1 or keep x0 ^ x1!
        # Let's restore x0, x1 so x is untouched outside:
        builder.cx(x0, x1)
        builder.cx(x1, x0)

        # 2. n3 = N6_0 * (N6_0 ^ N6_1) ^ N6_1^2  (6 Toffolis)
        gf2_3_sq_xor(builder, N6_1, n3)
        builder.cx(N6_0, N6_1)  # N6_1 now holds N6_0 ^ N6_1
        gf2_3_mul_keep(builder, N6_0, N6_1, n3, m_n3, inverse=False)
        builder.cx(N6_0, N6_1)

        # 3. inv_n3 = n3^{-1}  (3 Toffolis, +3 if gated by gate_ctl)
        gf2_3_inv_keep(builder, n3, inv_n3, p_n3, inverse=False)
        if gate_ctl is not None:
            builder.init_and(gate_ctl, inv_n3, inv_n3_use)

        # 4. d6_0 = inv_n3_use * (N6_0 ^ N6_1), d6_1 = inv_n3_use * N6_1  (12 Toffolis)
        builder.cx(N6_0, N6_1)  # N6_1 holds N6_0 ^ N6_1
        gf2_3_mul_keep(builder, inv_n3_use, N6_1, d6_0, m_d6[:6], inverse=False)
        builder.cx(N6_0, N6_1)
        gf2_3_mul_keep(builder, inv_n3_use, N6_1, d6_1, m_d6[6:12], inverse=False)

        # 5. tgt0 ^= d6 * (x0 ^ x1), tgt1 ^= d6 * x1  (36 Toffolis)
        builder.cx(x1, x0)  # x0 holds x0 ^ x1
        gf2_6_mul_keep(builder, d6, x0, tgt0, m_tgt[:18], inverse=False)
        builder.cx(x1, x0)
        gf2_6_mul_keep(builder, d6, x1, tgt1, m_tgt[18:36], inverse=False)
    else:
        # Uncompute in exact reverse order using 0 Toffolis!
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
        gf2_3_inv_keep(builder, n3, inv_n3, p_n3, inverse=True)

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


def gf2_12_mul_xor(
    builder: km.CircuitBuilder,
    x: km.array,
    y: km.array,
    tgt: km.array,
    control: km.q | None = None,
) -> None:
    """Computes `tgt ^= x * y` (or `tgt ^= control * (x * y)`) in GF(2^12) Tower Basis."""
    work6 = builder.alloc_qubits(6)
    if control is not None:
        x_use = builder.alloc_qubits(12)
        builder.init_and(control, x, x_use)
    else:
        x_use = x

    x0, x1 = x_use[:6], x_use[6:12]
    y0, y1 = y[:6], y[6:12]
    out0, out1 = tgt[:6], tgt[6:12]

    builder.cx(out0, out1)
    gf2_6_mul_xor(builder, x0, y0, out0, work6)
    builder.cx(out0, out1)

    builder.cx(x0, x1)
    builder.cx(y0, y1)
    gf2_6_mul_xor(builder, x1, y1, out1, work6)
    builder.cx(y0, y1)
    builder.cx(x0, x1)

    builder.cx(out0[:3], out0[3:6])
    out0_swapped = km.array([*out0[3:6], *out0[:3]])
    gf2_6_mul_xor(builder, x1, y1, out0_swapped, work6)
    builder.cx(out0[:3], out0[3:6])

    if control is not None:
        builder.del_and(control, x, x_use)
        builder.free(x_use)
    builder.free(work6)


def main() -> None:
    gf = km.GF2Field(12)
    T_std_to_twr, T_twr_to_std, exp_tab, log_tab = build_tower_isomorphism(gf)

    # Test gf2_12_inv_scaffold forward and reverse on all 4096 elements of GF(2^12)!
    b = km.CircuitBuilder()
    reg_x = b.create_quantum_register(12, name="x")
    reg_inv = b.create_quantum_register(12, name="inv_x")
    scaf = b.alloc_qubits(93)
    gf2_12_inv_scaffold(b, reg_x, reg_inv, scaf, inverse=False)
    b.free(scaf)
    # Wait: test forward value first, then test forward+reverse uncomputation!
    b2 = km.CircuitBuilder()
    rx2 = b2.create_quantum_register(12, name="x")
    rout2 = b2.create_quantum_register(12, name="out")
    inv_tmp = b2.alloc_qubits(12)
    scaf2 = b2.alloc_qubits(93)
    gf2_12_inv_scaffold(b2, rx2, inv_tmp, scaf2, inverse=False)
    b2.cx(inv_tmp, rout2)
    gf2_12_inv_scaffold(b2, rx2, inv_tmp, scaf2, inverse=True)
    b2.free(scaf2)
    b2.free(inv_tmp)
    circ_inv = b2.finish_circuit()
    print(f"gf2_12_inv_scaffold (fwd+rev): {circ_inv.operation_statistics().num_toffoli_like} Toffolis, {circ_inv.num_qubits} qubits.")

    xs_std = np.arange(4096, dtype=np.int64)

    def int_to_bits(arr: np.ndarray) -> np.ndarray:
        return ((arr[None, :] >> np.arange(12)[:, None]) & 1).astype(bool)

    def bits_to_int(bits: np.ndarray) -> np.ndarray:
        return (1 << np.arange(12)) @ bits.astype(np.int64)

    x_bits_std = int_to_bits(xs_std)
    x_bits_twr = ((T_std_to_twr @ x_bits_std.astype(np.int64)) % 2).astype(bool)

    sim = km.Simulator(batch_size=4096, seed=123)
    sim.use_same_registers_as(circ_inv)
    sim.write_across_shots(rx2, x_bits_twr)
    sim.do(circ_inv)
    assert not np.any(sim.read_phase_flipped_across_shots()), "Phase flipped during del_and!"
    out_bits_twr = sim.read_across_shots(rout2)
    out_bits_std = ((T_twr_to_std @ out_bits_twr.astype(np.int64)) % 2).astype(bool)
    inv_got = bits_to_int(out_bits_std)
    inv_exp = np.array([0 if x == 0 else exp_tab[(4095 - log_tab[int(x)]) % 4095] for x in xs_std])
    mismatches = int(np.sum(inv_got != inv_exp))
    print(f"Inverter verification on all 4096 elements: mismatches = {mismatches}/4096")
    assert mismatches == 0


if __name__ == "__main__":
    main()
