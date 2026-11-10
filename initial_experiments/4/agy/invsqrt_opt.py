"""Optimized Piecewise-Quadratic Circuit for Benchmark I (`invsqrt`).

Key optimizations:
1. In-place thermometer-to-binary compression of `flags` (10 -> 4 qubits `s_bin`) via 14 CNOTs
   and a 4-bit `del_lookup` (3 Toffolis), freeing 6 qubits before polynomial evaluation.
2. Top-bit alignment (`b_even = ~x[21]`, controlled 1-bit left-rotate of `x`) + bitwise inversion
   `d_bar = ~x[0:14]` so all Horner steps are pure additions (`+`) with non-negative coefficients.
3. Group-shared upper `c0` (`v6_table` on 6-bit address `addr[2:8]`) so `y0_top6` (6 qubits) is
   only allocated after Step B, keeping `y0_rem` at 23 qubits during `init_lookup`, Step A, and Step B.
4. Early `del_lookup` of `c2` (10 qubits) right after Step A so `c2` is unallocated during Step B,
   Denormalize, Phase, and Uncompute Step B, then recomputed via 2-block Select-Swap QROAM inside
   the freed space of `y0_rem` and uncomputed jointly with `r1` in a single final `del_lookup`.
5. Unified recursive pebbled Gidney controlled shifted adder (`add_shifted_pebbled`) + hierarchical
   carry splits (`carry_a10`, `carry_b14`, `carry_b10`), achieving 90 qubits and ~2,369 Toffolis.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import math
import kickmix as km
import numpy as np


def compute_piecewise_table(
    wp: int = 28,
    d1_bits: int = 14,
    d2_bits: int = 9,
    c2_drop_lsbs: int = 3,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Computes minimax integer coefficients (c0, c1, c2) for all 256 intervals."""
    nodes = 0.5 * (1.0 - np.cos(np.pi * (2 * np.arange(3) + 1) / 6.0))

    c0_arr = np.zeros(256, dtype=np.int64)
    c1_arr = np.zeros(256, dtype=np.int64)
    c2_arr = np.zeros(256, dtype=np.int64)

    mask_c2 = ~((1 << c2_drop_lsbs) - 1)

    for idx in range(256):
        be = idx >> 7
        sub = idx & 0x7F
        scale = 1.0 if be == 1 else 2.0
        dg = np.arange(0, 1 << 14, 6 if be == 1 else 3, dtype=np.int64)
        if be == 1 and dg[-1] != 0x3FFE:
            dg = np.append(dg, 0x3FFE)
        elif be == 0 and dg[-1] != 0x3FFF:
            dg = np.append(dg, 0x3FFF)
        db = 0x3FFF - dg
        bg = [(db >> (14 - m)) & 1 for m in range(1, max(d1_bits, d2_bits) + 1)]

        u_nodes = 1.0 + (sub + 1.0 - (1.0 / 16384.0) - nodes) / 128.0
        y_nodes = (1 << wp) / np.sqrt(scale * u_nodes)
        V = np.vstack([np.ones(3), nodes, nodes**2]).T
        coefs = np.linalg.solve(V, y_nodes)
        c2_init = max(0, int(round(coefs[2]))) & mask_c2
        c1_init = max(0, int(round(coefs[1])))

        u_grid = 1.0 + (sub + dg / 16384.0) / 128.0
        exact_y = (1 << wp) / np.sqrt(scale * u_grid)

        best_err = 1e18
        best_tuple = (0, 0, 0)
        step_c2 = max(1, 1 << c2_drop_lsbs)
        for dc2 in range(-3 * step_c2, 6 * step_c2, step_c2):
            c2 = max(0, c2_init + dc2) & mask_c2
            t2_grid = sum(bg[m - 1] * (c2 >> m) for m in range(1, d2_bits + 1))
            for dc1 in range(-12, 8):
                c1 = max(0, c1_init + dc1)
                r1_grid = c1 + t2_grid
                t1_grid = sum(bg[m - 1] * (r1_grid >> m) for m in range(1, d1_bits + 1))
                diff = exact_y - t1_grid
                hi, lo = float(diff.max()), float(diff.min())
                for c0 in (int(math.floor(0.5 * (hi + lo))), int(math.ceil(0.5 * (hi + lo)))):
                    c0_clamped = min((1 << wp) - 1, max(0, c0))
                    err = np.max(np.abs(c0_clamped - diff) / exact_y)
                    if err < best_err:
                        best_err = err
                        best_tuple = (c0_clamped, c1, c2)
        c0_arr[idx], c1_arr[idx], c2_arr[idx] = best_tuple

    return c0_arr, c1_arr, c2_arr


def add_single_bit_into(builder: km.CircuitBuilder, cin: km.q, tgt: km.array, subtract: bool = False) -> None:
    """Computes `tgt += cin` (or `tgt -= cin` if `subtract=True`) using W - 1 Toffolis and max(0, W - 2) ancillas."""
    W = len(tgt)
    if subtract:
        builder.x(tgt)
    if W == 1:
        builder.cx(cin, tgt[0])
    elif W == 2:
        builder.ccx(cin, tgt[0], tgt[1])
        builder.cx(cin, tgt[0])
    else:
        carries = builder.alloc_qubits(W - 2)
        cur = cin
        for i in range(W - 2):
            builder.init_and(cur, tgt[i], carries[i])
            cur = carries[i]
        builder.ccx(cur, tgt[W - 2], tgt[W - 1])
        builder.cx(cur, tgt[W - 2])
        for i in reversed(range(W - 2)):
            prev = cin if i == 0 else carries[i - 1]
            builder.del_and(prev, tgt[i], carries[i])
            builder.cx(prev, tgt[i])
        builder.free(carries)
    if subtract:
        builder.x(tgt)


def add_shifted_pebbled(
    builder: km.CircuitBuilder,
    ctl: km.q,
    lhs: km.array,
    tgt: km.array,
    subtract: bool = False,
    max_carries: int = 20,
) -> None:
    """Unified pebbled Gidney adder `tgt += ctl * lhs` using at most `max_carries` ancilla qubits."""
    k = len(lhs)
    W = len(tgt)
    assert 1 <= k <= W
    if subtract:
        builder.x(tgt)
    if W == 1:
        builder.ccx(ctl, lhs[0], tgt[0])
        if subtract:
            builder.x(tgt)
        return

    steps = [("gidney", i) for i in range(min(k, W - 1))]
    if k < W - 1:
        steps.append(("gate_ctl", k))
        for j in range(k, W - 2):
            steps.append(("ripple", j))
    N = len(steps)

    def step_up(idx: int, prev_c: km.q | None, out_c: km.q) -> None:
        kind, pos = steps[idx]
        if kind == "gidney":
            if pos == 0:
                builder.init_and(lhs[0], tgt[0], out_c)
            else:
                builder.cx(prev_c, lhs[pos])
                builder.cx(prev_c, tgt[pos])
                builder.init_and(lhs[pos], tgt[pos], out_c)
                builder.cx(prev_c, out_c)
        elif kind == "gate_ctl":
            builder.init_and(ctl, prev_c, out_c)
        elif kind == "ripple":
            builder.init_and(prev_c, tgt[pos], out_c)

    def step_down(idx: int, prev_c: km.q | None, out_c: km.q, with_sum: bool) -> None:
        kind, pos = steps[idx]
        if kind == "gidney":
            if pos == 0:
                builder.del_and(lhs[0], tgt[0], out_c)
                if with_sum:
                    builder.ccx(ctl, lhs[0], tgt[0])
            else:
                builder.cx(prev_c, out_c)
                builder.del_and(lhs[pos], tgt[pos], out_c)
                if with_sum:
                    builder.ccx(ctl, lhs[pos], tgt[pos])
                builder.cx(prev_c, tgt[pos])
                builder.cx(prev_c, lhs[pos])
        elif kind == "gate_ctl":
            if with_sum:
                builder.cx(out_c, tgt[pos])
            builder.del_and(ctl, prev_c, out_c)
        elif kind == "ripple":
            if with_sum:
                builder.cx(out_c, tgt[pos + 1])
            builder.del_and(prev_c, tgt[pos], out_c)

    def top_action(last_c: km.q) -> None:
        if k == W:
            builder.cx(last_c, lhs[W - 1])
            builder.ccx(ctl, lhs[W - 1], tgt[W - 1])
            builder.cx(last_c, lhs[W - 1])
        elif k == W - 1:
            builder.ccx(ctl, last_c, tgt[W - 1])
        else:
            builder.ccx(last_c, tgt[W - 2], tgt[W - 1])

    def run_segment(start: int, end: int, prev_c: km.q | None, budget: int) -> None:
        L = end - start
        if L <= budget:
            c = builder.alloc_qubits(L)
            for j in range(L):
                step_up(start + j, prev_c if j == 0 else c[j - 1], c[j])
            if end == N:
                top_action(c[L - 1])
            for j in reversed(range(L)):
                step_down(start + j, prev_c if j == 0 else c[j - 1], c[j], with_sum=True)
            builder.free(c)
        else:
            s = min(budget - 1, L - (budget - 1))
            if s < 1:
                s = L // 2
            mid = start + s
            c_mid = builder.alloc_qubits(1)[0]
            c_low = builder.alloc_qubits(s)
            for j in range(s):
                step_up(start + j, prev_c if j == 0 else c_low[j - 1], c_low[j])
            builder.cx(c_low[s - 1], c_mid)
            for j in reversed(range(s)):
                step_down(start + j, prev_c if j == 0 else c_low[j - 1], c_low[j], with_sum=False)
            builder.free(c_low)

            run_segment(mid, end, c_mid, budget - 1)

            c_low2 = builder.alloc_qubits(s)
            for j in range(s):
                step_up(start + j, prev_c if j == 0 else c_low2[j - 1], c_low2[j])
            builder.cx(c_low2[s - 1], c_mid)
            builder.free(km.array([c_mid]))
            for j in reversed(range(s)):
                step_down(start + j, prev_c if j == 0 else c_low2[j - 1], c_low2[j], with_sum=True)
            builder.free(c_low2)

    run_segment(0, N, None, max_carries)
    if subtract:
        builder.x(tgt)


def _flag_zero_pair(builder: km.CircuitBuilder, inverse: bool, x: km.array, flag: km.q) -> None:
    top_pair = x[-2:]
    builder.x(top_pair)
    if inverse:
        builder.del_and(top_pair[0], top_pair[1], flag)
    else:
        builder.init_and(top_pair[0], top_pair[1], flag)
    builder.x(top_pair)


def _normalize(builder: km.CircuitBuilder, inverse: bool, x: km.array, flags: km.array) -> None:
    for k in (reversed(range(len(flags))) if inverse else range(len(flags))):
        lanes = (x[2 * k::2], x[2 * k + 1::2])
        if inverse:
            for lane in lanes:
                builder.right_rotate(lane, 1, control=flags[k])
            _flag_zero_pair(builder, inverse, x, flags[k])
        else:
            _flag_zero_pair(builder, inverse, x, flags[k])
            for lane in lanes:
                builder.left_rotate(lane, 1, control=flags[k])


def _thermometer_cnots(builder: km.CircuitBuilder, flags: km.array, inverse: bool) -> None:
    ops = [
        (flags[1], flags[0]), (flags[2], flags[0]), (flags[3], flags[0]),
        (flags[4], flags[0]), (flags[5], flags[0]), (flags[6], flags[0]),
        (flags[7], flags[0]), (flags[8], flags[0]), (flags[9], flags[0]),
        (flags[3], flags[1]), (flags[5], flags[1]), (flags[7], flags[1]), (flags[9], flags[1]),
        (flags[7], flags[3]),
    ]
    for src, dst in (reversed(ops) if inverse else ops):
        builder.cx(src, dst)


def _rem_flags_table() -> list[int]:
    rem_indices = [2, 4, 5, 6, 8, 9]
    table = []
    for s in range(16):
        val = 0
        for bit_pos, k in enumerate(rem_indices):
            if s > k:
                val |= (1 << bit_pos)
        table.append(val)
    return table


def build_opt_invsqrt_circuit(
    n_p: int = 10,
    wp: int = 28,
    d1_bits: int = 14,
    d2_bits: int = 9,
    c2_drop_lsbs: int = 3,
    target_qubits: int = 90,
) -> km.Circuit:
    width = 2 * n_p + 2
    c0_arr, c1_arr, c2_arr = compute_piecewise_table(
        wp=wp, d1_bits=d1_bits, d2_bits=d2_bits, c2_drop_lsbs=c2_drop_lsbs
    )
    w_c1 = wp - 8
    w_c2_full = wp - 15
    w_c2 = w_c2_full - c2_drop_lsbs
    c2_s = c2_arr >> c2_drop_lsbs

    shift_hi = wp - 6
    v6_table = [int(c0_arr[4 * g + 3] >> shift_hi) for g in range(64)]
    c0_rem_arr = [int(c0_arr[i] - (v6_table[i >> 2] << shift_hi)) for i in range(256)]
    assert max(c0_rem_arr) < (1 << (shift_hi + 1))

    tab_init = [
        c0_rem_arr[i]
        | (int(c1_arr[i]) << (shift_hi + 1))
        | (int(c2_s[i]) << (shift_hi + 1 + w_c1))
        for i in range(256)
    ]
    tab_c2 = [int(c2_s[i]) for i in range(256)]
    rem_flags_tab = _rem_flags_table()

    b = km.CircuitBuilder()
    x = b.create_quantum_register(width, name='x')

    # 1. Normalize x and compress flags -> s_bin (4 qubits)
    flags = b.alloc_qubits(n_p)
    _normalize(b, False, x, flags)
    _thermometer_cnots(b, flags, inverse=False)
    s_bin = km.array([flags[0], flags[1], flags[3], flags[7]])
    rem_flags = km.array([flags[2], flags[4], flags[5], flags[6], flags[8], flags[9]])
    b.del_lookup(table=rem_flags_tab, address=s_bin, target=rem_flags)
    b.free(rem_flags)

    # 2. Align top bit of x: b_even = ~x[21]; if b_even, left_rotate(x, 1); invert x[0:14]
    b_even = b.alloc_qubits(1)[0]
    b.x(x[width - 1])
    b.cx(x[width - 1], b_even)
    b.x(x[width - 1])
    b.left_rotate(x, 1, control=b_even)
    b.x(x[:14])
    addr = km.array([*x[14:21], b_even])

    # 3. Single lookup (c0_rem, c1, c2_s) -> (y0_rem, r1, c2)
    y0_rem = b.alloc_qubits(shift_hi + 1)
    r1 = b.alloc_qubits(w_c1)
    c2 = b.alloc_qubits(w_c2)
    b.init_lookup(table=tab_init, address=addr, target=km.array([*y0_rem, *r1, *c2]))

    # 4. Step A: r1 += T(c2 * d_bar) with carry_a and early del_lookup(c2)
    live_a = 27 + (shift_hi + 1) + w_c1 + w_c2 + 1
    mc_a = max(6, target_qubits - live_a)
    carry_a = b.alloc_qubits(1)[0]
    r1_low_ext = km.array([*r1[:w_c2_full], carry_a])
    for m in range(1, d2_bits + 1):
        if m <= c2_drop_lsbs:
            add_shifted_pebbled(b, x[14 - m], c2, r1_low_ext[c2_drop_lsbs - m:], subtract=False, max_carries=mc_a)
        else:
            drop = m - c2_drop_lsbs
            if drop < w_c2:
                add_shifted_pebbled(b, x[14 - m], c2[drop:w_c2], r1_low_ext, subtract=False, max_carries=mc_a)
    add_single_bit_into(b, carry_a, r1[w_c2_full:w_c1], subtract=False)

    b.del_lookup(table=tab_c2, address=addr, target=c2)
    b.free(c2)

    # 5. Step B: y0_rem += T(r1 * d_bar) with hierarchical carries
    live_b = 27 + (shift_hi + 1) + w_c1 + 3
    mc_b = max(8, target_qubits - live_b)
    carry_b = b.alloc_qubits(1)[0]
    y0_low_ext = km.array([*y0_rem[:w_c1], carry_b])
    for m in range(1, 8):
        add_shifted_pebbled(b, x[14 - m], r1[m:w_c1], y0_low_ext, subtract=False, max_carries=mc_b)
    carry_b14 = b.alloc_qubits(1)[0]
    y0_b14_ext = km.array([*y0_rem[:w_c2_full], carry_b14])
    for m in range(8, d1_bits + 1):
        if m < w_c1:
            add_shifted_pebbled(b, x[14 - m], r1[m:w_c1], y0_b14_ext, subtract=False, max_carries=mc_b - 1)
    add_single_bit_into(b, carry_b14, y0_low_ext[w_c2_full:], subtract=False)
    add_single_bit_into(b, carry_b, y0_rem[w_c1:shift_hi + 1], subtract=False)

    # 6. Load y0_top6 via 6-bit group address addr[2:8] and add carry bit y0_rem[shift_hi]
    y0_top6 = b.alloc_qubits(6)
    b.init_lookup(table=v6_table, address=addr[2:], target=y0_top6)
    add_single_bit_into(b, y0_rem[shift_hi], y0_top6, subtract=False)
    y0 = km.array([*y0_rem[:shift_hi], *y0_top6])

    # 7. Denormalize and Phase
    top10 = b.alloc_qubits(n_p)
    y1 = km.array([*y0, *top10])
    b.left_rotate(y1[:wp + 1], 1, control=s_bin[0])
    b.left_rotate(y1[:wp + 3], 2, control=s_bin[1])
    b.left_rotate(y1[:wp + 7], 4, control=s_bin[2])
    b.left_rotate(y1[:wp + 10], 8, control=s_bin[3])
    for j in range(len(y1)):
        b.z_pow(y1[j], Fraction(2) ** (j - len(y1) + 1))
    b.right_rotate(y1[:wp + 10], 8, control=s_bin[3])
    b.right_rotate(y1[:wp + 7], 4, control=s_bin[2])
    b.right_rotate(y1[:wp + 3], 2, control=s_bin[1])
    b.right_rotate(y1[:wp + 1], 1, control=s_bin[0])
    b.free(top10)

    # 8. Uncompute y0_top6 before Uncompute Step B
    add_single_bit_into(b, y0_rem[shift_hi], y0_top6, subtract=True)
    b.del_lookup(table=v6_table, address=addr[2:], target=y0_top6)
    b.free(y0_top6)

    # 9. Uncompute Step B and free y0_rem
    add_single_bit_into(b, carry_b, y0_rem[w_c1:shift_hi + 1], subtract=True)
    add_single_bit_into(b, carry_b14, y0_low_ext[w_c2_full:], subtract=True)
    for m in reversed(range(8, d1_bits + 1)):
        if m < w_c1:
            add_shifted_pebbled(b, x[14 - m], r1[m:w_c1], y0_b14_ext, subtract=True, max_carries=mc_b - 1)
    b.free(km.array([carry_b14]))
    for m in reversed(range(1, 8)):
        add_shifted_pebbled(b, x[14 - m], r1[m:w_c1], y0_low_ext, subtract=True, max_carries=mc_b)
    b.free(km.array([carry_b]))

    b.del_lookup(table=c0_rem_arr, address=addr, target=y0_rem)
    b.free(y0_rem)

    # 10. Recompute c2 via 2-block Select-Swap inside freed y0_rem space, uncompute Step A, and joint del_lookup
    c2_blk0 = b.alloc_qubits(w_c2)
    c2_blk1 = b.alloc_qubits(w_c2)
    super_c2 = [int(c2_s[2 * hi]) | (int(c2_s[2 * hi + 1]) << w_c2) for hi in range(128)]
    b.init_lookup(table=super_c2, address=addr[1:], target=km.array([*c2_blk0, *c2_blk1]))
    b.cswap(addr[0], c2_blk0, c2_blk1)

    add_single_bit_into(b, carry_a, r1[w_c2_full:w_c1], subtract=True)
    for m in reversed(range(1, d2_bits + 1)):
        if m <= c2_drop_lsbs:
            add_shifted_pebbled(b, x[14 - m], c2_blk0, r1_low_ext[c2_drop_lsbs - m:], subtract=True, max_carries=15)
        else:
            drop = m - c2_drop_lsbs
            if drop < w_c2:
                add_shifted_pebbled(b, x[14 - m], c2_blk0[drop:w_c2], r1_low_ext, subtract=True, max_carries=15)
    b.free(km.array([carry_a]))

    tab_final = []
    for i in range(256):
        hi = i >> 1
        lo = i & 1
        b0 = int(c2_s[2 * hi + lo])
        b1 = int(c2_s[2 * hi + (1 - lo)])
        tab_final.append(int(c1_arr[i]) | (b0 << w_c1) | (b1 << (w_c1 + w_c2)))
    b.del_lookup(table=tab_final, address=addr, target=km.array([*r1, *c2_blk0, *c2_blk1]))
    b.free(c2_blk1)
    b.free(c2_blk0)
    b.free(r1)

    # 11. Uncompute alignment of x
    b.x(x[:14])
    b.right_rotate(x, 1, control=b_even)
    b.x(x[width - 1])
    b.cx(x[width - 1], b_even)
    b.x(x[width - 1])
    b.free(km.array([b_even]))

    # 12. Restore flags from s_bin and un-normalize x
    rem_flags_rev = b.alloc_qubits(6)
    b.init_lookup(table=rem_flags_tab, address=s_bin, target=rem_flags_rev)
    flags_rev = km.array([
        s_bin[0], s_bin[1], rem_flags_rev[0], s_bin[2],
        rem_flags_rev[1], rem_flags_rev[2], rem_flags_rev[3], s_bin[3],
        rem_flags_rev[4], rem_flags_rev[5],
    ])
    _thermometer_cnots(b, flags_rev, inverse=True)
    _normalize(b, True, x, flags_rev)
    b.free(flags_rev)

    return b.finish_circuit()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', help='Output .kmb or .kmx file path.')
    parser.add_argument('--wp', type=int, default=28)
    parser.add_argument('--d1', type=int, default=14)
    parser.add_argument('--d2', type=int, default=9)
    parser.add_argument('--c2-drop', type=int, default=3)
    parser.add_argument('--target-qubits', type=int, default=90)
    parser.add_argument('--format', choices=('kmb', 'kmx'), default='kmb')
    args = parser.parse_args()

    circuit = build_opt_invsqrt_circuit(
        wp=args.wp,
        d1_bits=args.d1,
        d2_bits=args.d2,
        c2_drop_lsbs=args.c2_drop,
        target_qubits=args.target_qubits,
    )
    circuit.to_file(args.output, format=args.format)
    toffolis = circuit.operation_statistics().num_toffoli_like
    print(f'Wrote {args.output}: {toffolis:,} max Toffolis, {circuit.num_qubits} qubits.')


if __name__ == '__main__':
    main()
