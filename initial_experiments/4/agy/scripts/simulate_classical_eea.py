"""Exact simulation and worst-case iteration analysis of the Proos-Zalka EEA."""

from __future__ import annotations

import kickmix as km
import numpy as np


def run_exact_eea_trace(syndromes_list: list[int], gf: km.GF2Field, n_iters: int = 222) -> dict:
    """Runs the exact step-by-step logic of Proos-Zalka EEA on a syndrome vector."""
    n = len(syndromes_list)
    m = n + 3
    t = n // 2

    reg_a = [0] * m
    reg_b = [0] * m
    reg_a[0] = 1
    for i in range(n):
        reg_a[3 + i] = syndromes_list[i]
    reg_b[m - 1] = 1

    n_A = n
    n_B = n + 1
    n_a = 1
    n_b = 0
    n_q = 0
    flag = 1
    counter = 0
    halting_counter = 0

    # Initial shift loop
    for step in range(n - 1):
        ctrl = 1 if reg_a[m - 1] == 0 else 0
        if ctrl:
            # left rotate reg_a[3 + step : m] by 1
            sub = reg_a[3 + step : m]
            reg_a[3 + step : m] = sub[1:] + sub[:1]
            n_A -= 1

    halt_it = None
    active_c0_nA = []
    active_c3_na = []

    for it in range(n_iters):
        is_c3 = 1 if counter == 3 else 0
        is_hc_zero = 1 if halting_counter == 0 else 0

        if is_hc_zero:
            c0 = 1 if counter == 0 else 0
            c1 = 1 if counter == 1 else 0
            c2 = 1 if counter == 2 else 0
            c3 = 1 if counter == 3 else 0

            # 1. _update_flag
            if c0:
                if n_A < n_B and n_q == 0:
                    flag ^= 1
                if n_A == n_B:
                    flag ^= 1
            elif c1:
                if n_A == n_B:
                    flag ^= 1
                flag ^= 1
                aux = (1 if n_B == 1 else 0) ^ 1
                if aux and reg_b[m - 2] == 0:
                    flag ^= 1
            elif c2:
                flag ^= 1
                aux = (1 if n_b == 0 else 0) ^ 1
                if aux and reg_b[0] == 0:
                    flag ^= 1
                if n_a == n_b + 1:
                    flag ^= 1
            elif c3:
                if n_a == n_b:
                    flag ^= 1
                if n_q == 1:
                    flag ^= 1

            # 2. _multiply_and_add_or_sub
            q_i = 0
            if c0:
                active_c0_nA.append((it, n_A, n_B, n_a, n_b, n_q))
                if reg_a[m - 1] != 0:
                    q_i = gf.div(reg_b[m - 1], reg_a[m - 1])
            nq_tmp = n_q + n_a - (1 if c3 else 0)
            for i in range(m - 1, -1, -1):
                mul = (1 if (c3 and i + 1 <= n_a) else 0) ^ (1 if (c0 and m - i <= n_A) else 0)
                if c0 or c3:
                    if i == nq_tmp:
                        q_i, reg_a[i] = reg_a[i], q_i
                if mul and q_i != 0 and reg_a[i] != 0:
                    reg_b[i] ^= gf.mul(q_i, reg_a[i])
            if c3:
                active_c3_na.append((it, n_A, n_B, n_a, n_b, n_q))
                if reg_a[0] != 0:
                    q_i ^= gf.div(reg_b[0], reg_a[0])

            # 3. _shift_b_array_right
            c0_nf = 1 if (c0 and flag == 0) else 0
            c3_nf = 1 if (c3 and flag == 0) else 0
            ctrl_nB = c0_nf ^ c1
            for i in range(m - 1, 0, -1):
                do_shift = (
                    (1 if (ctrl_nB and n_B >= m - i + 1) else 0)
                    ^ (1 if (c2 and n_b >= i) else 0)
                    ^ (1 if (c3_nf and n_B < m - i) else 0)
                )
                if do_shift:
                    reg_b[i - 1], reg_b[i] = reg_b[i], reg_b[i - 1]

            # 4. _update_metadata
            if c0:
                n_q += 1
            if c0_nf:
                n_B -= 1
            if c1:
                n_B -= 1
            if c2:
                n_b += 1
            if c3_nf:
                n_b += 1
            if c3:
                n_q -= 1

            # 5. _swap_a_and_b
            if c3 and flag == 1:
                n_A, n_B = n_B, n_A
                n_a, n_b = n_b, n_a
                reg_a, reg_b = reg_b, reg_a

        # _update_halting_counter
        if is_c3 and flag == 1:
            busy = 1 if n_b < n_B else 0
            running = busy ^ (1 if (busy and n_A == 0) else 0)
            if not running:
                if halt_it is None:
                    halt_it = it + 1
                halting_counter += 1

        # _update_counter
        if flag == 1:
            counter = (counter + 1) % 4

    return {
        "halt_it": halt_it,
        "final_n_A": n_A,
        "final_n_B": n_B,
        "final_n_a": n_a,
        "final_n_b": n_b,
        "final_n_q": n_q,
        "final_counter": counter,
        "final_flag": flag,
        "reg_b": reg_b,
        "reg_b_36": reg_b[36],
        "active_c0_nA": active_c0_nA,
        "active_c3_na": active_c3_na,
    }


def main() -> None:
    gf = km.GF2Field(12)
    prim = gf.primitive_element
    antilog = np.array([gf.pow(prim, k) for k in range(4095)], dtype=np.int64)
    rng = np.random.default_rng(20261022)

    max_halt_it = 0
    for weight in [0, 1, 2, 5, 10, 20, 30, 34, 35]:
        w_max_it = 0
        for trial in range(50 if weight > 0 else 1):
            if weight == 0:
                synd = [0] * 70
            else:
                pos = rng.choice(4095, size=weight, replace=False)
                vals = rng.integers(1, 4095, size=weight, endpoint=True)
                synd = []
                for i in range(70):
                    s = 0
                    for p, v in zip(pos, vals):
                        s ^= gf.mul(int(v), int(antilog[((i + 1) * int(p)) % 4095]))
                    synd.append(s)
            res = run_exact_eea_trace(synd, gf)
            assert res["final_n_q"] == 0, f"n_q={res['final_n_q']}"
            assert res["final_flag"] == 1, f"flag={res['final_flag']}"
            assert res["reg_b_36"] == 0, f"reg_b[36]={res['reg_b_36']}"
            w_max_it = max(w_max_it, res["halt_it"])
            max_halt_it = max(max_halt_it, res["halt_it"])
        print(f"weight={weight:2d}: max halt_it = {w_max_it}")

    print(f"Overall max halt_it across random trials: {max_halt_it}")


if __name__ == "__main__":
    main()
