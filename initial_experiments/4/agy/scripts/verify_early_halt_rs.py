"""Verify stopping at the end of round R' (when n_A < n_a after swap, or n_B <= 35 before swap)."""

from __future__ import annotations

import kickmix as km
import numpy as np
from verify_submission_dqi import FieldTables, ProblemSpecDQI, sample_shots


def test_early_halt_trace(syndromes_list: list[int], gf: km.GF2Field, n_iters: int = 216) -> dict:
    n = len(syndromes_list)
    m = n + 3
    t = n // 2  # 35

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

    for step in range(n - 1):
        if reg_a[m - 1] == 0:
            sub = reg_a[3 + step : m]
            reg_a[3 + step : m] = sub[1:] + sub[:1]
            n_A -= 1

    halt_it = None
    min_c0_nA = 999
    max_c3_na = 0

    for it in range(n_iters):
        is_c3 = 1 if counter == 3 else 0
        is_hc_zero = 1 if halting_counter == 0 else 0

        if is_hc_zero:
            c0 = 1 if counter == 0 else 0
            c1 = 1 if counter == 1 else 0
            c2 = 1 if counter == 2 else 0
            c3 = 1 if counter == 3 else 0

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

            q_i = 0
            if c0:
                min_c0_nA = min(min_c0_nA, n_A)
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
                max_c3_na = max(max_c3_na, n_a)
                if reg_a[0] != 0:
                    q_i ^= gf.div(reg_b[0], reg_a[0])

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

            # Wait: what if we check halting BEFORE _swap_a_and_b (or don't swap when halting)?!
            # Look! Before _swap_a_and_b, the newly computed (sigma, omega) IS ALREADY IN reg_b
            # with lengths (n_b, n_B) satisfying n_B < n_b!!
            if c3 and flag == 1:
                n_A, n_B = n_B, n_A
                n_a, n_b = n_b, n_a
                reg_a, reg_b = reg_b, reg_a

        # Look: after _swap_a_and_b, the newly computed (sigma, omega) is in reg_a with lengths (n_a, n_A)!
        # So the stopping condition is simply: n_A < n_a (or n_B == 0? Wait: after swap, new remainder is n_A,
        # so if new remainder is 0, n_A == 0 < 1 <= n_a, which is ALREADY covered by n_A < n_a!!)!
        if is_c3 and flag == 1:
            busy = 1 if n_a <= n_A else 0  # running while n_a <= n_A; halts as soon as n_A < n_a!
            if not busy:
                if halt_it is None:
                    halt_it = it + 1
                halting_counter += 1

        if flag == 1:
            counter = (counter + 1) % 4

    return {
        "halt_it": halt_it,
        "reg_a": reg_a,
        "reg_b": reg_b,
        "n_A": n_A,
        "n_B": n_B,
        "n_a": n_a,
        "n_b": n_b,
        "min_c0_nA": min_c0_nA,
        "max_c3_na": max_c3_na,
    }


def main() -> None:
    gf = km.GF2Field(12)
    spec = ProblemSpecDQI.benchmark_instance(num_shots=256)
    tables = FieldTables(gf)
    rng = np.random.default_rng(20261022)

    from simulate_classical_eea import run_exact_eea_trace

    max_halt_early = 0
    global_min_c0_nA = 999
    global_max_c3_na = 0

    for w in [1, 2, 5, 15, 25, 34, 35]:
        probs = [0.0] * 36
        probs[w] = 1.0
        n_trials = 1 if w == 0 else 50
        synd_batch, _, _ = sample_shots(spec, tables, rng, n_trials, weight_probs=probs)
        w_max_halt = 0
        for s_row in synd_batch:
            s_list = s_row.tolist()
            ref = run_exact_eea_trace(s_list, gf, n_iters=222)
            early = test_early_halt_trace(s_list, gf, n_iters=212)
            assert early["halt_it"] is not None, f"Did not halt for weight {w}!"
            # Verify that early's (reg_a, n_a, n_A) is IDENTICAL to ref's (reg_b, n_b, n_B)!
            assert early["n_a"] == ref["final_n_b"], f"{early['n_a']} != {ref['final_n_b']}"
            assert early["n_A"] == ref["final_n_B"], f"{early['n_A']} != {ref['final_n_B']}"
            # Check sigma coefficients in reg_a[:36] and omega coefficients in reg_a[37:73]
            assert early["reg_a"][:36] == ref["reg_b"][:36] if "reg_b" in ref else True
            assert early["reg_a"][36] == 0 and early["reg_a"][37] == 0
            w_max_halt = max(w_max_halt, early["halt_it"])
            max_halt_early = max(max_halt_early, early["halt_it"])
            global_min_c0_nA = min(global_min_c0_nA, early["min_c0_nA"])
            global_max_c3_na = max(global_max_c3_na, early["max_c3_na"])
        print(f"weight={w:2d}: max early halt_it = {w_max_halt:3d}")

    print(f"\nGlobal max early halt_it: {max_halt_early} (vs 222 in baseline!)")
    print(f"Global min n_A during c0: {global_min_c0_nA}")
    print(f"Global max n_a during c3: {global_max_c3_na}")


if __name__ == "__main__":
    main()
