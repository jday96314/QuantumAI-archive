"""Verify running-controlled swap logic in classical simulation."""

import kickmix as km
import numpy as np
from simulate_classical_eea import run_exact_eea_trace
from verify_submission_dqi import FieldTables, ProblemSpecDQI, sample_shots

def run_running_swap_eea(syndromes_list: list[int], gf: km.GF2Field, n_iters: int = 210) -> dict:
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

    for step in range(n - 1):
        if reg_a[m - 1] == 0:
            sub = reg_a[3 + step : m]
            reg_a[3 + step : m] = sub[1:] + sub[:1]
            n_A -= 1

    halt_it = None
    for it in range(n_iters):
        is_c3 = 1 if counter == 3 else 0
        is_hc_zero = 1 if halting_counter == 0 else 0

        if is_hc_zero:
            c0 = 1 if counter == 0 else 0
            c1 = 1 if counter == 1 else 0
            c2 = 1 if counter == 2 else 0
            c3 = 1 if counter == 3 else 0

            # 1. flag
            if c0:
                if n_A < n_B and n_q == 0: flag ^= 1
                if n_A == n_B: flag ^= 1
            elif c1:
                if n_A == n_B: flag ^= 1
                flag ^= 1
                aux = (1 if n_B == 1 else 0) ^ 1
                if aux and reg_b[m - 2] == 0: flag ^= 1
            elif c2:
                flag ^= 1
                aux = (1 if n_b == 0 else 0) ^ 1
                if aux and reg_b[0] == 0: flag ^= 1
                if n_a == n_b + 1: flag ^= 1
            elif c3:
                if n_a == n_b: flag ^= 1
                if n_q == 1: flag ^= 1

            # 2. mul/add
            lim_a = min(t, it // 6 + 1)
            q_i = 0
            if c0 and reg_a[m - 1] != 0:
                q_i = gf.div(reg_b[m - 1], reg_a[m - 1])
            nq_tmp = n_q + n_a - (1 if c3 else 0)
            for i in range(m - 1, -1, -1):
                mul = (1 if (c3 and i + 1 <= n_a and i < lim_a) else 0) ^ (1 if (c0 and m - i <= n_A) else 0)
                if c0 or c3:
                    if i == nq_tmp:
                        q_i, reg_a[i] = reg_a[i], q_i
                if mul and q_i != 0 and reg_a[i] != 0:
                    reg_b[i] ^= gf.mul(q_i, reg_a[i])
            if c3 and reg_a[0] != 0:
                q_i ^= gf.div(reg_b[0], reg_a[0])

            # 3. shift (single unary iteration)
            c0_nf = 1 if (c0 and flag == 0) else 0
            c3_nf = 1 if (c3 and flag == 0) else 0
            ctrl_nB = c0_nf ^ c1
            ctrl_nb = c2 ^ c3_nf
            lim_b = min(t, it // 2 + 1)
            for i in range(m - 1, 0, -1):
                do_shift = ((1 if (ctrl_nB and n_B >= m - i + 1) else 0)
                            ^ (1 if (ctrl_nb and i <= lim_b and n_b >= i) else 0))
                if do_shift:
                    reg_b[i - 1], reg_b[i] = reg_b[i], reg_b[i - 1]

            # 4. metadata
            if c0: n_q += 1
            if c0_nf: n_B -= 1
            if c1: n_B -= 1
            if c2: n_b += 1
            if c3_nf: n_b += 1
            if c3: n_q -= 1

            # 5. Swap conditioned on running
            running = 0
            if is_c3 and flag == 1:
                busy = 1 if n_b < n_B else 0
                running = busy ^ (1 if (busy and n_A == 0) else 0)
                if running:
                    n_A, n_B = n_B, n_A
                    n_a, n_b = n_b, n_a
                    reg_a, reg_b = reg_b, reg_a
                else:
                    if halt_it is None:
                        halt_it = it + 1
                    halting_counter += 1

        if flag == 1:
            counter = (counter + 1) % 4

    return {
        "halt_it": halt_it,
        "reg_b": reg_b,
        "n_b": n_b,
        "n_B": n_B,
    }

def main():
    gf = km.GF2Field(12)
    spec = ProblemSpecDQI.benchmark_instance(num_shots=100)
    tables = FieldTables(gf)
    rng = np.random.default_rng(20261023)

    for w in [1, 2, 5, 10, 20, 28, 32, 34, 35]:
        probs = [0.0] * 36
        probs[w] = 1.0
        synd_batch, _, _ = sample_shots(spec, tables, rng, 20, weight_probs=probs)
        for s_row in synd_batch:
            s_list = s_row.tolist()
            ref = run_exact_eea_trace(s_list, gf, n_iters=222)
            opt = run_running_swap_eea(s_list, gf, n_iters=210)
            assert opt["halt_it"] is not None and opt["halt_it"] <= 210
            assert opt["n_b"] == ref["final_n_b"]
            assert opt["n_B"] == ref["final_n_B"]
            assert opt["reg_b"] == ref["reg_b"], f"Mismatch on weight {w}!"
        print(f"weight={w:2d}: 20/20 shots EXACT match on reg_b with n_iters=210!")

if __name__ == "__main__":
    main()
