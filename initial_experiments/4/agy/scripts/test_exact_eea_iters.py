"""Print round-by-round trace of Proos-Zalka EEA on weight=1, 2, 35 errors."""

from __future__ import annotations

import kickmix as km
import numpy as np
from verify_submission_dqi import FieldTables, ProblemSpecDQI, sample_shots


def trace_rounds(syndromes_list: list[int], gf: km.GF2Field, verbose: bool = False) -> tuple[int, list[tuple]]:
    n = len(syndromes_list)
    m = n + 3
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

    rounds = []
    halt_it = None
    for it in range(222):
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
            if c0 and reg_a[m - 1] != 0:
                q_i = gf.div(reg_b[m - 1], reg_a[m - 1])
            nq_tmp = n_q + n_a - (1 if c3 else 0)
            for i in range(m - 1, -1, -1):
                mul = (1 if (c3 and i + 1 <= n_a) else 0) ^ (1 if (c0 and m - i <= n_A) else 0)
                if c0 or c3:
                    if i == nq_tmp:
                        q_i, reg_a[i] = reg_a[i], q_i
                if mul and q_i != 0 and reg_a[i] != 0:
                    reg_b[i] ^= gf.mul(q_i, reg_a[i])
            if c3 and reg_a[0] != 0:
                q_i ^= gf.div(reg_b[0], reg_a[0])
            assert q_i == 0, f"it={it} c={counter}: q_i={q_i} not cleared!"

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

            if c3 and flag == 1:
                n_A, n_B = n_B, n_A
                n_a, n_b = n_b, n_a
                reg_a, reg_b = reg_b, reg_a
                rounds.append((it + 1, n_A, n_B, n_a, n_b))
                if verbose:
                    print(f"  End of round at it+1={it+1:3d}: n_A={n_A:2d}, n_B={n_B:2d}, n_a={n_a:2d}, n_b={n_b:2d}")

        if is_c3 and flag == 1:
            busy = 1 if n_b < n_B else 0
            running = busy ^ (1 if (busy and n_A == 0) else 0)
            if not running:
                if halt_it is None:
                    halt_it = it + 1
                halting_counter += 1

        if flag == 1:
            counter = (counter + 1) % 4

    return halt_it, rounds


def main() -> None:
    gf = km.GF2Field(12)
    spec = ProblemSpecDQI.benchmark_instance(num_shots=256)
    tables = FieldTables(gf)
    rng = np.random.default_rng(20261022)
    syndromes, _, _ = sample_shots(spec, tables, rng, 16)

    print("Sampled DQI shot 0:")
    halt_it, rounds = trace_rounds(syndromes[0].tolist(), gf, verbose=True)
    print(f"Shot 0 halt_it = {halt_it}")

    # Also test weight 1 error using sample_shots with weight_probs=[0, 1, 0, ...]
    w1 = [0.0] * 36
    w1[1] = 1.0
    synd_w1, _, _ = sample_shots(spec, tables, rng, 4, weight_probs=w1)
    print("\nWeight 1 shot 0:")
    halt_it_w1, _ = trace_rounds(synd_w1[0].tolist(), gf, verbose=True)
    print(f"Weight 1 shot 0 halt_it = {halt_it_w1}")

    w2 = [0.0] * 36
    w2[2] = 1.0
    synd_w2, _, _ = sample_shots(spec, tables, rng, 4, weight_probs=w2)
    print("\nWeight 2 shot 0:")
    halt_it_w2, _ = trace_rounds(synd_w2[0].tolist(), gf, verbose=True)
    print(f"Weight 2 shot 0 halt_it = {halt_it_w2}")


if __name__ == "__main__":
    main()
