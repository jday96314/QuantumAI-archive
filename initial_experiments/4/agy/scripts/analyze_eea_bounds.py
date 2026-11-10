"""Exact reachability analysis of the Proos-Zalka Synchronized EEA state machine.

Computes exact per-iteration bounds on (n_A, n_B, n_a, n_b, n_q, counter, flag)
while halting_counter == 0.
"""

from __future__ import annotations


def analyze_eea(n_syndromes: int = 70, n_iters: int = 222) -> None:
    n = n_syndromes
    t = n // 2  # 35

    # State: (n_A, n_B, n_a, n_b, n_q, counter, flag, hc_zero)
    # Initially after _initial_shift_loop:
    # n_A can be 0..n (if all syndromes 0, n_A = 1? Wait: shift_step runs n-1=69 times,
    # each time if reg_a[m-1]==0 it decrements n_A by 1, so n_A in {1..70}).
    # n_B = n + 1 = 71, n_a = 1, n_b = 0, n_q = 0, counter = 0, flag = 1, hc_zero = True.
    reachable: set[tuple[int, int, int, int, int, int, int]] = set()
    for n_A_init in range(1, n + 1):
        reachable.add((n_A_init, n + 1, 1, 0, 0, 0, 1))

    c0_nA_min = []
    c0_nA_max = []
    c3_na_max = []
    q_idx_min = []
    q_idx_max = []

    total_c03_positions = 0
    total_baseline_positions = 0

    for it in range(n_iters):
        next_reachable: set[tuple[int, int, int, int, int, int, int]] = set()

        it_c0_nA = []
        it_c3_na = []
        it_q_idx = []

        for (n_A, n_B, n_a, n_b, n_q, counter, flag) in reachable:
            # Step 1: _update_flag
            # Depending on counter:
            possible_flags: list[int] = []
            if counter == 0:
                f = flag
                if n_A < n_B and n_q == 0:
                    f ^= 1
                if n_A == n_B:
                    f ^= 1
                possible_flags = [f]
            elif counter == 1:
                f_base = flag
                if n_A == n_B:
                    f_base ^= 1
                f_base ^= 1
                if n_B == 1:
                    # aux is False, so reg_b[m-2] is not checked
                    possible_flags = [f_base]
                else:
                    # reg_b[m-2] == 0 can be True or False
                    possible_flags = [f_base, f_base ^ 1]
            elif counter == 2:
                f_base = flag ^ 1
                if n_a == n_b + 1:
                    f_base ^= 1
                if n_b == 0:
                    possible_flags = [f_base]
                else:
                    possible_flags = [f_base, f_base ^ 1]
            elif counter == 3:
                f = flag
                if n_a == n_b:
                    f ^= 1
                if n_q == 1:
                    f ^= 1
                possible_flags = [f]

            for f in possible_flags:
                # Record active bounds during _multiply_and_add_or_sub
                if counter == 0:
                    it_c0_nA.append(n_A)
                    it_q_idx.append(n_q + n_a)
                elif counter == 3:
                    it_c3_na.append(n_a)
                    it_q_idx.append(n_q + n_a - 1)

                # Step 4: _update_metadata
                nA2, nB2, na2, nb2, nq2 = n_A, n_B, n_a, n_b, n_q
                if counter == 0:
                    nq2 += 1
                    if f == 0:
                        nB2 -= 1
                elif counter == 1:
                    nB2 -= 1
                elif counter == 2:
                    nb2 += 1
                elif counter == 3:
                    if f == 0:
                        nb2 += 1
                    nq2 -= 1

                # Step 5: _swap_a_and_b
                if counter == 3 and f == 1:
                    nA2, nB2 = nB2, nA2
                    na2, nb2 = nb2, na2

                # Restrict to valid RS decoding trajectories where locator degree <= t + 1 = 36
                if na2 > t + 1 or nb2 > t + 1 or nA2 < 0 or nB2 < 0:
                    continue

                # Step 6: _update_halting_counter
                # Halts if counter == 3 and f == 1 and not (nb2 < nB2 and nA2 != 0)
                if counter == 3 and f == 1:
                    busy = nb2 < nB2
                    running = busy and (nA2 != 0)
                    if not running:
                        # Halts! halting_counter becomes 1, so hc_zero is False on subsequent steps
                        continue

                # Step 7: _update_counter
                c2 = (counter + 1) % 4 if f == 1 else counter
                next_reachable.add((nA2, nB2, na2, nb2, nq2, c2, f))

        min_nA = min(it_c0_nA) if it_c0_nA else 0
        max_nA = max(it_c0_nA) if it_c0_nA else 0
        max_na = max(it_c3_na) if it_c3_na else 0
        min_q = min(it_q_idx) if it_q_idx else 0
        max_q = max(it_q_idx) if it_q_idx else 0

        c0_nA_min.append(min_nA)
        c0_nA_max.append(max_nA)
        c3_na_max.append(max_na)
        q_idx_min.append(min_q)
        q_idx_max.append(max_q)

        # Count how many positions i in 0..72 need multiplication at iteration it:
        # Need multiplication if (c3 reachable and i < max_na) or (c0 reachable and i >= 73 - max_nA)
        m = n + 3
        active_pos = 0
        for i in range(m):
            cond_a = (len(it_c3_na) > 0) and (i < max_na)
            cond_A = (len(it_c0_nA) > 0) and (i >= m - max_nA)
            if cond_a or cond_A:
                active_pos += 1
        total_c03_positions += active_pos
        total_baseline_positions += m

        if it < 15 or it % 20 == 0 or it == n_iters - 1:
            print(
                f"it={it:3d}: states={len(reachable):6d} | "
                f"c0 n_A in [{min_nA:2d}, {max_nA:2d}] | "
                f"c3 n_a <= {max_na:2d} | "
                f"q_idx in [{min_q:2d}, {max_q:2d}] | "
                f"active_mul_pos={active_pos:2d}/{m}"
            )

        reachable = next_reachable

    print(f"\nTotal active mul positions across {n_iters} iters: {total_c03_positions:,} vs {total_baseline_positions:,}")


if __name__ == "__main__":
    analyze_eea()
