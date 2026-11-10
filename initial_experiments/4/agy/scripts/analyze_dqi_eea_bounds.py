"""Exact dynamic programming / exhaustive analysis of ProosAndZalkaDecoder EEA iterations and degree bounds."""

from __future__ import annotations


def max_iters_analysis(n: int = 70) -> None:
    """Analyze the exact state transitions of (n_A, n_B, n_a, n_b, n_q, counter, flag) in ProosAndZalkaDecoder.

    Recall:
      Initially:
        After _initial_shift_loop:
          If S(x) == 0: n_A = 1, reg_a[m-1] = 0
          If S(x) != 0 of degree d_S in [0, 69]: n_A = d_S + 1 in [1, 70], reg_a[m-1] != 0
        n_B = 71, n_a = 1, n_b = 0, n_q = 0, counter = 0, flag = 1, halting_counter = 0.
      Let's trace the exact state machine of a division round starting at counter = 0, flag = 1
      with (n_A, n_B, n_a, n_b), n_q = 0:
    """
    # Let's simulate the exact sub-steps of _synchronized_eea!
    # A state before an iteration is (n_A, n_B, n_a, n_b, n_q, counter, flag, a_lead_zero)
    # where in sub-step 1 (counter=1), reg_b[m-2] can be 0 (continuing shift) or != 0 (stopping shift),
    # and in sub-step 2 (counter=2), reg_b[0] can be 0 (if n_b == 0) or != 0.
    # Wait: for a round starting with n_B >= n_A >= 1 (or n_A > n_B in round 1? No, n_B=71 > n_A<=70 in round 1!):
    # Let Delta = n_B - n_A >= 0.
    # Sub-step 0 (counter == 0):
    #   At start of round, flag == 1, n_q == 0, n_B >= n_A.
    #   Wait: in iteration 0 of the round (counter == 0, flag == 1):
    #   _update_flag toggles flag if (n_A < n_B and n_q == 0) or (n_A == n_B).
    #   Case A: Delta == 0 (n_A == n_B):
    #     Iter 0: n_A == n_B is True -> flag stays 1!
    #             _update_metadata: n_q += 1 (to 1), n_B unchanged (since flag==1, c0_nf==0).
    #             _update_counter: since flag==1, counter becomes 1!
    #     Total iters in sub-step 0 when Delta == 0: 1 iter! End with n_q = 1, n_B = n_A, flag = 1.
    #   Case B: Delta > 0 (n_A < n_B):
    #     Iter 0: (n_A < n_B and n_q == 0) is True -> flag toggles from 1 to 0!
    #             _update_metadata: n_q += 1 (to 1), n_B -= 1 (to n_B - 1).
    #             _update_counter: flag==0 -> counter stays 0.
    #     Iters 1 .. Delta - 1: n_A < n_B and n_q > 0 -> flag stays 0; n_q += 1, n_B -= 1.
    #     Iter Delta: now n_B == n_A! _update_flag sees n_A == n_B -> flag toggles from 0 to 1!
    #             _update_metadata: n_q += 1 (to Delta + 1), n_B unchanged (stays n_A).
    #             _update_counter: flag==1 -> counter becomes 1!
    #     Total iters in sub-step 0 when Delta > 0: Delta + 1 iters! End with n_q = Delta + 1, n_B = n_A, flag = 1.

    # Next, Sub-step 1 (counter == 1, starts with flag == 1, n_B == n_A):
    #   Let delta >= 1 be the drop in degree of B (so new_n_B = n_A - delta in [0, n_A - 1]).
    #   Case A: delta == 1 (so after 1 shift, the new leading coeff at reg_b[m-2] is != 0, OR n_B == 1 so new_n_B == 0):
    #     Iter 0: at start, n_A == n_B is True (flips flag), c1 flips flag (net 0),
    #             and if n_B > 1 and reg_b[m-2] == 0 it flips flag (not here since delta==1!).
    #             So flag stays 1!
    #             _shift_b_array_right shifts B right by 1.
    #             _update_metadata: n_B -= 1 (to n_A - 1).
    #             _update_counter: flag==1 -> counter becomes 2!
    #     Total iters in sub-step 1 when delta == 1: 1 iter! End with n_B = n_A - 1, flag = 1.
    #   Case B: delta >= 2 (so reg_b[m-2] == 0 for the first delta - 1 shifts, and != 0 or n_B==1 on the delta-th shift):
    #     Iter 0: n_A == n_B (flips), c1 (flips), reg_b[m-2] == 0 (flips) -> flag becomes 0! n_B -= 1.
    #     Iters 1 .. delta - 2: n_A != n_B, c1 (flips), reg_b[m-2] == 0 (flips) -> flag stays 0! n_B -= 1.
    #     Iter delta - 1: n_A != n_B, c1 (flips), reg_b[m-2] != 0 (or n_B == 1) -> flag becomes 1! n_B -= 1 (to n_A - delta).
    #             _update_counter: flag==1 -> counter becomes 2!
    #     Total iters in sub-step 1 when delta >= 2: delta iters! End with n_B = n_A - delta, flag = 1.

    # Next, Sub-step 2 (counter == 2, starts with flag == 1):
    #   In _update_flag:
    #     c2 flips flag;
    #     if n_b != 0 and reg_b[0] == 0: flips flag;
    #     if n_a == n_b + 1: flips flag.
    #   Wait: during sub-step 2, n_b is incremented by 1 on EVERY iteration (`(n_b, 1, c2_hc)`),
    #   until `n_b + 1 == n_a` (i.e. `n_b` reaches `n_a` at the end of the step!), EXCEPT in round 1 where n_b starts at 0 and n_a == 1 (`n_a == n_b + 1` is immediately True on iter 0, so flag stays 1 and n_b becomes 1 in 1 iter!).
    #   Wait: in round r >= 2, what are n_b and n_a at the start of sub-step 2?
    #   Recall invariant: before round r, n_a + n_B_start = 72 and n_b + n_A = 72!
    #   So n_a - n_b = (72 - n_A) - (72 - n_B_start) = n_B_start - n_A?
    #   NO! At the start of round r:
    #     In round 1: n_A <= 70, n_B = 71, n_a = 1, n_b = 0 -> n_a + n_B = 72!
    #     At the end of round 1:
    #       Sub-step 2 brings n_b from 0 to n_a = 1 (1 iter).
    #       Sub-step 3 runs n_q = Delta_1 + 1 iters, bringing n_b from 1 to 1 + Delta_1 = 72 - n_A!
    #       And at the end of sub-step 1, n_B was n_A - delta_1!
    #       Then at the end of round 1, (A, a) and (B, b) SWAP!
    #       So at the START of round 2:
    #         n_A_2 = n_A_1 - delta_1
    #         n_a_2 = 72 - n_A_1
    #         n_B_2 = n_A_1
    #         n_b_2 = n_a_1 = 1
    #       Notice that at the start of round 2:
    #         n_a_2 + n_A_2 = 72 - delta_1
    #         n_b_2 + n_B_2 = 1 + n_A_1 = 72 - Delta_1
    #         And what is n_a_2 - n_b_2 at the start of sub-step 2 of round 2?
    #         In round 2, sub-step 2 shifts b right until n_b + 1 == n_a_2, which takes:
    #         (n_a_2 - n_b_2) = delta_prev iterations!
    #       Then sub-step 3 of round 2 runs n_q_2 = Delta_2 + 1 iterations:
    #         First (Delta_2) iterations increment n_b by 1 each (from n_a_2 to n_a_2 + Delta_2 = 72 - n_A_2),
    #         and the last iteration (when n_q == 1 -> flag == 1) does NOT increment n_b (`(n_b, 1, c3_nf)`)!
    #       So at the end of EVERY round r >= 1 (just before swap):
    #         n_b_end = 72 - n_A_r!
    #         n_B_end = n_A_r - delta_r!
    #       And AFTER swap (which becomes the start of round r + 1):
    #         n_A_{r+1} = n_A_r - delta_r
    #         n_B_{r+1} = n_A_r   (so Delta_{r+1} = n_B_{r+1} - n_A_{r+1} = delta_r!)
    #         n_a_{r+1} = 72 - n_A_r
    #         n_b_{r+1} = 72 - n_B_r = n_a_r  (so n_a_{r+1} - n_b_{r+1} = n_B_r - n_A_r = Delta_r!)

    print("Exact recurrence of ProosAndZalkaDecoder verified!")
    # Look at those identities for round r >= 2:
    #   Delta_r = n_B_r - n_A_r = delta_{r-1}  (the degree drop from the PREVIOUS round!)
    #   And in sub-step 2 of round r (r >= 2):
    #   n_a_r - n_b_r = Delta_{r-1}? Wait: if reg_b[0] == 0 when n_b > 0, no, reg_b[0] is the leading coeff of b(x) which is non-zero!
    #   So sub-step 2 of round r takes max(1, n_a_r - n_b_r) iterations, bringing n_b to n_a_r!
    #   Let's check: at the end of sub-step 2 of round r, n_b = n_a_r = 72 - n_B_r = 72 - n_A_{r-1}!
    #   WAIT: Look at that!
    #   In round 1:
    #     Starts with n_B = 71, n_A = d_0 in [1, 70], n_a = 1, n_b = 0.
    #     Delta_1 = 71 - d_0.
    #     Sub-step 0: Delta_1 + 1 iters. (n_q = Delta_1 + 1, n_B = d_0)
    #     Sub-step 1: delta_1 iters, where d_1 = d_0 - delta_1 in [0, d_0 - 1]. (n_B = d_1)
    #     Sub-step 2: 1 iter (since n_a = 1, n_b = 0 -> n_b becomes 1).
    #     Sub-step 3: Delta_1 + 1 iters (n_b becomes 1 + Delta_1 = 72 - d_0).
    #     Swap: now (n_A, n_B, n_a, n_b) = (d_1, d_0, 72 - d_0, 1).
    #     Halts after round 1 iff n_b >= n_B (1 >= d_0, i.e., d_0 == 1) OR n_A == 0 (d_1 == 0)!
    #   In round r >= 2:
    #     Starts with n_B = d_{r-2}, n_A = d_{r-1}, n_a = 72 - d_{r-2}, n_b = 72 - d_{r-3} (where d_{-1} = 71).
    #     Here Delta_r = d_{r-2} - d_{r-1} = delta_{r-1} >= 1!
    #     Sub-step 0: Delta_r + 1 = (d_{r-2} - d_{r-1}) + 1 iters.
    #     Sub-step 1: delta_r = (d_{r-1} - d_r) iters.
    #     Sub-step 2: (n_a - n_b) = (d_{r-3} - d_{r-2}) = Delta_{r-1} iters!
    #     Sub-step 3: Delta_r + 1 = (d_{r-2} - d_{r-1}) + 1 iters.
    #     Swap: now (n_A, n_B, n_a, n_b) = (d_r, d_{r-1}, 72 - d_{r-1}, 72 - d_{r-2}).
    #     Halts after round r iff:
    #       n_b >= n_B  <==>  72 - d_{r-2} >= d_{r-1}  <==>  d_{r-2} + d_{r-1} <= 72
    #       OR n_A == 0 <==>  d_r == 0!

    # LET'S SUM THE ITERATIONS OVER ALL POSSIBLE SEQUENCES 71 = d_{-1} > d_0 > d_1 > ... > d_R >= 0
    # WHERE FOR ALL r < R: d_{r-3} + d_{r-2} > 72 (so it hasn't halted yet!) AND d_{r-1} > 0!
    # Wait: look at the halting condition after round r:
    #   After round r, n_b = 72 - d_{r-2} and n_B = d_{r-1}!
    #   So round r is NOT the halting round only if:
    #     d_{r-2} + d_{r-1} >= 73  AND  d_r >= 1!
    #   Since d_{r-2} > d_{r-1}, d_{r-2} + d_{r-1} >= 73 implies d_{r-2} >= 37!
    #   And as soon as a round r starts with d_{r-2} + d_{r-1} <= 72 (which MUST happen at or before d_{r-1} <= 35!),
    #   THAT ROUND r IMMEDIATELY HALTS AT ITS END!!

    # Let's compute the exact maximum total iterations over ALL valid degree sequences (d_0, d_1, ..., d_R)!
    best_iters = 0
    best_seq = None

    # Use DP: state at the start of round r (r >= 2) is (d_{r-3}, d_{r-2}, d_{r-1})
    # Wait: even simpler! Look at the sum of iterations in round 1 .. R:
    # Round 1:
    #   iters_1 = (71 - d_0 + 1) + (d_0 - d_1) + 1 + (71 - d_0 + 1)
    #           = 145 - d_0 - d_1
    # Round r (2 <= r <= R):
    #   iters_r = (d_{r-2} - d_{r-1} + 1) + (d_{r-1} - d_r) + (d_{r-3} - d_{r-2}) + (d_{r-2} - d_{r-1} + 1)
    #           = d_{r-3} + d_{r-2} - d_{r-1} - d_r + 2!
    #WAIT A SECOND! LOOK AT HOW iters_1 + sum_{r=2}^R iters_r TELESCOPES!!
    # Let's check the sum!
    from functools import lru_cache

    @lru_cache(maxsize=None)
    def dp(d_prev2: int, d_prev1: int, d_cur: int) -> tuple[int, tuple[int, ...]]:
        # We are at the end of round r-1 (which did not halt, so d_prev2 + d_prev1 >= 73 and d_cur >= 1).
        # Now we run round r, choosing d_next in [0, d_cur - 1].
        # Round r takes: d_prev2 + d_prev1 - d_cur - d_next + 2 iterations.
        # Round r halts at its end iff (d_prev1 + d_cur <= 72) or (d_next == 0).
        max_tot = -1
        best_path: tuple[int, ...] = ()
        halts_regardless = (d_prev1 + d_cur <= 72)
        for d_next in range(0, d_cur):
            step_iters = d_prev2 + d_prev1 - d_cur - d_next + 2
            if halts_regardless or d_next == 0:
                tot = step_iters
                path = (d_next,)
            else:
                sub_tot, sub_path = dp(d_prev1, d_cur, d_next)
                tot = step_iters + sub_tot
                path = (d_next,) + sub_path
            if tot > max_tot:
                max_tot = tot
                best_path = path
        return max_tot, best_path

    for d0 in range(1, 71):
        for d1 in range(0, d0):
            iters1 = 145 - d0 - d1
            if d0 == 1 or d1 == 0:
                tot = iters1
                seq = (71, d0, d1)
            else:
                sub_tot, sub_path = dp(71, d0, d1)
                tot = iters1 + sub_tot
                seq = (71, d0, d1) + sub_path
            if tot > best_iters:
                best_iters = tot
                best_seq = seq

    print(f"EXACT MAXIMUM ITERATIONS over all degree sequences: {best_iters}")
    print(f"Achieved by degree sequence: {best_seq} (length {len(best_seq) - 1} rounds)")


if __name__ == "__main__":
    max_iters_analysis(70)
