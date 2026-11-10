# Experiment 2: new DQI decoder, chemistry refit, invsqrt tweaks

Local overall score **4.5267** (experiment 1: 2.2545). Upload `submission.zip`.
Per-task verifier outputs are listed in `results/summary.json`. Every number below comes from
the unchanged host verifiers at default sizes and shot counts.

| Task | Experiment 1 | Now | Improvement vs. baseline | Verifier output |
| --- | ---: | ---: | ---: | --- |
| invsqrt | 329,225 | 311,390 | 6.67× | `results/inv-flagdel.json` |
| dqi | 15.549 T | **3.268 T** | **7.48×** | `results/dqi-lean2.json` |
| dfthcsos | 2.903 B | **1.804 B** | **1.86×** | `results/chem35r.json` |

## DQI: fixed-schedule ("generic") Euclid decoder (`generic_dqi.py` → `toom_dqi.py` → `lean_dqi.py` → `lean2_dqi.py`)

The host decoder runs 222 synchronized Proos–Zalka iterations. Each pays for a full
multiply sweep, an array shift and a full array swap. That is 3.7M Toffolis.

**Key observation.** For random errors, every quotient in the Euclid run on (x^70, S(x)) has degree exactly 1.
With w errors the run stops exactly after w steps. So step i (1..35) can use **classically known
slot positions**:
- Two 71-slot arrays hold [remainder | cofactor]. A starts as (x^70, 0) and B as (S, 1); odd steps update A and even steps update B.
- The two remainder slots that step i zeroes become the cofactor's two new top slots
  (cofactor degree d sits at slot 70−d). No degree registers, unary iteration, shifts or swaps are needed.
- q = q1·x + q0 is computed from the remainders' top coefficients (one gated GF(4096) inverse).
  After the update it is erased using the cofactor relation s_i = s_{i−2} + q·s_{i−1}, so every step is reversible.
- A 6-bit counter of active steps gates each step. Its low bit selects which array holds (σ, ω),
  and one controlled swap moves them into B for the search.
- Shots fail only when a leading coefficient vanishes by accident (≈ w/4096 per shot): 359/102,400
  phase errors, 0 value errors. The repetition factor rises to 1.68 (from 1.19), which the much
  lower cost more than pays for.

Further steps, each fully verified:

| Version | Mean Toffolis | Qubits | DQI score |
| --- | ---: | ---: | ---: |
| Experiment 1 (Proos–Zalka + CRT/tower arithmetic) | 3,768,745 | 1,861 | 15.549 T |
| + cheaper Forney kickback (`forney_dqi.py`) | 3,715,522 | 1,861 | 15.329 T |
| Generic Euclid (`generic_dqi.py`) | 819,571 | 1,832 | 4.643 T |
| + evaluation/interpolation update (`toom_dqi.py`) | 616,710 | 1,832 | 3.490 T |
| + counter flags, in-place sums and evaluations (`lean_dqi.py`) | 627,702 | 1,780 | 3.333 T |
| + measurement-erased quotient, reordered cleanup (`lean2_dqi.py`) | ~621,900 | 1,766 | **3.268 T** |

- **Evaluation/interpolation update.** (q1·x + q0)·y with deg y = n−1 needs only n+1 field
  products: q1·y_{n−1}, plus q(a)·y(a) at n points. Evaluation, interpolation and multiplication by
  constants are CNOT networks over GF(4096). Each product is XORed into all of its target slots and
  deleted by X-measurement. That is 73 products per step instead of 142.
- **Per-position search** (89 Toffolis): GF(4096) inverse 66 (tower), MBUC parity into one qubit
  12 (linear map selected by the classical mask bits), and zero test with kickback 11 (one
  multi-controlled Z). The parity is erased by measurement plus a Toffoli-free `gf2_phase_by_product`.

**Next lever (not built): repair the rare failures** (repetition factor 1.68 → ~1.2, about 1.4× on DQI).
A vanishing leading coefficient at step i is detectable (slot T−i is zero but slot T−i−1 is not).
A degree-2 quotient there leaves the arrays in an "offset" layout. In that layout the only change
to later updates is where the single product q1·Y[T−j] lands, and the q-relations read positions
shifted by one slot. Both can be routed with a few controlled swaps of 12-bit registers per step.
Extra factors of x in σ are harmless, because the search never evaluates at 0.

## Chemistry: GPU refit with fewer vectors (`chem/`)

- `chem/fit.py` fits the DF-THC vectors and weights to the two-body tensor on the GPU. It minimizes
  error_2 as the 5776×5776 residual of Σ X⊗X − G, in fp32 with TF32 disabled.
  Using C = B (instead of the earlier C = 29 compression) restores full per-ρ capacity, because
  Σ_c X⊗X = Σ_{bb′} (wwᵀ)_{bb′} P_b⊗P_{b′}.
- `chem/finalize.py` encodes vectors as 12-bit Givens words (exact inverse of the decoder) and refits w with
  the quantized vectors fixed. It then optimizes W alone for lambda_eff (error_2 does not depend on W)
  and encodes the eigenvectors of h1′ as the one-body rows.
- `chem/rotate.py` applies per-ρ orthogonal mixing of the C columns of (w, W). This leaves error_2, h1′,
  E_SOS and the angle words exactly unchanged and lowers Λ (319 → 289 at R=55).
- `chem/build_chem.py` writes the nine-word angle circuits (experiment 1's `chemistry_qrom.py`).

| (R, B, C) | error_2 (mHa) | lambda_eff | Toffolis | Score |
| --- | ---: | ---: | ---: | ---: |
| Experiment 1 (65,60,29) | 19.40 | 48.10 | 81,670 | 2.903 B |
| (55,60,60) + rotation | 14.92 | 48.43 | 73,958 | 2.654 B |
| (50,60,60) + rotation | 15.02 | 48.35 | 67,844 | 2.431 B |
| (40,60,60) + rotation | 18.50 | 49.04 | 55,538 | 2.018 B |
| **(35,60,60) + rotation** | 19.60 | 49.29 | 49,404 | **1.804 B** |

Rejected along the way: (30,60,60) and (35,55,55) end at error_2 20.2 and 21.1 mHa after quantization. B ≤ 50 configurations plateau at error_2 ≥ 56 mHa, and angle precision
below 12 bits breaks error_2 (51.7 mHa at 11 bits, `chem_brot.py`).

## invsqrt (`polynomial_invsqrt.py`)

- **`cshift`:** zero-aware controlled shifts replace the cyclic rotations (−126 Toffolis).
- **`del_flag_if_zero`:** the normalization flags are erased by X-measurement plus a conditional
  phase-if-zero fixup. The verifier counts executed Toffolis, so this saves about 11 on average.

Score 329,225 → 311,390.

## Rebuild

From the repo root run `source initial_experiments/2/env.sh`. Then:
- **DQI:** `$EXP_PYTHON -c "import sys; sys.path.insert(0,'initial_experiments/2'); import lean2_dqi; lean2_dqi.build_circuit().to_file('dqi.kmb', format='kmb')"`.
- **Chemistry:** `$EXP_PYTHON initial_experiments/2/chem/build_chem.py initial_experiments/2/candidates/chem35r_full.npz <dir>`.
- **invsqrt:** `$EXP_PYTHON initial_experiments/2/build_submission.py --task invsqrt --output <dir>`.
- **Packaging:** `package.py` assembles and scores the zip from the verified candidates.
- **Rotation note:** `chem/wopt.py` (re-optimizing W after rotation) gained nothing (lambda_eff 49.290 → 49.275).
