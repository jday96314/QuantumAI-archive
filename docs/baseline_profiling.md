# Host baseline scores and profiling

Measured 2026-10-22 using the unchanged circuits in
`from_host/sample_submission` and the unchanged host verifiers. All three
baselines pass full verification at the default problem sizes and shot counts.

## Scores

These are **raw task costs; lower is better**, not leaderboard improvement
ratios. Normalizing these baselines against themselves gives an improvement
ratio of 1 for each task and a geometric mean of 1.

| Task | Raw score | Scored qubits | Mean Toffolis per shot | Verified shots |
| --- | ---: | ---: | ---: | ---: |
| `invsqrt` | 2,076,774.233715282 | 260 | 7,735 | 200,540,928 |
| `dqi` | 24,439,463,583,672.93 | 1,893 | 5,725,108.382109375 | 102,400 |
| `dfthcsos` | 3,355,599,588.3155375 | 730 | 93,237 | 254,464 |

Full JSON results, including circuit hashes, verification metrics, and resource
counts: [invsqrt](profiling/invsqrt/result.json),
[DQI](profiling/dqi/result.json), [chemistry](profiling/dfthcsos/result.json).
Seeds use the host's submission hashes, an empty salt, and no explicit override.

### Inverse square root: arithmetic dominates cost

The score is `Q * (Toffolis + total_T / 4)`. The 7,735 Toffolis plus
rotation synthesis give **7,987.593206597238 Toffoli equivalents**. There is one
explicit T gate and 37 non-Clifford rotations; their combined equivalent cost
is 252.5932, or 3.16% of the total. The allowed synthesis error per rotation
is `1.5289254939242122e-12`. The worst simulated error/tolerance ratio is
0.47691161, at `x = 2,588,411`.

| Circuit stage, including computation and uncomputation | Toffolis | Share of raw Toffolis |
| --- | ---: | ---: |
| Horner polynomial evaluation | 2,142 | 27.69% |
| Multiply normalized input by squared estimate | 2,086 | 26.97% |
| Multiply estimate by Newton correction | 1,734 | 22.42% |
| Square initial estimate | 884 | 11.43% |
| Denormalizing shifts | 600 | 7.76% |
| Normalizing shifts | 230 | 2.97% |
| Correction and interval selection | 59 | 0.76% |

The first four rows account for **88.51%** of Toffolis. Improving reversible
multiplication, truncation/precision choices, and polynomial evaluation has
more scope than optimizing the final phase gates alone. All 260 circuit
qubits count in the score; 22 hold the input and 238 are workspace.

Counts come from a seeded batch spanning the input domain, with each
annotated stage simulated separately while preserving state. The stage sum
matches the full verifier's Toffoli count exactly. See
[per-stage measurements](profiling/invsqrt-costs.json) and the host's
[reference circuit](../from_host/verifiers/invsqrt/nr_reference_circuit.py).

### DQI: Euclidean decoding dominates cost

The score is `R(eps_plus) * Q**2 * mean_Toffolis`. There were **zero value or
phase failures**, but the confidence bound is nonzero:
`eps_plus = 0.0008660587301657926` at `log2_delta = -128`.
This gives a repetition multiplier of **1.1912604648115832**.

| Stage | Mean Toffolis in attribution batch | Share |
| --- | ---: | ---: |
| Forward extended Euclidean algorithm (EEA) | 2,109,450 | 36.85% |
| Reverse EEA | 2,109,450 | 36.85% |
| Chien search / Forney phase application | 1,444,394.2266 | 25.23% |
| Initial alignment and its reverse | 61,824 | 1.08% |

The forward/reverse EEA pair contributes **73.69%**. Inside each EEA pass,
`multiply_and_add_or_sub` accounts for 1,613,222 static Toffolis (76.48% of
the pass), mostly its `sweep_positions` work (1,436,954). This makes the
finite-field multiplication/position sweep the largest arithmetic target.

The attribution batch contains 256 inputs drawn by the host sampler, with
input RNG and measurement seed 20261022. Its total is 5,725,118.2265625
Toffolis per shot; the full 102,400-shot score uses 5,725,108.382109375.
Conditional gates explain the small difference. The static upper bound is
5,725,461; it should not replace the executed mean in scoring.

The circuit uses 840 syndrome qubits and 1,053 scratch qubits. Because width
is squared in the score, a 10% width reduction would reduce cost by 19% if
the other terms stayed fixed. See [attribution data](profiling/dqi-costs.json)
and the host's [decoder](../from_host/verifiers/dqi/dqi_reference_circuit.py).

### Chemistry: angle lookups dominate Toffolis; state preparation sets width

The score is `lambda_eff * T_total * Q_total`, with
`lambda_eff = 49.30136925057439`, `T_total = 93,237`, and `Q_total = 730`.

| Component | Toffolis per block encoding | Share |
| --- | ---: | ---: |
| ANGLE table lookups, forward and reverse | 79,480 | 85.25% |
| Inner state preparation | 9,461 | 10.15% |
| Givens rotations | 3,450 | 3.70% |
| Outer state preparation | 846 | 0.91% |

There are ten ANGLE circuits, each with 3,974 Toffolis and 119 qubits.
Their two traversals account for the 79,480 lookup Toffolis. The reference
factorization is `R=65, B=60, C=60`, with 12-bit angles and batches of eight.

The scored width is **186 persistent + 544 shared temporary qubits**.
The 186 comprise 96 angle-word qubits and 90 state-preparation qubits.
Inner state preparation sets the 544-qubit workspace maximum; ANGLE needs
only 11 temporary qubits. Reducing ANGLE workspace alone therefore does not
reduce the scored width. Lookup optimization targets Toffolis; improving
inner state preparation or the factorization targets the width bottleneck.

The recovered words give one-body error **3.1852 mHa** and two-body error
**18.6396 mHa**, both below 20 mHa. `Lambda = 308.45593973309445` and
`E_gap = 3.965477244291833` satisfy the strict acceptance floors. The
two-body error has only about 1.36 mHa of headroom, so lowering angle
precision needs re-verification. The host's
[cost model](../from_host/verifiers/dfthcsos/verify.py) supplies the full
resource breakdown in the JSON result.

## Where evaluator runtime goes

Runtime is distinct from quantum circuit cost: cheap Clifford gates still
take simulator time, and input packing, parsing, and correctness checks do
not enter the quantum score.

Measurements use explicit `perf_counter` timers around simulator methods
and verifier stages, plus `cProfile` for Python internals. In this Python
3.13 / bundled native-extension combination, `cProfile` omitted substantial
native-boundary time; using its visible function totals as total evaluation
time would misidentify the bottleneck. The timing wrapper measures those
boundaries directly. Worker times are summed across tasks and are **not wall
time**; nested function timers are inclusive and must not be added together.

| Full verification run | Workers | Observed wall time |
| --- | ---: | ---: |
| Inverse square root, no profiling | 8 | 125.59 s |
| DQI, explicit timing | 4 | 86.31 s |
| Chemistry, explicit timing | 1 | 37.66 s |

These measurements include circuit/input loading and final score computation
inside the harness, but exclude interpreter startup and imports. The initial
inverse-square-root run overlapped other work; the final DQI and chemistry
timing runs were sequential. Do not add summed worker times to wall time.

### Inverse square root: scalar phase extraction is a major overhead

The [timing sample](profiling/invsqrt-timing-repeat/result.json) checks
262,144 shots (`x=1..4096`, 64 repeats), in 1,024 simulator batches of 256.
It took 1.335 s overall, with 1.152 s inside the 16 worker tasks:

| Worker operation | Seconds | Share of worker-task time |
| --- | ---: | ---: |
| `read_shot_phase`, including Python Fraction construction | 0.468 | 40.57% |
| Native `Simulator.do` | 0.395 | 34.26% |
| Reading input/scratch bits for validation | 0.062 | 5.37% |
| Packing/unpacking input bits | 0.027 | 2.31% |
| Remaining checks, conversion, transfers, and timing overhead | 0.202 | 17.49% |

The verifier calls `read_shot_phase` once per shot, constructs Python rational
numbers, and converts each to fixed-point units. Python profiling confirms
`Fraction.__new__`, `math.gcd`, and `half_turns_to_units` as the visible Python
hot spots. The full verifier makes **200.5 million** such phase reads.
A bulk fixed-point phase-read API could remove substantial Python/native
crossing and rational-conversion work. This is an optimization candidate,
not a measured speedup; the small-input sample need not have identical phase
conversion cost to the whole domain.

### DQI: native simulation dominates, followed by input packing

The [full timing run](profiling/dqi/result.json) contains all 400 batches of
256 shots. Summed worker-task elapsed time is **267.848 s** across four
processes:

| Worker operation | Summed seconds | Share of worker-task time |
| --- | ---: | ---: |
| Native `Simulator.do` | 229.553 | 85.70% |
| Boolean input packing (`to_bits`) | 18.692 | 6.98% |
| Writing simulator inputs | 13.759 | 5.14% |
| Sampling syndromes and expected phases | 5.313 | 1.98% |
| Other work | 0.531 | 0.20% |

Outside the workers, loading the 846 MiB circuit takes **7.874 s** and
canonical serialization/hashing takes **2.985 s**. The full run is 86.31 s
wall time; the remaining parent work includes static operation statistics,
reaction depth, process coordination, and final summaries.

The cost-attribution batch also shows why runtime and Toffoli cost differ:
Chien/Forney takes about **49.8% of the segmented native simulation time**,
despite only 25.2% of the Toffolis. It executes about 34.3 million CX gates
per shot, versus about 19.5 million across the two EEA passes. This single
batch is diagnostic rather than a full-distribution runtime estimate.

Native simulator/circuit improvements are the primary runtime opportunity.
Avoiding repeated large integer-array broadcasts when packing the 4,095
12-bit measurement registers, and caching unchanged circuit loading/hashes
in a local iteration loop, are secondary opportunities. Skipping checks or
reducing shots would change the verification workload and confidence.

### Chemistry: simulation first, then reference parsing and tensor metrics

The [full timing run](profiling/dfthcsos/result.json) uses 254,464 shots at
once and executes ten circuits forward and backward, for 20 native calls:

| Nonoverlapping work | Seconds | Share of wall time |
| --- | ---: | ---: |
| Native `Simulator.do` | 20.058 | 53.26% |
| Reference-integral parsing and expansion (`load_reference`) | 6.842 | 18.17% |
| Decoding words / tensor accuracy and normalization metrics | 5.580 | 14.82% |
| FCIDUMP checksum check | 1.120 | 2.97% |
| Other loading, simulator transfers, validation, and summaries | 4.060 | 10.78% |

The encompassing word-recovery stage is 23.752 s (63.07%); this **includes**
the native simulation above. Loading all inputs is 7.995 s (21.23%),
including the reference parse and checksum. The analytical quantum resource
model itself takes only **0.00061 s**.

Python profiling identifies about 4.285 million FCIDUMP line splits and
8.564 million pair-index calculations; tensor contractions in
`evaluate_dfthc_vectors` / `reconstruct_h2` dominate the mathematical stage.
For repeated evaluations, caching the validated reference tensors in a
long-lived process could avoid reparsing the 128 MB FCIDUMP. Preserve its
checksum validation and float32-rounding semantics. Simulator execution
remains the largest single stage; neither caching nor faster analytical
cost arithmetic eliminates that work.

## Reproduce

Use the environment setup in the [top-level README](../README.md). For the
measurements here, the venv is `/tmp/quantumai-venv.wOrZGC`, created from
`/mnt/data01/miniconda3/envs/torch2.13/bin/python`. Its C++ runtime is supplied
by `LD_LIBRARY_PATH=/mnt/data01/miniconda3/envs/torch2.13/lib`.
Versions: Python 3.13.15, NumPy 2.5.4, SciPy 1.18.1. Hardware: AMD Ryzen 7
5800X3D, eight physical cores / sixteen hardware threads, 62 GiB reported
RAM. `OPENBLAS_NUM_THREADS=1` was set to avoid nested BLAS threading.

After activating a compatible environment and setting `PYTHONPATH`:

```bash
# Choose fresh output directories; the harness refuses to mix runs.
python scripts/profile_baselines.py invsqrt --max-workers 8 --output /tmp/invsqrt-score
python scripts/profile_baselines.py dqi --timings --max-workers 4 --output /tmp/dqi-timing
python scripts/profile_baselines.py dfthcsos --timings --output /tmp/chemistry-timing

# Short inverse-square-root timing sample; NOT a full-domain score.
python scripts/profile_baselines.py invsqrt --timings --max-workers 1 --invsqrt-max-x 4096 --output /tmp/invsqrt-timing

# Python profiles include worker processes, not just the waiting parent.
python scripts/profile_baselines.py dqi --profile --max-workers 4 --output /tmp/dqi-cprofile
python scripts/profile_baselines.py dfthcsos --profile --output /tmp/chemistry-cprofile

# Gate attribution; checks outputs, phases, and operation counts against
# the original circuit. DQI temporarily exports a large text circuit.
python scripts/profile_circuit_costs.py invsqrt --output /tmp/invsqrt-costs.json
python scripts/profile_circuit_costs.py dqi --output /tmp/dqi-costs.json
```

The scoring harness calls the host verifier functions without changing their
acceptance rules, circuit files, or mathematical cost formulas. Profiling
instrumentation lives in `scripts/`, outside `from_host/`. Each result JSON
records versions, settings, wall time, and the full verifier response; timed
runs also aggregate per-worker timing counters. Raw `.prof` files and
individual worker timing files remain available locally but are gitignored.

Some initial runs overlapped on the same machine, and the data live on a
shared mounted filesystem. These are observed timings, not portable speed
guarantees or confidence intervals. Profiling introduces overhead; the
reduced inverse-square-root run samples only `x=1..4096` and is used for
runtime attribution, never as the reported baseline score. DQI and chemistry
timing results use the full verification workloads.
