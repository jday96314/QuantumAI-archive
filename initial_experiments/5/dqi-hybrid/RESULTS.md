# Round 5 / dqi-hybrid: register-design stop, no replacement

**Outcome:** no complete hybrid below 1,200 live qubits was found in the
register-liveness/design pass. Work stops at the requested scope limit, before
rewriting the decoder. There is **no verified task-cost win** and no replacement
submission. The round-4 Codex-selected archive is preserved unchanged:
combined **9.367722153134604**, inverse **182,956.5271606425**, DQI
**1,598,204,229,384.9265**, chemistry **708,547,388.9553173**.

The deliverables are exact selected-builder replays, complete per-stage live
register timelines and native gate breakdowns, a working reversible overwrite
component, its native circuits and independent tests, and the obstruction/design
analysis below. The audit passes: [results/audit.json](results/audit.json).
All new files are inside this experiment. No other round-5 directory was read,
searched, or copied; prior rounds and `from_host/` were only read.

## Exact controls and gate attribution

The starting point is **round-4 Codex-selected**, not an inherited round-3
submission. Its actual builder is
`../../4/codex/build_cubic_compact_decoder.py`, which selects
`compact_product_decoder.build(70,3,22,4095, linear_padding=True,
virtual_zero=True,virtual_lead=True,low_inverse=True,lazy_flags=True)` and
`cubic_inverse_product_kernel.apply`. The selected source manifest and all
underlying pinned dependencies were rehashed successfully.

For comparison only, the exact round-3 artifact comes from `../../3/stream_capture.py`
with prefix `stream10_offset3_select_noscale_compact`, `build_stream.py`,
`stream_forney.phase_kernel`, and lifetime compaction. Both decoder-core tapes
were rebuilt from their actual selected sources. **All 15,032,464 round-3 rows
and 68,957,627 round-4 rows match the selected tapes byte for byte.** Bookkeeping
instrumentation changes no instructions. The round-4 replay stops immediately
after core comparison, before the builder writes a wrapper or circuit.

| Exact selected circuit | Forward core T | Reverse core T | Complete scan T, static | Complete T, static | Complete mean T, native | Q | R | Native cost |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Round 3 comparator | 159,664 | 133,839 | 352,252 | 645,755 | 645,755 | 1,429 | 1.2230139976 | 1,612,737,419,463.0198 |
| Round 4 starting control | 324,901 | 324,901 | 370,152 | 1,019,954 | 1,019,865.856484375 | 961 | 1.6968460733 | 1,598,204,229,384.9265 |

The scan columns include normalization/unpadding, evaluator preparation,
all 4,095 runtime-mask field points, and restoration; they are not just the
inner phase kernel. They were parsed by the native Kickmix parser from the
actual selected split KMX files. The core-plus-scan static sums exactly match
the selected full-circuit native statistics. No phase-kernel substitution was
made. The round-4 dynamic-mask mean is 88.143515625 below its static maximum;
the unchanged core contributes a fixed 649,802 T.

Evidence:

- [results/round3_liveness.json](results/round3_liveness.json),
  [results/round3_timeline.csv](results/round3_timeline.csv): all forward/reverse stages.
- [results/round4_liveness.json](results/round4_liveness.json),
  [results/round4_timeline.csv](results/round4_timeline.csv): all forward/reverse stages.
- [results/native_breakdown.json](results/native_breakdown.json): native scan statistics and exact sums.
- [results/scan_liveness.json](results/scan_liveness.json): allocation replay of the exact scan code.
- [results/control_round3_native.json](results/control_round3_native.json),
  [results/control_round4_native.json](results/control_round4_native.json): copies of the
  **inherited** full 102,400-shot host results, bound again to actual circuit bytes.
  These are not newly run verifier results. The counts are respectively
  0 value / 8 phase errors and 0 value / 380 phase errors.

The selected archive SHA-256 remains
`9168da25a9ac35f661885edad3e208db215bf67fa6525cca0ccc5d8831318b68`.
The selected DQI KMB SHA-256 remains
`4989da5cb633529e8851b9aa0c8c52abeb8990423a3a9463a70f673661a909c0`.

## Where storage actually lives

**Round 3 BM.** The 70 syndrome words always reserve 840 physical qubits.
`SIG` materializes 36 separate 12-bit locator slots. `BB` uses three private
words, then aliases syndrome words S0..S33, borrowing 408 input qubits. Its
fixed coefficient-window rotations are physical swaps. The final 35 input
words S35..S69 remain original syndrome data until the phase preparation
transforms them in place into the evaluator. There is no separate 420-bit
evaluator allocation. S34 is not included in that evaluator bank.

The late forward persistent reservation is exactly
`840 input + 432 locator + 36 private template + 19 stack + 52 fixed = 1379`.
The fixed 52 comprise discrepancy 12, degree counter 8, emitter pivots 3,
active/filled 2, stop index 9, parking index 6, and parked syndrome 12.
The stack has a 12-bit discrepancy and 7-bit tag. The journal is classical
builder metadata, not a bank of quantum history; its semantic information is
retained through the stack, flags, counter, and reversible polynomial state.
At round 68, 49 temporary qubits raise the core peak to **1,428**. The scan
raises the actual live peak to **1,429**, including its root, derivative,
evaluator accumulations, inverse tower, and measurement-phase repair scratch.

| BM forward round | Live peak | Persistent live at end |
| --- | ---: | ---: |
| 0 | 1,002 | 971 |
| 10 | 1,080 | 1,031 |
| 20 | 1,140 | 1,091 |
| 30 | 1,200 | 1,151 |
| 40 | 1,260 | 1,211 |
| 50 | 1,320 | 1,271 |
| 60 | 1,380 | 1,331 |
| 68 | 1,428 | 1,379 |

**Round 4 compact Euclid.** The same 840 input qubits carry *all* remainder,
cofactor, and inverse-scale words. With current cofactor degree k, the stored
lengths are `(k-1, k, 70-2k, 1)`, summing to 70 words; leading monic
coefficients are implicit. Completion flags are allocated lazily, reaching
35 bits. They preserve exceptional/completion decisions; no quotient-history
bank is retained outside the input array. The block-of-three transform holds
a temporary six-word head (72 bits). At k=19, input 840, completion flags 22,
and temporary/fixed workspace 99 produce the **961** peak. Serial updates
from k=22 need less scratch; their peak rises only to 914 as flags accumulate.
The final cofactor pair occupies 408 + 420 bits and the inverse scale occupies
12, still exactly the original 840.

| Euclid forward entry | Stored previous/current/remainder/scale words | Live peak |
| --- | --- | ---: |
| First transform | 0 / 1 / 68 / 1 | 874 |
| Block k=1 | 0 / 1 / 68 / 1 | 943 |
| Block k=10 | 9 / 10 / 50 / 1 | 952 |
| Block k=19 | 18 / 19 / 32 / 1 | 961 |
| Serial k=22 | 21 / 22 / 26 / 1 | 904 |
| Serial k=34 | 33 / 34 / 2 / 1 | 914 |

At the scan, all core-private scratch is clean and released; 840 input and
35 completion bits remain. The phase uses the two cofactors and derivative
through the inverse-product identity, avoiding a separate evaluator bank.
The existing selected scan reaches **959 actually live qubits**: 875 retained
plus 84 scratch. Its native wire frontier remains 961 because the wrapper
initially reserves the core frontier. The scan allocation experiment uses
two field points, replaying the unchanged source; point-dependent constants
and classical masks change no allocation/free sequence. Full gate totals use
the actual 4,095-point files, not extrapolation from two points.

The JSON timelines partition physical ownership; `input_array` remains 840
even when its logical role changes. The CSV separately records cofactor and
remainder lengths. Unclaimed temporary operand pivots are included in
`fixed_or_anonymous`; this label does not assert they are persistent state.
Event boundaries use allocator-event indices, avoiding ambiguity when a
release and a later allocation share the same gate index.

## Concrete overwrite component and proof

[consume_component.py](consume_component.py) emits native gates for

```
a = A[L]
b = sum(A[i] * S[i], i < L)
x = S[L]
(x, park=0) -> (a*x + b, [a=0]*x)
```

All field operations use the actual attributed Gu/shared-path arithmetic.
Let z=[a=0] and armed=a XOR z*1. The implementation:

1. Computes z reversibly from a; conditionally swaps x into the 12-bit park.
2. Temporarily replaces a=0 by 1, leaving every nonzero a unchanged.
3. Multiplies x in place by armed using Gu's multiplication and inverse-based
   measurement-phase recovery, then XORs the remaining discrepancy terms.
4. Restores a and erases the zero-test scratch.

The inverse subtracts the other terms (XOR), multiplies by armed inverse,
and reverses the controlled swap. Since armed is always nonzero, the inverse
exists on every tested operand, including the zero locator polynomial.
The same gate map is a permutation even for arbitrary initial park contents;
the displayed discrepancy formula assumes an initially clean park. The
component tests use that clean-park contract. A phase sandwich applies the
parity of the discrepancy under a **dynamic 12-bit classical mask**, then
restores all syndrome and locator operands and clears park and scratch.

| L | Public quantum bits, including park | Peak Q | Forward T | Inverse T | Phase-sandwich T |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1 | 60 | 104 | 209 | 209 | 418 |
| 17 | 444 | 488 | 881 | 881 | 1,762 |
| 34 | 852 | 896 | 1,595 | 1,595 | 3,190 |

These are **native component counts**, not full-decoder resources. Forward T
is `42*L + 167`. Private scratch is 44 bits beyond the public operands and
park. The 896-qubit L=34 component omits the rest of the full decoder's
syndrome/template/control state and cannot be priced as a competition entry.
Its park must stay live until the inverse; it cannot be counted as free space.

Native outputs are in `components/consume_l{1,17,34}_{forward,inverse,both,phase}.kmb`
and matching KMX/NPY, with builders and hashes in the result JSONs. Tests cover:

- **1,000,000 classical cases** with exact host-law syndrome marginals, uniform
  distinct error positions, nonzero values, and the host's optimal weight law.
  Locator operands are independently varied across every degree 0..34,
  including forced zero leading coefficients. Zero reconstruction failures.
  This is a representation test, **not one million complete BM decodes**.
- **102,400 native shots per emitted circuit**, checking the independent field
  reference, inverse, composition, scratch, syndrome restoration, locator
  preservation, phase, and classical mask preservation. L=34 has two
  independent salts; L=1 and L=17 have one each: **1,638,400 component shots**
  total, zero value/scratch/phase failures.
- All 4,096 pivot values, all 4,096 input values on the zero-pivot branch,
  each exceptional locator degree, and an all-zero locator control are
  explicitly included in every native suite.

See [results/law_law-independent-a.json](results/law_law-independent-a.json),
[results/consume_l34_independent-a.json](results/consume_l34_independent-a.json),
[results/consume_l34_independent-b.json](results/consume_l34_independent-b.json),
[results/consume_l17_independent-a.json](results/consume_l17_independent-a.json),
and [results/consume_l1_independent-a.json](results/consume_l1_independent-a.json).
These JSONs are generated by this experiment's native component checks;
none is presented as an official full DQI verifier result.

## Tested hypotheses, failures, and stop decision

**Direct syndrome recycling into the BM template:** already implemented by the
selected round-3 code, and reproduced exactly here. It leaves 36 locator words
outside the input array. Under this representation, even deleting every other
ancilla leaves `840 + 432 = 1272` reserved qubits. Making a locator coefficient
implicit saves only 12, still 1260. This is a storage-policy bound, not an
information-theoretic bound for arbitrary encodings.

**Overwrite a syndrome with its discrepancy without parking:** rejected by a
strict injectivity counterexample. Hold all other operands fixed, choose a=0,
b=7. Inputs x=1 and x=2 both become discrepancy 7. There are 4,096 preimages,
so at least 12 retained bits are needed elsewhere on this fiber. The emitted
component implements precisely that retention. A replacement by an armed
nonzero multiplier alone computes x+b on this branch and therefore does not
compute the actual discrepancy. Parking is a reversible solution, not a width
saving when counted honestly.

**Clear more original syndromes using the locator recurrence:** after processing
through zero-based round r with degree L, the established equations allow
recovering S[t] from later processed words only when `t+L <= r`, and only with
a nonzero relevant leading coefficient. On the generic late path r=68, L=35,
this reaches t<=33: exactly the 34-word prefix already occupied by the template.
The next candidate t=34 requires a relation at index 69, not yet established.
At the final round one additional word may become eligible; that cannot remove
the earlier width peak. The recurrence concerns original syndromes, not the
template data already overwriting those slots. Erasing a template coefficient
requires another recovery identity. The exceptional zero-pivot case also needs
the retained information proved above.

**Compact cofactors plus a remainder:** credible and already implemented at
Q=961 in the selected round-4 control, but its exact forward/reverse cost is
649,802 T. Copying its state layout does not preserve BM's 293,503 T core.
The missing operation is a cheap reversible coordinate change/reconstruction
between BM's original-syndrome/template state and the compressed remainder/
cofactor state. The existing dense triangular transforms pay for that change.
No new implementation with lower native cost was obtained here.

**Early BM-to-compact handoff or quotient-only checkpointing:** assessed as
unresolved implementation/algebraic work, not emitted candidates. The selected
BM first reaches 1200 at round 30, so a handoff must occur earlier. It must
encode the processed prefix, preserve the original syndrome inverse, recover
or erase the template and branch history, and produce the compact remainder
without transiently retaining both full representations. Replaying BM to
recover an omitted template, or retaining a quotient list while materializing
both cofactors, has no demonstrated schedule meeting the width and gate
budgets. No gate-only estimate from those proposals is claimed as evidence.

There is **no proven global information-theoretic obstruction** to the target.
For example, locator plus evaluator determine a processed syndrome prefix by
formal-series division when the locator constant is nonzero; the template is
then a deterministic function of that prefix and the BM convention. Thus
there is real redundancy that might be exploitable. Round 4 already proves
low-width reversible state storage is possible. What remains unproved and
unimplemented is inexpensive reversible recovery compatible with the shared
arithmetic path and exceptional branches. The local zero-pivot obstruction
must not be mistaken for a universal lower bound.

No credible **new shared-path hybrid** below approximately 1200 qubits emerged.
The requested stopping condition therefore applies. There was no new complete
DQI circuit to submit for the final host default-shot/salted acceptance step,
and no new task score was inferred from component checks. The existing phase
scan was kept for all comparisons; cheaper phase integration was not pursued.

## Objective and preservation

Every full-circuit comparison above uses the unchanged host objective
`R(ClopperPearsonUpper(errors,102400,-128)) * Q**2 * mean_Toffolis`.
For context, the target Q=1050, T=600000, R=1.23 costs 813,645,000,000,
about 1.964 times below the selected control. This is a target, not an entry.

Even granting a hypothetical Gu circuit 442,000 host Toffolis at Q=1837,
the host's **zero-error** 102,400-shot confidence bound implies minimum
R=1.1912604648 and cost approximately **1,776,835,860,881.6**. This optimistic
pricing example already exceeds the selected control; it is not a verified
conversion of Gu's published pipeline to the host convention.
No inverse, chemistry, scoring rule, verifier confidence, or canonical sample
selection was changed. No winning circuits or NPZ are claimed. The existing
complete selected circuits and chemistry NPZ remain in the protected round-4
archive, bound by its digest and retained native evidence.

## Reproduction

From the repository root:

```bash
source initial_experiments/5/dqi-hybrid/env.sh
"$EXP_PYTHON" initial_experiments/5/dqi-hybrid/profile_selected.py round3
"$EXP_PYTHON" initial_experiments/5/dqi-hybrid/profile_selected.py round4
"$EXP_PYTHON" initial_experiments/5/dqi-hybrid/profile_selected.py scan
"$EXP_PYTHON" initial_experiments/5/dqi-hybrid/profile_scan_storage.py
"$EXP_PYTHON" initial_experiments/5/dqi-hybrid/consume_component.py --l 34 --shots 102400 --salt independent-a
"$EXP_PYTHON" initial_experiments/5/dqi-hybrid/consume_component.py --l 34 --shots 102400 --salt independent-b
"$EXP_PYTHON" initial_experiments/5/dqi-hybrid/consume_component.py --l 17 --shots 102400 --salt independent-a
"$EXP_PYTHON" initial_experiments/5/dqi-hybrid/consume_component.py --l 1 --shots 102400 --salt independent-a
"$EXP_PYTHON" initial_experiments/5/dqi-hybrid/law_check.py --shots 1000000 --salt law-independent-a
"$EXP_PYTHON" initial_experiments/5/dqi-hybrid/audit_results.py
```

The round-4 core replay takes a few minutes and transiently materializes its
2.4-GB tape in memory for exact comparison, without writing another copy.
The component exporter invokes the read-only round-4 `translate_prefix`
executable and writes all intermediate/native artifacts inside this directory.
`env.sh` disables Python bytecode writes to inherited source directories.

Arithmetic and the upstream BM construction are adapted from
[andigu/uncompute](https://github.com/andigu/uncompute/tree/main/tasks/dqi-reed-solomon),
revision `64b5a13aa249f1360846657392e6c66817b26d0f`, under Apache-2.0.
The complete vendored source and license remain at `../../3/vendor/uncompute`.
This experiment imports them read-only and does not redistribute an altered
vendor tree. The host field conversion, inverse gadgets, and selected compact
Euclid code are also imported from the attributed prior-round implementations.

The host's dynamic classical masks select parity operations and permit cheap
classical field maps; they do not make the quantum locator degree a classical
allocation decision. Omitting high locator slots on a typical-degree argument
would change the implemented map on substantial input mass (the exact host law
has P(weight=35)=0.0276659241). Physical Q is the circuit's wire count, not a
per-shot mean. No such truncation is introduced here.

[results/artifact_manifest.json](results/artifact_manifest.json) binds every
local artifact except transient logs and the manifest itself. Re-running
`audit_results.py` regenerates the CSVs, checks, and manifest.
