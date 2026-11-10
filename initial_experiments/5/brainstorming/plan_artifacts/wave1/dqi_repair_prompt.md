# Round 5 / dqi-repair — diagnose and repair compact DQI exceptions

You own `initial_experiments/5/dqi-repair`. Work on the **round-4 Codex compact Euclidean DQI decoder** only. This is an architecture-level exception-repair experiment, not a general optimization campaign.

**Starting facts:** 961 qubits, approximately 1,019,866 mean Toffolis, 380 phase errors and zero value errors in 102,400 shots; repetition factor ~1.69685; DQI cost **1.598204229385e12**. A correct or near-correct decoder could move R toward 1.191–1.223, but the reason for the 380 errors has not been verified. The published Gu et al. DQI decoder uses a one-zero-pivot repair, but its algorithm is different: do not assume the failures have the same cause.

## Work plan

1. Locate the exact round-4 selected builder and full verifier output; reproduce its counts before editing. Read the round-4 consolidated report, the selected Codex report, and `from_host/verifiers/dqi/` for the syndrome/mask phase contract.

2. Construct a **fast, independent classical failure census** under the verifier's actual weight/position/value/mask sampling law. Target >=1,000,000 cases if practical. Compare the compact decoder's locator/evaluator/phase with an exact classical RS decoder; record the **first point of divergence** and causes (degree-gap/zero pivot, stopping/guard state, quotient update, locator root, Forney kernel, etc.). Separately compare a representative subset against native kickmix so the fast simulator cannot silently mis-model the circuit. A wrong decoded error does not always flip the phase for a random mask; distinguish event incidence from phase-error incidence.

3. Implement the **smallest reversible correction** for the dominant classified failure mode. Candidates include a one-exception realignment, an exceptional quotient update, or retaining enough position/history to unwind the exceptional branch. Exploit state reconstruction / already-consumed syndrome words where possible. Every quantum-controlled arithmetic operation costs gates even on rare inputs: do not price a fallback as though it runs only on failing shots.

4. Prove/test the correction on adversarial and random inputs including weight 35, zero syndromes, one-error syndromes, arbitrary incoming ancilla conditions where relevant, and independent random classical masks. Check complete syndrome restoration, clean scratch, exact phase, and inverse execution before native full scoring.

5. Report `R`, executed `T`, total `Q`, and end-to-end DQI cost for every retained candidate. Compare candidates using **matched explicit diagnostic seeds** first, then independent final native salts. Do not select on a lucky failure count.

## Success criteria and scope limit

- Strong result: expected <=8 phase failures/102,400 at approximately Q<=980 and T<=1.05M; cost gain ~1.3x or greater (measure rather than assume). Bigger width or gate overhead is acceptable **if the full cost still improves**.
- Stretch: use the corrected branch structure to eliminate enough duplicated update work to approach a >=1.5x DQI reduction.
- If the 1M-case failure census does not find a dominant, reproducible exception, report that result and one concrete alternative, rather than spending the remaining run on unrelated optimizations.

**External reference:** https://github.com/andigu/uncompute/tree/main/tasks/dqi-reed-solomon/circuits/discovered

## Shared constraints

- Write only inside your assigned experiment directory. **Do not read, search, or copy any other round-5 agent directory.** Prior rounds, shared baseline assets, and `from_host/` are read-only. Do not modify the verifier or official scoring rules. Discover actual repository paths rather than assuming links in summary reports are present verbatim.
- Start with the **round-4 Codex-selected**, not round-3 or another agent's inherited, artifacts. The local combined reference is **9.3677221531**: inverse cost **182,956.5272**, DQI **1,598,204,229,385**, chemistry **708,547,388.9553**. Preserve those candidates until a replacement passes.
- Use the **actual objective** rather than a gate-only proxy. For DQI: `R(epsilon_plus) * Q**2 * mean_Toffolis` using the host's Clopper–Pearson definition; for chemistry: `lambda_eff * T_total * Q_total`; for inverse: `Q * (mean_Toffolis + total_T/4)`.
- Maintain `RESULTS.md` with tested hypotheses, implementation changes, checks, exact/verifier-measured vs estimated numbers, failed candidates, and clear reproduction commands. Deliver any winning circuits/NPZ, builders, and native verifier JSONs. A faster component or a floating-point fit is **not** a verified task-cost win.
- For final DQI acceptance, use all **102,400** shots and independent validation salts; for chemistry, complete the native ANGLE/oracle and physical validation (not just the heuristic resource model). Do not exploit canonical-hash/sample variation, alter verifier confidence, or claim a win from shortened tests. Keep independent controls and attribute reused Apache-2.0 public code.