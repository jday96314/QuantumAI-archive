You own `initial_experiments/5/dqi-census`. Goal: find decoder architectures that could reach **≥2.5× on DQI** (cost ≤~0.64e12). Four rounds and four recent agents found the compact-Euclid family locally optimal:
- The core has only 3.9% mechanical slack.
- The scan costs 88 T per point because the end state lacks the evaluator.
- The division-free scans as built need ≥1,383 Q.

The leaderboard top (18.9) implies a better design exists in at least one task.

**Read first:** RESULTS.md in `initial_experiments/5/{dqi-repair,dqi-phase,dqi-hybrid,dqi-core,dqi-dispatch}`, the round-4 report, Khattar et al. arXiv:2510.10967, and https://github.com/andigu/uncompute/blob/main/tasks/dqi-reed-solomon/circuits/discovered/README.md.

**Numbers:**
- Objective R·Q²·T. Incumbent Q 961, T 1,019,866 (core 649,802 + scan 370,152), R 1.697.
- Zero failures gives R 1.191. Q ≥ 840 syndrome qubits.
- For ≥2.5× at R=1.191: T ≤693k at 880 Q, ≤582k at 960 Q, ≤487k at 1,050 Q, ≤373k at 1,200 Q.
- Components: field multiply 42 T, inverse 66 T. Gu BM 227k forward+reverse at ~1,700–1,837 Q; BM storage floor 1,272 Q under its representation.
- Scan bodies: 88 (inverse-product), 80 (Forney, needs ω), 23 (needs ω and Ψ), 11 (needs F). Ψ lifetime 88.9k T at 1,383 Q.
- dqi-dispatch's classical local-transfer model handles every gap pattern.

**Phase 1, ideation (no building):** propose ≥6 genuinely distinct architectures (different algorithm, state representation or phase extraction, not schedule variants). For each give: the state through decoding, scan and reversal; a register table at the peak; a Toffoli estimate from the counts above; the failure model and R; the objective; and the single biggest uncertainty. Include at least:
- (a) an end state linear for the scan: replace A by F′=(A·B′)⁻¹ mod B in place, folding β^70 into the classical masks. It pays only below ~157k T per direction.
- (b) lazily allocated BM with coefficient windows for early rounds, handing off to compact storage.
- (c) a decoder holding (σ, ω) in syndrome storage throughout.
- (d) something not derived from Euclid or BM, e.g. structured Hankel/PGZ or a transform-domain method.
- (e) designs that spend width to buy large T cuts, priced honestly under Q².

**Phase 2:** rank by expected objective. For the top two, build the minimal native component that tests the biggest uncertainty, with exact counts, and update the estimate.
**Phase 3:** if a design projects ≥2× total on DQI, write a concrete build plan.

**Stopping:** never stop because a first draft loses to the incumbent. Drop a candidate only when its biggest uncertainty is resolved against it, and record why.
dqi-dispatch is changing the decoder concurrently: read its RESULTS.md for status, but don't build on its in-progress files. CPU cap 2.

## Shared constraints
- Write only inside your assigned directory. Reading other directories is allowed but read-only; treat other agents' in-progress files as unstable. `from_host/` and prior artifacts are read-only; never modify the verifier, confidence or scoring.
- Current local reference: combined 9.3677 (inverse 182,956.5272; DQI 1,598,204,229,385; chemistry 708,156,073.7876 banked in chem-audit). Preserve selections until a replacement passes.
- The public leaderboard top is 18.9 (≈8.2× our product of ratios). A failed approach is evidence about that design only.
- Objectives: DQI `R(epsilon_plus) * Q**2 * mean_Toffolis` (host Clopper–Pearson); chemistry `lambda_eff * T_total * Q_total`; inverse `Q * (mean_Toffolis + total_T/4)`.
- Accept only on full default native verification (DQI: all 102,400 shots plus an independent salt; chemistry: native ANGLE/oracle and physical checks; inverse: all 200,540,928 shots). No hash or salt grinding.
- The machine is shared (8 cores, 1 GPU): respect your CPU cap and wrap full verifier runs in `flock /tmp/r6-fullverify.lock`.
- Report a numeric first checkpoint within a few hours. Keep RESULTS.md current: hypotheses, changes, checks, measured versus estimated numbers, failures and reproduction commands. Attribute reused Apache-2.0 code.