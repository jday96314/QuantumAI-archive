# Follow-up: dispatcher v2, width-first

Your local-transfer model is the right algorithm (all 102,400 host shots correct). The matrix/Bennett realization failed mainly on **width** (≥1,107 Q). My +90k Toffoli kill threshold was too tight. With ~0 failures, R ≈ 1.191, so at 961 Q the decoder breaks even up to +433k mean Toffolis, keeps 1.3× up to +98k (+124k with dqi-core's trims), and keeps 1.2× up to +191k. Width binds: each 10 qubits ≈ 2.1%.

1. **Base:** adopt dqi-core's trimmed middle products and budgeted inverse retention (`initial_experiments/5/dqi-core/components/trimmed_middle.py`, `budget_inverse.py`; −25,596 T at 961 Q) plus your completion counter. Use dqi-core's full trimmed builder only once its RESULTS.md marks it ACCEPTED.
2. **Design:** keep the incumbent's fixed boundaries (block3 at k=1..19, serial k=22..34) and regular heads.
   - Handle (2,1), (1,2), (3), deferrals, immediate completion and the final two slots by **routing inside the shared cofactor passes**. Passes accept multipliers of degree ≤2; quotient coefficients come from your moment recurrence; a deferred quotient is a carried pending-degree state plus a controlled one-word bank realignment.
   - Never materialize explicit 2×2 matrix entries or a Bennett bridge.
3. **Fixed-schedule classical model first:** convert `local_transfer.cpp` to the exact fixed schedule and carried bits you will emit. Require ≤10 failures on the 102,400 host shots and ≤~20 per million. Assert every intermediate state against your variable model and EEA.
4. **Width plan on paper before emission:** a live-qubit table per stage, target peak ≤965 (stretch ≤950).
5. Emit components, then the full decoder. Run controls (zero, one-error, weight-35, ordinary, independent masks), then the full 102,400-shot run plus an independent salt under the flock.

**Success:** ≤10 failures, peak ≤975, mean T ≤1.145M (≥1.23×). Stretch: T ≤1.07M at ≤950 Q (≈1.4×).
**Kill:** peak >990 Q, projected T >+300k, or the fixed-schedule model exceeds 20 failures on the host shots. Before stopping, report the best width/T point reached and what blocks it.
Read anything (read-only); write only here. CPU cap 3.