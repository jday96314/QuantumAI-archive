# Coordinated quantized dictionary search

These are classical chemistry fits, excluded from the verified submission.

`qubo_anneal.cpp` and `lattice_qubo.py` solve small binary quadratic proposal
models using deterministic-seeded annealing. The solver matched exhaustive
optima on 50 random 12-bit instances. The paired-binary encoding of ternary
increments matched 30 exhaustive 8-coordinate cases. These checks concern the
proposal model; chemistry acceptance still requires its native verifier.

`anneal_coordinate.py` extends experiment 3's variable-projection coordinate
search with coordinated annealed proposals. Its ternary variant permits an
individual uphill move within a jointly improving group. Actual projected
physical loss, rather than the quadratic proposal, decides acceptance.

All following 700-vector fits use eight-bit words, C=350, ridge=1e-6, and the
same inherited starting snapshot. Results are saved in
`results/lattice_anneal_comparison.json`.

| Search | Saved error_2 (mHa) |
| --- | ---: |
| Matched greedy control, 120 steps | 20.328942 |
| Binary annealing, 200 steps | 20.302388 |
| Ternary annealing, 120 steps | 20.294246 |
| Wider ternary proposals, 160 steps | 20.253110 |

Every one misses the 20 mHa limit. The different step/proposal budgets prevent
attributing their final differences entirely to the solver; matched log points
are preserved for that reason.

`augment_lattice.py` adds residual-directed existing quantized atoms and
reprojects the weights to rank350. Starting with the wider fit, B=714 gives
20.003783mHa, B=715 gives19.981015mHa, and B=716 gives19.965766mHa. These are
two-body model errors, not full acceptance. Physical alignment toward the
checkpoint and a 240-step constrained LBFGS refit did not produce a cheaper
candidate. Its valid saved start has effective normalization72.1486, versus
48.8952 for the checkpoint; the final logged iterate still exceeds20mHa.

`block_anneal.py` applies residualized feature tangents and joint ternary
proposals to four 256-vector ensembles. Each actual trial reprojects the
weights before acceptance. In160steps the model improves from32.044843mHa to
approximately30.8mHa, still far outside the20mHa limit. Logs and snapshots are
retained. No circuits were exported from any fit in this document.

All reference data and inherited initializations were read from experiment3.
All outputs were written under the current experiment.

`results/rank_error_lower_bounds.json` also records the unrestricted positive
semidefinite rank-approximation floor from the reference pair-space spectrum.
At rank250 it is36.86mHa, rank27525.28mHa, rank30016.70mHa, and rank3507.55mHa.
Thus reducing a single ensemble to250 or275 columns cannot meet the20mHa bound,
even with an unrestricted dictionary. This bound concerns the two-body tensor;
it is not a complete lower bound on submission cost.
