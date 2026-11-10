"""Update final report from audited native numbers, never from estimates."""
from census_common import *
a=json.loads((HERE/'results/audit.json').read_text());s=json.loads((HERE/'selection.json').read_text())
assert a['status']=='PASS'
rows=[]
for r in a['research_evaluations']:
    mt=r['metrics'];rc=r['resources']
    rows.append(f"| {r['candidate']} | {r['status']} / {r['verified_shots']} | {1000*mt['error_1']:.6f} | {1000*mt['error_2']:.6f} | {mt['lambda_eff']:.6f} | {rc['toffoli']:.0f} | {rc['qubits']} | {r['conditional_objective']/1e6:.6f} |")
mt=s['metrics'];rc=s['resources']
rows.append(f"| structured banked | PASS twice / {s['shots_per_run']} each | {1000*mt['error_1']:.6f} | {1000*mt['error_2']:.6f} | {mt['lambda_eff']:.6f} | {rc['toffoli']:.0f} | {rc['qubits']} | {s['score']/1e6:.6f} |")
bounds=[]
for b in a['independent_bounds']:
    bounds.append(f"| {b['file']} | {1000*b['span_bound']:.6f} | {1000*b['rank_psd_bound']:.6f} | {1000*b['saved_error2']:.6f} | {'yes' if b['coefficient_rescue_excluded'] else 'no'} |")
extended=json.loads((HERE/'results/repair_B703_extended.json').read_text())
growth=json.loads((HERE/'results/growth511_sequential.json').read_text())
last=[r for r in growth['rows'] if r['B']==511]
growth_text='; '.join(f"{r['mode']} {1000*r['refit']['error2']:.6f} mHa (span bound {1000*r['span_bound']:.6f} mHa)" for r in last)
global_rank=a['global_rank'];grc=global_rank['conditional_R1_B719_minimum_C']
final=f'''
<!-- audited-final-start -->
## Final measured result

**Verified improvement: {s['score']:,.4f}, {s['improvement']:.9f}x over the
banked 708,156,073.7876. The 1.25x target ({TARGET:,.4f}) was not reached.**
Local selected candidate: `candidates/structured`; prior external selections
remain untouched. `selection.json` binds the exact candidate and both native
records. Native submission digest: `{s['submission_sha256']}`.

| Candidate | Native status / shots | error1 mHa | error2 mHa | lambda_eff | T | Q | Product M |
|---|---|---:|---:|---:|---:|---:|---:|
{chr(10).join(rows)}

Products for rejected rows are conditional arithmetic, not valid scores.
All rows use complete unchanged native ANGLE recovery, default 64 shots per
address, fixed physical shifts, native oracle semantics and every physical
gate. Full runs hold `/tmp/r6-fullverify.lock`. The selected candidate passes
two predeclared independent salts with identical digest and score. No salt
or hash search occurs.

## Top-route uncertainties and stopping decision

**Integer continuation:** the 240-round unrestricted repair reaches
{1000*extended['final']:.6f} mHa at B703, demonstrating that changing integer
words can cross the first 16-row boundary. It does not establish the
low-normalization path to B511. The separate physical-generator transfer,
6,000-step constrained weight/offset fit, finite offset optimization and native
physical test quantify the normalization cost. A failed early 60-round
endpoint was extended rather than discarded. This family is not ruled out
globally; it remains an open route to deeper thresholds.

**Integer growth:** both residual-mode insertion and incumbent reinsertion
use exact integer atoms before all fitting, matched 8-row batches and rank-C
solves. A second 479->511 test repairs the poor one-shot seed by sequential
deletion first. Its endpoints are {growth_text}. At 639->671, residual growth
is 34.984890 mHa versus 24.534255 mHa for reinsertion. Independent span bounds
exclude all real coefficient-only rescue of these particular final growth
dictionaries. This resolves the finite candidate pool's capacity uncertainty
against it. Different seeds, refreshed pools, word moves and broader integer
growth remain unresolved, not globally excluded.

Independent CPU SciPy QR/eigensolves use exact native-decoded vectors and
the fresh signed target. Bounds have a 1e-8 Ha exclusion slack; these are
floating-point certificates, not interval-arithmetic proofs.

| Frozen research dictionary | Span bound mHa | PSD rank-C bound mHa | Saved residual mHa | Coefficient rescue excluded |
|---|---:|---:|---:|---|
{chr(10).join(bounds)}

**Dictionary-independent rank audit:** every native two-body model is a PSD
matrix of rank at most R*C in the orthonormal symmetric pair coordinates.
The full fresh native spectrum gives necessary **R*C >= {global_rank['minimum_RC_for_20mHa']}**
at 20 mHa, independent of dictionary, integer grid or factor gauge. This is
the PSD/rank approximation theorem evaluated in float64, not an interval
proof. It excludes standalone C256/C224 strategies even with different words
at R1. At B719 and this minimum C, the existing J7 compiler still gives
conditional **{grc['cost_at_banked_lambda']/1e6:.6f}M** at banked normalization.
`results/global_rank.json` retains the adjacent threshold margins and full
spectral-tail numbers. This does not exclude combined B/rank/oracle changes.

**Other census routes:** the unchanged incumbent dictionary cannot reduce C
to 320 (21.363025 mHa optimum), far short of the C reduction needed for the
target. Uniform one-body precision reduction fails at 123.108293 mHa;
one-body and two-body bounds are independent, so unused one-body tolerance
cannot be reassigned to error2. The exact host map shows a split width cliff
at B513->512->511. The structured ANGLE fallback preserves the physical
factorization, increases total T but reduces persistent word storage by 24Q,
and improves the actual product by {(1-1/s['improvement'])*100:.6f}%.
The compiler architecture is reused from chem-qrom and freshly fitted to the
banked words; this is not evidence for a new low-error physical basin.

These are finite, scoped candidate decisions. The census stops after the
declared continuation extension, constrained physical refit, matched growth
controls, fixed-dictionary capacity audits and complete native tests. It does
not claim that any shape, all integer growth, all ANGLE designs, or chemistry
as a whole is infeasible below the target. No experiment job is left running.

## Reproduction and artifact audit

From this directory:

```bash
source env.sh
"$CENSUS_PY" probe.py > results/probe.log 2>&1
bash run_searches.sh
bash run_followup.sh
bash run_physical.sh
"$CENSUS_PY" final_audit.py > results/audit.log
"$CENSUS_PY" report.py
```

`results/structured_components.json` checks all valid addresses, arbitrary
incoming words, both directions, clean scratch and phase under two independent
seeds. `results/audit.json` verifies all native digests/default workloads,
exact physical permutation, preserved banked NPZ and local file hashes.
`results/source_manifest.json` attributes and hashes read-only reused source;
`artifact_manifest.json` binds every numerical checkpoint and delivered
candidate. `NOTICE.md` attributes native Apache-2.0 and prior research code.
Fresh stochastic trajectories can differ numerically; delivered exact words,
weights, circuits and native digests bind the measured results.
<!-- audited-final-end -->
'''
path=HERE/'RESULTS.md';text=path.read_text().split('<!-- audited-final-start -->')[0].rstrip()
text=text.replace('Status: Phase 1 complete; Phase 2 numerical prototypes complete, extended\ncontinuation and full native checks in progress. No selection is replaced.',
    f'Status: complete. Local verified compiler improvement selected: {s["score"]:,.4f};\n1.25x target not reached. External selections remain untouched.')
path.write_text(text+'\n\n'+final)
print(json.dumps(dict(score=s['score'],improvement=s['improvement'],target_achieved=s['target_achieved'])))
