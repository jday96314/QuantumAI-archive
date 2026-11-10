"""Record the follow-up probe and complete archive, with measured score labels."""
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
probe=json.loads((HERE/'results/subtree_probe.json').read_text())
assembly=json.loads((HERE/'results/package_assembly.json').read_text())
metric=json.loads((HERE/'results/package_metric.json').read_text())
audit=json.loads((HERE/'results/package_independent_audit.json').read_text())
assert metric['status']==audit['status']=='PASS'
current=next(r for r in probe['permutations'] if r['name']=='current')
optimized=next(r for r in probe['permutations'] if r['name']=='hamming_nn_2opt')
lines=['## Follow-up: exact subtree probe and complete package','',
    '**Probe outcome: zero projected Toffoli saving, including optimized address order.**',
    'Every one of the 795 full seven-word XOR-difference payloads is distinct',
    'in every oracle. Every one of the 315,615 distinct-input pair XORs is also',
    'distinct in every oracle. Thus a row permutation cannot create a constant',
    'subtree containing two valid payloads, nor repeated pair-difference leaves.',
    'The selected emitter already skips empty ranges and collapses constant',
    'subtrees; the final artificial padding row is already handled.','',
    'The exact selected-emitter count is 396 lookup + 49 low-bit selection +',
    '37 ANF erasure = **482 Toffolis per oracle per pass**. The fixed 38',
    'one-body pairs activate all 49 bank bits, so permuting dictionary rows',
    'does not reduce that component either. The counter equals the loaded',
    'accepted KMB operation statistics for all eleven oracles.','',
    '| Oracle | Current T/pass | Current round-trip T | Current-order saving | Optimized-order saving |',
    '|---|---:|---:|---:|---:|']
for k in range(11):
    r=next(r for r in probe['per_oracle'] if r['permutation']=='current' and r['k']==k)
    o=next(r for r in probe['per_oracle'] if r['permutation']=='hamming_nn_2opt' and r['k']==k)
    lines.append(f'| {k} | {r["current_forward_T"]} | {2*r["current_forward_T"]} | {r["forward_reverse_saving"]} | {o["forward_reverse_saving"]} |')
lines+=['| **Total** | **5302** | **10604** | **0** | **0** |','',
    'Savings columns include both forward and reverse passes. All fourteen',
    'tested common permutations give the same count: current order, eleven',
    'per-oracle lexical orders, a combined lexical order, and a common',
    'nearest-neighbour/2-opt order. One-body rows stay fixed; dictionary U,',
    'angle-word and weight rows travel together, with exact inverse-permutation',
    'array checks. The equality-collapse objective is identically zero, so the',
    'zero optimum is certified independently of the heuristic search.','',
    f'The secondary Hamming path objective improves from {current["hamming_path_cost"]:,} to {optimized["hamming_path_cost"]:,}',
    f'({100*(1-optimized["hamming_path_cost"]/current["hamming_path_cost"]):.4f}%); this does not save Toffolis in this emitter.',
    'Chemistry T remains **30,504**. The 3% build threshold is 915.12 T, hence',
    'at least **916 integer Toffolis**. Projected saving is **0 T / 0%**, so',
    '**no pruned family was built or native-verified**. This certificate is',
    'specific to exact whole-payload equality under permutations. It does not',
    'bound bitwise factoring, alternate gauges, linear transformations or',
    'other oracle architectures.','',
    '### Packaged accepted artifacts','',
    '- Chemistry: unchanged `candidates/banked`, canonical digest',
    f'  `{assembly["canonical_chemistry_digest"]}`.',
    '- Inverse: round-4 Codex `candidates/zero_head_affine.kmb`.',
    '- DQI: `initial_experiments/5/dqi-core/candidates/trimmed_decoder.kmb`.',
    '  The acceptance marker appeared while packaging was being prepared.',
    '  The first assembly guard stopped before writing an archive; the updated',
    '  default/independent/paired full native reports were then checked.',
    '  All have 102,400 shots. Paired runs have identical failure counts, R',
    '  and Q, with exactly 25,596 fewer mean Toffolis. The default accepted',
    '  report is the reference; the lower independent-salt cost was not selected.','',
    f'Archive: **[submission.zip](submission.zip)**, **{assembly["bytes"]:,} bytes**, exactly 14 members.',
    f'SHA-256: `{assembly["sha256"]}`.','',
    '**Full unmodified host-metric validation PASS.** `from_host/metric.py`',
    'read/validated every circuit and required NPZ array, ran all three native',
    'verifiers with full defaults and verification_shot_divisor=1, then computed',
    '`metric.score`. CPU affinity was restricted to two CPUs; native worker',
    'defaults therefore used at most two workers. Only the local input-directory',
    'configuration was set; no verifier, confidence or scoring function changed.','',
    '| Task | Verified shots | Measured cost in full metric run |',
    '|---|---:|---:|']
for task in ['invsqrt','dqi','dfthcsos']:
    lines.append(f'| {task} | {metric["verified_shots"][task]:,} | {metric["task_details"][task]["cost"]:.12g} |')
lines+=['',
    f'**Measured local combined score: {metric["local_combined_score"]:.13f}.**',
    f'Accepted-report reference combined score: **{assembly["accepted_reference_score"]:.13f}**.',
    'These are deliberately separate: DQI sampling failures and the host',
    'Clopper–Pearson repetition factor vary with the salt. The metric run used',
    'one predeclared fixed local salt, with no seed/hash/salt search:',
    f'`{metric["configuration"]["salt"]}`.','',
    'The archive was independently reopened after verification. Every member',
    'was read to EOF for CRC/content SHA checks; inverse/DQI circuit hashes,',
    'chemistry canonical digest, full shot counts, task costs and combined',
    'score all match the fresh native reports. The original accepted artifacts',
    'and all round-4 inputs remain untouched.','',
    'Packaging layout/hash binding adapts read-only round-4 Codex',
    '`package_submission.py` and `audit_submission.py`; the metric and native',
    'Apache-2.0 source stay unchanged under `from_host/`. DQI attribution and',
    'its accepted full-run evidence are retained in',
    '`results/dqi_core_selection_evidence.txt` and `package_trimmed_dqi_proofs.json`.','',
    '```bash',
    'source env.sh',
    '"$AUDIT_PY" subtree_probe.py > results/subtree_probe.log 2>&1',
    '"$AUDIT_PY" assemble_archive.py > results/package_assembly.log 2>&1',
    'OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 flock /tmp/r6-fullverify.lock "$AUDIT_PY" package_metric.py > results/package_metric.log 2>&1',
    '"$AUDIT_PY" package_audit.py > results/package_independent_audit.log 2>&1',
    '"$AUDIT_PY" record_package_results.py',
    '```','',
    'Evidence: [probe JSON](results/subtree_probe.json), [all permutation/oracle counts](results/subtree_probe.csv),',
    '[assembly](results/package_assembly.json), [fresh metric result](results/package_metric.json),',
    '[independent ZIP audit](results/package_independent_audit.json), [complete summary](submission.summary.json).','']
block='<!-- package-followup-start -->\n'+'\n'.join(lines)+'\n<!-- package-followup-end -->'
p=HERE/'RESULTS.md';text=p.read_text()
if '<!-- package-followup-start -->' in text:
    start=text.index('<!-- package-followup-start -->');end=text.index('<!-- package-followup-end -->')+len('<!-- package-followup-end -->')
    text=text[:start]+block+text[end:]
else:text+='\n\n'+block+'\n'
banner=f'**Complete archive:** [submission.zip](submission.zip). Full host-metric PASS; local combined score **{metric["local_combined_score"]:.10f}**. Exact subtree-pruning saving is **0 T**, so chemistry remains banked. DQI uses the newly accepted trimmed decoder.\n\n'
if banner not in text:
    text=text.replace('# Round 5 chemistry assumption audit\n\n','# Round 5 chemistry assumption audit\n\n'+banner,1)
p.write_text(text)
