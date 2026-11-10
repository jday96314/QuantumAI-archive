"""Finalize the fixed-schedule continuation's evidence and honest disposition."""
import hashlib,json
from pathlib import Path

ROOT=Path(__file__).resolve().parent
def read(name):return json.loads((ROOT/'results'/name).read_text())
cases=[('Exact host',102400,2026102405,'host'),('Law stream A',1000000,2026102405,'law'),('Law stream B',1000000,2026102406,'law'),('All-weight controls',4608,2026102423,'controls')]
data=[read(f'fixed_{mode}_{n}_{seed}_g3.json') for _,n,seed,mode in cases]
for r in data:assert r['fixed_schedule'] and r['wrong_decodes']==0 and r['unsupported']==0
names=['routed_cofactor_d19_g1.json','routed_cofactor_d19_g2.json','routed_cofactor_d19_g3.json','shared_slot_n0.json','shared_slot_n18.json','local_schur.json','moment_probe_k19_n8.json','moment_probe_k19_n9.json','routing_merge_r19.json','completion_counter_b3.json','completion_counter_b2.json']
components={n:read(n) for n in names}
for name,r in components.items():
    assert r['status']=='PASS',name
    assert all(x['status']=='PASS' for x in r['circuits']),name
screen=read('fixed_route_resource_screen.json');assert screen['status']=='KILL_SEPARATE_CHART_HEADS' and screen['partial_projected_added_T']>300000
assert screen.get('base_accepted'), 'Run adopt_accepted_base.py and wait for profile_fixed_route.py to finish before auditing.'
embedded=read('embedded_components.json');assert embedded['max_Q']==973
base=read('core_components_adoption.json');assert base['forward_T']==base['reverse_T']==312103 and base['peak_qubits']==961
baseline=ROOT/'../../4/codex/candidates/cubic_trace_product_n70_m4095.kmb'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
baseline_hash=sha(baseline);assert baseline_hash=='4989da5cb633529e8851b9aa0c8c52abeb8990423a3a9463a70f673661a909c0'
adopted=json.loads((ROOT/'parameters/base_adoption.json').read_text());assert adopted['status']=='ADOPTED_READ_ONLY_ACCEPTED_BASE'
assert sha(Path(adopted['circuit']))==adopted['sha256']
checks=sum(r['fixed_statistics']['intermediate_checks'] for r in data)
eea=sum(r['full_eea_state_audits'] for r in data)
table='\n'.join(f'| {label} | {n:,} | {seed} | {r["wrong_decodes"]} | {r["independent_mask_phase_errors"]} | {r["fixed_statistics"]["intermediate_checks"]:,} |' for (label,n,seed,_),r in zip(cases,data))
native_table='\n'.join(f'| {label} | {components[name]["circuits"][0]["Q"]} | {components[name]["circuits"][0]["T"]:,} |' for label,name in [('Shared cofactor, d19/g1','routed_cofactor_d19_g1.json'),('Shared cofactor, d19/g2','routed_cofactor_d19_g2.json'),('Shared cofactor, d19/g3','routed_cofactor_d19_g3.json'),('Carry-finishing field slot, 18 words','shared_slot_n18.json'),('Six-word routed Schur head','local_schur.json'),('Conditional 8-moment head','moment_probe_k19_n8.json'),('Conditional 9-moment head','moment_probe_k19_n9.json'),('Merged prefix realignment, k19','routing_merge_r19.json')])
embed_table='\n'.join(f'| {r["name"]} | {r["Q"]} | {r["T"]:,} |' for r in embedded['rows'])
text=f'''# Fixed-schedule compact DQI dispatch

**Disposition: fixed-schedule classical model passes; native shared components
pass; the separate carry-chart moment realization crosses the revised +300k T
kill limit. No exception-correct replacement emitted or accepted here.**

The best compiled arithmetic embedding reaches **973 Q**. The incomplete
schedule's resource screen is **planned Q=973, partial projected mean T=
{screen['partial_projected_mean_T']:,.3f}**, including dqi-core's 25,596-T saving.
The width is not an integrated-decoder measurement, and the T projection still
omits required work. Its added T is **{screen['partial_projected_added_T']:,}**,
already above the revised kill limit. The blocker is separate moment arithmetic
for the carry charts, followed by canonical bank routing. A common moment map
that shares those products remains open; this result does not rule it out.

The **accepted dqi-core trimmed decoder is now adopted read-only as the base**:
Q=961, static T=994,358, default mean T={adopted['default_mean_T']:,.6f},
382 phase / 0 value failures in its native default run, cost=
{adopted['default_cost']:,.3f}. Its SHA256 is `{adopted['sha256']}`.
The acceptance includes a full paired-seed comparison (exactly 25,596 fewer mean
T with identical decoding), a full default run, and a full independent salt.
`parameters/base_adoption.json` records the verified artifact and builder.
This experiment did not rerun or modify them.

The protected round-4 incumbent also remains intact: Q=961, mean T=1,019,865.856484375,
380 phase / 0 value failures, R=1.6968460732519504, DQI cost=
1,598,204,229,384.9265. Its SHA256 was rechecked:
`{baseline_hash}`. This experiment reports no exception-correct replacement score
or new combined score. Prior assets, verifiers, confidence calculations, and
scoring are read-only.

## Exact fixed-schedule model

`local_transfer.cpp ... fixed` uses the original initial slot k=0, block3 at
k=1,4,7,10,13,16,19, and individual serial slots k=22..34. It does not select new
block boundaries. Updates use direct quotient/cofactor recurrences, not matrix
entries. The previous variable model is used only as an assertion oracle.

| Input | Cases | Seed | Wrong decodes | Independent-mask phase failures | Intermediate checks |
| --- | ---: | ---: | ---: | ---: | ---: |
{table}

The controls cover every weight 0..35, 128 cases each. The host data is the exact
original 102,400-shot file independently checked by dqi-repair. Every decoded
locator matches Berlekamp–Massey and every true error magnitude matches, excluding
extra roots as well. Thus every original host mask has the correct phase; no
particular mask or salt is selected. Both independent law streams are reported.

There are **{checks:,} intermediate checks**, including **{eea:,} full-EEA
quotient comparisons**. Every internal degree-slot cut and every outer boundary
packs the public state, reconstructs it solely from that state and the declared
metadata, and checks against the current EEA state. Completed quotients also
match the previous variable-transfer model. Neither oracle chooses a branch.

The host run has 403 deferred outer boundaries and 746 missing-degree units.
The fixed model repairs the serial degree-three failure of the variable model.
It also covers (2,1), (1,2), (3), multiple gaps, immediate completion and the final
two serial slots. Logical prefixes reach ten moments, but at most six fresh
words are read at a block; the rest belong to the carried prefix.

## Exact carried-state contract

Persistent metadata is six completion-count bits, two six-bit sparse gap-stage
tags, and three mode bits. A cubic occupies two gap tags. More than two missing
degree units is an explicit unsupported extension, counted as failure; none
occurred in the reported samples. This is not a universal no-failure claim.
Mode is also independently derived from gap tags and completion count in every
chart assertion. No hidden degree counter is used.

At nominal k the cofactor prefix has 2k-1 field words, the retained tail has
70-2k words, and pinv has its fixed word. Zero ghost words pad the prefix. The
decoder retains 70 field words throughout, using the incumbent's physical frame
handles; the tail is not shifted when the cofactor chart changes.

A quadratic's first slot scales Bprev by beta and pinv by inverse(beta), then
defers both quotient MACs to its second slot. The omitted high remainder is
(0,r), with r=inverse(pinv). A cubic's first slot does nothing; its second
prepares A=beta*Bprev+q2*x^2*B; its last adds q1*x*B+q0*B and the free x^3
monomial shift. The two-slot cubic omits (0,0,r,s), where
s=(q2+B_sublead)*r. Each arithmetic multiplier has degree at most two.
Completion pads both cofactors by x and increments the binary count.

## Base adoption and width plan

Initially only `../dqi-core/components/trimmed_middle.py` and `budget_inverse.py`
were adopted through `core_components.py`. Replaying them reproduces **312,103 T
per direction, 624,206 core T, 961 Q**, a **25,596-T** saving. During final checks
dqi-core's full report changed to ACCEPTED. Its verified full builder and KMB are
now the read-only base reference. No unaccepted full builder was used. New conditional contexts use the low
inverse rather than assuming the original retention-body allowance is sufficient.

[WIDTH_PLAN.md](WIDTH_PLAN.md) was written before emission. Its target allocation
is 840 state + 22 metadata/control + 72 copied-head + 39 arithmetic bits =
**973 Q**. The scan estimate is <=952. A cubic pending prefix must reuse a public
carry coefficient and the one cache; overlapping another cache raises the plan
to 985. That complete prefix handoff is not emitted or width-verified.

Compiled native arithmetic embeddings, with live-operand padding:

| Embedding | Q | Static T |
| --- | ---: | ---: |
{embed_table}

These are **building-block embeddings**, not a complete state transform, scan,
controller, or decoder. They exclude pending-prefix operand construction.
Consequently there is no measured full-decoder Q/T point to promote.

## Native components and checks

| Component | Q including operands | Forward static T |
| --- | ---: | ---: |
{native_table}

The shared cofactors test all previous/current quotient degrees 1..3, including
zero exceptional coordinates, native-field polynomial references, inverse and
both composition orders, clean scratch and phase. The carry-finishing slot
selects either X <- cX or X <- X+cY, including c=0 on the MAC branch; public
c and Y are preserved. Its inverse costs 228 T more, included in the screen.

The six-word Schur head uses three regular inversions, conditional word rotations
and triangular corrections; it does not reconstruct moments in a Bennett bank.
It supplies the fixed unit layouts: 21 -> (beta,0,q1,q0,gamma,alpha),
12 -> (beta,alpha,gamma,0,q1,q0), and 3 -> (1,0,beta,q2,q1,q0).
All four quotient patterns have 1,024 native cases in each direction, with
public source preservation, scratch and phase checks. Classification flags are
supplied in these component tests; the global classifier is not emitted.

Conditional series products mask one evaluated operand into a twelve-bit cache
and erase it after the ordinary product. Conditional scaling chooses a nonzero
diagonal; it does not condition a demolition measurement on a quantum bit.
The 8/9-moment probes pass 512 cases each per direction and both composition
orders. They are normalized arithmetic probes, not complete carried heads.

The prefix merge uses one rotation for modes 1/4/5 because their trailing zero
ghosts make the extended maps equivalent. It passes 4,096 native cases per
direction, all seven modes, source flags and phase. This saves **30,168 T** from
the earlier routing screen. Empty field slots are pruned as well.

An early Schur-head version had correct values and scratch but random phase
errors: a constant one was removed after measurement erasure rather than before.
It is rejected and retained as `results/local_schur_rejected_phase_v1.json`.
The corrected version and final slot layout pass all current native checks.

## Revised resource kill

The screen in `profile_fixed_route.py` charges all emitted paths on every input;
there is no frequency discount for rare carries. It includes one controlled
ordinary moment map and additional one-slot/two-slot moment charts, with shared
field slots and routed Schur arithmetic. Controlled bank maps are counted once
per unit/pass; return rotations are deliberately omitted.

| Contribution relative to incumbent, both directions | Delta T |
| --- | ---: |
| dqi-core trims and retained inverses | -25,596 |
| Conditional moment charts | {screen['conditional_moment_delta']:,} |
| Shared scaling/MAC slot | {screen['shared_slot_delta']:,} |
| Routed six-word Schur heads | {screen['routed_schur_delta']:,} |
| Merged bank routing, one way only | {screen['one_way_bank_routing_delta']:,} |
| **Partial added T** | **{screen['partial_projected_added_T']:,}** |

The partial mean projection is **{screen['partial_projected_mean_T']:,.3f} T**.
It excludes carry normalization, chart-return rotations, upper cofactor grafts,
mode derivation, sparse-history updates/cleanup, completion detection/padding,
serial corrections, and pending-head Schur changes. Relative to the accepted
trimmed base, added T is **{screen['partial_added_T_over_accepted_base']:,}**.
These are **estimates from
measured components and exact native-emitter counts**, not a native full decoder.
The +300k kill condition is already met; full assembly and verification stop.

The best concrete change from the old experiment is width-feasible arithmetic,
but separate moment charts consume the T budget. A single shared moment path
across all carry charts, and direct transitions between working bank layouts,
could change this result. Neither is implemented here. This failure rejects the
screened separate-chart realization only, not the correct fixed algorithm.

No full decoder, native zero/one/weight-35/ordinary phase controls, full 102,400
native run, or independent native salt is claimed. No candidate is accepted on
classical failures or component counts. The actual objective remains
R(epsilon_plus)*Q^2*mean_Toffolis with host Clopper–Pearson; no R or replacement
cost is inferred from these components. The accepted trimmed base is retained;
the original round-4 asset remains protected.

## Reproduction and attribution

Run here with at most three CPU processes and single-thread BLAS/OMP:

```bash
source env.sh
g++ -O3 -std=c++17 local_transfer.cpp -o local_transfer
./local_transfer 102400 2026102405 host 3 fixed > results/fixed_host_102400_2026102405_g3.json
./local_transfer 1000000 2026102405 law 3 fixed > results/fixed_law_1000000_2026102405_g3.json
./local_transfer 1000000 2026102406 law 3 fixed > results/fixed_law_1000000_2026102406_g3.json
./local_transfer 4608 2026102423 controls 3 fixed > results/fixed_controls_4608_2026102423_g3.json
"$EXP_PYTHON" core_components.py
"$EXP_PYTHON" adopt_accepted_base.py
"$EXP_PYTHON" routed_cofactor.py --degree 19 --max-quotient 3 --shots 4096
"$EXP_PYTHON" routed_cofactor.py --degree 19 --max-quotient 2 --shots 4096
"$EXP_PYTHON" routed_cofactor.py --degree 19 --max-quotient 1 --shots 4096
"$EXP_PYTHON" shared_slot.py --n 18 --shots 4096
"$EXP_PYTHON" shared_slot.py --n 0 --shots 512
"$EXP_PYTHON" local_schur.py --shots 4096
"$EXP_PYTHON" moment_probe.py --k 19 --n 8 --shots 512
"$EXP_PYTHON" moment_probe.py --k 19 --n 9 --shots 512
"$EXP_PYTHON" routing_merge.py
"$EXP_PYTHON" embed_components.py
"$EXP_PYTHON" profile_fixed_route.py
"$EXP_PYTHON" fixed_final_audit.py
```

Any future full default native verifier and independent salt must run under
`flock /tmp/r6-fullverify.lock`; this resource kill precedes that stage.
All new files stay in this directory. Inputs and all reused source are read-only.
Round-4 Codex supplies compact identities, regular arithmetic, compiler and field
conversion. dqi-repair supplies host-checked cases/field tables and prior diagnosis.
dqi-core supplies the two checked trim components. Field/gate arithmetic adapts
Apache-2.0 andigu/uncompute revision 64b5a13aa249f1360846657392e6c66817b26d0f;
the source and license remain in `../../3/vendor/uncompute`. The prior matrix
experiment is retained in [RESULTS.matrix-history.md](RESULTS.matrix-history.md).
'''
(ROOT/'RESULTS.md').write_text(text)
sources=list(ROOT.glob('*.py'))+list(ROOT.glob('*.cpp'))+list(ROOT.glob('*.hpp'))+[ROOT/'env.sh',ROOT/'RESULTS.md',ROOT/'WIDTH_PLAN.md']
refs=[ROOT/'../dqi-core/components'/n for n in ['trimmed_middle.py','budget_inverse.py']]+[ROOT/'../../4/codex'/n for n in ['build_cubic_compact_decoder.py','compact_product_decoder.py']]
refs += [ROOT/'../dqi-core'/n for n in ['build_accepted_decoder.py','components/trimmed_decoder.py','results/accepted_decoder.json']]
report=dict(status='PASS',disposition='KILL_SEPARATE_CARRY_MOMENT_REALIZATION',baseline_sha256=baseline_hash,
            adopted_base_circuit=adopted['circuit'],adopted_base_sha256=adopted['sha256'],
            full_decoder_emitted=False,full_native_verified=False,replacement_accepted=False,
            classical_cases=[dict(label=c[0],cases=c[1],seed=c[2],wrong=r['wrong_decodes']) for c,r in zip(cases,data)],
            native_components=names,max_compiled_arithmetic_embedding_Q=embedded['max_Q'],
            partial_projected_mean_T=screen['partial_projected_mean_T'],partial_projected_added_T=screen['partial_projected_added_T'],
            source_sha256={str(p.relative_to(ROOT)):sha(p) for p in sources},readonly_dependency_sha256={str(p.relative_to(ROOT)):sha(p) for p in refs},
            native_artifact_sha256={str(p.relative_to(ROOT)):sha(p) for p in (ROOT/'candidates').glob('*.kmb') if p.name.startswith(('routed_','routing_merge','local_schur','shared_slot','moment_probe'))},
            evidence_sha256={str(p.relative_to(ROOT)):sha(p) for p in (ROOT/'results').glob('*.json') if p.name!='fixed_final_audit.json'})
(ROOT/'results'/'fixed_final_audit.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if not k.endswith('sha256')}))
