# Content-aware chemistry QROM compilation

**A. Verified improvement:** [`candidates/final`](candidates/final) passes the
unchanged full native chemistry verifier with two independent salts. Cost is
**697,327,965.8381622**, versus **708,547,388.9553173**: a **1.5834% cost
reduction**, or **1.016089× improvement**. The 1.10× and 1.20× goals were not reached.

The accepted R=1/B=719/C=350 physical factorization is retained. The selected
NPZ only reorders dictionary rows and corresponding U/W rows. One-body rows
remain in their own sector. No physical fitting, tolerance, scoring-rule,
reference-tensor, or host-verifier change was made.

All project writes are confined to this experiment. Read-only references were
the accepted round-4 Codex implementation, completed `chem-qfit` builders/report,
completed `inv-dyadic` report/streaming eraser, earlier compiler dependencies,
and immutable host infrastructure. No concurrently active agent directory was
consulted. No dependencies were installed and no Git state changed.

## Final native resources

| Resource | Accepted baseline | Frozen final |
|---|---:|---:|
| lambda_eff | 48.59417497035283 | 48.59417497035283 |
| Batch J / oracle count | 7 / 11 | 4 / 19 |
| ANGLE forward + reverse Toffolis | 10,604 | 11,708 |
| Inner / outer preparation Toffolis | 17,333 / 317 | 17,333 / 317 |
| Givens Toffolis | 2,250 | 2,250 |
| Total Toffolis | 30,504 | 31,608 |
| Persistent host qubits | 141 | 117 |
| Maximum temporary host qubits | 337 | 337 |
| Total host qubits | 478 | 454 |
| Maximum individual oracle qubits | 123 | 128 |
| Chemistry cost | 708,547,388.9553173 | 697,327,965.8381622 |

The score win comes from reducing persistent angle-word storage by 24 qubits.
ANGLE and total Toffolis increase relative to J=7. Inner preparation still sets
the temporary maximum. Oracle 0 uses **292 T /124 Q**; each of the other eighteen
uses **309 T /128 Q**.

Final evidence:

- [`final_verify_A.json`](results/final_verify_A.json) and
  [`final_verify_B.json`](results/final_verify_B.json): unchanged full verifier,
  **50,880 shots each, 64 per valid address**, identical scored resources.
- [`final_check.json`](results/final_check.json): every oracle and valid address,
  eight arbitrary incoming target samples/address, both directions, unchanged
  selection, clean scratch and no phase; also complete-family recovery,
  including signed words and final word padding.
- [`final_refit_rebuild.json`](results/final_refit_rebuild.json): feature search
  and quadratic fits rerun from source, then all nineteen circuits and the NPZ
  rebuild **byte-for-byte**. A separate frozen-model rebuild passes in
  [`final_rebuild.json`](results/final_rebuild.json).
- [`final_freeze.json`](results/final_freeze.json) binds exact bytes and the two
  subsequent independent salts. [`audit.json`](results/audit.json) checks those
  hashes, complete native workloads, cap, host width and unchanged accepted
  source bytes. [`source_manifest.json`](results/source_manifest.json) records
  compiler and read-only dependency hashes.

Final submission digest:
`1f51d105b22d3c1eb7715b3f6ba2303fd9fcf18d22aac3f5a27d9b31d57b80d0`.
Exact circuit/NPZ SHA256 values are in the freeze manifest.

Native physical metrics are error1 `0.011986315442796814`, error2
`0.019995304490553204`, Lambda `303.5283147431962`, gap `3.915157472365536`.
Differences from baseline are summation roundoff after permutation. The
unpermuted quadratic candidate has exactly identical physical metrics and a
byte-identical physical NPZ, and also improves cost (698,078,065.5230045).

## Baseline and actual-table analysis

[`baseline_reproduction.json`](results/baseline_reproduction.json) records all
eleven rebuilds: **482 T /123 Q each, exact byte matches** to
`initial_experiments/4/codex/candidates/norm_cycle2_signed`.
[`baseline_verify.json`](results/baseline_verify.json) reproduces the native
incumbent score. A local reference is in `candidates/baseline`; accepted source
files remain untouched.

The table has 795 distinct rows and 75 eight-bit angles. Every word obeys
bit7=bit6. Conjugating public bit7 by bit6 before and after lookup permits seven
payload bits while preserving the XOR interface on **arbitrary incoming words**,
including incoming values that are not sign-extended angles.

[`content_analysis.json`](results/content_analysis.json), reproduced by
`analyze.py`, records every J=7 delta table:

- Every paired-difference bank has **rank49**, all 49 bits active and 398 distinct
  differences. A global output XOR-basis change alone cannot shrink this bank
  in the incumbent row order.
- Delta-word entropy is **6.22–6.91 bits** out of seven. Zero-word fractions are
  0–4.87%; difference bit density is close to one half.
- Complete multi-output subtrees of sizes 2,4,8,16 are distinct. The only
  constant pair is the duplicated last row for address padding. This scopes
  the observed failure of repeated-subtree elimination to this representation.
- All 795 care-prefix vector ANF coefficients are nonzero. Shared address
  monomials still leave 784 nonlinear monomials.

## Compiler 1: select lookup and care-prefix streamed erasure

The accepted compiler loads the first row of each pair into public output and
the pairwise XOR difference into clean scratch. One controlled XOR per bank
coordinate selects the second row. Measurement and ANF phase repair erase the
bank without disturbing arbitrary incoming public values.

Baseline erasure zero-pads 398 bank rows to 512, creating monomials needed only
for invalid addresses. The Boolean Möbius transform is triangular in numeric
order: coefficient m uses only rows that are bit-subsets of m. Setting
coefficients numbered >=398 to zero therefore preserves every valid address.
This changes unspecified invalid-address behavior, **not** padding words at
valid addresses.

[`erase.py`](erase.py) prunes unused monomials and chooses a split. The 9-bit
paired address needs **30 ANDs rather than 37** at split 4/5. Keeping the small
bank and streaming the other needs **14 scratch wires rather than 30 resident
monomial wires**, at the same Toffoli count. This extends the completed inverse
experiment's streaming method with care-set completion and monomial pruning.
Conditioned CZs between address monomials repair measurement phase before freeing
workspace.

The complete `pruned7_l8` family passes native verification: 30,350 total T,
Q=478, cost **704,970,274.5473996**. A parity-constrained row pairing also gives
`parity0`: rank 44 instead of 49 in oracle 0, total 30,340 T and cost
**704,737,994.3910412**. Its order/inverse are in `parity0_transform.json`.

## Compiler 2: direct multi-output streamed ANF

[`algebraic.py`](algebraic.py) computes shared address monomials, XORs each
complete coefficient vector into public output, and deletes the monomial by
measurement. It uses no difference bank or lookup trie, with at most eight
monomial wires simultaneously live.

It costs **784 T/oracle**. J=7 uses 74 oracle qubits but costs 862,874,324.84
overall. J=13 fits six oracles at 122 Q each and lowers total T to **29,308**,
but raises host width to **526**, for cost **749,128,190.10**. Both families
pass all-address arbitrary-target and complete-family native component checks.
These costs use compiled circuits and unchanged host formulas; they are not
claimed as full native chemistry verifier runs. Neither beats the incumbent.

## Compiler 3: four-row lookup with nonlinear bank compression

For a=4h+l, expand four payloads in the low bits x,y:

```
D(h,x,y) = c0(h) XOR x*c1(h) XOR y*c2(h) XOR (x*y)*c3(h)
c1 = D0 XOR D1
c2 = D0 XOR D2
c3 = D0 XOR D1 XOR D2 XOR D3
```

At J=4, c0 is 28 bits and the three difference words need 84 clean bits. The
199-row high lookup costs 197 T, selection 84 plus one for x*y, and care-prefix
erasure 19: **301 T /133 Q**, violating the cap. Shrinking prefix recursion
incurs large equality-check overhead: **355 T /130 Q** at subtree depth 4.
All 66 component-grid circuits were compiled and checked on every valid address,
including transient over-cap research points. Over-cap points were not saved
as valid submissions. Grid Toffoli counts are static compiled counts; deep
prefix equality cleanup can execute fewer gates on some measurement branches.
All competitive score claims instead use the full native mean. For the final
depth-seven compiler, static and native Toffoli counts coincide.

The first working remedy, [`local4.py`](local4.py), omits five c1 columns.
On three high-address regions of 64,64,71 rows, each omitted bit is an exact
linear parity of the remaining 79 columns. Streaming those parities fits 128 Q.
The family costs **314 T/oracle** and fully verifies at total 31,832 T, Q=454,
cost **702,269,799.0559472**. Decoder matrices are retained in its build JSON.

The better [`quadratic.py`](quadratic.py) uses a global shared decoder:

```
missing(z) = A(z) XOR sum_k output_k * (L_k(z) AND R_k(z))
```

z is the retained79-bit bank; A is five affine parities; L/R are affine parities;
output_k is a five-bit vector. Fixing L_k/output_k makes the remaining fit
exactly linear over GF(2). There are **995 exact constraints** (199 rows times
five omitted bits). Eight shared products provide 1,040 scalar affine
coefficients before gauge dependencies. Exact elimination finds a solution
for every table; every fitted row is checked. There is no numerical fit or
Hamiltonian approximation.

Each product is XORed into its output vector with CNOT conjugation and one CCX.
Five decoded wires plus two parity wires use exactly the seven wires left after
public registers and the 79-bit bank. Applying decoded bits costs five controlled
XORs. Measuring them leaves an affine/quadratic phase in z; Z and CZ of recomputed
affine parities repair it with **zero further Toffolis**. No decoded nonlinear
value remains live during bank erasure.

Actual count is
`197 lookup +79 retained select +1 low monomial +8 products +5 decoded select +19 erasure =309`.
The unpermuted complete family passes full native verification at **31,642 T /
454 Q**, cost **698,078,065.5230045**.

## Address order, affine decoding and exact gauges

Greedy Hamming order reduces mean adjacent distance across the 525-bit table
from **259.90 to 223.72**, but actual J=7 cost remains 475 T/oracle after
care-prefix repair. Lexicographic order also leaves the count unchanged. Both
compiled families pass native component checks. Lower Hamming distance is not
used as evidence of lower Toffoli cost.

[`hybrid.py`](hybrid.py) groups rows in fours by three independent linear
features of oracle 0's 28-bit payload. Every feature-class count is divisible by
four **separately in the physical sectors**. The last two-body row is duplicated
only for invalid-address completion. An invertible 28-bit XOR basis puts these
features into output coordinates. Their low-bit differences vanish, eliminating
**nine bank columns**, 84 ->75. Oracle 0 then costs 292 T /124 Q; the other eighteen
use the 309 T quadratic compiler on the same permuted dictionary.

`U[:,new_b,:]=U[:,old_b,:]` and `W[:,new_b,:]=W[:,old_b,:]` use the identical
two-body permutation; H2_I stays fixed. One-body rows stay in their sector and
the native verifier matches their eigenvalues by overlap. Maximum reconstructed
X change is **9.99e-16**. Reversing the saved permutation recovers the original
words and physical arrays exactly. `final_transform.json` contains order,
inverse and output-basis matrix; `transform_inverses.json` contains its inverse
matrix and exact CNOT/SWAP decoding sequence.

[`gauges.py`](gauges.py) searches all **128 linear high-bit Givens gauges**.
At each reversed edge, either sign of the accumulated norm shifts that angle
by pi. Choose the sign so bit7 is any selected parity of bits 0..6; earlier
edges absorb the signed norm, leaving a global vector sign. All 128 tables
decode to the intended projectors, maximum vector error up to sign **8.88e-16**.
Re-encoding every table with the original signed convention restores **every
original word exactly** (`gauge_search.json`, `gauge_inverses.json`).

All 128 gauges retain rank49 for every J=7 difference bank. Five representative
families (masks0,1,63,64,127) were actually compiled and passed arbitrary-target
and complete-family native checks; all cost 475 T/oracle. A reversible Gray-code
payload transform also gives 475 T/oracle. An initial Gray implementation mixed
the redundant high bit into the payload and was rejected by the arbitrary-target
test. Masking to seven payload bits before shifting fixes it; the corrected
family passes. No rejected circuit or favorable random sample enters final claims.

## Complete measured comparison

Every row uses complete compiled circuits and native component checks. 'Full'
means the unchanged chemistry verifier also ran at 64 shots/address, including
physical reconstruction and caps. Other costs use unchanged host formulas and
held-fixed normalization.

| Family | J/K | ANGLE T, both passes | Host Q | Max oracle Q | Cost | Full |
|---|---:|---:|---:|---:|---:|---|
| Baseline |7/11|10,604|478|123|708,547,388.96|yes|
| Care-prefix two-row |7/11|10,450|478|123|704,970,274.55|yes|
| Parity-paired two-row |7/11|10,440|478|123|704,737,994.39|yes|
| Two-row J=6 |6/13|12,168|470|108|732,409,461.39|component|
| Two-row J=5 |5/15|13,830|462|93|757,255,663.05|component|
| Two-row J=4 |4/19|17,252|454|78|819,638,337.98|component|
| Four-row J=3 |3/25|14,000|446|104|734,714,769.05|component|
| Direct ANF J=7 |7/11|17,248|478|74|862,874,324.84|component|
| Direct ANF J=13 |13/6|9,408|526|122|749,128,190.10|component|
| Local linear four-row |4/19|11,932|454|128|702,269,799.06|yes|
| Quadratic four-row |4/19|11,742|454|128|698,078,065.52|yes|
| **Frozen hybrid** |**4/19**|**11,708**|**454**|**128**|**697,327,965.84**|**yes, twice**|

## Reproduction

Run from the project root. The environment uses existing compatible Python,
kickmix and C++ runtime installations, disables bytecode writes to read-only
imports, restricts numerical threads to one, and keeps caches/temporary files
inside this experiment.

```bash
source initial_experiments/5/chem-qrom/env.sh

# Re-run feature search and quadratic fits; require identical frozen bytes.
"$EXP_PYTHON" "$QROM_ROOT/reproduce.py" --refit
"$EXP_PYTHON" "$QROM_ROOT/reproduce.py" --rebuild --audit

# Arbitrary incoming words, every valid address, both directions.
"$EXP_PYTHON" "$QROM_ROOT/check.py" "$QROM_ROOT/candidates/final" \
  --repeats 8 --seed 2026102493

# Full unchanged native checks; salts were fixed after freezing bytes.
"$EXP_PYTHON" -m dfthcsos.verify "$QROM_ROOT/candidates/final" \
  --salt chem-qrom-final-independent-A-20261024
"$EXP_PYTHON" -m dfthcsos.verify "$QROM_ROOT/candidates/final" \
  --salt chem-qrom-final-independent-B-20261024

# Alternative architectures/batching, including all-address component tests.
"$EXP_PYTHON" "$QROM_ROOT/survey.py"
"$EXP_PYTHON" "$QROM_ROOT/algebraic.py" --batch 13
"$EXP_PYTHON" "$QROM_ROOT/check.py" "$QROM_ROOT/candidates/anf13"

# Gauge/permutation research and actual-table statistics.
"$EXP_PYTHON" "$QROM_ROOT/gauges.py"
"$EXP_PYTHON" "$QROM_ROOT/permutation.py" greedy
"$EXP_PYTHON" "$QROM_ROOT/permutation.py" lex
"$EXP_PYTHON" "$QROM_ROOT/permutation.py" parity
"$EXP_PYTHON" "$QROM_ROOT/analyze.py"

# Unpermuted and combined four-row builders; final/ is never overwritten.
"$EXP_PYTHON" "$QROM_ROOT/local4.py"
"$EXP_PYTHON" "$QROM_ROOT/quadratic.py"
"$EXP_PYTHON" "$QROM_ROOT/hybrid.py"
```

`*_build.json` includes exact resources and decoder coefficients.
`*_check.json` binds native component checks to hashes. The original low-depth
J=7 resource failure is retained as `pruned7`; it is not `pruned7_l8` or final.
No cap-violating family is presented as an accepted submission.

## Assessment and limits

**A. Verified improvement:** final beats the incumbent by 11,219,423.1171551
cost units, retains the physical factorization, and passes two independent
complete native runs and exact rebuilds.

**B. Validated component advances:** care-prefix streaming removes seven
Toffolis per two-row oracle. Shared quadratic interpolation and Clifford phase
repair reduce four-row peak width from 133 to 128 and make smaller batching
useful. These are preserved working compilers, not abstract gate estimates.

**C. Scoped negative results:** direct shared ANF loses after persistent
storage is priced. Greedy/lexicographic Hamming order and 128 affine Givens
gauges do not improve this compiler's count. Low-workspace prefix equality
incurs substantial gates. These results do not establish global incompressibility.

**D. Unresolved:** jointly choosing variable quadratic factors, output vectors,
omitted columns and address order may reduce decoder cost. Larger batches
need stronger bank compression under the cap; smaller batches need fewer passes.
At Q=454, 1.10× incumbent improvement requires about **2,411 further total
Toffoli savings** beyond final. Small constant decoder changes cannot supply
that scale. A stronger compiler discovery remains necessary; no unsupported
global lower bound or optimality claim is made.
