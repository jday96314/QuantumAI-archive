# Shared-tree chemistry representation research

**Complete; no replacement selected.** The incumbent remains the native-verified
R1/B719/C350 eight-bit model costing **708,547,388.9553173**. Fresh complete native
verification reproduces its original digest and cost.

Twenty tree graphs (17 distinct binary amplitude lattices) were screened on
each of four controls. Adaptive topology search improves B351's rounded fit,
and conditional integer repair reaches **54.288581 mHa**, versus **58.690849
mHa** under the incumbent tree. Both remain invalid. New-tree B479 repair
reaches **22.500832 mHa**, but its independently checked unlimited-rank PSD
lower bound is **20.087505 mHa**. The legal shifted-grid implementation also
fails accuracy and resource competitiveness. No 1.20x or 1.30x score gain is
claimed. This is a scoped negative result with modest representation gains.

All writes are confined to this directory. Completed `chem-qfit`, `chem-discrete`,
and `chem-psd` reports were read before implementation. Those experiments and
round-4 inputs are read-only. Their exact-dictionary exclusions do not exclude
other ordered trees or integer dictionaries.

## Host semantics established

The unchanged `dfthcsos/scoreinfra/dfthc_math.py` validates 75 integer edges,
vertices in [0,76), visited={0}, parent visited, child unvisited on every edge.
Thus every rooted spanning tree is allowed, but arbitrary edge permutations are
not. Root must remain physical orbital 0. Decoding starts at e0 and applies
`(u_p,u_c) <- (cos(theta)u_p-sin(theta)u_c,
sin(theta)u_p+cos(theta)u_c)` in listed order, theta=pi*word/128.
Trees are shared across every address, including all 76 one-body vectors.
Physical orbital labels are never permuted in the Hamiltonian.

Sibling rotations sharing a parent generally do not commute. Disjoint edge
swaps do commute and do not create a different representable set. Reverse-order
Givens elimination provides a continuous encoder for every permitted tree.
Sign gauges differing by integer multiples of pi have equivalent word lattices
at eight bits; this will receive an explicit independent test.

## Controls and acceptance

Controls include the actual incumbent, B639 float and integer basins,
R2/B351/C189 continuous basin, and saved exact R2/B479/C175 split. The search
uses balanced, heap, chain, star, hierarchical correlation/amplitude, host-tensor
clusters and sensitivity-selected Huffman trees. Exact decoded projector
vectors, the full signed native pair tensor, rank-constrained PSD weight refits
and normalization determine the conclusions. Exact-word and topology mutations
follow the initial screens; matched budgets and best states are preserved.
Only a complete unchanged native physical/oracle pass with lower actual
lambda_eff*T*Q may be reported as a verified improvement.

## Calibration and first matched screens

`results/encoding_tests.json` tests 16 tree graphs on random dense vectors and
all 76 basis vectors. Continuous inverse round trips, an independent full
rotation-matrix product, the unchanged host decoder and signed canonicalization
agree at floating-point roundoff. The incumbent re-encodes with **zero changed
words**. Native error2=19.9953044905532 mHa; direct native pair error differs by
less than 2e-13 mHa. Native Lambda/gap/lambda_eff reproduce the accepted values.
`results/search_gradient_control.json` independently checks the integer-search
Jacobian (max error 5.4e-12) and physical loss gradient (3.1e-11 absolute).

Twenty tree graphs per control, including two Huffman trees and two host-tensor
orbital clustering trees, use identical
12-cycle rank-C refit budgets. The tensor Huffman tree weights orbital masses
by ||Z_i||², Z_i=sum_c w_ic X_c. This minimizes the isotropic part of expected
tensor-Jacobian quantization sensitivity, rather than mean vector L2.

| Control | Source error mHa | Incumbent-tree refit mHa | Best alternative-tree refit mHa |
|---|---:|---:|---:|
| Accepted R1/B719/C350 | 19.9953 | 19.5330 | 94.3860 (tensor Huffman) |
| R1/B639/C350 continuous | 15.1118 | 101.5521 | 100.3060 (tensor Huffman) |
| R2/B351/C189 continuous | 19.9788 | 172.0871 | 170.2827 (power-balanced) |
| R2/B479/C175 saved exact | 22.6233 | 22.6200 | 94.6754 (correlation clusters) |

These are exact-word numerical tensor fits, not native candidate passes.
Single-block solves reach their analytic rank-C optima; two-block 12-cycle
values are matched-budget upper bounds, with later longer refits kept separate.
`results/screen_*.json`, `results/screen.log`, `results/screen_huffman.log`
and `candidates/<control>__<family>.npz` preserve all tables and exact words.

An important equivalence: reversing cluster representatives changes graph
topology/depth but leaves decoded quantization unchanged to roundoff. The
induced **binary amplitude partition**, reconstructed by reverse union of edge
supports, identifies these cases. On an even uniform angle grid, swapping the
two sides of a binary split uses theta -> pi/2-theta plus signs, which preserves
the lattice. `binary_signature` prevents counting these as independent basin
changes. Balanced, star, chain, and differing physical cluster partitions are
not collapsed into this equivalence. Gauge changes alone also preserve the
uniform-grid projector lattice.

## Implementation corrections retained

The preliminary canonicalizer changed +pi/2 to -pi/2 without propagating a
radius sign. Independent projector tests caught this before conclusions were
drawn. Those invalid preliminary outputs are quarantined in
`results/aborted_encoder/`; all reported screens were rerun with the corrected
reverse encoder. A joint-search JSON NumPy-bool serialization error terminated
an early run; `joint_b479_aborted_serialization.log` is retained and the run
was restarted. Synthetic width testing also rejects J7 for dense nine-bit
words (130 qubits including internal unlookup scratch); legal nine-bit
codebooks use J6, explicitly charged, instead of assuming the old batch width.

## Topology and codebook searches

`topology_search.py` runs orbital-label mutations, sibling order changes and
subtree reparenting with exact re-encoding/refitting. A 240-proposal B351
representation search improves 172.087 to 154.874 mHa at the 12-cycle budget;
300 additional cycles reach about **140.400 mHa**. This remains invalid.
The 12-cycle score is used throughout selection; polishing is not credited as
a topology gain. Joint searches start from actual exact dictionaries and
compare every topology proposal against a matched fixed-word refit.

`phase_grid.py` tests shared per-edge half-step offsets. These are implemented
as **nine-bit host words with constant low bits**, preserving the eight-bit
grid spacing and number of signed code points. There is no changed host
interpretation. Actual output widths, all ANGLE circuits, and nine-bit Givens
costs are charged. This is a different legal representation, not an eight-bit
native score claim.

## Frozen geometry diagnostics and complete physical/resource screens

The native pair-space union-span relaxation permits arbitrary cross-block
coefficients, without positivity or rank restrictions. QR reconstruction of the
orthogonal outside-span residual gives these fixed-dictionary bounds:

| B351 dictionary | Rank-C fit mHa | Union-span lower bound mHa |
|---|---:|---:|
| Continuous source, C189 | 19.9788 | 14.5572 |
| Incumbent-tree rounded words | 172.0871 (12 cycles) | 140.1961 |
| Optimized new tree | 140.3998 (300 more cycles) | 124.3681 |
| Optimized half-step phase codebook | 145.8616 (300 more cycles) | 129.9280 |

These bounds exclude those particular word tables even with unrestricted
coefficients. They do not exclude other exact words or other trees. The
continuous B351 feature blocks have condition numbers 1525/2210, weight
Frobenius norm 55.18 and tensor row-sensitivity norm 25.25. The incumbent has
condition number 172, weight norm 3.72 and sensitivity norm 4.40. These are
diagnostics of strong cancellation and quantization sensitivity, not universal
capacity limits. `results/geometry_audit.json` preserves values and hashes.

All following rows optimize every identity offset with native shifts, fit
one-body integer words, run eight Gram-preserving column-gauge sweeps, compile
every ANGLE circuit, and use unchanged host resource arithmetic. These finite
normalization optimizations are not normalization lower bounds.

| Exact dictionary | error1 / error2 mHa | lambda_eff | T | Q | Conditional product |
|---|---:|---:|---:|---:|---:|
| Topology-only B351/C189 | 14.9799 / 140.3998 | 121.7500 | 22516 | 461 | 1,263,749,931.82 |
| Phase-grid B351/C189, brot9/J6 | 14.9515 / 145.8616 | 120.2700 | 20688 | 778 | 1,935,777,129.54 |
| Tensor-Huffman B639/C350 | one-body passes / 100.3060 | 92.7869 | 27794 | 478 | 1,232,723,173.58 |

No row passes two-body physics. Nine-bit pricing has a material discontinuity:
native inner preparation chooses k2=5 instead of 4, increasing temporary width
from 321 to 640 for B351/C189. Its break-even lambda_eff is only **44.0221**,
versus **68.2616** at brot8. Nine-bit output bits cannot be treated as free
even when their low bit is shared and emitted with X gates.

`grid_oracle_controls.json` passes synthetic worst-case forward/reverse XOR,
arbitrary incoming targets, all addresses, clean scratch and phase under two
independent seeds. `actual_oracle_components.json` passes all actual circuits
for the three frozen screens. The shifted-grid candidate additionally follows
the full unchanged native recovery path with **49,792 shots**:
`grid9_b351_submission_native.json`, status **REJECT**, solely error2 above
20 mHa; all recovered words match exactly. This is a legal representation and
verified oracle component, not a verified chemistry improvement.

The preserved incumbent separately passes full unchanged native verification
with **50,880 shots**, independent salt, original digest and original score
708,547,388.9553173 (`results/incumbent_native.json`).

## Matched searches and integer repair endpoints

Both adaptive and random B351 searches start with the same 15 distinct binary
families, evaluate 240 additional trees with 12 two-block coordinate cycles
each (5,760 block minimizations), and polish their selected state for 300 cycles.
This comparison matches proposal and spectral-refit work. Construction/QR
overhead and shared-device scheduling differ; wall times are descriptive.

| B351 representation search | Selected 12-cycle error mHa | After 300 more cycles mHa |
|---|---:|---:|
| Adaptive population: orbital/sibling/reparent moves | 154.874226 | 140.399754 |
| Independent random balanced trees | 167.496453 | 152.594524 |

The subsequent exact-word solver changes 1/2/4/8 angles in each proposed row,
combines up to 128 rows, screens actual decoded projector increments using the
full residual and dictionary-orthogonal curvature, then accepts only after
rank-C physical refitting. A fixed-word refit precedes each proposal, so weight
convergence alone cannot count as an integer improvement.

| Conditional B351 repair | Search block updates | Common final polish updates | Final error mHa | Final union-span bound mHa |
|---|---:|---:|---:|---:|
| Selected new tree + word moves | 4004 | 600 | 54.288581 | 45.455235 |
| Incumbent tree + word moves | 4016 | 600 | 58.690849 | 49.720653 |
| Selected new tree, fixed words only | 4600 | included | 131.803563 | 124.368118 |
| Incumbent tree, fixed words only | 4600 | included | 146.962010 | 140.196071 |

Limits round to complete block cycles; the old-tree repair receives 12 extra
block updates. The new-tree repair is about 7.5% better **conditional on the
chosen tree**. Selecting that tree consumed an additional 5,760 proposal block
updates plus initial family evaluations. This is not an equal-total-work claim
for the entire new-tree pipeline. The fixed-word controls also receive slightly
less work than the repair-plus-polish totals (4600 versus 4604/4616). All four
are invalid; the exact-word span bounds exclude their saved dictionaries.

Joint annealing on B479 evaluates 320 tree mutations, accepts 22 different
topologies including 21 genuinely uphill states, and resets to the best state
every 40 steps. It finds **no better topology**. Its best checkpoint has the
original tree and zero changed source words; improvement to 22.510480 mHa is
entirely coefficient refitting, not a topology gain.

`local_nni.py` enumerates binary nearest-neighbor restructurings
`A,(B,C) -> (A,B),C` / `(A,C),B` on supports of at most eight orbitals. Only
two continuous angles per row change, up to exact lattice gauge maps. It tests
88 distinct neighbors on B479 and another 88 on the saved exact B639 basin.
No nearest re-encoding improves the matched fixed-word baseline. B639 retains
its exact original words at **29.843748 mHa**. The best rejected B479 topology
is separately repaired by integer moves, rather than being discarded after
nearest rounding:

| Conditional B479 word repair | Search / polish block updates | Error mHa | lambda_eff | Actual T / Q | Conditional product |
|---|---:|---:|---:|---:|---:|
| New local binary topology | 1604 / 600 | 22.500832 | 71.758326 | 27981 / 461 | 925,627,945.09 |
| Original topology | 1604 / 600 | 22.512244 | 71.476967 | 27981 / 461 | 921,998,629.65 |

The small accuracy advantage loses on normalization and conditional product.
NNI selection cost is additional to this conditional repair comparison.
This is an explicit example of reconstruction improvement failing the actual
resource objective.

After the same eight-sweep normalization screen, repaired new-tree B351 has
error1/error2 **14.955621/54.288581 mHa**, lambda_eff **88.331252**, T/Q
**22516/461**, conditional product **916,867,444.73**. The repaired incumbent-tree
B351 control has error2 **58.690849 mHa**, lambda_eff **100.490239**, product
**1,043,076,222.35**. Neither is a scored submission. For a 1.20x improvement,
B351/C189 requires lambda_eff <= **56.884702**, B639/C350 <= **44.443546**,
and B479/C175 <= **45.774488**, in addition to physical validity.

## Convex certificates for the new B479 words

`psd_capacity.py` fits unlimited-rank positive covariances in the correctly
whitened block metric. For residual D=H-M, sequentially subtracting each
positive block compression yields Q_r D Q_r.T <= 0. A further explicit scalar
safety shift makes the largest compressed eigenvalues smaller than -1e-11.
Then every admissible block-PSD model satisfies

    ||H-M_admissible||² >= 2<H,D> - ||D||².

The native signed target, exact words and full dual arrays are preserved.
An independent NumPy QR/eigensolve and the unchanged host vector decoder check
the certificates in `results/independent_dual_audit.json`.

| Frozen B479 word dictionary | Full-PSD feasible residual mHa | Feasible dual lower bound mHa |
|---|---:|---:|
| New NNI topology + integer repair | 20.160130 | **20.087505** |
| Original topology + integer repair | 20.149991 | **20.071715** |

Both exact dictionaries are excluded at every C, not merely C175. The duals
have roundoff safety slack but are not interval-arithmetic proofs. Their
0.072–0.088 mHa excluded margins are much larger than the calibrated native
pair discrepancy. The small rank-C accuracy gain of the new tree does not
improve its relaxed PSD capacity. Further coefficient-only work on these
exact words cannot solve the task; different words or a different geometry
are necessary.

## Why the tested continuous basins remain fragile

For feature p_i=sym(u_i u_i.T), F_c=sum_i w_ic p_i and Z_i=sum_c w_ic F_c,
an angle perturbation gives

    dM_i = dp_i Z_i.T + Z_i dp_i.T
    ||dM_i||² = 2||dp_i||²||Z_i||² + 2(dp_i.Z_i)².

For a tree split carrying probability mass m_i, ||dp_i||²=2m_i. Thus the
isotropic expected tensor error from independent uniform small angle errors
is proportional to sum_i ||Z_i||² times the sum of internal split masses.
Weighted Huffman coding minimizes this term exactly. `compare_scores.py`
checks that its tensor-Huffman tree attains the minimum across all screened
families. This does not include anisotropic sensitivities, coherent rounding,
nonlinear tensor terms or weight refitting.

Surrogate rankings are empirical, not universal: after excluding chain/star
and equivalent binary trees, B639 Spearman correlation with the rank-C error
is 0.727 for physical frozen-weight scoring versus 0.464 for vector RMS.
For B351, vector RMS performs better (0.776 versus 0.675 for frozen tensor
error); isotropic expectation is only 0.376. Incumbent retreeing is poorly
predicted by all these local surrogates. Exact decoded rank-C fits therefore
remain the acceptance criterion (`results/surrogate_comparison.json`).

`span_robustness.py` finds that B351 quantization perturbs feature matrices by
operator norms **27–42 times their smallest continuous singular values**.
Median span cosines remain about 0.99975, but the worst principal angles are
82–86 degrees. The optimized tree has similarly fragile directions (82–88
degrees). Small average vector errors can therefore erase a few critical
cancellation directions while leaving most of the span nearly unchanged.
This supports a basin-conditioning explanation in addition to lattice choice;
it is not a theorem that the B351 shape or all eight-bit trees are infeasible.

Different gauges, signed orbital conventions and representative changes that
preserve binary splits do not expand the even uniform eight-bit projector
lattice. An eight-bit nonuniform codebook emitted as ordinary eight-bit host
words is only a subset/relabeling of the already allowed words. More physical
angles require more legal output precision and must pay the host costs, as
the nine-bit experiment explicitly demonstrates.

## Verification, preserved artifacts and reproduction

Complete unchanged native records for frozen new submissions all reject only
the two-body accuracy gate; one-body accuracy, strict floors/domain, all ANGLE
words and per-oracle width pass:

| Record under results/ | Native status | Full shots | Error2 mHa |
|---|---|---:|---:|
| incumbent_native.json | PASS | 50880 | 19.995304 |
| grid9_b351_submission_native.json | REJECT | 49792 | 145.861573 |
| repair_topology_b351_submission_native.json | REJECT | 49792 | 54.288581 |
| repair_nni_b479_submission_native.json | REJECT | 66176 | 22.500832 |

The native digests match every delivered circuit/NPZ family. The incumbent's
digest remains `ebf8e357a1b693fc8589d5c48a5d6fa10343dfa5944724d8f7a77b8edc2157f7`.
`results/artifact_audit.json` independently checks **1183 exact-word checkpoints,
939893 decoded vectors, 934 tree graphs and 914 binary lattices**. Decoder
disagreement is zero; stored projector disagreement is <=1.34e-15.
Quarantined invalid preliminary encoder outputs are explicitly excluded.

`trajectories/*_dictionaries/` retain the best exact integer table per evaluated
ordered tree. Compact archives contain uint8 words/edges, brot, model shape
and refit error; decoding requires only `trees.decode`. Main `candidates/*.npz`
also retain explicit float64 weights, and physical NPZs retain identity offsets
and one-body words. `archive_search.py` reconstructs early trajectories from
saved decisions and verifies every tree signature; it does not search or select
new scores. Archives cover 255 adaptive, 251 random and 311 annealing trees,
plus the NNI and family screens. Source/host hashes, normalization records,
exact circuits, certificate arrays and all native JSONs are retained.

`artifact_manifest.json` binds delivered bytes; `selection.json` explicitly
selects no replacement. `NOTICE.md` attributes all completed read-only source
methods and plotting dependencies. No host assets, reference tensors, shared
Git state or external environments were changed. Resource-heavy scripts limit
CPU threads to two; some independent jobs overlapped on the RTX 4090. Timing
differences are not used as evidence of optimizer superiority.

From this directory:

```bash
source env.sh
bash verify_saved.sh                    # exact artifact and certificate audits
bash verify_saved.sh --native           # also repeat all complete native paths
bash reproduce_searches.sh              # full numerical research reproduction
"$TREE_PY" compare_scores.py
"$TREE_PY" span_robustness.py
"$TREE_PY" plot_results.py
```

Fresh search trajectories may differ in low float64 digits or wall time;
delivered exact words, explicit weights and native digests bind these results.
Standalone figures are `results/tree_screens.{png,pdf,svg}` and
`results/search_trajectories.{png,pdf,svg}`; machine-readable screens are in
`results/tree_screens.csv`. Offset gradients independently match finite
differences within 8.7e-9, and the offset Lambda/gap reproduce the native
incumbent (`results/offset_semantics.json`).

## Final assessment and remaining hypothesis

**A. Verified improvement:** none. The original accepted candidate is preserved
and fully reproduced; target costs 590,456,157.46 (1.20x) and 545,036,453.04
(1.30x) are not reached.

**B. Validated research progress:** legal arbitrary-tree encoders/independent
decoders, binary-lattice equivalence, an exact isotropic tensor-sensitivity
Huffman construction, and modest conditional topology/word gains are validated.
These are working research components, not a claim of a quantization breakthrough
or a verified chemistry improvement.

**C. Scoped negative result:** all tested rounded-tree, joint annealing, local
binary restructuring and shifted-grid pipelines are uncompetitive. Exact-word
span and PSD certificates explain why the preserved failing dictionaries
cannot be rescued by coefficient refitting. The resource/normalization screens
also fail independently of accuracy.

**D. Unresolved hypothesis:** broader shared-tree/exact-word cooptimization
remains plausible. The strongest next direction is to construct a well-
conditioned dictionary directly in a different binary lattice while jointly
penalizing actual LCU normalization, rather than preserving the fragile small
singular directions of these particular continuous basins. This experiment
does not establish a global precision floor, shape-wide infeasibility or
topology irrelevance. All numerical searches are terminal.
