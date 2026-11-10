# Final handoff — 9.3677221531, experiments stopped at user request

The user requested wrapping up experiments and writing a report. That report is
[EXPERIMENT_REPORT.md](EXPERIMENT_REPORT.md). All experiment jobs are terminal;
no further experiments should run without a new user request. The 14.0 target
remains unmet. The last selected change produced a verified score improvement.
Currentarchive submission.zip SHA256
9168da25a9ac35f661885edad3e208db215bf67fa6525cca0ccc5d8831318b68.
Both archive audits PASS all14members and full native default shots. Previous
9.332656488711933 archive/summary/audits are preserved in
revisions/cubic_trace_product/. DQI and chemistry selected bytes are unchanged.

NEW SELECTED INVERSE: candidates/zero_head_affine.kmb,
results/zero_head_affine_verify.json FULLPASS200540928,
cost182956.5271606425,Q90,meanT1782.000090874218,static1793,34rotations,
worsterror/tolerance0.9796714395489131at958102 (unchanged numerical phase).
KMB SHA27b14e544481db8eae7bf7f7e61bef93660ada06585cb5e398d3a18af9b50af1.
Independentfullbuilderreexecution ..._rebuilt.kmb is byteidentical; evidence
results/zero_head_affine_rebuild_check.json and
parameters/zero_head_affine_source_manifest.json. Rebuild:
source env.sh; "$EXP_PYTHON" zero_head_affine_invsqrt.py.
package_submission.py defaults now select this inverse and native JSON.

ZERO_HEAD_AFFINE_RESEARCH.md explains both changes. All3133452input mantissas
have bit26=1, range67108897..134217726. Restore the affine offset modulo2^26,
leaving that borrowed input bit fixed; combine constant3*2^23 into a one-AND
three-bit increment and parity into a CNOT. Saves12T. Clean private phase-
headroom moves use AND/CNOT forward, CNOT/HMR-AND cleanup backward, saving10T;
the first known-one move saves another1T. Normalization swaps are untouched.
Sourcechain zero_head_affine_invsqrt.py -> fixed_top_affine_invsqrt.py ->
unchanged affine_invsqrt.py. Earlier intermediate constant3_affine and
fixed_top_affine each also passed full200540928native shots but are not selected.

NEW RESEARCH: goldschmidt_row_model.py now has complete-domain numerical
models using actual Booth-row rounding and a small high-part square.13models
and component evidence in results/gold_row_checkpoint.json. Floored monomial
square f8/s10/r27/o27/sq18 passes maxrelative4.438247414e-8. No full seed/residual
quantum inverse exists. split_square_correction.py does the square using an
unchanged high source, low accumulator and signed3bit counter; all16384cases
perforward/inverse/bothfreshcompositions pass, scratch/phase too. With20bitR,
split9,shift11,37additional livequbits, uses71Q/425T per direction undercap90.
The square remains expensive before table replays/products/phase costs.

New GF4-linear tensor search gf4_tensor_seed.py/gf4_tensor_flip.cpp keeps the
CRT rank14 seed overGF4 and allows projective field scalings.64mproposals,
64chains,cap3,seed771415,restart10000: no rank13found. Exact36GF4basis-pair
checksPASS; binaryupperbound42 unchanged. Not a lowerbound.
results/gf4_tensor_flip_summary.json. Initial seed omitted CRT reduction and
was rejected before search; corrected final seed and field tensor check PASS.

FINAL COMPONENT RESULT: applying existing controlled_linear transvections to
direct_repad cuts N70 cofactor padding from 28069T to 7140T per direction.
transvection_repad.py now passes 4096 arbitrary-input cases each forward,
inverse, and both fresh composition orders, including flags, scratch and phase.
Evidence: results/transvection_repad_n70_check.json. Compiled Q875 comprises
840 public coefficient qubits and 35 public flags; the earlier resource-only
profile declared Q878 including three unused scratch qubits. This invalidates
the old cost assumption used to reject Wronskian construction. Still need H-division
adapter, cubic trace inverse-only kernel, and a TOTAL reversible Sylvester
map handling degree gaps; current generic Sylvester62,098T/dir rejects gaps.
Do not claim a whole-decoder win from these counts. Detailed older research
and selected DQI construction are preserved below.

# Prior authoritative continuation — 9.3326564887

Newest turn: PROGRESS in research components; no score change.
Completed a new24-input quadratic Boolean field-multiplication search,
which permits mixed left/right forms beyond the old bilinear tensor moves.
quadratic_field_flip.cpp:64m proposals,64chains,cap4,restart10000,seed1174819,
best42unchanged. All276quadraticcoefficients and8192random products PASS in
native kickmixGF12basis (not the upstream m.gmul towerbasis). Evidence
results/gf12_quadratic_flip_summary.json. Not an optimality proof; alljobsdone.

Added and tested
`lookahead_norm_reciprocal.py` (4096 cases/direction, Q101/T854),
`lookahead_six_quantum.py` (16384 cases/direction, Q126/T2094), and
`opening_euclid_lookahead.py` (N22,1024 cases/direction including fresh arbitrary
inverse/forward composition, Q390/T6078). All final tests PASS exact values,
clean scratch and phase, including forced singular cases. Six-moment transform
now exists as a total reversible quantum map, agreeing with the independent
Padé formulas on nonzero D2,p, including a zero first pivot.

The opening-block integration is NOT selected: initial four-degree prefix at
N70 costs28114T/Q970 vs selected20360T/Q943. It copies six words and repeats
the2094T head transform three times. Do not spend a full native decoder build
on this nonwinning prefix without first improving the representation. No full
global cofactor adapter for these new coordinates exists yet. Details and
equations in LOOKAHEAD_RESEARCH.md; evidence results/lookahead_six_checkpoint.json,
results/lookahead_six_quantum_check.json, results/opening_lookahead_check.json,
results/opening_lookahead_prefix_profile.json. All own jobs terminal.

Also added physical_factor_rank_bound.py, which proves necessary B>=10 for
ALL column gauges and dictionaries at the chemistry20mHa error threshold:
(d^T G d-.02*76)/(Tr G+.02*sqrt2926)=9.5685957243. This excludes separate
B2/B4/B8 physical factors even with joint gauge optimization, but is not a
general score bound. Larger separate factors lose the angle-table resource
screen. NORMALIZATION_BOUND.md and results/physical_factor_rank_bound.json.
Selected archive, selected circuit hashes and selected DQI source manifest
were rechecked unchanged. Goal remains active, with no external blocker.

The14.0goal is ACTIVE and NOTmet. This turn made progress: a genuinelydifferent
GF16cubic-trace DQI phase kernel is now fully native verified, packaged, audited
and independently rebuilt. All own research processes are TERMINAL. No external
blocker; no subagents authorized. Keep working toward14.0; do not mark complete.

**submission.zip** SHA256
81c3ea50eb2a3e11b80744cba015915e9db12ba324c01783817ecba153f62d20;
**overall_score9.332656488711933**. Both current archive audits PASS14members.
Previous9.310120488836658archive/summary/audits preserved in
revisions/compact_product_slim/. Inverse/chemistry unchanged, stillselected
candidates/affine_best.kmb andcandidates/norm_cycle2_signed with full evidence.

SelectedDQI candidates/cubic_trace_product_n70_m4095.kmb;
results/cubic_trace_product_n70_m4095_verify.json FULLPASS102400,
value_errors0,phase_errors380,Q961,meanT1019865.856484375,staticT1019954,
R1.6968460732519504,cost1598204229384.9265.
CanonicalSHA4989da5cb633529e8851b9aa0c8c52abeb8990423a3a9463a70f673661a909c0.
`package_submission.py` defaults now selectthisDQIanditsnativeJSON.
`build_cubic_compact_decoder.py` defaultfullyrebuildsselectedN70/M4095.
Independentfullreexecution to..._rebuilt.kmb hasSAMESHA;
results/cubic_trace_product_rebuild_check.json PASS. Source manifest
parameters/cubic_trace_product_n70_m4095_source_manifest.json references old
unchanged compact-coremanifest plusnewkernel/modules. CoreKMBtapeunchanged.

NEW CUBIC_TRACE_DQI_RESEARCH.md hasfullmath/tests/resources/reproduction.
Corefiles cubic_trace_field.py,cubic_inverse_product_kernel.py. GF4096as
GF16[t]/(t^3-2), Q=a^2+2bc,N=aQ+2b^3+4c^3,
Trace12(1/x)=Trace4(Q/N). Cubes4AND each, inverseGF16fiveAND (BoyarPeralta
innercore with independentlyfoundbasis maps), Norm/Q26AND,gatedparity15,
inverseproduct42+46=88T,Q78incl36public vsold92T/Q80.
Runtimeclassicalmask maps totracepairingk; classicallycomputeinvk andlinearly
rescalequantumdivisorby1/kusingclassicallycontrolledCNOTs, restoringallinputs.
ClassicalANDadapterusesBIT_STOREconditions, notquantumCCX. Zero maskskip.
Exhaustive61440kernelshotsPASSall4096divisors×12maskbasisvectors,all4096masks,
zero divisors/nonroots/source/mask/scratch/phase. FieldidentityPASS4096.
SamecoreN14/m63sameexplicitseed77291old/new both4phaseerrors,0values;
old45010T,new44756.324mean. No generic-degree-gap accuracyimprovementclaimed.

REJECTED narrowerjointvariant: joint_cubic_inverse_product_kernel.py combines
17quadraticfront-endANDsinto16 viaquadraticflip/plussearch. Frozen
parameters/cubic_joint_quadratic_16.txt, all256truthrowschecked, kernel87T/Q78,
full61440kernelPASS. candidates/joint_cubic_trace_product_n70_m4095.kmb
FULLnativePASS102400,0values,406phaseerrors,Q961,meanT1015772.10875,
R1.7223319602461098,cost1615697003965.9805 — WORSE selected. KMBhash
7783891a5e4eb1fd3a20cbfedcf8d488e1d0e00e9d69699a0107994a7ee048b7.
Nativeinputsdependoncircuithash; twolegitimategatecountvariantshavedifferent
finite-sampleR. Don'ttryresource-equivalenthash/seedtuning. Furthersearchdidn't
reach15AND. A32mproposalGF12tensorflip/plussearchdidn'treduce42multiplier;
results/gf12_tensor_flip_summary.json, retainedcoverPASS144basis pairs.

NEW LOOKAHEAD_RESEARCH.md: lookahead_pair_quantum.py totalreversibleextension
of(u,v)=((d*b+a^2)/D2,(e*b+d*a)/D2),D2=d^2+ea; validPademeaningD2nonzero,
handlese0usingdpivot. SingularD2and e=d=0 haveexplicitreversibleextension.
Forward/inverse/freshcomposition16384shotsPASS, Q102incl48public,T654/dir.
lookahead_moment_map.py sixmomentclosedformsP2/P3/pinv allowfirstpivotzero;
8994of12000casesvalid,minorsD2andD3nonzero, fullEuclid/inverse/determinantPASS.
NOTa completequantumhead/globalcofactoradapteryet; nostatedwholedecodersaving.
Terminal/directwhole-seriesblocksalsonew: terminal_euclid_block.pyN22b4
PASS512fwd/inv/arbitrarycompositionQ295/T8840; series_euclid_block.pyN22k7b2
PASS512Q293/T7024. N70all78latefullseriesprofilesLOSEselectedserialToffolis;
Q902butearly961Qpeakunchanged. No fullcandidatefromthesenonwinningcounts.

Earlierthisturn:R28B64C32Adam12000ends17.020772619mHaUNQUANTIZED;
quantizedfixedPSDC64refitsq7/q8remain229.739/131.728mHainvalid.
R24B64C32LBFGS1000ends22.206614165mHa,weightnorm625.47invalid;
fullcovarianceAdam24.434mHa invalid. AllTERMINAL.GROUPED_FACTOR_RESEARCH.md.
CorrectedcubicNUMERICinversephasecandidatecubic_phase_corrected_q94.kmb
PASS131072quantumshotsbutQ94/T2458LOSselected90Q/T1816; excluded.
PHASE_AND_GROUPED_RESEARCH.md recordsalloldBoolean/spooky/Walsh/Goldmodels,
newcubicfloatfitsandfullpairdictionarynormalizationcertificate62.32194.

Readthisentryfirst. Older entries below can contradictlatestselectedpaths and
livejobs. Preserveverifiedsources, frozenparametersandallarchives; continued
radicalresearchremainsnecessary. No hostverifierchanged; workonly4/codex.

# Latest continuation — goal active, no live jobs

Verified archive still **9.310120488836658**, SHA256
6b7eb993f3076695d772b9c769f49bc99abc044f03abf8a50a21a146024931ae.
Goal14.0notmet; noexternalblocker, no subagentsauthorized. All own research
jobs are TERMINAL. Previous turn made research progress but no accepted gain.
See **PHASE_AND_GROUPED_RESEARCH.md** and the appendedR24/R28 section in
**GROUPED_FACTOR_RESEARCH.md** for the authoritative new results.

R24/B64/C32 Adam12000best24.645975575mHa; L-BFGS1000finished
22.206614165mHa, invalid, weightnorm625.47. Covariance/dictionary alternating
4000finished24.433965819mHa. R28/B64/C32 Adam12000finished
**17.020772619mHa unquantized**, weightnorm44.68. However fixed quantized
R28dictionary PSDrefits at7/8bits give229.738643/131.728289mHa, invalid.
R24quantized7/8/9bits withfullC64PSDrefits give279.057878/171.904133/
116.963129mHa. Exact snapshots retained inparameters/grouped*_current_q*.
No grouped native candidate. Don't mistake17mHa float model for acceptance.

New cubic fittedfrac6/fine7table27/19/11/5modelpasses all3133452inputs,
worst5.64446021034e-8. A new compressed intermediate cleanup schedule yields
`candidates/cubic_phase_corrected_q94.kmb`, **quantumPASS131072shots**,
inputs/scratch/exactmodelphase, Q94/T2458. It LOSES selectedinverse90Q/T1816.
No fullnative run needed for that resource loser. Defaultbuilder/check now use
safe tableparameters/cubic6f7_w27_19_11_5_safe.txt. Earlyuncorrectedresource
screens are REJECTED: affineR0overflow. SafeR0usesfineq7and−2^17bias,
all-x range16318..8257535. `cubic_phase_invsqrt.py` undoeshighcarryintegration
before measuring a consumed source's lookup-only highbits, restores/reintegrates
forinverse; this reusable idea passesquantum butdoesn'tprovidecurrentcostwin.

DirectWalsh normalizedmantissa synthesis needs2048terms toreach3.7141e-8,
before parity/exponentmodulation; no winningcircuit. Exact2^21-pointprobe:
results/walsh_mantissa_21_probe.json.

New fixedfullsparsepairdictionary dualcertificate:
results/full_pair_dictionary_bound.json, sourcefixed_dictionary_bound.py.
All76onehots+2850pairs, anyrank/columngaugewithin20mHa ⇒ numericalanalytic
boundlambda_eff>=62.321940. Testedbrot3/C276,300,350heuristic-onlycosts
1.824b/1.976b/2.293b vsselected708.547m. This boundsTHATdictionaryonly,
notallshapesoroverallscore. Featureinverseindependentroundtrip1.11e-16.

Boolean/spookycompilerandGoldschmidtmodels from prior continuation are now
recorded inPHASE_AND_GROUPED_RESEARCH.md. Both spookycompiler unitchecksPASS,
butactual795ANDpolygraph needs too manyqubits/recomputations. Carry-save863AND
worse; ABC1236ANDworse. Goldschmidtnumericall-x modelsPASSseveralparameter
setsbutNOquantumresourcebuilder. Don't countwhole-productroundingasBoothrows.

Preserve the selected affineinverse/compact_product_slimDQI/norm_cycle2_signed
chemistry, their full native results, frozen parameters and both archive audits.
No selected files or hostverifiers changed. Older continuation entries below
are historical and can contradict the latest state; follow this entry first.

# Active objective

NEWEST AUTHORITATIVE CHECKPOINT: **9.310120488836658**, ZIP SHA256
6b7eb993f3076695d772b9c769f49bc99abc044f03abf8a50a21a146024931ae.
Selected DQI is now `candidates/compact_product_slim_n70_m4095.kmb`, full
nativePASS102400, zero value errors,365phase errors, Q961/T1036334,
repetition1.6820375683380908, score1609838144370.2598, hash
781e1c04fe1bca206adbf9f105ef537b31329dfbe53225c6bcbc4a07e3c2b2b4.
The new complete architecture barely beats inherited1612737419463.0198;
native phase pricing erased most raw resource savings. Both ZIP audits PASS.
Inverse and chemistry remain selected affine/norm_cycle2_signed unchanged.
package_submission.py defaults now select this DQI and its native result.
All jobs from this turn are TERMINAL; live-job notes below are historical.
Previous9.304538archive preserved inrevisions/norm_cycle2/.

Selected full DQI rebuild:
```bash
source env.sh
"$EXP_PYTHON" compact_product_decoder.py --n 70 --block 3 --serial-after 22 --m 4095 \
  --linear-padding --virtual-zero --virtual-lead --low-inverse --lazy-flags \
  --output candidates/compact_product_slim_n70_m4095.kmb
```
Full source manifest parameters/compact_product_slim_source_manifest.json.
Build uses upstream Apache2field primitives read-only via../../3/import_published;
extra local files controlled_linear.py andlow_workspace_inverse.py. The latter
requires E.J.recording=False, as in this standalone builder.

Important width discoveries: lowinverse removes norm-chain live ranges before
inverse output allocation, erases Y usingNorm(output), and measures output
before cleanup recomputation; it retains102T and passesall4096fieldvalues.
Full cofactor scaling still needs12inverse+12product wires, so bodypeak24,
instead30 oldpeak. Completionflags are now allocated lazily and released as
the inverse proceeds. During the scan, pinv itself is the priorcofactor's
explicit leading coefficient (no copy); B's monic leading coefficient is
virtual and evaluated by linear XORs of unary completion flags, saving12more.
The common-x rotations act onB's35storedcoefficients modulo35 (ratherthan36),
and the virtualmonic term supplies x^L. All matchedN14/m63 tests haveexactly
the same2phase errors on4096fixedshots, zero value errors.

The original fullversion scored1678643393596.0247 and was not selected.
The new slim source haschanged sinceoriginal; preserveoldcandidate+nativeJSON
as evidence ratherthan expecting currentbuilderdefaults to reproduceoriginal.
The selectedcandidate hasnow been independently rebuilt and byte-hash compared:
results/compact_product_slim_rebuild_check.json PASS samecanonicalSHA256.
All jobs areTERMINAL. The archive andbothaudits remaincurrent9.310120488836658.

NEW radical inverse research front end `boolean_polynomial_phase.py` builds a
combinationalXOR/AND graph for the selectedHorner, affinerestoration andphase
shift. Numericoutputexactlymatches independentinteger poly.model on32857
random/boundaryx. Input74bits=normalizedx22+flags5+coefficients23/17/7.795AND
nodes,maxANDdepth45,maxancestors668. This is NOTaquantumcandidate andhasno
resourcescore. ExistingNetwork.emit retainsparentsuntilchildrenerase andwill
needfar too manyqubits; a genuinelydifferent measurement/spooky pebbling or
globalBooleanoptimization would be needed. Do notassume795ANDs is a win:
lookup~510T andnormalization/rotation stillcost, coefficientinputs74qubits are
live. A wholegraphoriginalx22-input synthesis could shareconstantROMlogic but
ABC earlierROM mappingswerepoor(~4kAND), so not an automatic gain. No new
Booleanquantumemitter orABCwholearithmetic joblaunched. GoalACTIVE14notmet.
Continue larger algorithmresearch; nativepipeline andsourcebytesaresecured.
The older continuation entries below are historical.

LATEST WORK (2026-10-24): a completely assembled compact inverse-product DQI
decoder now has full native evidence. First version
`candidates/compact_product_n70_m4095.kmb` PASS102400, zero value errors,
317 phase errors, Q985/T1058846, repetition1.6340042594,
score1678643393596.0247. This is4.09% WORSE than selected DQI1612737419463.0198,
so the archive remains9.304538097994794. Native hash
07974777f89e866829ea51d7897eb61e8e886d27026231886180986de88fea0c.

LIVE FOLLOWUP: `compact_product_slim_n70_m4095` build then native verification
in one `set -e` shell, session65346, shell PID32193, builder PID32197. N70,
block3, serial_after22, m4095, flags --linear-padding --virtual-zero
--virtual-lead --low-inverse --lazy-flags. Logs/results use that stem under
results; candidate undercandidates. Build takes~5min, then native102400 shots,
4workers. CHECK ACTUAL PROCESS AND LOGS before treating a tool observation
timeout as terminal. It was serializing after scanposition3840 at~4min.
Smaller N14/m63 matched4096-shot tests all pass zero value errors and exactly
the same2 phase errors as first version; latest stemcompact_product_lazy_n14_m63.

The new architecture consumes syndromes into compact cofactors, scales prior
cofactor bypinv, unpadscommon x factors by4812T coefficient rotations and uses
`error=x^N/(pinv*Bprev(x)*B'(x))` for nativeS1-based syndromes. No Wronskian
or squarefree padding in complete decoder. The inverse-product kernel92T,
4096 field testsPASS. Full builder `compact_product_decoder.py`; all local.

NEW reusable gadgets: `controlled_linear.py` factors an invertible affine
GF2 map into <=n rank-one transvections, each oneToffoli after a free Clifford
basis change. It combines fake-undo withmonomialpadding; random controlled
maps/inverse/composition PASS512 atdimensions2..64. `low_workspace_inverse.py`
uses Y=Norm(x^-1) to clean thetowerintermediates early and measures inverse
output before cleanup-chain allocation. Same102T; exhaustive4096fieldvalues,
randomXORtargets,scratch andphasePASS atQ45 including24publicinput/outputbits.
It assertsE.J.recording isfalse; do NOT drop it into journaledbaselinecode.
`streamed_unlookup.py` stores onlyone ANFbank andDFSstreams theother: n11
needs31workspaceinstead83,same83T; exhaustiveaddress testsn2..11PASS.
All these research tests are outside currentarchive.

Mixed monomial core N70block4serialfrom21 PASS512 forward/inverse/arbitrary
composition Q1004/T345260 including direct squarefree repadding. Direct+
suffix-only core413318T. Earlier thin=True cofactorband pass is REJECTED:
arbitrary-input composition changes values andphase (results/thin rejected).
SeeDEFERRED_COMPLETION_RESEARCH.md. Exact stageDP in
results/adaptive_euclid_profiles_n70.json: b3thenserial core362780T including
28069T repad; b4thenserial345260T. Full decoder omits that repad.

Grouped chemistry fits are ALL TERMINAL andinvalid; seeGROUPED_FACTOR_RESEARCH.md.
R10B64C35 Adam10000 best130.352mHa; R16B64C24 best63.840mHa, L-BFGS1000
62.566mHa; R16B64C64 fullrankprojection+L-BFGS1000 best55.778mHa.
No angle quantization or native chemistry build from these fits. Simple sparse
pair dictionary719coordinates error505.9mHa;2500coordinates38.81mHa, invalid.
guard28 affineinverse resource screen199530 bestestimated, loses185026; no
native inverse run. Earlier frozen production inverse/chemistry untouched.

`package_submission.py` now accepts --dqi and --dqi-result, default stillbaseline.
If slim fullnative beatsselected DQI, package it and independentlyaudit; otherwise
preservearchive. Current9.304538 archive newly backedup inrevisions/norm_cycle2/.
GoalACTIVE14.0notmet, noexternalblocker orsubagentsauthorized. Continue working.
The older continuation entries below are historical.

LATEST AUTHORITATIVE CHECKPOINT: **9.304538097994794**, ZIP SHA256
51a94509b0b3b35b191b1ca9f74ecdc89e208a69044743065ccd37e51703d024.
Chemistry`candidates/norm_cycle2_signed`, full nativePASS50880,
cost708547388.9553173,lambda48.59417497035283,Q478,T30504. Canonicaldigest
ebf8e357a1b693fc8589d5c48a5d6fa10343dfa5944724d8f7a77b8edc2157f7.
All11oracles alsoPASS2048 random-targetshots each; independent exporter rebuild
and both ZIPaudits PASS14members. Previous9.2870497 saved inrevisions/norm_cycle1/.
package_submission.py defaults to this candidate. Inverse/DQI unchanged.

NEW **DEFERRED_COMPLETION_RESEARCH.md** records a corrected decoder completion
architecture: monomial padding during the core, squarefree common-factor
padding at the end. N22 forward/inverse/arbitrary-input composition PASS512,
Q404/T51675perdirection. N70 classicalPASS3486cases/62665roots,14generic gaps
skipped. FullN70block4 forwardcount1004Q/451699T, downfrom1006Q/608137T;
still excludesHmap,Hdivision,phase scan,repetition. No fullDQI gain claimed.
Initial Fibonacci-undo attempt was REJECTED because completion insideblocks
revivedremainders; see rejectedJSON. Do not reuse that incorrect version.

NEW **NORMALIZATION_BOUND.md** proves a global native lambda_eff lowerbound
46.78617959 from fixedshifts and thegapfloor. Selected48.59417 has under3.87%
normalization headroom. Larger gains need gate/width or factorization resources.
This is not a bound on the overall score. Upper-gap-only convex probe was worse.

All experiment jobs from this turn are TERMINAL. Error-only700 process28022
finished20.076379370mHa (stillinvalid). Error-only500 process29414/session70548
finished38.114424354mHa (invalid); frozen seedparameters/lll500_q8_seed.npz SHA
910039f64a8e4eedb93c73278b9ba8f0e44fb12069b68d15162f165815782757.
Normalizationcycle2 process29022 was deliberately stopped afterstep79 to refit
factors and offsets using its improved physicalerrorroom. Frozenrefitseed SHA
d230c24f7b51d73e558aff5f76d25a94f0e5990cdc5b7ebbf1da2fa97dca13b7.
The8000steprefit,onebodyquantization andfullnativepipelinefinishedsuccessfully.
Frozenfinalfit`parameters/norm_cycle2_final_fit.npz` SHA
10688688f0660fd1dc9580ba3541e7db2b6278622cfce48ed18a3c501b11ced6.

Additionalinversefits are researchonly: affine8_fitted_invsqrt.py shrinks
curvature widths with frac8/fine9tables, all-x modelsPASS but resource screens
lose. affine8f9_17_6 best86Q/2002T estimate192812. Expandedkeep3..5/product17orNone
screenbestaffine8_17_7 88Q/1962T estimate193776. Neither beats185026verified.
Fine9fits(27,17,5),(27,17,6),(27,18,5),(27,16,6)allmodelpass; nofullnativejob.
guard28 fitsw1=18,w2=6 or7 modelpass(391/294leaves), not circuit-screened.
Selectedinversebuilderandarchiveunchanged;generatedaffine_build_source.py is scratch.

Goal ACTIVE,14.0 notmet, noexternalblocker orsubagents authorized. Keep working.
All older checkpoint/job observations below are historical and superseded.

NEWEST CHECKPOINT: **9.287049724889583**, archive SHA256
d9409bfab8dc74725bf9963de29d127626c3490732601a815b00af06faae3709.
Selected chemistry`candidates/norm_lattice719_signed`, nativePASS50880,
cost712557711.9823189,lambda48.869214215291805,Q478,T30504. Both archive
audits PASS14members. Previous9.2854358 preserved inrevisions/signed_gauge/.
Read **NORM_DICTIONARY_RESEARCH.md** for the new normalization-aware dictionary
search, frozen parameters, evidence, failed column-exchange probes and live jobs.
Current new jobPID29022/session91305 (`norm_lattice719_cycle2`) usesall719rows,
maxstep32,512trials,16moves/iteration,200steps and gap projection. Atstep3 best
model-valid lambda48.8135560, not nativeverified. Error-only700 searchPID28022/
session31184 remains live, step75 best20.079770mHa, stillinvalid. All other new
normalization refits and first-cycle searches are terminal. Goal ACTIVE below14.
The older checkpoint/job observations below are historical and superseded.

CURRENT CHECKPOINT supersedes older numbers below: **9.285435820940434**,
submission.zip SHA256
5de660862765b7b4a2cefab452e1a731180f6c7a49e105fa30440feaca72e7ab.
Chemistry selected`candidates/signed_gauge7_prefix8`, nativePASS50880,
cost712929326.0737408, lambda48.89470055602427,Q478,T30504. Canonicaldigest
9db8ca09bffcd4ea9293e4b2eb3ef0b9c1962f26050300e4fd625d4e8f276883.
All11oracles also pass2048random-target shots each and independently rebuild
byte-for-byte; all NPZ arrays rebuild exactly. Archive audit and independent
archive audit both PASS all14members. Old9.2439070508 saved inrevisions/jacobi/.
Inverse and DQI unchanged. package_submission.py now defaults to this chemistry.
See **ANGLE_GAUGE_RESEARCH.md** for the exact encoding gauge, reconstruction
commands and failed inverse precision screens. None of the new inverse fits
beats the selected affine circuit; no full native inverse job was launched.

Chemistry jobs: process26361/session10923 was intentionally stopped after
step172(best20.119658mHa) to use a better completed branch. Paired escape run
session29690 is FINISHED after26steps at **20.114304597mHa**. It used frozen
parameters/lll700_pair_escape_seed.npz SHA
2cbc6bd9b0461f5374b76a3cd997e885ab869d92e982aca56c2c95f381348384.
It searched96single rows plus24paired150-angle blocks, maxstep32,padding.0001,
128trials,8accepted moves periteration,100steps. Output/log are
candidates/lll700_pair_escape_q8_c350.npz andresults/lll700_pair_escape_q8_c350.log.
The current live continuation is **process28022/session31184**, all700rows,
maxstep32,padding.0001,256trials,32moves/iteration,300steps, initialized from
the terminal paired best. Output/log stem`lll700_pair_polish_q8_c350`.
Inspect its actual process/log before deciding terminal. All are research fits,
not valid new submissions. Do not stop/restart solely on an observation timeout.

Latest continuation: **BLOCKED_EUCLID_RESEARCH.md** contains the new decoder
architecture and exact evidence. Quantum tests now pass for generic single
steps, quotient blocks, a full armed N20 sequence, and a completion-padded N22
sequence. The quadratic padding's classical N70 check passes34850cases and
626265roots, zero artificial collisions,150generic degree gaps skipped.
The padded N22 sequence uses406Q/65605T per direction, with arbitrary-gap-input
composition restoring every input, flag, scratch and phase. However the whole
N70 profile is1006Q/608137T **per direction before H and scan**, so the current
completion schedule is too expensive. All results are research, not a new
submission. Do not build a huge whole decoder on the unpadded arithmetic count.

The fixed-degree H map can reuse the existing Sylvester prototype after
reversing the monic B/Bprev pair: C=E²+zO²,D=U²+zV² gives Hrev=EU+zOV,
E0=U0=1. Details and zero-pivot audit requirements are in that handoff. The
source psi_forward already has epsilon degree-defect repair; inspect it before
assuming every gap requires a new generator. Payload/scaling modifications
still need non-generic value-restoration tests. No H adapter has been built.

The optimized normalized reciprocal passes512shots atn35,Q423/T11166 andn70,
Q843/T32622. The Möbius map and continued-fraction quotient extraction are
implemented. SERIES_RESEARCH.md and BLOCKED_EUCLID_RESEARCH.md supersede their
older resource numbers in the historical text below.

Chemistry search is still live: process26361, session10923,
`results/lll700_multi_q8_c350.log`, best20.144238mHa atstep44 when last recorded.
It accepts up to32individually nonlinear-checked moves per iteration; saves
best actual-reference errors atomically. Input was the intentionally stopped
long trust8 continuation. Output`candidates/lll700_multi_q8_c350.npz` is a
research fit, still failing20mHa. Inspect process/log for newer status. Previous
trust8-100 run finished20.214431mHa. No other research jobs remained live after
the whole-profile and larger classical checks completed.

Current archive is unchanged at9.243907050802967. Goal remains ACTIVE, no
subagents authorized, no blocker. Keep working; do not pause/complete below14.

Build a submission with overall score >=14.0. Work only in
initial_experiments/4/codex. Do not modify sibling agents' directories or prior
experiments. No subagents have been authorized. Goal remains ACTIVE; do not mark
complete at the current score.

The latest user explicitly set a fresh ACTIVE goal: continue until14.0 is met
or interrupted, investigate radical changes, and do not give up. They report
someone has reached16.6. get_goal confirms ACTIVE. An optional async question
asks for a public notebook or per-task breakdown; work does not depend on it.
This turn made progress: verified an exact rotation solver, checked fixed-Gram
bounds, packaged a small full-native chemistry improvement, and implemented a
new LLL/Babai integer-angle search that accepts coordinated jumps beyond ternary
moves. See RADICAL_SEARCH.md for live-job handles and the next research steps.
Do not treat old resource estimates as impossibility proofs. No blocker.

# Previous checkpoint details (superseded by signed-gauge archive above)

submission.summary.json: overall9.243907050802967 (target not reached).
submission.zip SHA256
46b101eb757311e7a878682c8428bc1c511c7b77652a063f89b47f4647bdc091.
Archive audit results/archive_audit.json PASS binds all14members and native hashes.
The previous9.2438728734 archive/summary/audit are in revisions/affine/,
zipSHAc01c862078ad658999cd37cc17f4f3d5e2d6f2fc93049b88e853d5fb0737e13b.

Inverse candidates/affine_best.kmb, results/affine_verify.json PASS200540928,
score185026.55749858925, Q90, meanT1805.0004279625155, maxT1816, rotations34,
worst_relative_tolerance .97967143955, circuitSHA
89dc5547e068941b138fab5a1d68d33877ab43230a1b180043a84282852d34c1.
Build with build_affine_invsqrt.py --output candidates/affine_best.kmb.
The rebuild candidates/affine_rebuilt.kmb matched byte-for-byte before packaging.
The previous9.2155702694 archive/summary/audit are in revisions/delayed/,
zipSHA712774c8d4ff564d4764d6d7c344e6b3bd0fcb20eb465fea1331b4b22230b2e6.
Table candidates/margin_curve8_t5.9e-08.npy copied
from experiment3 unchanged. Sign-integrate-only precursor also fullPASS,
score190335.92375, results/sign_integrate_verify.json.

The prior 9.1631467535 archive, summary, and audit are preserved under
revisions/sign_anf/ (zip SHA68bb0f6b...a2e07c6). build_invsqrt.py still rebuilds
that older sign_anf circuit byte-for-byte: default-false research options did
not alter it. Its rebuild hash was explicitly checked after the source edits.
The new delayed candidate was also rebuilt byte-for-byte before packaging.

Latest inverse change: AFFINE_RESEARCH.md. Store C0+center/4-A(parity) in23bits;
A=(961/1024,1217/1024). Restore the baseline at the phase boundary, borrowing
the normalized input's top bit as the phase word's bit26. A cheap7bit address
subtraction plus two constant increments and one parity increment restores it.
Delay two residual bits21..22, keep4curvebits, Qlimit90, defer3headroombits until
after curvature release. affine_invsqrt.py uses guarded source replacements in
the unchanged binary_centered_invsqrt.build; all generated native gates are
covered by the full host verification. Source changes fail fast if hooks move.

The preceding delayed-only inverse loaded bits22..25 of C0 after the outer
Booth accumulation, immediately before its carry integration. Borrow those
zero bits together with the ten output headroom bits during arithmetic. After
reverse integration, HMR the four C0 bits and reuse them as scratch. Append
their classical outcomes to the final ANF unlookup, whose table includes their
original coefficient values. This has an independent 1024-shot/random-table
deferred-erasure check and the full unchanged host verifier PASS above.

Other new research is excluded from the archive:
- COMPACT_EUCLID_RESEARCH.md: classical nonmonic and monic Euclid state maps
  pass full-Euclid and inverse comparisons. Monic uses70fieldwords, but serial
  field-product estimate411516T/direction is too high. First quantum monic step
  PASS512 atN70,Q875,T3012/direction, including armed zero pivots. Ungated
  inverse-mask phase primitive PASS57344,all4096sources and12maskbasisvectors,
  Q34,T36. These are research primitives, not a new decoder.
- New fine-coordinate quadratic tables(27,17,7)/(27,17,8) both pass all-x
  classical models. Coarse-coordinate(27,19,10)/(27,19,11) fits also pass their
  C++ model. Actual compacted resource screen finds all worse than checkpoint;
  results/quadratic_fit_resources.json is authoritative, not the earlier raw
  wire-index fine8 screen. No such fits entered the archive.
- SERIES_RESEARCH.md: 3-scratch-qubit normalized reciprocal (length35 T12804,
  length70 T37494), each quantum PASS512; serial-CF inversion sum517818 per
  direction is too expensive. Generic in-place Bezout and Sylvester transforms
  also PASS512, but degree-gap cases are excluded and their cost does not win.
- LATTICE_RESEARCH.md: new C++ QUBO annealer and coordinated hard-angle search.
  Best700-vector modelerror20.253mHa still fails20mHa. Adding15atoms reaches
  19.981mHa but factor/offset refit has much worse normalization; no win.
- phase_headroom.py: input-bit borrowing for phase shifts PASS2048 across all22
  normalization shifts. space_increment.py and signed_counter.py have clean
  native checks. Resource searches produced no valid improvement. Small-carry
  candidates still overflow; do not select by gate counts alone.
- monotone_lookup.py: binary-search preparation PASS1024, but delayed-C0 variant
  has1860staticT vs selected1835, same90Q. Its smoke PASS is not a full result.
- Boolean synthesis: vendor/abc revisionb63f68252df1bb0e37d52c76a28c1882d7a67703,
  license retained. Mapped circuits are much more expensive than the existing
  lookup. ESOP all-output minimization hit an upstream assertion; log retained.
- Cubic classical models with widths(28,20,12,5), (28,20,13,6), (29,21,13,6)
  pass the relative-error threshold in a full classical all-x model, but their
  preliminary circuit resources do not fit100Q. No complete winning cubic.

The experiments listed above completed in earlier turns. Current LLL research
may still be running; inspect the actual handles recorded in RADICAL_SEARCH.md.
Goal remains ACTIVE, with no blocker.

Chemistry candidates/jacobi_l1, results/jacobi_l1_verify.json
PASS50880, score722581182.7790365, lambda48.89470055602427, Q478,
T30916.972543238993, beta8,J7,K11,maxANGLEwidth128. All angle circuits match
selectswap7_anf byte-for-byte. The new NPZ rotates native W and H2_I, preserving
the combined coefficient Gram matrix to7.2e-13. All arrays independently rebuilt
bit-for-bit as parameters/jacobi_native_rebuilt.npz. Hostdigest
0a04a6801d388eeec3078d410e8cd0c4f60759f4cd2995bd5e5ac8e15bf88d50.
Rebuild: build_selectswap_chem.py --out candidates/selectswap7_anf --monomial;
jacobi_l1.py --passes 20 --pairs 16000, then copy
the original angle circuits and parameters/jacobi_native.npz to jacobi_l1/.
The preceding selectswap7_anf result remains valid at722589197.6207536.
All11oracles additionally pass random-target XOR/phase/cleanliness checks in
results/selectswap_random_targets.json (1024shots each). Precursor with builtin
unlookup candidates/selectswap7 fullPASS cost724645921.1053.

DQI unchanged from baseline zip, inherited fullPASS102400 result
baseline/dqi_verify.json; cost1612737419463.0198, Q1429,T645755,8phaseerrors
priced through repetition_factor1.223013998. The archive member hash is checked
against its circuit_sha2566976734728a4263d2233d1d7f1cab3fb3a61bab452b734ca067e7c47aed835ac.
Do not mistake this accepted decoder for a zero-error one.

# Implemented changes

1. Signed carry ±counter can be added by complementing high accumulator bits
under coordinate sign, incrementing by counter, and undoing the complement.
This removes a temporary sign-extension AND in each of4integrations.
2. monomial_unlookup.py: HMRtable, ANF coefficient transform, banks of all
monomials in address halves, classically controlled CZ between banks. For9bits,
37ANDs vs builtin41. Banks reuse the now-zero measured target wires. Unitchecks
n2..10 alladdresses PASS. Inverse and chemistry use this.
3. Chemistry two-entry lookup: load F(2h) into output and D(h)=F(2h)^F(2h+1)
into scratch; CCX(addresslow,D,output) selects the XOR payload for arbitrary z.
EraseD usinghiaddress lookup. Ordinary QROM onh useslowbits6/highprefix3 tofit
J7's56+56wordbits+10selection+6workspace=128. Persistent storage56vs104 saves48
scored qubits though total Toffolis increase. The later column rotation changes
the factorization arrays while preserving the represented Hamiltonian.

# Failed / bounded explorations

optimize_workspace.py broad q78..100, curvekeep0/2/4/6, productlimit8/12/16/20,
splits(17,8),(18,9),(19,10),(20,11). Predictedbest187327 fromsplit17,8 is INVALID:
results/workspace_best_verify.json FAIL x2phaseerror .001953. Carry wraps more
than one bit can represent. Preserve split18,9 unless a new exact arithmetic
model verifies smaller accumulators. Search script estimates do not imply PASS.

parity_phase optional experimental branch in binary_centered_invsqrt.py computes
singleparity coefficient tables and adds sqrt(1/2)-1 phase under flag0. Linear
f9/f10 cost profiles did not beat selected quadratic circuit (f9 also too
inaccurate without more lookup refinement); no claim of a passing finalvariant.
Do not accidentally enable it in final builder.

# Larger remaining opportunity / evidence

Current score still needs product-of-cost reduction(14/9.1631)^3 ~=3.57.
Microoptimizations alone are insufficient. DQI current T645755 comprises
293503 Toffolis in syndrome-consuming BM tape plus ~352k in omega/Forney work.
Forney recursive inverse-parity kernel costs80/point over4095points. Existing
experiment3 vendor/uncompute already implements the half-degree Euclidean
reciprocal polynomial route with23/point, but its width~1837 and full adapted
score~1.829e12 were worse thanstream1429's1.6127e12. Main obstacle: original
route zeroes35syndrome words after omega to hold Bezout coefficients; stream
route borrows those35input words as live BMbackup, leaving only49free wires.
Do not naïvely add432cofactor wires and assume lowerTmeanslowerQ²Tscore.

Read-only primary source consulted:
https://andigu.github.io/projects/uncompute/dqi-reed-solomon.html
and local vendor circuit.py build() around1910..2098. It explains sigma=E²+zO²,
AE+BO=1 => inverse sigma' =B²+zA², reusing the zero omega tail for cofactors.
Could investigate a reversible in-place interpolation/magnitude map or a way to
release the BMbackup before the low-T scan. No such improvement implemented yet.

Other prospective work: fused or differently structured inverse arithmetic to
beat ~1770T/Q94; a substantially smaller accurate8bit chemistry dictionary
(B<=512) would lower both inner preparation and angle costs. Prior experiment3
contains extensive failed angle/lattice fits; read its report before repeating.
No new chemistry optimizer was launched in this turn.

# Environment / reproducibility

source env.sh from thisdirectory sets EXP_PYTHON=/tmp/quantumai-venv.wOrZGC/bin/python,
LD_LIBRARY_PATH torch2.13/lib, PYTHONPATH local+host+legacy3/2. Allnewwriteslocal.
Native scipy/numpy/kickmix available; torch not in thisvenv. CUDA4090 wasfree at
inspection, separate /mnt/data01/miniconda3/envs/torch2.13/bin/python providesGPU.
No external messages sent. RootREADME/gitignore changes predated this work;
do not revert them. Siblingagy/claude directories permissiondenied intentionally.

# Latest DQI research (after packaging)

See WRONSKIAN_RESEARCH.md. Newly implemented/tested research primitives:
- inverse_mask_kernel.py: root-gated inverse parity with runtime CLASSICAL mask,
  **50T**,56Qincluding24quantuminputs; all4096divisors PASS, clean phase/scratch.
- template_permutation.py: recursive rearrangeable37-word network, selectorL6bits.
  **2738T,24extraQ,11layers**, quantum512shots PASS;630random classicalpermutations
  N2..64 PASS. Fits existing49free wires. Inverse isreverse layer order.
- capture_template.py reran the existing BM generator hook-only, recording
  template37words/counter8bits/unborrowedS34. Every OLD hook field exactlymatches
  experiment3 stream10_offset3_select_noscale_compact hook, including allwireids,
  5885230sweepindex,15032464tapeops,293503totalT,1428peak. Newhook underparameters.
  No tapewritten; reuse prior3tape read-only. templatefirst3words904..915,
  935..946,947..958;slots3..36borrowS0..33. Counter852..859. S34unborrowed408..419.
- check_template_prefix.py +translate_prefix.cpp compile actualfirst5885230
  tapeops into candidates/template_prefix.kmx (~200MB). **PASS1024shots** comparing
  ALL70mutatedSwords,36sigmawords,37templatewords,andcounteragainstclassical
  probe_stream_template.model. This provesmodelvalueagreement ontestedshots,
  NOTphaseorwhole-decodercorrectness. PrefixQ1428,T159664.

Wronskian identities: qformlocator sigma, correctly normalizedtemplate tau:
  e(x)*sigma'(x)*tau(x)=x^(2L) atactualerrorroots;
  H²=(sigma*tau)'=sigma'*tau+sigma*tau', degreeH<=L;
  e(x)=1/(H(x)*x^(-L))².
H*x^(35-L) fits36coeffs; fold x^-35 andfield-square intoCLASSICALmask atscan,
so the50Tkernel applies. ComputeWusingtwohalfdegreepolyproducts (~2*36GFmul)
ifoutofplace, but allocating432extraQwouldlose Q²Tscore. MainUNRESOLVEDpiece:
**in-place structured Sylvester map fromtemplate→H**, pluspayloadpadding/garbage
parking whilepreservingallinputs evenexceptional paths. In-placeHtransforms
need roughly<50kToffoliper direction torepay~30*4095scan-gatesaving afterlayout.

Important correction: simple37cyclicunrotationd works ONLYunconsumedpublicBM.
Streamingbank rotatesPREFIXofgrowingsize. probe_stream_template.normalize_bank
correctlyreverses prefixrotationsr=69downthrough2L+1 inclusive forL<=34;
L35needsoneoppositerotation. A degree-selectednetwork realizes thispermutation.
Afterthat validtemplate degree<=L+1 andconstant0. Highgarbage k>=L+2 canpotentially
park intozeroomega[k-2] slots (omega highcoeff>=Lzero). Thisparking isNOTyet
implementedquantummechanically. Inclassicmodeltrimmingthosehighslots works.

**100000-stream-model diagnostic** results/stream_template_normalized_100000.json:
6failures, ALLalreadyhaveincorrectlocator sigma. NoaddedWronskian failures among
99994correctlocators. Stops94063done,379secondpark,5553none,4stack,1bank.
Thisisuniformweight0..35, NOTthehostweightdistributionorfullnativeverification.
Originalunconsumedmodel100k had28failures; doNOTcallitPASS. Bothdiagnosticsretain
fullfailurelists. Initialstreammodelwithoutprefixnormalization failed; retained
results/stream_template_1000.json is that known-invalid exploration.

Counter convention: probe models returncounterproxy d=actualfinaldW+1 (theykeep
round68incrementforconvenientdegreeformula; itdoesn'taffectthatlastupdate).
Oncorrectpaths L=(69-actual_signed_dW)//2. Using7-bitmodarith69-dW thenbits1..6
cancomputeL reversibly inexistingcounter. Nativeprefixcheckconfirmsactual
counter==(modelproxy-1)mod256. DoNOTusemodelproxy directlyashostcounter.

Potentialomega preprocessing optimization: needonlyzeroomegahighgarbageparking,
notnormalizedphysicalomega. Use nonmonicreciprocal convolution. Existingomega()
recursion supportsnonunitdiagonal byscaling eachleafwordbysigma_recip[0], instead
ofnormalizingsigmacoefficients. Restoreσorientationafteromega beforeWronskian.
Could includeunborrowedS34 as36thsyndrome (tailS69..34); itsfinalomega35 isalways
zero oncorrectpaths, giving12cleanborrowedqubits with~36extraGFproductwork.
Needexplicitborrowing ofpublicS34 ratherthanassumingallocatorcanfreeregister.
Noimplementationyet.

At the end of that earlier research stage the archive was9.1631467535,
hash68bb0f6b...a2e07c6. The authoritative checkpoint at the top of this file
supersedes it. No experimental DQI component has entered the submission.
