# FOLLOW-UP COMPLETE — sign-mixing experiment and report

User asked to test Claude's 3b/chem/fit_mix32.py penalty on our verified eight-bit
719/C350 fit, rebuild if useful, then write a report for all three tasks.
This follow-up is complete. Current fully verified score **9.00962290835304**;
archive SHA-256 6ebc780a1c8027f26080e47818ea5db63dde58ce9ffbf90f46c4e7411b904800. Current submission.zip is byte-identical to
revisions/mixing_mu1e7/submission.zip. Milestones5–9 remain unchanged; milestone9
continues to hold the original9.008540688 checkpoint.

Best chemistry: candidates/mixing_mu1e7_submission,
results/mixing_mu1e7_verify.json PASS50880shots, score757163489.4986277,
lambda_eff48.89518731780483, T29440,Q526,6circuits,8bit,J13.
Inverse compactdefer_k4_q94 and DQI stream10_offset3_select_noscale_compact
unchanged, with existing full native PASS results bound to archive source bytes.
results/final_archive_audit.json checks all9members, chemistry source digest,
archive checksum, and unchanged milestone5–9 checksums.

Matched fixed-U refits:16000Adam steps2e-6,cap.019995,whitened,wBoffset.
No penalty control native lambda48.89575173155617. Weak physical-error-scaled
mu1e-7 native48.89518731780483, mixing79.77083. Mu1e-6 lowestmix77.26907 but
lambda49.01626544689377, gap3.93929; mu3e-6 lowestmix73.60242 butlambda49.41868,
gap4.01665. Both stronger runs retained baseline as best-costcheckpoint.
Literal added coefficient1e-6 native48.89566422098973, worse thanweakwinner.
Matched8000step QR rotations: control48.90337850324619 vs regularized48.91210486116479.
Seven new candidates all FULL nativePASS50880shots. All retain identicalU/tree/
two-body angle word hashes. Regenerated one-body words for refits; exactrotation
preserves h1 and could reuse originalwords.

Penalty conventions are explicit in constrained_factor --mix-normalization:
physical_error adds(mu/cap**2)*M; objective adds literalmu*M. Defaultmu0
unchanged. New exports atomic, stripped staleonebody, can retain_lowmix candidates.
rotate_qr optionalmix coefficient defaults0. measure_mixing.py records exact
physicalmetrics andhashes. results/mixing_experiment_summary.json containsfull
comparison, costs, parameters andnativepaths. BaselineM81.08235,l1sq603.02953;
Claude's lowmix511_off is11bitB511C375,M41.67816,l1sq563.63622, a differentdictionary.
Weakwinner saves0.036031%chemistry vsoldsubmission, only0.001154% vsmatchedcontrol;
overall improves0.00108222, notpredicted0.03. No globaloptimality claim.

EXPERIMENT_REPORT.md is the finished per-task report with worked/failed ideas,
verified progression tables, mixingablation, numericalpitfalls andprovenance.
README updated. No optimization jobs remain from this follow-up. Experiment3b
was read-only; no externalmessages sent or subagents used. Older notes below
are historical.

# GOAL ACHIEVED — 2026-10-23 13:44 EDT

Fully verified overall score **9.008540688342578 >=9.0**. All task verifiers PASS
at default problem sizes and shot counts. Canonical submission.zip and immutable
milestones/9/submission.zip have SHA-256 ded3c54c685f583d3f6cbc9e6f76324078119117e99e7b95b5da693871d1b8eb.
Milestones5–8 preserved. The final source files are:
- inverse candidates/compactdefer_k4_q94.kmb; results/compactdefer_k4_q94_verify.json;
  score190711.9554548813,Q94,T1778.000379mean/max1789,200540928shots.
- DQI candidates/stream10_offset3_select_noscale_compact_full.kmb;
  results/stream10_offset3_select_noscale_compact_verify.json; score1612737419463.0198,
  102400shots (hostPASS, phase errors included in scored repetition factor).
- chemistry candidates/residual719_lbfgs2_submission;
  results/residual719_lbfgs2_verify.json; score757436402.4877431,
  lambda_eff48.91281116246895,T29440,Q526,50880shots.

The final inverse build is reproduced by build_compact_invsqrt.py; rebuilt KMB
SHA4e1f02b42cf2843a7bf56507b64af9f2fa72fdab757f3be79f8ced4d52cdcbaf
matches the verified circuit. Native milestone package gates passed. Final archive
audit results/final_archive_audit.json binds archive bytes to their verified sources
and checks all milestone checksums. README has been rewritten for the final result.
Own remaining optimizers2236388,2236690,2236691 were stopped after the goal passed;
independent experiment3b processes were left alone. No further research is needed
for the requested9.0 goal. Older LIVE sections below are historical/stale.

# LIVE UPDATE 2026-10-23 13:40 EDT

Continue until fully verified overall >=9.0 or user interrupts; do not give a partial-progress final. Only edit experiment 3. No agents. Do not touch/kill independent experiment 3b jobs.

Canonical archive now **8.964068556176516**, SHA 0b182d780e0303985e7cb88a36ae2e2f8f743de261284243b9e1c49a06c124ac. Inverse tailproducts_l17_q97 PASS200540928 cost193403.45305897677; DQI unchanged stream10_offset3_select_noscale_compact_full; chemistry residual719_highlr_submission PASS50880 cost758067179.5737915 lambda48.95354472, T29440,Q526. Results same tags +_verify.json. Milestones5-8 preserved, 9 still pending.

PROMISING full inverse verification ACTIVE: candidates/compactdefer_k4_q94.kmb, results/compactdefer_k4_q94_verify.json/.log, session58101, parentPID2237072 workers2237083..85. Predicted cost190711.9 with94Q/T1789max; **should take total score ABOVE9 with already verified chemistry** if PASS. Full native200M shots --max-workers3 running. Another full verifier for defercurve_k4_q96 parent2236864/workers2236868..70 started2min earlier; ignore if94Q passes. No result yet as of13:40. If94Q PASS package canonical and immutable milestone9, update docs, inspect finalarchive, stop own optimizationjobs then finish task. Do NOT stop before verification/package.

Breakthrough code: binary_centered_invsqrt --defer_curve_phase (Python flag) measures unused curvature high bits after firstHorner stage and postpones phase repair until curve is reloaded for inverse stage, when ordinary Z conditional on stored measurement bits repairs it at zeroToffoli cost. Restoring uses original address (temporarily unflips centered x sign). --borrow_integration allows idle y-high bits as clean carry workspace for the two integration adders. product_limit17 in booth_adder defers top partial product after low products are erased. Existing default behavior unchanged. Final physical allocator holes require compact_wires with22 public input bits pinned; this reduces96 allocated wire IDs to94 live wires. build_compact_invsqrt.py reproduces final candidate. Source fixed table candidates/margin_curve8_t5.9e-08.npy. Final settings frac7,widths(27,18,8),lookup_frac8,split[18,9,None],carry1,rounded/Booth,reuse_normalized_top,output_extra10,defer_workspace,qubit_limit94,curve_keep4,defer_curve_phase,product_limit17,borrow_integration. Grids results/compactdefer_grid.json,borrowdefer_grid.json,defercurve_grid.json.

Active own chemistry jobs: residual719_native_long PID2236690 and partialwb_long2236691 (120ksteps8e-6, cap.019997, whitened, offsetLRfactor.3, nativeW / wb_scale.1); best valid~48.936 as of10k, outputs respective candidates. LBFGS2 PID2236388 at~200outer best48.9131458 but slowing/plateau; source highlr48.953. These need onebody8 --target.012 then build_global_chem and full native verification before any package; use safe bytes snapshots for live npz. All may be stopped once milestone9 achieved. Bad Adam16e-6 long was stoppedPID2236383 after diverging toLambda368. constrained_factor new --wb-scale and --offset-lr-factor default preserve native math, exports actualW.

Completed/failed this continuation: inverse variant oldtables grid no gain. Positive integer128 finished .0225887; positive bound L-BFGS-B719 .0220635 and grown988 .0216618 fail20mHa. positive_grow.py gradient-guided adds atoms but tiny initialweights/Adam largely retain22.59; larger dictionary minprojection does improve. Positive grown740/775AL endlambda52.13/50.29 (noncompetitive), physicalminimum rank350 19.71/19.39mHa. Original positive719integer physicalminimum rank35019.897. Joint angle24k hard20.5449mHa fail; hard snapshot refit impossible rank350 minimum20.215mHa; stoppedrefit2234676. Block projected8bit job stopped2222529 after2h and~31mHa. Positive projectedinteger stopped2232077 after~22.7mHa. C340whitened end52.18, gap4.03, noncompetitive. QR and500kexactJacobi currentw350 only tinygain (296.4H2). Truncating finalphysicalfactor to349 columns yields20.0095mHa fail, 34020.330mHa; no trivialrankreduction. ColumnclusteringR2B511C175 10seeds rawprojection error.18-.58Ha fails badly. DavioBoolean curvaturehighbits search41AND too costly vs lookup; no integration used. Most prior experimentjobs have ended/stopped; inspectps before changing anything.

# VERIFIED ARCHIVE UPDATE 2026-10-23 12:35 EDT

Canonicalnow8.711653852079, SHA762bf55a8a6734f2356a0e5a5e3e535fbe00c8f93c6c9d9e505abc8432232e93.
Currentchemresidual719_q8_initial_submission /same_verify.json823241449.1849973;
invbudgetunusedbooth_q98.kmb194025.3005777288;DQIunchangedstream10... .
C350factorrefinementbestlambda~51.8at6400,stillrunning16k. NEWsnapshotfullnative
pipelineCURRENTLYRUNNING fromparameters/residual719_q8_c350_factor_snapshot.npz,
outputsresidual719_q8_c350_onebody.npz, _submission, _verify.json/.log.
IfPASSandlowercostpackagecanonical, keepgoingto9. ExactneededlambdaC350rough48.24.
C340ALfinished12kbestlambda53.03,worsethanC350current,skipitsnativefornow.
Positive719anneal24k at7800soft.0186/hard.0795best.0342(sourcepositivefixedangle
.02861wasbetter;firstwstepjumped,don'tmistakebestforimprovement). Willsettlelater.
New--positive-weightsclampsw>=0whilejointθ/wfit. FixedUpositive71920kfinished
.02837(notaccurate), unshiftedLCU130.462, CPcheapnormbutgeometryneedswork.
No goalreachedyet; NO FINAL, continueuntil>=9oruserinterrupt.

# LIVE UPDATE 2026-10-23 12:31 EDT

NO PARTIALFINAL. Continue until fullyverified >=9.0 or userinterrupt.
Workonlyexperiment3,noagents,nevertouch3b/rootuserfiles. Readsubmission.summary.json
forcanonicalexactscore; newestpackageinflightfromresidual719_q8_initial_verify.json.
Fullyverified8bitchem PASS score823241449.1849973,lambda54.01359098,T28976,Q526,
brot8,J13,6circuits. circuits candidates/residual719_q8_initial_submission.
Inverse nowbudgetunusedbooth_q98.kmb PASS200540928 score194025.3005777288,Q98,
T1728.999795mean,max1740,34rotations. DQIunchangedstream10...1612737419463.0198.
Beforethis8bitchem updatecanonical8.6572499489 SHAa6d7803972206211f52cd5acc8280a61e40c88deec3282cddc0a36e218c5035c.
Milestones5-8preserved,9stilloutstanding. Threshold chemistry746623289.53withcurrent
inv/DQI. ForB719C3408bitT28976Q526 needlambda<=48.9866 toreach9.

ACTIVEIMPORTANT:
1) residual719_q8_factor.log constrained_factor12000 lr5e-6 cap.01998, at8400
   currenterror.020023violates butBESTVALIDlambda53.1702. Inputresidual719_q8_nudge.
   Outputcandidates/residual719_q8_factor.npz. NeedONEBODY8+build+nativeverifywhenready.
2) residual719_q8_c350_factor.log16000 lr5e-6 cap.01998, at2800bestlambda53.0867,
   error.01998267gap3.915949. Morephysicalerrorroommightallowfurthernormreduction.
   PipelineC350: project_factor.py fromsame719C340snapshot->offsets->align_columns
   --physicalreferencepositive650floatC350->QR6000->nudgegapfloor->AL16000.
   Outputresidual719_q8_c350_factor.npz needsONEBODY8+build+nativefullverify.
3) residual719_positive_q8_anneal.log NEW24kregularized_angle_fit --positive-weights,
   sourceparameters/residual719_positive_angle_init.npz (fixedanglespositivefit.02861),
   θlr.02,wlr2e-5,gamma0,startstrength1e-4,end100,unwhitened. Keepsw>=0 andglobal
   tracefixed, unshiftedLCU130.462. Bestchancefornormsmallenough9ifaccuracy<20mHa.
   Outputresidual719_positive_q8_anneal.npz (hardU/w/words coherent) and_soft.
4) residual719_q8_c340.log global_coordinate512moves/pairs1024/energy/padding.003,
   fromgrowth719ofresidualcandidates. At~240 physical.01995609 plateau; outputK
   whilerunning. Snapshotparameters/residual719_q8_c340_snapshot.npz .01996812.
   NEWgrowthfrom700reseedC340.020679added19atomsresidualpool, crossed20mHa~step70.
5) chem4b256c128_q8_projected.log15000steps, --whitened --project-every5cycles3,
   now~6000sof.022/hard.04875. Slow~.8sec/step, likelyfails30mHafinalbutnotfinished.
6) residual719_q8_positive.log fixed-vectors positive_global_fit20000step lr.0005,
   fixedtraces,at12kbest.02861014 (fails), endssoon. NonnegativeinitNNLS5000steps
   nonnegative_weights.py usedreferencephysicalpositive650X, Xerror.05039. Joint
   theta/wpositiveannealaboveisnextstep.
7) positive700_c335_oldtree_quadratic.log fromfloated700C335(.01250) roundingoldtree,
   at~150 .0245fails. Newtreeversionstoppedplateau.0238(PID2227001).
Otherolderstalledjobsmostlyendedorstopped. Checkps, ownjobsonly.

CODEANDRECENTFINDINGS:
- short_adder.add prefix_bits Cuccaro lowprefix+measurementuppercarry. qubit_limit
  picksminprefixbasedonliveQ. Boothusesunusedcleanyhighbits+diff ascarryworkspace.
  binary_centered defer_workspace delaysyhigh/countersuntilafterlookup, freesbefore
  deletion. budget_booth_unused_grid.json best98Q/T1740. Curvepartialfreegridbest
 195825(worse), don'trepeat. Adderexhaustive8192casesvariedm/n/prefix/sub/carrypassed,
  bestcandidateFULL200M nativepassedandpackaged.
- align_columns --physical NEW alignsrepresentedX matricesacrossdifferentB,
  referenceC>=targetC (selecttopphysicalnormCcolumns), rotatesw/Wexactly. For719C340
  rawH21412->aligned513->QR354.403. Offsets gap4.082914, gramLCUlowerbound277.04.
  Nudgewithtight.01999capcouldonlylowergapto3.9931(lambda54.0136), becauseinitial
  modelerror.019968 leftlittleerrorbudget. C350projectionerror~.01957 thennudgefloor
  .0196471,gap3.916,lambda54.002. ALc350bettererrorroom.
- global_coordinate/global_fit FINALwexport nowtrimsrequested--columns instead
  ofexportingfullB withzerocolumns. ActiveolderprocessesmaystillwritefullB atfinish;
  trimmanually/project_factor.py asneeded. WhileactiveKcheckpoints convertusing
  eigenKtopC. project_factor.py NEWfixedUprojectiongpu, outputscoherentwC.
- quantize_onebody.py --target .0098 NEWstopwhenpasswithmargin; baseline8bit330steps
  .0097988nativepassed. Alwaysregenerateonebodyafterw/Wchanges.
- constrained_factor.py NOWsavesvalidINITIALcandidate beforefirstAdamstep;
  olderactiverunsdon'thavethiscodebuttheyhavebestvalidsnapshotsnow.
- positive_global_fit.py --fixed-vectors freezesU toexactquantizeddictionary,
  --fixed-tracesglobalphysicaltrace. regularized_angle_fit.py --positive-weights
  clampsraw w>=0 afterstep, incompatiblewithwhitened; defaultsunchanged.

RESEARCHBRANCHESCOMPLETED:
- Unionprunepositive700+oldgrow725: unionfullerror~.010butbackward700->.02107;
  biasedprune3e-7/1e-6retainsallfirst700->.019716fullrank. LowrankC340 plusQP
  plateau.021019. target-rank340selectionultimatelysametrajectory; bothstopped
  PIDs2227592/2227778. Noexactduplicates! I initiallymisreadremovedrowcountsas
  duplicatesandcorrectedtheuserimmediatelyafterdirectcheckfound700unique700.
- reseed_dictionary.py removes16weaksourceatoms andaddsresidualeigenmatrix
  directioncandidatesquant8. Fullerror.020077butQP C340->.020679(betterthanold
  .0214). Sourceparamspositive700_reseeded_q8.npz, residualpool
  parameters/positive700_residual_pool_source_pool.npz. Thiswasgrowth719seed.
- swap_dictionary.py monotoneoneswapatatimepreservesactualrankfullloss;64swaps
  accepted2tinyimprovements only. q340plateau.02103, stoppedPID2228195(alreadyended).
- optimize_atoms.py continuouscandidateoptimizationor--latticeSTEangleoptimization,
  maximizesprojectedtargetcapturegain, quantizes/proposes12candidates andchecks
  actualrank340errorbeforeaccept.8cyclesfloataccepted4tinyimprovements(.0206808),
  8cycleslatticeaccepted1(.02067977), notenough. Logspositive700_optimized_atoms_q8
  andpositive700_lattice_atoms_q8. Optionalpoolfromresidualscriptabove.
- old700reseedQPc340stoppedPID2227840at.02067937; snapshotpreserved.
- global500_q8_c350_quadratic failedplateau.03835, ended.
- global500_lbfgs_fast2STOPPEDPID2217635after440 .01681255; validcoherentfinalsnapshot
  parameters/global500_lbfgs_final_snapshot.npz preservedbeforeSIGINT.
- positive650_varproj_q8finished12k .03152699fails; blockwhitened18kfinished.0315317
  fails. Directcoordinatebasisfeaturedictionaryscreenerror500basis.73Ha,1200.259,
  needs>2600bases20mHa; simpleaxis/pairstructuredQROMnotpromising.

# NEWEST STATUS 2026-10-23 ~11:55 EDT

Keep trying until fully verified >=9.0 or user interrupts. NO PARTIAL FINAL.
Only edit experiment3; no agents or touching3b/rootfiles. Canonical NOW
8.456865172711836, sha268a42bae877b7f2e51e18ba154e0ffeecb1fe0e5569af3b2dc806a904af4a9b.
- invsqrt candidates/margin_curve8_t5.9e-08.kmb, same_verify.json PASS200540928,
  cost208146.9321285,Q114,T1575mean,max1586,eps1.00156e-13.
- DQI stream10_offset3_select_noscale_compact unchanged1612737419463.0198.
- chem lbfgs500_q10c375_submission /same_verify.json PASS838859431.3087063,
  lambda53.21747,Q833,T18923,brot10,7circuits. Milestones5-8 untouched.

NEW INVERSE BREAKTHROUGH: short_adder.add --prefix_bits (Pythonparameter) hybrid
Cuccaro low prefix + temporary-AND upper carries, tradesextraT forlowerQ.
--qubit_limit automatically choosesminimum prefix percall for livequbitbudget
usingbuilder.num_allocated_qubits. Booth uses unusedborrowedyhigh bits ascarry
workspace, plusdiffaux. binary_centered --defer_workspace delayshighy/counters
untilafterlookup, freesbeforelookupdelete. Defaultsunchanged; newcandidatefixed
coeff tablemargin_curve8_t5.9e-08.npy. Results:
- segmented_adder_test.log exhaustive8192casesm/n1/1,1/6,3/6,6/6,random20/28,
  prefixes1,3,6,10,sub±,carry±,borrow1: all value/phase/junkPASS.
- budget_booth_unused_grid.json cap98bestQ98,T1740max, predictedcost194025.32
  (vs208146). candidates/budgetunusedbooth_q98.kmb **FULL200MVERIFIERNOWRUNNING**
  results/budgetunusedbooth_q98_verify.json/.log, --max-workers3. OncePASSpackage
  iflowercost; shouldyieldoverall~8.66 andraiseallowablechemistrycostfor9.
  Conventionalprefixgridbest15,0Q100,T1746 worse. unusedworkspacecap98~minimal.

ACTIVE chemistry jobs (inspectps/logs; someoldjobsstopped):
- positive700_q8_energy: global_coordinate --moves512 --pair-candidates1024
  --candidate-energy --curvature-padding.003, from.020931, now~.020815plateau.
  NEW global_coordinate opt-in quadraticcoupledproposals fromapproximate
  variable-projectionHessian, binarygreedyICM; exactobjectivecheckedbeforeaccept.
  Canchangemanyangleswithinonerow. Gradientenergyranking +reducedpaddinghelps.
- positive700_c335_relax: global_fit --columns335 --steps2000 --lr.0001 --psd
  --spectrum-rank512 fromenergyquant700. Floatat1600 .012543. Shellautomatically
  followswithoptimize_tree6000(newcommon8bittree), thenglobal_coordinate3000
  moves512,pairs1024,energy,padding.003 topositive700_c335_newtree_quadratic.npz.
  WATCHnextlogspositive700_c335_tree and...newtree_quadratic.
- global500_q8_c350_quadratic: newbestfloat500.016865 nearest8->.1325,
  quadraticmoves512 descends.03835after240steps, stalls(fails), mayprune.
- global500_lbfgs_fast2 still1000outer(~380,float.01685),slow~20sec/outer.
- chem4b256c128_q8_whitened18k: at~12k .0335hard/.027soft (fails).
- chem4b256c128_q8_projected15k: --whitened --project-every5 --project-cycles3,
  weightlr0. At~3200 soft.0166/hard.10; ~.8sec/step,likelyhours.
- positive650_varproj_q8 12k: at~11k hard.03154, noncompetitive(fails).
- positive700_norm_positive rotate_positive3000: at2500H2~431, worseQR372.51.
  Canpruneafterfinish/stopownjobifneeded. rotate_positive.py NEW triesrotation
  towardnonnegativew, preservesK. NotyetbetterthanstandardQR.

Chemnormprobe fromQPlattice700(snapshoterror.02093fails): offsets gap4.083,
rawLambda1489,Kgramlowerbound279.72. Align_columns withreferencepositive650w
PAD50ZERO rows. Jointw/wB alignmentwithWref0wasbadH2~1052. NEW --weights-only
alignmentH2531.79 thenQR6000->H2372.5097 (positive700_norm_qr.npz). Guesslambda
~55aftergapfloor, stillneedmorecostreductionorbetterlowLCUfit. Nativecan'tuse
thisuntilphysicalerror<.02. NMFpositivecolumnrotationisworse431.

Stoppedduplicativeoriginal700C350/C320integerjobs, pair128job, quadratic32,
quadratic128C350PID2224608,quad512PID2224762,quad1024PID2225096. C320quad128
PID2224612 stillrunning~.0226; canpruneorupgrade. Oldblock_integerstopped.
TenbitintegerC375finished2000 .019539, outputwfullB500atend; snapshotcanonical
usedearliererror.019744, couldrefinefuturebutnotpriority.

# VERIFIED UPDATE 2026-10-23 ~11:25 EDT

CANONICAL NOW 8.44000931221495, submission.zip sha256
16161df74d91068e9dadf9cae6eac7d6f9941e67616af3b5dd34ddbf891faab4.
Same inverse209396.5186 andDQI1612737419463.0198; chemistry now
candidates/lbfgs500_q10c375_submission, results/lbfgs500_q10c375_verify.json:
PASS36864shots, score838859431.3087063, lambda53.21746717,Q833,T18923,
brot10,J11,7circuits. Fullpipeline source snapshot then offsets/align/Jacobi/QR,
nudge/constrained_factor8000/onebody. Onebody was stopped after~400steps because
error .0023 already far belowthreshold; shell built and fullverified afterward.
READMEscoreupdated; oldmilestones5-8 preserved. Need chemistry<=691815743 ifothersfixed.
No finaluntil9oruserinterrupt. Keepgoingautonomously.

NEW ongoing jobs sincepreviousLIVEsection:
- positive700_grown_q8_c350 (PID2222292): 6508bitangleanneal .03457 grown700from
  positive725C500float candidate pool, projectioninitial.02827, integer4000steps;
  at~500 .02151 (fails), Kgramlowerbound139.5 (promisingnorm).
- positive700_grown_q8_c320 (PID2222619): sameinit rank320, at~400 .0233fails.
- positive700_q8_pair(PID2223483): fromC350step300snapshot, pairsearchonlywhen
  singlemovesfail, at180 .021524. Could stopoldduplicativeC350/pairjobsifneeded.
- positive700_q8_quadratic (PID2224293) NEW coupledmove search, fromlatest.021513.
  global_coordinate --pair-candidates192 --pair-trials12 now buildsvariable-
  projection approximate off-diagonal Hessian, greedybinaryquadraticproposal
  up32angles (multipleperrowallowed), triesjointproposalFIRST; laterpairtrials.
  Defaultglobal_coordinate behaviorunchanged unlesspair-candidates>0.
- chem4b256c128_q8_projected (PIDlookupps): NEWwhitenedjointanglefit withPSDbloch
  weightprojections every5steps3cycles, weightlr0,15000steps,startstrength.003.
  At1400 soft.0148/hard.152. Slow~.75sec/step, likely2hours.
  regularized_angle_fit --project-every --project-cycles optional, preserves
  projectedAcolumnorientation bybatchedSVDalignment, resetsweightmomentum.
- chem4b256c128_q8_whitened18000: at6400soft.0185/hard.0717, gamma0,projectoff.
- positive650_varproj_q8: at6400soft.01974/hard.03397,12ksteps, threads1,spectral512.
- anneal_lbfgs500_q10_c375_integer: nearfinishing2000, ~.01954; bestsnapshotmay
  improvefutureconstrainedfitbutcanonicalalreadyverifiedfromearlier.019744snapshot.
- global500_lbfgs_fast2: at~300 fullfloaterror.0169, slow1min/5outersteps.

positive700_ridge_grid.log testedridge1e-5..1e-8 onbothC320/C350geometries; lowering
ridgeimprovedactualerroronly5microHa andKboundslightlyincreased, not useful.
block_integer oldjobPID2221290stoppedafter~.032 (notcompetitive).
Inverse coefficientbitsGF2rank testfoundnoextraexactlinearrelationsbesidestopC0=1;
noqubitsavingthere. DQIzeroerrorswouldonlyslightlyreduceRconfidencepenalty, not18%.
UpstreamPsi methodrequiresclean35cofactorwords; streambankconsumedsyndromesoccupy
these, addingworkspaceundoesTgain. NoimplementedbetterDQIalgorithmyet.

# LIVE STATUS 2026-10-23 ~11:10 EDT (newest)

Continue until FULLY VERIFIED overall >=9.0 or user interrupts. Do not send a
partial final. No agents; only edit experiment3, never touch user3b jobs/files.
Canonical remains 8.2142104791 (same sources as older status below).

Current own active experiments (inspect ps/logs before restarting):
- anneal_lbfgs500_q10_c375_integer: global_coordinate from finished q10 anneal,
  C375 rank projection. PASSED MODELED .01963 at step850, continues2000steps.
  Snapshot parameters/lbfgs500_q10c375_snapshot.npz .01974479 converted K to w375.
  offsets -> align_columns to canonical -> exactJacobi100k -> QR5000 -> nudge:
  candidates/lbfgs500_q10c375_nudge.npz modeled error .01981529 gap3.91600075,
  lambda_eff54.34513. Now constrained_factor8000 lr1e-5 cap.019995 running,
  out lbfgs500_q10c375_factor.npz. NEED ONEBODY and BUILD + FULL VERIFIER.
- positive650_varproj_q8: anneal_angles12000, C350, spectral512, startstrength.003,
  threads1. At2100 soft.01523/hard.0677. Full-target checkpoints correct.
- chem4b256c128_q8_whitened: NEW regularized_angle_fit --whitened18000,
  source best block float .01406785, gamma0, weightlr.0001,startstrength.001.
  Whitened factor A=L^T w parameterization with differentiable fp64 Cholesky
  and fp32 physical loss; snapshots store actual decoded U and w. At1400
  soft.01415/hard.178. Slow .12sec/step but promising better conditioning.
- positive700_grown_q8: NEW grow_quantized from positive650_q8_pure(.03457) using
  positive725_c500 float target. Adds50 neighboring rank-one features. Shell
  automatically launches global_coordinate4000 C350 moves32 afterward. Could
  lower quantization error while retaining a smaller normalization.
- chem4b256c128_q8_integer: NEW block_integer.py alternating PSD block refits
  (12cycles) and discrete moves with actual fixed-w loss checks. At310 .03212,
  improves too slowly, may prune later. Out snapshot U/w/words coherent.
- global500_lbfgs_fast2: spectral512 correct descending eigenvalue sorting,
  still running1000 outersteps, around200 actual full error .0170.

Finished:
- chem4b256c128_fit 30000 steps normalize-trace: float .01406785 (strongestblock).
- chem4b256c128_q8 first gamma1e-6 anneal: .034672.
- chem4b256c128_q8_pure gamma0/betterfloat: .033582, still fails.
- positive650_continue float .0189023; positive725_c500 float .01722058.
- positive650_q8_pure anneal gamma1e-6: hard .034572, fails.
- chem4b256c128_float_offsets lambda85.68; exactJacobi4*100k ->lambda60.2368,
  gramLCU lowerbound337.27, gap4.0809. QR may lower further if quant8everpasses.
- chem4b256c90_splitfit ~.0231, chem8b128c45_splitfit ~.0404, worse than otherbasin.

New code: block_integer.py, align_columns.py; regularized_angle_fit generalizedR
and added --whitened/--start-strength. anneal_angles optional --spectrum-rank
uses largest512 eigenvectors, full-target hard-error recomputation each100steps,
--threads. Defaults unchanged. README firstline fixed8.214 andmilestone8status.

# LIVE UPDATE 2026-10-23 ~10:18 EDT

No partialfinal; keeptryinguntilverified9.0oruserinterrupt. Milestones5-8saved.
CanonicalNOW8.2142104791: inverse adaptivebooth_f7_m8_w27_18_8.kmb andsame_verify.json
cost209396.5186399,Q114,T1574mean,34rot,PASS200540928shots. DQIstream10unchanged,
chemglobal500_factor_polishunchanged. Alternativeadaptivecurve9nativePASS209783(worse).
Readnextstatussectionforpaths. READMEfirstlinecurrently8.188stale(updatefromsummary).

Newactiveexperiments:
- positive725_continue.log50kAdamfrompositive725_c35025kfinalerror.02482, at20kcontinued
  ~.0212. Allw>=0. UnshiftedLCUrow bound130.465 (tracefixedphysically). FullSOSoffset
  invariantlowerbound~270: tracesnorm521.94,offnorm1081.04,traces.off=-748.58,
  soactualnormalizedgainlimitedbutmayhelptoreach9withB600-6508bits.
- positive600_c350.log30kAdamfromgrow600_float projected16bitC350; ~.031after19k.
- positive725_trace.log30kAdam .0003 frompositivecontinue~.021, --fixed-traces actually
  normalizes GLOBAL sum_c(traceXc)^2 to referenceI G I, notpercolumn. Aimsstability.
- regularized725_trace_q8.log16k jointangle/w annealing, gamma1e-5,rowLCU penalty,
  weightlr1e-5,θlr.02, globaltrace-normalizedsignedw. regularized_angle_fit.py NEW.
  Initialregularized725_q8 (withouttracefix,weightlr.0002) wasunstableandstopped.
- anneal_lbfgs500_q10_c350.log12k, ~soft.0186/hard.028after3500. Active.
- global500_lbfgs_fast.logNEW1000outersteps, lr1, spectrum-rank512, fromsnapshot
  global500_lbfgs_polish2 ~.0171. global_fit.py --spectrum-rank speedsgradientusing
  F=PZ,H=FF^T. Every5steps ERRORRECOMPUTEDAGAINSTFULLTARGET; every10snapshotusesfull
  targeterrorforbest. Verifiersunaltered. Oldphase2jobPID2215446stoppedafter~150;
  snapshotparameters/global500_lbfgs_fast_init.npz preservedvalidinput.
- grow700_q8_c350 old~.0205 andpruned700_q8_c350~.0204 stillrunningbutfailing.
- grow725_q8_c350~.0196 active. ItsinitialfullpipelineVALIDbutcost999160863 worse
  thancanonical909956199: lambda63.598,Q526,T29868,brot8,6circuits51264shotsPASS.
  NormafterJacobi/QR~505 vs355forcanonical500. Couldtrystrongernormobjective/pruning.
- Positive725_lbfgs (snapshotmidAdam) STOPPED~25outer .0266; finalAdamwasbetter. Direct
  positivefitwAdam25k takes~6min; PyTorchLBFGS1000wouldtakehours, notworthstaleinit.

Newfields inpositive_global_fit.py --optimizer{adam,lbfgs}, --fixed-traces; lbfgsuses
squaredweights. Adamclampsweightspositive; explicitfp32matrixresidualstable near.02.
PPTfeasibilitytest reference5776partialtranspose: eigenmin .00288572 positive. Doesnot
provecompletelypositivebutdoesn'truleitout. CanonicalXonly28/375PSDcolumns. Lowrank
screenofXshoweddiscardedsquaredmass1.967atB12, .48atB32, .176atB48, solowrankC1
ensembleideaunprovenandnotlaunched. Browsersearchedbilinearfieldalgorithms; no usable
improvementoverour42ANDGF4096/15ANDGF64found. Noexternalcodeexecuted.

Onebody8optimizedhostfunctionPASS.0094006717456 exactlyvsmodel.00940067174562.
725full8quantizedonebodyoptimized.0097577nativePASS. quantumonebodyconstraintsolved.

# NEWEST STATUS 2026-10-23 ~10:05 EDT

Continue experiments until FULLY VERIFIED overall >=9.0 or user interrupts.
NO PARTIAL FINAL. Workonly /initial_experiments/3; noagents, don'ttouch3b/rootfiles.
Milestone8 NOW SAVED at8.0920837321 in milestones/8/submission.zip. Milestones5-8saved,
9outstanding. Canonical archive now8.1881338541, exact sources:
- invsqrt candidates/booth_borrow_f7_w27_18_9_e1.kmb
  results/booth_borrow_f7_w27_18_9_e1_verify.json PASS200540928shots,
  cost211403.482332,Q116,T1561.99939mean,max1573,35rotations.
- DQI stream10_offset3_select_noscale_compact_full.kmb / corresponding_verify.json,
  cost1612737419463.0198,Q1429,T645755,8phaseerrors,0valueerrors102400shots.
- chemistry global500_factor_polish_submission / global500_factor_polish_verify.json,
  cost909956199.7708666,lambda53.61396,Q833,T20375,brot11,C375,B500,36864shotsPASS.
Need chemistry<=~685M ifothercostsstayfixed toreach9.

ACTIVE JOBS (inspectps/logsbeforerestart, ownprocessesonly):
- positive725_c350.log positive_global_fit.py25kAdam maynearlyfinished: error.03ish.
  NEWidea nonnegativew givesPSDphysicalX, no coefficientcancellation. Referencepartial
  transpose test foundsmallestpositive.0028857 (PPTdoesnotruleoutPSD X decomposition).
  New positive_global_fit.py supports --optimizer lbfgs (w=square(raw), float32explicit
  residual), defaultAdamclampsw>=0. positive725_lbfgs.log1000outerstepslr1 fromsnapshot.
  Needaccuracy<.02; normalizationcouldbemuchsmaller, potentiallybestchancefor9.
- global500_lbfgs_polish2.log400steps continuation ofphase1300steps. Phase1finished
  .017284vsoldfloat.01851! Phase2~.01718after70. Source global500_lbfgs_polish.npz;
  treebaselinefallback. L-BFGS slow~8secouter butimproving. Couldenablebrot10/C350.
- anneal_lbfgs500_q10_c350.log12ksteps, newfloat.017284, startstrength.003,lr.02.
  anneal_angles.py now --start-strength (default1e-4unchanged). Soft~.01835,hard.037early.
- grow725_q8_c350.log6000integersteps, ~.0197after900 (validmodeled). B7258,C350.
- grow700_q8_c350.logoriginal700from650,~.0206after1000fails.
- pruned700_q8_c350.log725backwardprunedto700,~.0205after400failsbutbetter700basin.
- adaptivebooth_f7_m8_w27_18_8_verify.json full200Mshotverifier running.
  BuiltadaptiveBooth table268leaves,Q114,T1585,worstrel5.943e-8(nearlyboundary).
  Estimates~208620butrotationmargincouldraiseit. Alternate adaptivebooth...18_9:
  257leaves,Q115,T1582,rel5.872e-8,notyetnativeverified. Bothcouldimprove211403slightly.

7258bit FULL PIPELINE VERIFIED PASS BUT NONCOMPETITIVE999160863 (notpackaged):
grow725_q8_c350_submission / results/grow725_q8_c350_verify.json,
brot8,Q526,T29868,lambda63.598,6circuits,51264shots (checkexactshotsinjson).
Projected725frominteger~.01991 -> offsets gap4.082, rawLambda1649 -> Jacobi500k H2~832
->QR6000H2~505 ->nudgegap3.916,error.019930 ->onebody8optimized ->nativePASS.
Norm istooexpensive;positivefactors/pruningaimtolowerit. Canonicalchemstill11bits.

ONEBODY8 SOLVED: quantize_onebody.py optimizesintegerrowsjointlyforfixedH/eigenvalues.
Analyticexactsingle-anglecandidates using firstandsecondderivatives, weightedmatrixloss,
rowgroupmoveswithbacktracking. Probe canonical11physicswith8bitonebody: initial~.05,
optimized~.009374. Unchangedhosth1_thc_reconstruction_error reproduced.0094006717456
vsmodel.00940067174562(midrun). 7258actualonebody~.009757passes.
build_global_chem.py now accepts optional onebody_angle_words+onebody_brot in fit;
assertbitsmatch. Generateonlyafterallw/Wchanges. Theseonebodyrowsaren'torthogonal,
whichisallowed; hostHungarianassignmentcheckedwithitsunchangedfunction.

BOOTH INVERSE BREAKTHROUGH:
booth_adder.py radix4signedproductaccumulation; eachpairlocalbits sharesoneadder.
Magnitudecontrolsodd=b0 XOR bprev,two=(NOTodd)&(b1 XORb0),sign=b1. Productrow has
odd*v +2*two*v, twoANDsperbit, butonlyonecarrychain. ManualHMRproductXORclaims
correctphase. Exhaustive4096shotunitcases bothsigns,5/6coordbits andtruncatedshiftsPASS.
short_adder.py optional carry_workspace letsusecleanborrowedcarrywires, defaultunchanged.
Boothclearslineardiffbeforeadder,reusesitforcarry,thenrecomputesforANDcleanup: saves1Q
withoutToffolichange. binary_centered uses normalizedx21known1 asconstantC0topbitAND
auxiliary: Xto0duringBooth, restore1beforeintegrate. Caveatx21onlyknown1betweentarget
highintegrations; codeorderscorrectly. Allcurrentbestbinariesnativeverified.
binary_centered model/build now booth=True or perstageboollist, extra_frac works.
C++ fit_coefficients modes booth or boothlast(stage0only); fractionN/globalphase flags
compose. fit_adaptive.cpp accepts argv8 booth; compiledfit_adaptive_booth.

Verified inverse progression (milestone8usesfirst):
booth_f8_w27_17_7 cost219021,Q112,T1712mean (PASS200M).
booth_aux_f8_w27_17_7 cost217065,Q111, sameT (PASS200M).
booth_f7_w27_18_9_e1 cost213226,Q117,T1562mean (PASS200M).
booth_borrow_f7_w27_18_9_e1 current211403,Q116,sameT (PASS200M).
booth_more_fit.json /booth_more_estimates.json holdsvariouspassesbutbestabove.
Mixedboothlastw16failed. Globalphase48casesafterBoothallfailed. AdaptiveBoothgrid
results/adaptive_booth_fit.json holdsmodels/nativecandidates forfullverification.

Failed/prunedexperiments:
- STE6008Adam12k stopped~2500,errortripledfrom.0357->.11. PreconditionedSTE+jitter
  stoppedaftererror.15within100. global_quant_fit.py now supportsPSD,rankconstraint,
  source tree, optionalpreconditionedNewtonSGD andjitter; monotoneintegerfitbetter.
- anneal6008finished20k .0367566(C320); fullrank350 .035718. coordinate6008C350
  stopped~800 .02662; grow650neighbors+integerstopped~735 .02293. OutputsKsaved.
- grow700from650selected50neighbors went.0212then.0206, activeabove.
- grow725from700selected25neighbors passes~.0199; backwardpruned700initno-ridge
  error.019754 butrank350ridge.0206afterfit. Cantryiterativegrow/prunetoreachB6508.
- grow51111C350modeledpass.01995butnative982968212lambda58.62 (worse), stopped~3k.
- global500ridge4Adam6k finishedbestearly.01901, later.023 butGramLB138 vs154;
  quant11 projections~.0256fail. Nothelpful.
- optimize_tree.py weighted8bitHuffman/greedy+30kpermutations improvedproxy.355->.295,
  trueprojection6008still.1002(initialold.116); notpassing. Savedtreeparameters.
- roundingmasksglobalgrid~400nonepass; perintervalmaskcoverage4patternsneeded
  forf8width26/17butnotcheap, f7width27/17still33intervalsuncovered. Tieeven/odd120
  fitsnonepass. Coarsecurvetablesandpartialdeletecurve_keep_bits allnoncompetitive.

Readoldcontextbelowforbaselines,packager,environment,DQIimplementationdetails.

# Latest continuation update (2026-10-23 ~08:55 EDT)

KEEP WORKING UNTIL FULLY VERIFIED SCORE >=9.0 or user interrupts; no partial final.
Current archive is7.9461644218 with inverse rounded_reuse1_extra10 (230949.5054),
DQI stream8_offset3_select_compact (1615254361201.19), chemistry global500_factor_polish
(909956199.77). All defaults full shot counts passed. Milestones5,6,7saved;8/9open.
A later verified DQI stream10_offset3_select_noscale_compact improves slightly again
(1612737419463.02,Q1429,T645755,8phaseerrors,0valueerrors). Not packaged yet.

Active jobs: inspect ps/logs before restarting.
- rounding_mask_fit.log: fit_rounding_masks.py grid with mixed per-row rounding.
  New optional rounding_masks supported in binary_centered model/build, C++ roundmask
  mode reads masks from argv8 comma-separated. Hoping f7 widths27/17 or f8 widths26/16
  fit when some rows truncate instead of round. Native verify any cheaper candidate.
- grow600_optimized_tree.log: optimize_tree.py brot8 30k proposals. Tests Huffman,
  unbalanced greedy quantization, leaf permutation annealing, weighted vector error.
  Output parameters/grow600_optimized_tree{,_quantized}.npz. Needs true PSD projection
  then lattice fit; tree proxy alone is not physical accuracy.
- grow511_q11_c350.log: canonical11bit dictionary +11 selected neighboring vectors
  (grow_quantized.py) then rank350 integercoordinate6000steps. ~.02008 after700.
  B511 costs same inner QROAM blockcount as B500, fewer C may help. Must reach.02.
- anneal500_q10_coordinate.log: 10bit rank350 annealing FAILED .0227078; full projection
  .021899 then integercoordinate6000. ~.02025 after2800, still failing. Let finish.
- anneal600_q8_c320.log: 8bit600 annealing20k nearingfinish (~.03675 hard), FAILS. Float
  grow600_float_fit finished .01177. Need a better quantization tree/more dictionary.
All published continuous fits finished; none beat old500 significantly.

New inverse workspace optimization: build option reuse_normalized_top=True uses input
x[21] (known one after normalization) for the constant top C0 bit. output_extra_bits=10
also removes unused integer-phase bit. rounded_reuse1_extra10 verified200540928shots,
cost230949.5054 Q108,T1892mean. Source binary_centered supports defaults unchanged.
Adaptive lookup experiment (fit_adaptive.cpp, results/adaptive_interpolation_fit.json)
found valid f7 tables but costs242k+, no improvement. Coarse-curve and fraction searches
found no cheaper passes. Optional curve_keep_bits implements deleting high curvature
bits after C1, restoring before inverseC1. It works at construction but all estimates
worse (best keep5 Q106,T1971). del_lookup tiny-target needs local workspace: include
known-zero unused high y bits in target, no extra allocation. No native verify needed
for these noncompetitive variants.

Chemistry polish: rotate_qr lowmix .02,lr.0002,smooth.00005,12k fromfactor_restart3
->global500_factor_polish.npz (lambda_h2=355.368 vs356.368), full native PASS909956199.
Overcomplete425 nearidentity QR gives355.259 but C increase too costly. Large-mix10
QR failed badly, kept original. offsets_score onpolish barely improves53.61396->53.61295,
not worth rebuilding. Fixed canonical U rank errors: C340 .020754,C350 .020407,C360
.020160,C365 .020062,C370 .019979,C375 .0199055. Could pipeline C370, slight costgain.

DQI major new changes (stream_capture.py selects options by output tag):
1. stream_forney.py keeps original sigma0 scalar, implicitly uses monic coefficient1.
   Chien acc_s now sigma2, sums other even coefficients + derivative accumulator,
   flips low bit for the implicit1. Removes oldzero12 storage. Independently validated
   test_implicit_monic.py:256shots71positions, values/phase/junk PASS.
2. Bank offset4 originally had1426phaseerrors! Classical model consumed alloddrounds,
   while circuit stoppedconsuming at63. Actualbank guard68 stoppedvalid35errorcases
   because borrowedS32 not consumed. Fix extends consumption to65 foroffset4. Model
   still offset4correct. Gating discrepancies was a failed diagnosis (sameerrors,+32kT),
   only used if tag contains gated; never use forbest.
3. offset3 extends consumption through67 and defers odd-round bank materialization/
   rotation untilafter consuming the newly free input. Stack guard before discrepancy,
   bank guard skipped only when trailing slot is newly materialized (conceptually0).
   stream7_offset3_compact PASS102400 Q1433,T652283,6phaseerrors,cost1629523997345.
4. select quotient: use six-bit gated inverse norm (ctl*yd) for final inverse product.
   mq=(active XOR ctl)*pv + ctl*inv(pv). ctl impliesactive. No dinv/mraw12 registers;
   keeps18chain+6selectednorm throughyield. Saves5712T vscompactold withsameQ.
   stream8_offset3_select_compact PASS102400 Q1433,T646571,6phaseerrors,cost1615254361201.
5. noscale: in odd discrepancy, use armedSIG[l] directly as scale instead of storing
   active?SIG[l]:1 in extra12qubits. On correct inactive-completed paths SIG[l]=0 so
   armedscale=1 anyway. Rareabortedpath remains unitarywithnonzeroscale. stream10
   offset3_select_noscale_compact PASS102400 Q1429,T645755,8phaseerrors,cost1612737419463.
   BM peak may nowbelowphasepeak; inspect hook. Keep compactname to prematerializeSIG
   before same-round scratch isfreed; omitting risks journal allocator replay conflicts.

Debug tooling diagnose_stream_failures.py now --tag --result --limit options. Original
stream3 default result is stream3_verify.json but actual legacyfile stream3_full_verify.json;
pass --result explicitly ifneeded. Offset4 fixtures saved stream_offset4_failures.npz,
forwarddebug log showed127/128 stoppedcode273(bank68), confirming consumption bug.

Older saved context follows (some canonical info superseded):

# Continue until verified overall score >=9.0 or user interrupts

User explicitly requires sustained experiments, no partial-progress final, and
submission.zip checkpoints at integer thresholds5,6,7,8,9. Milestones5,6,7saved.
Noagents. Don'ttouchuser'sindependent3bjobs. RootREADME/.gitignorepreexisting.
Goaltoolstatusoldblockedcannotresume; create_goalfailsunfinished; keepworking.

Canonical submission.zip / submission.summary.json is now7.8101966011:
- invsqrt candidates/rounded_f8_w27_16_8.kmb, results/rounded_f8_w27_16_8_verify.json
  PASS200540928shots, cost233305.8854,Q109,T1894mean,34rotations,worsterror.9646?
- dqi candidates/stream4_full.kmb, results/stream4_full_verify.json
  PASS102400shots,4phaseerrors,0valueerrors,Q1470,T643427,R1.209464,cost1681616324033.48.
- chemistry candidates/global500_factor_restart3_submission,
  results/global500_factor_restart3_verify.json, cost911197629.032,lambda_eff53.687102,
  Q833,T20375,11bits,J10,8circuits,36864shotsPASS.
Milestone7 is preservedat7.0239098891 (earlierinv298066,dqi1772681819628,
chemglobal500_coord11_c375_nudge cost930185660.98). Don'toverwriteoldmilestones.

Environment: source initial_experiments/3/env.sh fromroot. EXP_PYTHONkickmix313.
GPU Python /mnt/data01/miniconda3/envs/torch2.13/bin/python. RTX4090,16CPU,62GB.

## Most recent active experiments (inspectps/logsbeforererun)
- anneal500_q10_c350.log: anneal_angles.py16000steps,10bits,rank350,lr.02.
- grow600_float_fit.log:5000stepglobalPSDfitfromgrow600_q8_init, thenautomatic
  anneal600_q8_c320.log20000steps8bits,rank320. Theseaimat8bitmemorycostjump:
  beta8 forcesinnerQROAMk2=4,Q~530vs833;B600C320ifaccuratecouldapproach9.
- published450_fit andpublished500_fit18000steps, nearingcompletion:450~.0247
  (fails),500~.01893(float,notyetquantized),normbound151.
Allotherlongrunsbelowarefinishedorstopped.

## Inverse sqrt breakthrough
binary_centered_invsqrt.py now supports accum_bits list,carry_bits=1,rounded=True.
Horner accumulation runs on lowbits plusoneoverflowbit, thenintegrates signedcarry
intohighwordonce. Sign isknownfromcenteredlocalcoordinatefirstbit. Split[18,7,None]
withwidths28,17,7verifiedcost246037.74,Q112,T1944. Thenround-to-nearestpartialrows
throughcarry_inputallowedwidths27,16,8 withsplit[17,6,None],verifiedcost233305.8854.
Added /3/short_adder.py (copied /2 thenextendedoptionalcarry_in); defaultbehavior
unchanged. Exhaustivecarry-inaddertestsPASSallsmalloperandsandbothsigns.
fit_coefficients.cpp supports rounded,reachable,deep,fractionN,globalphase modes.
Reachablemantissa: oddk onlyevenm;evenk,m>max_xonlymultiples4. C++fittersfilterthese.
split_interpolation_fit.py grid322configs; rounded_interpolation_fit.json234configs
hasseveralpasses. Bestactualisrounded_f8_w27_16_8. C++deepsearchdidnotimprove.
Interception-onlyfraction1/2 (binary_centered extra_frac) grid235configs: nopasses
below233k. Globalphaseoffsetscan312nearpasscases: nopasses. Noextra_frac/globalphase
candidateinarchive. Codeparametersdefaultsoff,canonicalbinaryalreadyverified.
PotentialnewideaNOTimplemented: adaptive dyadic coefficienttable witharithmetic
coordinatefbase6/7 andfmax8address. Fitmergedwindowswherepossible,pruneQROMtrie;
couldreducetablecost512to~128withonly1-2morearithmeticrows. Uniformf7lowwidths
missedaccuracy, butfinewindows couldfitsharedcoarsecoordinate andmergeothers.

## DQI bugfix and diagnostics
Canonicalstream4 fixesunnecessaryround69guards in stream_capture.round_head:
if r andr<69. Round68doesn'tincrementdW(speciallastrounddoesn'tneedcounter),
butouroldstackguard69usedstalesignandabortedvalidweight35shots! Lastpushis
unconditionalswapandcanpreserveoldstackinDw; neitherstacknorbankguardneeded69.
Reproducedstream3's29failures via diagnose_stream_failures.py (unchangedhostfunctions).
results/stream3_failures.npz contains29S,mask,expectedphase fixtures.20None-reason
caseshadexactclassicalSIGbutactualcircuitabortedcode278at69. Remaining6stack
and3bank67 truncationsarelegitimateapproximations. Fixfullverifier4phaseerrors.
stream_model.py decode(trace=True) nowreturns3tuple(sig,why,history), default2tuple.
Failurediagnosticslogs stream3_failure_bm_tests/round_trace/sigma_values/stop_values.
stream_forney.prefix_zero nowhierarchicalchunks8 (extra2412T,lowerphaseworkspace).
stream4 isNONcompact,Q1470. Compact2 BM1452 butphase1465; +4896T meansmostlytie.
compact_wires.cpp lifetimecompactiononlysaved1Q. Proxrotationdidnotimprove.

## Chemistry best and work
Bestfloatshared500dictionary chem_global500_ridge_fit .0185,rank1,B500.
Project_global solvesPSDpositiveG,actualGusedforerror. Referencehasnegativenorm.00029765.
Pipeline: project_global --brot --columns --ridge -> offsets_cpu (~1s) ->
rotate_global500k -> rotate_qr -> nudge_gap -> constrained_factor -> build_global_chem
-> unchangedhostverify. packagervariablecircuitcountsupported.
Successfulquantized11coordinate global500_coordinate11 (2000steps,.019487full).
C375 trunc .0199055passed. QR2lambda_h2372.65. nudge_gap.py useserrorbudgetrank1PSD
updatealongprojectedoccupiedmatrixtolowerminSOSgap4.08->3.916,keepsC375,err.0199709.
Then constrained_factor.py fixedU,w/Woptimization loweredlambda54.805->53.945;
verifiedglobal500_coord11_c375_factor_submission915586030.75. Seed3QRrestart
loweredlambda_h2359.91->356.37 andverifiedcanonical911197629.03.
rotate_qr.py now--seed,--mix,--smooth; seeds1-3,mix1,smooth.002,8000stepssimilar356.
rotate_prox.py3000stepsdidnotbeatbaseline (candidatefallbacksourceexact).
constrained_factor secondrunlr3e-7cap.019985neverfoundfeasibleartifact, bestInfinity.

Global500coord10 6000steps ended .021666(fails). tree_screen commonGivens trees:
bestgeometrytreesinparameters/global500_tree_q10_0_tree.npy,etc. Tree10coordinate
fromrank330float ended .020218(fails). Rank350coord11 stopped~.02027(fails).
Floatrank330globalfit6000steps .01999, butrounding11 .0252; needangleawarefit.
anneal_angles.py newmethod: startcontinuousangles, exactKprojection (optionalrank),
graduallyincreaseperiodicpenaltyweightedbycoordinatecurvature. SavesbestHARD
candidate mainnpzwithK+U+words; _soft.npz Uiscontinuous,notvaliduntilquantized.
Bothhardnpzsneedproject_global forcolumns/factorizationbeforebuild.

Pruning/globalshapes: prune_global.pybackwardFrobselection somewhatbetterthan
tracegreedy, but450init .18. Global450backwardfit10000 ended .02463(fails).
Gradual500->495->...450 with2500stepsper5drop:495.01907,490.01967,485.02019,
480.02084,450.02565 (only495/490floatpass). Noquantizedtestyet,maybe13bits495.
Multi-ensemble block_fit30000steps R4B128C90 ended .0799,R2B256C180 .0449;
exactPSDblockprojection1500cycles .07816/.04349. Notcompetitive.

## Published THC data fetched (attribution retained)
remote_thc.py HTTP-range indexed9.2GBtar WITHOUT downloadingwholearchive.
Sourcehttps://zenodo.org/records/4248322 Lee etal PRXQuantum2,030305(2021).
Downloadedonly12MB: Li76orbitalM350,400,450,500,600.h5 +LICENSE/README/metadata,
in vendor/published_thc. CC-BY-4.0metadata,originalLICENSEretained,SOURCE.mdwritten.
Noexternalcodeexecuted. h5pyinstalled--no-depsin/tmp/quantumai-h5py forGPU Python.
UsePYTHONPATH=/tmp/quantumai-h5py:$PYTHONPATH toreadh5files. etaPp normalizedrows
isTHCU,MPQcoefficientkernel,U76onebodyeigenbasis. Preparedparameters/published_thcB_init.
Theirrawspansdon'tfitthebenchmark'sfixedsymmetryshift well: q12fullerrorsB350.085,
400.0726,450.0639,500.0602,600.0540. Adding76eigenvectorsShelped450to.0373,
notpassing. Lower-ridgefloatscreenalsofailed. FreshglobalPSDfitsstarted450/500
frompublishedUasnewbasins;logsactiveabove.500isapproachingoldbest,450worse.

## New decoder

Apache2 upstream vendor/uncompute revision64b5a13aa249f1360846657392e6c66817b26d0f.
import_published.py captures upstream BM+Psi, published_dqi.py/build_published.py
replace fixed masks with runtime masks and convert GF4331/public <->4179/host.
Published candidate verified 1.829622334585T, Q1843,T444030,5 phase errors.
published2 Q1837 had8 errors and slightly worse score. Canonical superseded it.

stream_capture.py: adapt public BM to consume old syndromes into its discrepancy,
then reuse the resulting zeros as backup bank storage. Uses prefix bank rotations,
active stop flag, one parked input for a vanished leading coefficient, and stack
overflow guard. stream_model.py classical model passed5000 sampled locators.
stream3 capture: Q1470,T291269 for forward/reverse BM. Tape and hook in parameters.
stream_forney.py converts public locator+retained last35 syndromes to host field,
reverses locator, forms reverse-sequence omega in place, and scans forward points
with transformed masks. Reverse error value e'=e*alpha^(71j), so masks compose
with alpha^(-71j). RecursiveForney kernel80T/point; transposed Toom middle product
forms omega35 in6174T, independently tested at n2,3,6,17,35.
build_stream.py --tape-prefix stream3 --phase --m4095 produces canonical decoder.
Full verifier PASS102400shots,0valueerrors,29phaseerrors,Q1473,T641109.

Important bugs fixed: freeze zip iterators in Em.swap_fn for journal replay;
unique9bit stop codes per guard; C++translator cat(empty file) must check EOF
before streambuf insertion, otherwise ostream fails and reverse pass is omitted.
Older stream/stream2 candidates are failed experiments. stream3_identity and
stream3_sigma tests now both PASS1024shots with no value/junk/phase errors.
debug_stream older logs were generated before empty-file fix, not current failures.

## Inverse sqrt

Centered and binary-centered polynomial interpolation: centered_invsqrt.py,
binary_invsqrt.py,binary_centered_invsqrt.py. fit_coefficients.cpp searches nearby
integer coefficients exhaustively; fit_interpolation.py records results.
Two centered candidates passed200540928shots ~301799 and302355 cost.
Binarycentered f8 widths(28,17,7) passed full verifier, score298066.97754,Q119.

## Chemistry (current promising direction)

Previous R30 Adam fits plateaued .0207; lbfgs_chem.py doubleprecision polish reached
.0192317557 (float) after300 outer steps/~47min. Quantization still needed.
block_chem.py fixed-vector PSD block coordinate descent R30 plateaued .022178.
constrained_chem.py augmented Lagrangian on actual lambda_eff and error/gap constraints.
R35C45 final verified modest improvement, already canonical. R40 worse resources.
Merge experiments R12B128C90 reached .01846 float but normalization huge; constrained
fit still Lambda610,gap30,poor. Alternating global spectral/local dictionary projections
for81x25x4 failed .24 error (including weighted variant).

reference_spectrum.py cached full symmetric2926 basis eigensystem. Rank lower bounds:
292=>.019206,300=>.016697,350=>.007552,400=>.004125,500=>.001300.
Reference has tiny negative eigenvalues with total norm~.00029765.

global_dictionary.py greedily picks a global rank-one dictionary from existing
vectors; global_fit.py eliminates the optimal linear coefficient matrix each
step and optimizes U only. R1 B500: init error.34577 -> .01886 after5000steps
ridge1e-7; Gram normalization bound grew to~403. Ridge1e-5 continuation improves
conditioning and keeps bound~159 while error~.0192, still running10000steps.
R1 B400 greedy-pruned from500 started .3888, currently ~.037 with12000steps planned.
global_quant_fit.py STE Givens angles+same exact coefficient solve; brot12 run8000steps
currently best~.0214, not yet valid.

CRITICAL CURRENT INVESTIGATION: global_fit/global_quant_fit optimize against the
reference's slightly indefinite G, yielding indefinite K. Naively clipping K in
Euclidean coordinates distorts physical error (.019 -> .068). Correct PSD projection
is in the feature Gram metric, or simply fit G_positive=(V*max(eigenvalues,0))V.T,
whose difference from G is only .00029765. project_global.py has just been changed
in its --ridge branch to solve using G_positive, then evaluate against actual G.
Three screens currently running in results/global500_psd_q12/13/14.log using
parameters/global500_ridge_snapshot.npz. Inspect these before modifying long fits.
The ongoing global_fit and global_quant_fit processes still use original indefiniteG;
future runs should use positiveG if screens confirm this repairs physical accuracy.
Keep actual-G checks for saved candidates. Do not claim modeled fits are verified.

Next pipeline once a quantized global dictionary passes:
project_global.py --fit... --brot... --ridge1e-5 --columns... exports a fit.
optimize_offsets.py minimizes convex SOS gap for fixed X using L-BFGS, and moves to
accepted gap3.916 if necessary. It is written but not yet run.
rotate_global.py + rotate_columns.cpp sample Jacobi column rotations to minimize
LCU norm, preserving XX and one-body correction. Written/compiled, not yet run.
Then export submission.npz with U shape(1,B,76), W=lowercase fitw,
H2_I=offset+sum_b w, old ignored H1_I/Hsym_lower, tree_edges; encode onebody
eigenvectors of h1'=h0-sum(offset*X), concatenate two-body angle_words.
Build ANGLE circuits with chemistry_qrom.build(delta,brot), batch size suitable
for128Q cap (brot12=>9,11=>10,13=>8). Existing build_chem.py hardcodes12/9,
so write a generic builder. package.py currently asserts12files (=9ANGLE); relax
only to verifier-derived ceil(75/batch_size) if using a new batch size.

Global dictionary resource possibility: A=76+B only~576 vs2176. If normalization
can stay near160-250 and accuracy passes, chemistry cost could fall substantially.
Need actual host verification, onebody quantization, and normalization optimization.

Inspect ps/logs before restarting jobs: server restarts have occurred, but long
processes usually finish on disk. Periodic user updates <=60sec; keep working.
