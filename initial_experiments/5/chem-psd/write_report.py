"""Bind the completed numerical investigation to its preserved artifacts."""
import csv,hashlib,json,platform,sys
from pathlib import Path
import numpy as np
import scipy,torch

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
def read(name):return json.loads((HERE/f'results/{name}.json').read_text())
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()
def main():
    native=read('incumbent_native');resources=read('resources')['rows'];A=read('experiment_a')
    spectra=read('span_diagnostics');actual=read('actual_resources')['records']
    phys=[]
    for p in sorted((HERE/'results').glob('screen_*_physical.json')):
        z=json.loads(p.read_text());phys.append((p.name.removesuffix('_physical.json'),z))
    assert len(phys)==8 and all(not z['passes_physical_gates'] for _,z in phys)
    assert len(actual)==len(phys) and native['status']=='PASS'
    assert native['score']==708547388.9553173 and native['submission_sha256']=='ebf8e357a1b693fc8589d5c48a5d6fa10343dfa5944724d8f7a77b8edc2157f7'
    decoded=[p.name for p in (HERE/'results').glob('screen_*_decoded_psd.json')]
    assert len(decoded)==8
    for name in decoded:assert json.loads((HERE/'results'/name).read_text())['diagnostics']['dual_error_lower_bound']>.02
    baselinesource=ROOT/'initial_experiments/4/codex/candidates/norm_cycle2_signed'
    threshold=read('claude_r2b351_threshold_capacity')['ranks'];demonstrated=min((r for r in threshold if r['direct_error_2']<=.02),key=lambda r:r['C'])
    row189=next(r for r in resources if (r['R'],r['B'],r['C'],r['brot'])==(2,351,189,8))
    selection=dict(status='COMPLETE_NO_REPLACEMENT',selected_source=str(baselinesource),native_cost=native['score'],
        native_T=native['toffoli_count'],native_Q=native['num_qubits'],native_lambda_eff=native['lambda_eff'],chemistry_improvement=1.,
        target_1p3_reached=False,smallest_demonstrated_continuous_C=demonstrated['C'],
        continuous_error_2_at_that_C=demonstrated['direct_error_2'],continuous_is_not_a_submission=True,
        all_eight_new_native_word_screens_fail=True,no_new_full_native_submission_PASS=True)
    (HERE/'selection.json').write_text(json.dumps(selection,indent=2)+'\n')
    lines=['# Round 5 continued: block-PSD chemistry factorization','',
      '**Complete; no replacement. The 1.3x target was not reached.** The selected round-4 chemistry remains R1/B719/C350, eight bits, T=30,504, Q=478, lambda_eff=48.59417497035283, cost **708,547,388.9553173**. Fresh unchanged-host verification passes 50,880 shots with an independent salt. No numerical search or experiment process remains running.','',
      'The main positive result is a continuous R2/B351 factorization at **C189, 19.978781 mHa**, with conditional T/Q **22,516/461** (1.4047x at incumbent normalization). This is not an integer or scored submission. The exact rounded dictionary has a **140.196114 mHa joint-span lower bound** and a roughly 145.48 mHa block-PSD optimum. Increased columns help the continuous representation but cannot rescue those decoded vectors.','',
      'The best newly constructed exact eight-bit dictionary is the residual-augmented **R2/B375** pivot: full-PSD error **51.459643 mHa**, explicit dual lower bound above 51 mHa. Its actual rank-C192 native word screen is **53.181669 mHa**, also invalid. These are fixed-dictionary failures, not global limits on these shapes or on chemistry factorization.','',
      '## Scope and calibration','',
      'All writes are confined to this directory. Inputs are the unchanged host, the round-4 accepted candidate, completed `chem-blocks/RESULTS.md`, `chem-blocks/jad_search/RESULTS.md`, both resource-frontier implementations, and completed `chem-qfit/RESULTS.md` and its referenced round-4 checkpoints. No concurrently active round-5 agent directory was used.','',
      '`core.py` constructs the reference directly through native `load_inputs`: float32 FCIDUMP integrals, native float32 fixed S, and fixed beta_1. It uses `P[b,pq]=u[b,p]*u[b,q]*fac[pq]`, with fac=sqrt(2) off diagonal. The signed target is never clipped to PSD. Its negative spectral component is 0.297649784 mHa. Direct four-index, pair-space and native incumbent errors agree to about 2e-16 Ha.','',
      '`results/calibration.json` records the native control and the analytic single-block PSD control; `envelope_gradient_control.json`, `hessian_control.json`, `offset_controls.json` record finite-difference checks. The final Hessian implementation uses symmetric-pair reshuffling without a dense N^4 intermediate. Its directional error is about 8.7e-12. Native gap/Lambda formulas match the incumbent within 1e-8; fixed shifts are never fitted.','',
      '## Fixed-dictionary convex optimization','',
      'For each block, QR of P_r.T gives an orthonormal row dictionary Q_r and triangular R_r. Optimize S_r>=0 in `Hhat=sum Q_r.T S_r Q_r`; atom covariance is `K_r=R_r^-1 S_r R_r^-T`. This avoids projecting PSD in an incorrect atom metric. The solver caches `Q_r Q_s.T` and `Q_r H Q_r.T`. FISTA uses the safe Lipschitz bound R, adaptive restart, and full-capacity eigenvalue projection. Rank-C comparisons use exact cyclic block minimization by positive top-C spectral projection.','',
      'Loss contractions are checked against direct pair-space tensors and explicit factor weights. All generator modes can activate; extra zero factor columns do not remain artificially inactive. Rank-constrained fits are nonconvex feasible upper bounds, with multiple warm starts and rank-aware projected-gradient diagnostics. They are not global rank-optimality claims.','',
      'For E=Hhat-H, the dual tensor `D=E+delta*Pi_union` is feasible when every `Q_r D Q_r.T>=0`. It gives squared-error lower bound `-||D||^2-2<D,H>`. During optimization a conservative cached formula uses an upper bound on union rank. `certify.py` independently reconstructs dense final dual tensors, checks positive feasibility slack and direct weight reconstruction. Certificates use floating-point safety slack, not interval arithmetic. The reported excluded margins are far larger than numerical reconstruction errors.','',
      '### Experiment A: prior dictionaries','',
      '| Dictionary | C88 | C128 | C160 | C175 | C192 | C224 | Full PSD feasible / dual lower (mHa) |',
      '|---|---:|---:|---:|---:|---:|---:|---:|']
    for rec in A:
        curve={r['C']:r['direct_error_2']*1000 for r in rec['rank_curve']};d=rec['full_PSD']
        cells=[f'{curve[c]:.3f}' if c in curve else '—' for c in [88,128,160,175,192,224]]
        lines.append('| '+rec['name']+' | '+' | '.join(cells)+f" | {1000*d['direct_error_2']:.3f} / {1000*d['dual_error_lower_bound']:.3f} |")
    lines += ['',
      'Extra C cannot rescue any of these exact prior dictionaries. The common-span relaxation further distinguishes their obstructions:','',
      '| Dictionary | Unrestricted joint-span PSD lower bound (mHa) | Interpretation |','|---|---:|---|']
    meanings=['Union span already excludes 20 mHa','Joint span permits <20; within-block PSD geometry excludes it','Joint span permits <20; within-block PSD geometry excludes it','Continuous full block-PSD fit passes','Quantized union span itself excludes 20 mHa','Integer moves improve the span but do not pass','Residual additions improve it further but do not pass']
    for r,meaning in zip(spectra['dictionaries'],meanings):
        lines.append(f"| {Path(r['input']).stem} | {1000*r['union_span_PSD_lower_bound']:.3f} | {meaning} |")
    lines += ['',
      'The joint-span relaxation permits cross-block coefficients and therefore is only a lower bound, not an achievable block representation. SVD dimensions and smallest retained singular values are recorded in `span_diagnostics.json`. Native unrestricted PSD rank must be at least 290; the resource CSV explicitly excludes shapes failing this necessary test.','',
      '### Stronger continuous basin and column threshold','',
      'The completed qfit report identifies `4/claude/chem/L_r2b351c175.npz`, whose continuous error is about 21.115 mHa. Its fixed-dictionary full-PSD interval is approximately **18.428–18.538 mHa**. Cold truncation of the relaxed covariance sometimes underperformed the original rank-C weights; both starts were retained and refined.','',
      '| C | Best demonstrated continuous error (mHa) |','|---:|---:|']
    best={}
    for tag in ['claude_r2b351_cont','claude_r2b351_warm_cont','claude_r2b351_threshold']:
        for r in read(tag+'_capacity')['ranks']:best[r['C']]=min(best.get(r['C'],float('inf')),1000*r['direct_error_2'])
    for C in [160,175,184,188,189,190,191,192,200,208,224]:lines.append(f'| {C} | {best[C]:.6f} |')
    lines += ['',
      '**C189 is the smallest demonstrated passing continuous rank.** C188 reached 20.037224 mHa in the tested refinement; that is not a certificate that C188 or a smaller C cannot pass with another rank optimizer. No physically viable C exists for the exact decoded dictionary, because its full-rank relaxation and common-span bound both fail. `claude_r2b351_*_capacity.json` preserves complete curves and warm/cold evidence.','',
      '![Rank-versus-error curves](results/rank_curves.png)','',
      '## Experiment B: different geometry and exact lattice moves','',
      'Two new non-JAD initializations were tested: (1) an accepted-covariance graph with high-coupling vectors shared across ensembles; (2) individual target-generator eigenvectors, spatial-participation grouping, and nonorthogonal simultaneous projector pivots. The graph R4/B224 start shares 59 vectors and covers all 719 incumbent atoms in its union. It nevertheless has a very poor block-PSD optimum, showing that coverage alone is insufficient.','',
      '`geometry_trust.py` improves unit vectors by PSD variable projection with curvature history and a physical preconditioner. A feasible-covariance Armijo test guarantees physical descent before the next covariance solve. An ordinary inexact-inner L-BFGS control terminated prematurely and is retained as a numerical optimization failure. The trust searches use about 100 accepted geometry steps, thousands of convex covariance updates, and independent long final PSD solves. GPU warm runs were continued on CPU after a timing control; exact saved checkpoints bind the reported results.','',
      '| Geometry | Initial full-PSD error (mHa) | Final feasible / dual lower (mHa) | Accepted steps |','|---|---:|---:|---:|']
    for tag in ['trust_r4b224_cpu','trust_prior_r2b351_cpu','trust_eigen_r2b351','trust_graph_r4b224']:
        z=read(tag);d=z['final'];lines.append(f"| {tag} | {1000*z['initial']['direct_error_2']:.3f} | {1000*d['direct_error_2']:.3f} / {1000*d['dual_error_lower_bound']:.3f} | {z['outer_steps']} |")
    lines += ['',
      'The first row continues the preserved GPU warm geometry (the original R4/B224 control was 27.564 mHa). Geometry budgets and exact stop reasons are in `geometry_budgets.json`. These terminal fixed-dictionary certificates do not certify global geometry minima.','',
      '`integer_trust.py` then changes coupled Givens words across whole vectors, screens exact row tensor increments, tests combined decoded moves, and refits unrestricted PSD covariances. A second run uses validated per-vector physical Newton curvature. A bounded projected-trial lookahead pilot is also retained; its sparse logged history is explicitly marked as such. Every accepted move uses the true signed physical residual.','',
      '| Exact eight-bit geometry search | Initial / final full-PSD error (mHa) |','|---|---:|']
    for tag in ['integer_claude_r2b351','integer_newton_r2b351','integer_projected_r2b351']:
        z=read(tag);lines.append(f"| {tag} | {1000*z['initial']['direct_error_2']:.3f} / {1000*z['final']['direct_error_2']:.3f} |")
    lines += ['',
      'One additional resource-screened family, R2/B375, appends 24 exactly encoded atoms per block to an integer checkpoint. Four additions are shared high-leverage residual vectors and twenty are block-specific. The pool combines incumbent vectors and individual eigenvectors of positive residual modes; the actual optimization still retains the full signed target. Adding these decoded nonorthogonal vectors lowers full-PSD error from 55.831479 to 51.459643 mHa. A convincing dual bound still excludes 20 mHa for this exact augmented dictionary.','',
      '![Geometry convergence](results/geometry_curves.png)','',
      '## Experiment C: decoded physical and normalization checks','',
      'Each row below restores rank C, decodes exact signed eight-bit words, refits weights, optimizes the actual SOS gap with all identity offsets, preserves native symmetry shifts, and evaluates native physics. Selected cases also receive column-gauge L1 optimization that preserves the combined weight/identity Gram matrix. Since every two-body dictionary fails even its relaxed PSD gate, further one-body polishing and constrained task-cost optimization cannot rescue these exact geometries. One-body results below are measured screens, not one-body impossibility claims.','',
      '| Screen | C | Continuous / decoded error2 (mHa) | Native error1 (mHa) | Lambda | Gap (Ha) | lambda_eff |','|---|---:|---:|---:|---:|---:|---:|']
    for tag,z in phys:
        m=z['native_metrics'];lines.append(f"| {tag.removeprefix('screen_')} | {z['rank']} | {1000*z['continuous']['direct_error_2']:.3f} / {1000*m['error_2']:.3f} | {1000*m['error_1']:.3f} | {m['Lambda']:.3f} | {m['E_gap']:.6f} | {m['lambda_eff']:.3f} |")
    lines += ['',
      'Every new row fails both error gates; all pass the strict raw-Lambda/gap floors and gap-domain test. Independent decoded full-PSD refits in `decoded_psd_summary.json` have lower bounds above 20 mHa for all eight screened dictionaries. Rank-C physical fits can be worse than those relaxations because their optimizer is nonconvex. This difference is retained rather than mislabelled as a quantization lower bound.','',
      '## Complete host resource accounting','',
      f'`resources.csv` contains {len(resources)} shape/precision points, with unmodified host chi10, automatic preparation workspace choices, explicit signed-ANGLE bounds and necessary signed-spectrum rank screens. Host preparation/Givens arithmetic is exact. ANGLE bounds are constructive worst cases; synthetic dense tables attain them and pass arbitrary-target XOR/phase/scratch tests under two independent seeds.','',
      '| Shape (eight bits) | T | Q | Conditional improvement at incumbent lambda | lambda_eff ceiling for 1.3x |','|---|---:|---:|---:|---:|']
    for shape in [(4,180,88),(4,180,128),(4,180,160),(4,224,128),(2,351,175),(2,351,189),(2,351,200),(2,351,224),(2,375,192)]:
        r=next(r for r in resources if (r['R'],r['B'],r['C'],r['brot'])==(*shape,8));lines.append(f"| R{shape[0]}/B{shape[1]}/C{shape[2]} | {r['T_total_bound']:,} | {r['Q_total_bound']} | {r['conditional_improvement']:.4f}x | {r['lambda_max_for_1p3']:.4f} |")
    lines += ['',
      '`resource_actual.py` also compiles **all eleven actual ANGLE circuits for each of the eight physical screens**, retained under `resource_circuits/`, and applies the exact host persistent-plus-maximum-temporary resource formula. Actual compiled T/Q agree with the bounds for these tables. The best two new physical products are still above 1.06 billion at their measured normalization, and every model is physically invalid. These products are not accepted task scores. No full native PASS is claimed for an invalid family.','',
      'Only the incumbent has complete native chemistry verification: `results/incumbent_native.json`, canonical digest `ebf8e357a1b693fc8589d5c48a5d6fa10343dfa5944724d8f7a77b8edc2157f7`, 50,880 shots. Its original twelve files remain read-only. There is no new viable candidate requiring a new accepted native submission.','',
      '## Artifacts, limits, and reproduction','',
      '`selection.json` records no replacement. `artifact_manifest.json` hashes all research code, exact NPZs, JSON evidence, plots, actual circuit families, source inputs and baseline files. `parameters/` preserves continuous covariances, explicit rank-C weights, integer words, decoded fits, offsets and normalization gauge inputs/outputs. CSV, PNG and PDF curves are standalone artifacts.','',
      'The numerical conclusions concern the saved dictionaries and finite geometry searches. C189 is an empirical constructive threshold, not a proven minimum rank. Floating trajectories can vary by device and scheduling; physical metrics and certificates are recomputed from the exact hashed checkpoints. The failed initializations, premature L-BFGS control, finite lookahead pilot and CPU restarts are preserved. None establishes a global multiblock, eight-bit, normalization or chemistry limit.','',
      'Run commands from this directory after `source env.sh`; set `PSD_DEVICE=cpu` for CPU reproduction.','',
      '```bash',
      '$PSD_PY calibrate.py',
      '$PSD_PY resources.py',
      '$PSD_PY resource_controls.py',
      '$PSD_PY experiment_a.py',
      '$PSD_PY dictionaries.py --R 4 --B 224 --shared 59 --strategy graph --tag overlap_graph_r4b224',
      '$PSD_PY dictionaries.py --R 2 --B 351 --shared 32 --strategy eigenatoms --tag eigen_r2b351',
      '$PSD_PY solve.py --input parameters/eigen_r2b351.npz --out parameters/eigen_r2b351_psd.npz --steps 2400',
      '$PSD_PY geometry_trust.py --input parameters/eigen_r2b351_psd.npz --tag trust_eigen_r2b351 --steps 100 --inner 30 --radius .05',
      '$PSD_PY capacity_scan.py --input parameters/claude_r2b351_psd.npz --tag reproduce_rank_curve --ranks 175 184 188 189 192 200 224 --full-steps 0',
      '$PSD_PY integer_trust.py --input parameters/claude_r2b351_q8_psd.npz --tag reproduce_integer --steps 100 --inner 30',
      '$PSD_PY integer_trust.py --input parameters/integer_newton_seed.npz --tag reproduce_newton --steps 100 --inner 30 --newton',
      '$PSD_PY residual_pivot.py',
      '$PSD_PY physical.py --input parameters/residual_r2b375_psd.npz --tag reproduce_residual --rank 192 --steps 1000 --jacobi 8 --build',
      '$PSD_PY quantized_certificates.py',
      '$PSD_PY resource_actual.py',
      'PYTHONPATH="$PSD_DIR/.plotdeps:$PYTHONPATH" $PSD_PY plot_results.py',
      '```','',
      'Host DFTHC imports retain the Apache-2.0 license at `from_host/verifiers/LICENSE.txt`. Signed encoding/oracles reuse read-only round-4 `angle_gauge_probe`, `build_selectswap_chem` and `monomial_unlookup`, plus completed chem-blocks builders. Exact integer-row increments and offset/gauge semantics are attributed in the local source. The plotting dependencies were installed only into `.plotdeps/`, with licenses included there; no host or global environment was changed.']
    (HERE/'RESULTS.md').write_text('\n'.join(lines)+'\n')
    # Exclude vendor packages/caches and this self-referential manifest.
    files={str(p.relative_to(HERE)):sha(p) for p in sorted(HERE.rglob('*')) if p.is_file() and not any(part.startswith('.') or part=='__pycache__' for part in p.relative_to(HERE).parts) and p.name!='artifact_manifest.json'}
    inputs=[ROOT/p for p in ['from_host/verifiers/dfthcsos/verify.py','from_host/verifiers/dfthcsos/scoreinfra/dfthc_math.py',
        'from_host/verifiers/dfthcsos/scoreinfra/costs_prepare.py','from_host/verifiers/dfthcsos/scoreinfra/reference.py',
        'from_host/verifiers/dfthcsos/FeMoCo76o113e.fcidump','from_host/verifiers/dfthcsos/FeMoCo76o113e.fixed_shifts.npz',
        'initial_experiments/4/claude/chem/L_r2b351c175.npz','initial_experiments/5/chem-blocks/RESULTS.md',
        'initial_experiments/5/chem-blocks/jad_search/RESULTS.md','initial_experiments/5/chem-qfit/RESULTS.md',
        'initial_experiments/5/chem-blocks/frontier.py','initial_experiments/5/chem-blocks/jad_search/resources.py',
        'initial_experiments/5/chem-blocks/signed_oracle.py','initial_experiments/5/chem-blocks/oracle_controls.py',
        'initial_experiments/4/codex/angle_gauge_probe.py','initial_experiments/4/codex/build_selectswap_chem.py',
        'initial_experiments/4/codex/monomial_unlookup.py','initial_experiments/4/codex/jacobi_l1.py',
        'initial_experiments/4/codex/jacobi_l1',
        'initial_experiments/5/chem-blocks/parameters/r4b180c88_pca_lbfgs.npz',
        'initial_experiments/5/chem-blocks/parameters/r2b351c175_lbfgs.npz',
        'initial_experiments/5/chem-blocks/jad_search/candidates/r4_b224_c88_float.npz']]
    inputs += [p for p in baselinesource.iterdir() if p.is_file()]
    manifest=dict(status='COMPLETE_NO_REPLACEMENT',runtime=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,torch=torch.__version__),
        baseline_canonical_digest=native['submission_sha256'],artifacts=files,readonly_input_sha256={str(p.relative_to(ROOT)):sha(p) for p in inputs})
    (HERE/'artifact_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(artifacts=len(files),physical_screens=len(phys),resource_rows=len(resources),selection=selection)),flush=True)
if __name__=='__main__':main()
