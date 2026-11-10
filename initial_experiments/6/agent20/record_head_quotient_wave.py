"""Bind this completed research wave without promoting a submission."""
from pathlib import Path
import hashlib,json
HERE=Path(__file__).resolve().parent
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def record(p):
 p=HERE/p;assert p.is_file(),p
 return dict(path=str(p.relative_to(HERE)),sha256=digest(p),bytes=p.stat().st_size)

sources=['CONTINUATION.md','RESULTS.md','record_head_quotient_wave.py',
 'shared_gap_quotient_peel.py','check_shared_gap_quotient_peel.py','results/code/shared_gap_quotient_peel_v1.py',
 'continuant_frame_affine.py','continuant_head_frame.py','check_continuant_head_frame.py','check_continuant_frame_oracle.py','profile_continuant_frame.py','profile_direct_head_parameters.py',
 'results/code/continuant_head_frame_tree_v1.py','results/code/continuant_head_frame_cached_resultant_v2.py',
 'compact_lookahead_parameters.py','check_compact_lookahead.py','counter_quotient_jacobi.py','check_counter_quotient_jacobi.py',
 'high_copy_jacobi.py','check_high_copy_jacobi.py','profile_high_copy_jacobi.py',
 'prefix_high_copy_jacobi.py','check_prefix_high_copy.py','profile_prefix_high_copy.py','profile_prefix_counter.py','prefix_counter_model.cpp','results/code/prefix_counter_model_v2.cpp','prefix_counter_oracle.py','prefix_counter_reference.py','check_prefix_counter_oracle.py',
 'mixed_quotient_jacobi.py','check_mixed_quotient_jacobi.py','mixed_prefix_oracle.py','mixed_prefix_reference.py','check_mixed_prefix_oracle.py','mixed_prefix_model.cpp','check_mixed_prefix_model.py',
 'routed_quantum_control.py','profile_routed_H.py','quadratic_multiplier_search.cpp','results/code/quadratic_multiplier_search_v1.cpp',
 'chem_positive_integer_refine.py','chem_native_factor_admm.py','results/code/chem_native_factor_admm_v1.py','results/code/chem_native_factor_admm_v2_weighted.py',
 'triple_gap_frame.py','check_triple_gap_frame.py']
reports=['results/shared_gap_quotient_checks.json','results/routed_shared_gap_quotient_checks.json',
 'results/continuant_frame_affine_formulas.json','results/continuant_head_frame_tree_v1_checks.json','results/continuant_head_frame_cached_resultant_v2_checks.json','results/continuant_head_frame_checks.json','results/continuant_frame_direct_block_profiles.json','results/direct_head_parameter_profiles.json',
 'continuant_frame_b2/results/cfraction_n14_m127_agent20-continuant-frame-small-a.json',
 'results/compact_lookahead_checks.json','results/compact_lookahead_full_profile.json','compact_lookahead/results/cfraction_n14_m127_agent20-compact-lookahead-small-a.json',
 'results/counter_quotient_jacobi_checks.json','results/high_copy_jacobi_b3_checks.json','results/high_copy_jacobi_block_profiles.json','results/prefix_high_copy_b3_checks.json','results/prefix_high_copy_block_profiles.json',
 'results/prefix_counter_profiles.json','results/prefix_counter_limit935.json','results/prefix_counter_model_b6_t18_v2.json','results/prefix_counter_model_b6_t20_s3.json','results/prefix_counter_small_b3_t3_checks.json',
 'results/mixed_quotient_jacobi_odd_checks.json','results/mixed_quotient_jacobi_even_checks.json','results/mixed_quotient_jacobi_p18_checks.json','results/mixed_prefix_unrouted_small_checks.json','results/routed_mixed_prefix_small_checks.json','results/mixed_prefix_b6_t18_l935_profile.json','results/routed_mixed_prefix_profile.json','results/mixed_prefix_model_b6_t18.json','results/mixed_prefix_model_agreement.json','results/mixed_prefix_cost_screen.json','results/routed_H_core_profile.json',
 'results/positive_native_integer_R4B240.json','results/native_factor_admm_R4B240.json','results/native_factor_admm_weighted_R4B240.json','results/native_global_atom_admm_R4B240.json','results/triple_gap_frame_checks.json']

# Native bytes in older reports remain tied to their archived source versions.
patterns=['shared_gap_quotient_n*.kmb','routed_shared_gap_quotient_n*.kmb','continuant_head_frame_b*.kmb','continuant_head_frame_tree_v1_b*.kmb','continuant_head_frame_cached_resultant_v2_b*.kmb','compact_lookahead_n*.kmb','counter_quotient_jacobi_n*.kmb','high_copy_jacobi_b3_n*.kmb','prefix_high_copy_b3_n*.kmb','mixed_quotient_jacobi_n*.kmb','prefix_counter_n14_m127_b3_t3_s2_l935.kmb','mixed_prefix_n14_m127_b2_t2_l935.kmb','routed_mixed_prefix_n14_m127_b2_t2_l935.kmb','triple_gap_frame.kmb']
candidates=sorted({p for pattern in patterns for p in (HERE/'candidates').glob(pattern)})
accepted=digest(HERE/'best/submission.zip');milestone=digest(HERE/'milestones/10/submission.zip')
assert accepted=='df518f8d20a0b4d4ea72700280f689bc17d0685471e0d1ddb2cc998b1c60a52e'
assert milestone=='59bf461be714e0a4786fb0cc825d8b105622fec35ce193acd44ab1298004c2a0'
old=json.loads((HERE/'results/repaired_H_build.json').read_text())
assert all(digest(Path(p))==sha for p,sha in old['source_sha256'].items())
for name,source in [('native_factor_admm_R4B240','results/code/chem_native_factor_admm_v1.py'),('native_factor_admm_weighted_R4B240','results/code/chem_native_factor_admm_v2_weighted.py'),('native_global_atom_admm_R4B240','chem_native_factor_admm.py')]:
 r=json.loads((HERE/'results'/f'{name}.json').read_text());assert r['finished'] and r['source_sha256']==digest(HERE/source),(name,'source/terminal state')
agreement=json.loads((HERE/'results/mixed_prefix_model_agreement.json').read_text());assert all(digest(HERE/p)==sha for p,sha in agreement['source_sha256'].items())
parameters=['parameters/positive_native_integer_R4B240.npz','parameters/native_factor_admm_R4B240_best.npz','parameters/native_factor_admm_R4B240_last.npz','parameters/native_factor_admm_weighted_R4B240_best.npz','parameters/native_factor_admm_weighted_R4B240_last.npz','parameters/native_global_atom_admm_R4B240_best.npz','parameters/native_global_atom_admm_R4B240_last.npz','parameters/quadratic_multiplier_best.txt','parameters/forney_norm_pair_quadratic_best.txt']
out=dict(status='COMPLETED_RESEARCH_NO_PROMOTION',accepted_score=10.128965274112238,independent_combination_score=10.144446072819372,accepted_archive_sha256=accepted,first_crossed_10_archive_sha256=milestone,goal20_complete=False,new_integer_milestones=[],jobs_live=False,sources=[record(p) for p in sources],reports=[record(p) for p in reports],native_artifacts=[record(p.relative_to(HERE)) for p in candidates],parameters=[record(p) for p in parameters],finite_search_logs=[record('results/quadratic_multiplier_search.log'),record('results/quadratic_forney_search.log')],exceptions=['prefix_counter_model_b6_t18.json is obsolete; corrected completion model report is _v2','direct_head_parameter_profiles.json is partial: b4/cache4 explicitly aborted','component results and complete small oracles are not full-host benchmark acceptance','new global-atom chemistry control completed150 steps without improving56.729335mHa init'])
path=HERE/'results/head_quotient_research_manifest.json';path.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(status=out['status'],sources=len(sources),reports=len(reports),native_artifacts=len(candidates),accepted_archive_sha256=accepted,jobs_live=False)),flush=True)
