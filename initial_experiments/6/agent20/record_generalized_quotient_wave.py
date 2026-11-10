"""Hash the closed generalized-quotient research wave, without promotion."""
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def record(path):
    assert path.is_file(), path
    return dict(path=str(path.relative_to(HERE)), bytes=path.stat().st_size, sha256=digest(path))


sources = ['CONTINUATION.md', 'RESULTS.md', Path(__file__).name,
           'triple_gap_head.py', 'check_triple_gap_head.py', 'conditional_ccz_builder.py',
           'triple_gap_frame.py', 'check_triple_gap_frame.py',
           'triple_gap_block.py', 'check_triple_gap_block.py', 'profile_triple_gap_block.py', 'profile_triple_extraction.py',
           'generalized_gap_quotient.py', 'check_generalized_gap_quotient.py', 'profile_generalized_gap_quotient.py',
           'quotient_cofactor_frame.py', 'check_quotient_cofactor_frame.py',
           'prefix_lookahead_quotient_block.py', 'check_prefix_lookahead_quotient.py', 'profile_prefix_lookahead_quotient.py',
           'chem_native_factor_admm.py', 'chem_positive_dictionary_fit.py', 'chem_fixed_dictionary_union_floor.py']
sources += [str(p.relative_to(HERE)) for p in (HERE/'results/code').glob('*triple*')]
sources += [str(p.relative_to(HERE)) for p in (HERE/'results/code').glob('*generalized_gap_quotient*')]
sources += ['results/code/chem_native_factor_admm_v3_global.py', 'results/code/chem_positive_dictionary_fit_v1.py']
reports = ['triple_gap_frame_checks', 'triple_gap_head_checks', 'triple_gap_head_selected_v1_checks',
           'triple_gap_head_shared_encoded_v2_checks', 'triple_gap_block_checks', 'triple_gap_block_profiles',
           'triple_extraction_profile', 'general_gap_quotient_n3_7_checks', 'general_gap_quotient_n17_35_checks',
           'general_gap_quotient2_n2_3_7_checks', 'general_gap_quotient3_n2_3_7_checks',
           'general_gap_quotient3_sparse_n3_7_checks', 'general_gap_quotient3_profiles',
           'quotient_cofactor_checks', 'prefix_lookahead_quotient_b6_checks', 'prefix_lookahead_quotient_profiles',
           'native_global_atom_admm_R4B240', 'native_global_atom_weighted_inner20',
           'positive_R11B64C32', 'positive_R11B64C32_union_floor', 'generator_rank_trace_invariant']
bindings = [('native_global_atom_admm_R4B240', 'results/code/chem_native_factor_admm_v3_global.py'),
            ('native_global_atom_weighted_inner20', 'chem_native_factor_admm.py'),
            ('positive_R11B64C32', 'chem_positive_dictionary_fit.py')]
for report, source in bindings:
    r = json.loads((HERE/'results'/f'{report}.json').read_text())
    assert r['finished'] and r['source_sha256'] == digest(HERE/source), (report, source)

patterns = ['triple_gap_*.kmb', 'general_gap_quotient*.kmb', 'quotient_cofactor_n*.kmb',
            'prefix_lookahead_quotient_*.kmb']
native = sorted({p for pattern in patterns for p in (HERE/'candidates').glob(pattern)})
params = sorted({p for pattern in ['native_global_atom_weighted_inner20*.npz', 'positive_R11B64C32*.npz']
                 for p in (HERE/'parameters').glob(pattern)})
accepted = digest(HERE/'best/submission.zip')
milestone = digest(HERE/'milestones/10/submission.zip')
assert accepted == 'df518f8d20a0b4d4ea72700280f689bc17d0685471e0d1ddb2cc998b1c60a52e'
assert milestone == '59bf461be714e0a4786fb0cc825d8b105622fec35ce193acd44ab1298004c2a0'
out = dict(status='COMPLETED_RESEARCH_NO_PROMOTION', accepted_score=10.128965274112238,
           accepted_archive_sha256=accepted, first_crossed10_archive_sha256=milestone,
           new_integer_milestones=[], goal20_complete=False,
           sources=[record(HERE/p) for p in dict.fromkeys(sources)],
           reports=[record(HERE/'results'/f'{p}.json') for p in reports],
           native_artifacts=[record(p) for p in native], parameters=[record(p) for p in params],
           source_report_bindings=[dict(report=f'results/{r}.json', source=s) for r,s in bindings],
           limits=['Component checks and small oracles are not full benchmark acceptance.',
                   'Large v3 generalized reconstruction is count-only; large v1 passed native checks.',
                   'Prefix serial_matches compares incompatible head encodings; it is not a tail correctness test.',
                   'Chemistry union floor excludes only the saved dictionary, not all shapes.',
                   'Earlier phase probe records a subsequently fixed dirty-flag cleanup bug.',
                   'Separate public-component full verification is ongoing, outside this closed research wave.'])
path = HERE/'results/generalized_quotient_research_manifest.json'
path.write_text(json.dumps(out, indent=2)+'\n')
print(json.dumps(dict(status=out['status'], sources=len(out['sources']), reports=len(reports), native_artifacts=len(native), accepted_archive_sha256=accepted)))
