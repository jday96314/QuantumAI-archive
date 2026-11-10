"""Run the precommitted native checks and accept strictly on paired behavior/T.

The host verifier is invoked unchanged. Both paired runs use the same explicit
seed; defaults and independent salt are reported without selecting between them.
Run only after the incumbent paired run and the builder have finished.
"""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parents[2]
VERIFIER=PROJECT/'from_host/verifiers/dqi/verify_submission_dqi.py'
R4=PROJECT/'initial_experiments/4/codex'

def read(path):return json.loads(Path(path).read_text())

def run(name,extra):
    path=HERE/'results'/f'{name}.json';log=path.with_suffix('.log')
    cmd=['flock','/tmp/r6-fullverify.lock',sys.executable,str(VERIFIER),
         str(HERE/'candidates/trimmed_decoder.kmb'),'--max-workers','2',*extra]
    meta=dict(command=cmd,shots=102400,start_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
    (HERE/'results'/f'{name}_command.json').write_text(json.dumps(meta,indent=2)+'\n')
    print(json.dumps(dict(stage='start',name=name,command=cmd)),flush=True)
    with path.open('w') as out,log.open('w') as err:
        subprocess.run(cmd,check=True,stdout=out,stderr=err)
    report=read(path)
    assert report['status']=='PASS' and report['verified_shots']==102400,report
    assert report['num_qubits']==961 and report['max_toffoli_count']==994358,report
    print(json.dumps(dict(stage='complete',name=name,R=report['repetition_factor'],
        mean_T=report['toffoli_count'],cost=report['score'],value_errors=report['value_errors'],phase_errors=report['phase_errors'])),flush=True)
    return report

def main():
    os.chdir(HERE);os.environ['OMP_NUM_THREADS']='1';os.environ['OPENBLAS_NUM_THREADS']='1'
    plan=read(HERE/'parameters/verification_plan.json')
    base=read(HERE/'results/incumbent_seed_2026102405.json')
    assert base['status']=='PASS' and base['verified_shots']==102400
    assert base['circuit_sha256']=='4989da5cb633529e8851b9aa0c8c52abeb8990423a3a9463a70f673661a909c0'
    build=read(HERE/'results/trimmed_decoder_build.json')
    assert build['Q']==961 and build['static_T']==994358
    paired=run('trimmed_seed_2026102405',['--seed',str(plan['paired_explicit_seed'])])
    match={key:(base[key]==paired[key]) for key in ['value_errors','phase_errors','verified_shots','num_qubits','effective_error_rate','repetition_factor']}
    comparison=dict(seed=plan['paired_explicit_seed'],matches=match,
        incumbent={key:base[key] for key in ['value_errors','phase_errors','num_qubits','toffoli_count','repetition_factor','score']},
        trimmed={key:paired[key] for key in ['value_errors','phase_errors','num_qubits','toffoli_count','repetition_factor','score']},
        mean_T_reduction=base['toffoli_count']-paired['toffoli_count'])
    (HERE/'results/paired_comparison.json').write_text(json.dumps(comparison,indent=2)+'\n')
    assert all(match.values()),f'Paired behavior mismatch: {comparison}'
    assert abs(comparison['mean_T_reduction']-25596)<1e-8,comparison
    default=run('trimmed_default_verify',[])
    independent=run('trimmed_independent_verify',['--salt',plan['independent_salt']])
    for report in [paired,default,independent]:
        assert report['circuit_sha256']==build['circuit_sha256']
        assert report['value_errors']==0
    accepted=dict(status='ACCEPTED',name='trimmed decoder',basis='lower T at equal Q and exactly matched decoding on the precommitted explicit seed; no selection between default and salt outcomes',
        Q=961,static_T=994358,core_T=624206,scan_T=370152,static_T_reduction=25596,
        circuit='candidates/trimmed_decoder.kmb',builder='build_accepted_decoder.py',
        circuit_sha256=build['circuit_sha256'],paired_comparison=comparison,
        default={key:default[key] for key in ['verified_shots','value_errors','phase_errors','repetition_factor','toffoli_count','score']},
        independent={key:independent[key] for key in ['verified_shots','value_errors','phase_errors','repetition_factor','toffoli_count','score']},
        independent_salt=plan['independent_salt'],native_verifier_sha256=hashlib.sha256(VERIFIER.read_bytes()).hexdigest())
    (HERE/'results/accepted_decoder.json').write_text(json.dumps(accepted,indent=2)+'\n')
    results=HERE/'RESULTS.md';history=results.read_text()
    section='ACCEPTED: trimmed decoder\n\n'
    section+='Accepted on lower T at equal Q and exact paired-seed decoding, not on failure-count luck. '
    section+='Full native Q=961, static T=994358 (25596 fewer than incumbent). '
    section+='The schedule and cubic scan are unchanged.\n\n'
    section+='Builder: `build_accepted_decoder.py`; KMB: `candidates/trimmed_decoder.kmb`; '
    section+=f'SHA-256 `{build["circuit_sha256"]}`.\n\n'
    section+='| Full 102400-shot run | Value / phase failures | R | Mean T | Cost |\n'
    section+='| --- | ---: | ---: | ---: | ---: |\n'
    for name,report in [('Incumbent, paired seed 2026102405',base),('Trimmed, paired seed 2026102405',paired),('Trimmed, default',default),('Trimmed, independent salt',independent)]:
        section+=f'| {name} | {report["value_errors"]} / {report["phase_errors"]} | {report["repetition_factor"]:.12g} | {report["toffoli_count"]:.12f} | {report["score"]:.12g} |\n'
    section+='\nAll four runs used the unchanged full-default host verifier under '
    section+='`flock /tmp/r6-fullverify.lock`, with at most two CPU workers. '
    section+='The paired runs have identical failure counts, R and Q, and exactly 25596 fewer mean T. '
    section+='JSONs: `results/incumbent_seed_2026102405.json`, `results/trimmed_seed_2026102405.json`, '
    section+='`results/trimmed_default_verify.json`, `results/trimmed_independent_verify.json`; '
    section+='acceptance: `results/accepted_decoder.json`. Independent salt: '
    section+=f'`{plan["independent_salt"]}`.\n\n'
    section+='Rebuild from repository root: `source initial_experiments/5/dqi-core/env.sh`, then '
    section+='`"$EXP_PYTHON" initial_experiments/5/dqi-core/build_accepted_decoder.py`. '
    section+='Re-run verification after the incumbent paired report exists with '
    section+='`"$EXP_PYTHON" initial_experiments/5/dqi-core/verify_and_accept.py`. '
    section+='Commands and fixed sampling plan are recorded in `parameters/verification_plan.json` and `results/*_command.json`.\n\n'
    section+='The earlier research-only stopping result below is superseded by this full-circuit acceptance. '
    section+='The protected round-4 incumbent remains unchanged. Work stops here.\n\n---\n\n'
    results.write_text(section+history)
    memo=HERE/'DESIGN_MEMO.md'
    memo.write_text('Full-circuit follow-up: the trimmed decoder is now ACCEPTED; see the acceptance block at the top of RESULTS.md. The research-only acceptance status below is historical.\n\n'+memo.read_text())
    print(json.dumps(accepted,indent=2),flush=True)

if __name__=='__main__':main()
