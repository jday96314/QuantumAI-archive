"""Run unmodified from_host/metric.py on the complete extracted archive.

Only EVAL_INPUTS_DIR is configured for local input files. Default shot counts,
confidence, scoring and verifier functions are unchanged. CPU affinity limits
the native process pools to two CPUs. A passive Python profile observer records
native returns and parent Future completions; it never changes their results.
"""
from pathlib import Path
import os,sys,json,hashlib,time,math,zipfile,importlib.util

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
VALID=HERE/'package_validation';INPUTS=VALID/'eval_inputs';EXTRACT=VALID/'extracted'
SALT='chem-audit-complete-archive-validation-2026-10-24-fixed-v1'
BASE=dict(invsqrt=2076774.233715,dqi=24439463583672.929688,dfthcsos=3355599588.315537)
EXPECTED=dict(invsqrt=200540928,dqi=102400,dfthcsos=50880)
CPUS=sorted(os.sched_getaffinity(0))[:2];os.sched_setaffinity(0,CPUS)
assert os.process_cpu_count()==2
for p in [VALID,INPUTS,EXTRACT]:p.mkdir(parents=True,exist_ok=True)
import numpy as np
import pandas as pd

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(8<<20),b''):h.update(chunk)
    return h.hexdigest()
def default(o):
    if isinstance(o,np.generic):return o.item()
    if isinstance(o,np.ndarray):return o.tolist()
    if isinstance(o,Path):return str(o)
    raise TypeError(type(o).__name__)
def dump(path,data):Path(path).write_text(json.dumps(data,indent=2,default=default)+'\n')

assembly=json.loads((HERE/'results/package_assembly.json').read_text());archive=HERE/'submission.zip'
assert sha(archive)==assembly['sha256']
with zipfile.ZipFile(archive) as z:
    assert set(z.namelist())==set(assembly['members'])
    for name in z.namelist():
        assert not Path(name).is_absolute() and '..' not in Path(name).parts
    z.extractall(EXTRACT)
for name,digest in assembly['members'].items():assert sha(EXTRACT/name)==digest
solution=pd.DataFrame([dict(row_id=k,benchmark=k,Usage='Public',seed_salt=SALT,baseline_cost=BASE[k]) for k in BASE])
solution.to_csv(INPUTS/'solution.csv',index=False)
dump(INPUTS/'metadata.json',dict(Hyperparameters=dict(verification_shot_divisor=1)))
metric_path=ROOT/'from_host/metric.py';metric_hash=sha(metric_path)
spec=importlib.util.spec_from_file_location('chem_audit_unmodified_host_metric',metric_path)
metric=importlib.util.module_from_spec(spec);spec.loader.exec_module(metric)
metric.EVAL_INPUTS_DIR=INPUTS
config=dict(archive_sha256=assembly['sha256'],salt=SALT,verification_shot_divisor=1,cpu_affinity=CPUS,process_cpu_count=os.process_cpu_count(),
    metric_source_sha256=metric_hash,solution_sha256=sha(INPUTS/'solution.csv'),metadata_sha256=sha(INPUTS/'metadata.json'))
dump(HERE/'results/package_metric_config.json',config)

progress=HERE/'results/package_metric_progress.jsonl';progress.write_text('')
root_pid=os.getpid();active=None;native={};completed={k:0 for k in BASE};chunk_shots={};observer_errors=[]
modules={'invsqrt.verify_submission_invsqrt':'invsqrt','dqi.verify_submission_dqi':'dqi','dfthcsos.verify':'dfthcsos'}
def log(record):
    record.update(time=time.time(),pid=os.getpid())
    with progress.open('a') as f:f.write(json.dumps(record,default=default)+'\n')
    if record['event']!='progress':print(json.dumps(record,default=default),flush=True)
def observer(frame,event,arg):
    global active
    if os.getpid()!=root_pid or event not in ['call','return']:return
    try:
        module=frame.f_globals.get('__name__');name=frame.f_code.co_name
        if module in modules and name=='verify':
            benchmark=modules[module]
            if event=='call':
                active=benchmark
                if benchmark=='invsqrt':chunk_shots[benchmark]=frame.f_globals['BATCHES_PER_CHUNK']*frame.f_locals['batch_size']
                elif benchmark=='dqi':chunk_shots[benchmark]=frame.f_locals['batch_size']
                log(dict(event='native_begin',benchmark=benchmark,expected_shots=EXPECTED[benchmark]))
            elif isinstance(arg,dict):
                native[benchmark]=arg;dump(HERE/f'results/package_native_{benchmark}.json',arg)
                log(dict(event='native_end',benchmark=benchmark,status=arg.get('status'),verified_shots=arg.get('verified_shots'),score=arg.get('score')))
                active=None
        elif module=='concurrent.futures._base' and name=='result' and event=='return' and active in chunk_shots:
            item=arg[0] if isinstance(arg,list) and len(arg)==1 else arg
            size=5 if active=='invsqrt' else 4
            if isinstance(item,tuple) and len(item)==size:
                completed[active]+=1
                log(dict(event='progress',benchmark=active,completed_chunks=completed[active],
                    completed_shots_lower_bound=min(EXPECTED[active],completed[active]*chunk_shots[active])))
    except Exception as error:observer_errors.append(repr(error))

# Child simulators need no profiler; their test inputs and computations are
# untouched. Observe only the parent's returned Future results for progress.
os.register_at_fork(after_in_child=lambda:sys.setprofile(None))
os.chdir(VALID);started=time.time()
log(dict(event='metric_read_and_validate_begin',archive=str(archive)))
sys.setprofile(observer)
try:metric.generate_standard_submission(str(EXTRACT))
finally:sys.setprofile(None)
assert not observer_errors,observer_errors
assert set(native)==set(BASE),native.keys()
for k,r in native.items():assert r['status']=='PASS' and r['verified_shots']==EXPECTED[k]
assert native['dqi']['log2_delta']==-128.0
assert native['invsqrt']['circuit_sha256']==assembly['members']['invsqrt.kmb']
assert native['dqi']['circuit_sha256']==assembly['members']['dqi.kmb']
assert native['dfthcsos']['submission_sha256']==assembly['canonical_chemistry_digest']
assert sha(metric_path)==metric_hash
scored=pd.read_csv(VALID/'submission.csv');solution_read=pd.read_csv(INPUTS/'solution.csv')
score,details=metric.score(solution_read,scored,'row_id',verification_shot_divisor=1)
for k,r in native.items():assert math.isclose(details[k]['cost'],r['score'],rel_tol=1e-13)
record=dict(status='PASS',archive=str(archive),archive_sha256=sha(archive),local_combined_score=score,
    accepted_reference_score=assembly['accepted_reference_score'],task_details=details,
    verified_shots={k:r['verified_shots'] for k,r in native.items()},configuration=config,
    native_reports={k:str(HERE/f'results/package_native_{k}.json') for k in BASE},
    generated_submission_csv=str(VALID/'submission.csv'),seconds=time.time()-started,
    scope='Unmodified host metric readers, full-default native verifications and metric.score; one fixed predeclared local salt. Salt-dependent DQI cost is reported as measured, never optimized over salts.')
dump(HERE/'results/package_metric.json',record)
summary=dict(assembly,status='PACKAGE_FULL_METRIC_PASS',overall_score=score,local_metric_score=score,
    task_scores={k:r['score'] for k,r in native.items()},sources=record['native_reports'],metric_result=str(HERE/'results/package_metric.json'))
dump(HERE/'submission.summary.json',summary)
log(dict(event='metric_complete',status='PASS',local_combined_score=score,seconds=record['seconds']))
print(json.dumps(record,indent=2,default=default),flush=True)
