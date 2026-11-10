"""Reproduce the numerical models, components, and two frozen full candidates.

Run with this experiment's env.sh sourced. --full also reruns the unchanged
host verifier, with all default inputs and all 64 repetitions per input.
"""
import argparse,json,subprocess,shutil
from pathlib import Path
from build import HERE,build,save
from certify import certificate

# name, mode, starting depth, output/slope/curvature fractional bits,
# relative fitting target, depth cap, slope search radius.
SPECS=[
 ['local_q28','local-quadratic',5,28,19,9,5.6e-8,9,18],
 ['shared_q28','shared-quadratic',6,28,19,9,5.6e-8,9,18],
 ['local_linear29','local-linear',8,29,20,10,5.6e-8,13,18],
 ['local_q27','local-quadratic',5,27,18,8,5.6e-8,9,24],
 ['shared_q27','shared-quadratic',6,27,18,8,5.6e-8,9,24],
 ['local_mixed28','local-mixed',5,28,19,9,5.6e-8,12,18],
 ['shared7_q27','shared-quadratic',7,27,18,8,5.6e-8,9,24],
 ['local7_q27','local-quadratic',7,27,18,8,5.6e-8,9,24],
 ['local_q29','local-quadratic',5,29,19,9,5.6e-8,9,24],
 ['shared6_precise_q27','shared-quadratic',6,27,19,9,5.6e-8,9,24],
 ['shared7_narrow_q27','shared-quadratic',7,27,17,7,5.6e-8,10,32],
 ['local_narrow_q28','local-quadratic',6,28,18,8,5.6e-8,9,32],
 ['shared8_narrow_q27','shared-quadratic',8,27,17,7,5.6e-8,10,32],
 ['local_q26','local-quadratic',6,26,18,8,5.6e-8,10,32],
 ['shared7_loose_q27','shared-quadratic',7,27,18,8,5.9e-8,9,48],
 ['local_accurate_q28','local-quadratic',5,28,19,9,5e-8,9,24],
 ['shared7_precise_q26','shared-quadratic',7,26,19,9,5.9e-8,10,32],
 ['local7_q26','local-quadratic',7,26,18,8,5.9e-8,10,48],
 ['mixed_tail1_32_q28','local-mixed@0.96875',5,28,19,9,5.6e-8,12,18],
 ['mixed_tail1_8_q28','local-mixed@0.875',5,28,19,9,5.6e-8,12,18],
]

def run(args):
    subprocess.run(list(map(str,args)),check=True)

if __name__=='__main__':
 import sys
 p=argparse.ArgumentParser();p.add_argument('--models',action='store_true');p.add_argument('--components',action='store_true')
 p.add_argument('--circuits',action='store_true');p.add_argument('--full',action='store_true');a=p.parse_args()
 if not any(vars(a).values()):a.models=a.components=a.circuits=True
 if a.models:
    run(['g++','-O3','-std=c++17',HERE/'fit.cpp','-o',HERE/'fit'])
    for name,mode,base,w0,w1,w2,tol,cap,radius in SPECS:
        path=HERE/'models'/f'{name}.json'
        run([HERE/'fit',mode,base,w0,w1,w2,tol,cap,path,radius])
        cert=certificate(json.loads(path.read_text()))
        (HERE/'results'/f'{name}_certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
 if a.components:
    run([sys.executable,HERE/'recenter.py']);run([sys.executable,HERE/'stream_erase.py'])
 if a.circuits:
    for name,model_name,limit in [('final_shared','shared7_loose_q27',90),('final_local','local7_q27',92)]:
        model=json.loads((HERE/'models'/f'{model_name}.json').read_text())
        cert=certificate(model);model['affine_residual_range']=cert['affine_residual_range']
        c,profile=build(model,limit=limit,keep=4,delay=2,product_limit=17,deferred=3)
        # Check saved frozen artifact, rather than overwrite it during a verifier run.
        c,path,resources=save(c,name+'_rebuilt',profile['stage_live_peaks'])
        original=HERE/'candidates'/f'{name}.kmb'
        assert original.read_bytes()==path.read_bytes(),f'{name} rebuild differs'
        run([sys.executable,HERE/'check.py',path,HERE/'models'/f'{model_name}.json','--random','32768','--repeats','8'])
 if a.full:
    verifier=HERE.parents[2]/'from_host/verifiers/invsqrt/verify_submission_invsqrt.py'
    for name in ['final_shared','final_local']:
        with (HERE/'results'/f'{name}_reverified.json').open('w') as out:
            subprocess.run([sys.executable,str(verifier),str(HERE/'candidates'/f'{name}.kmb'),'--max-workers','4'],check=True,stdout=out)
