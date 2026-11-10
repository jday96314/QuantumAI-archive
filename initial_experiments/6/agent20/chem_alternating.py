"""Finite alternating exact-word/physical-weight search with preserved evidence."""
import argparse,hashlib,json,shutil,subprocess,sys,time
from pathlib import Path

HERE=Path(__file__).resolve().parent


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--init',type=Path,required=True);p.add_argument('--cycles',type=int,default=5)
    p.add_argument('--tag',default='alternating_B703');args=p.parse_args();current=args.init.resolve();records=[];best=float('inf');started=time.time()
    for cycle in range(args.cycles):
        for phase in ('words','weights'):
            tag=f'{args.tag}_{cycle}_{phase}'
            script=HERE/('chem_fixed_weight_repair.py' if phase=='words' else 'chem_tangent_refit.py')
            command=[sys.executable,str(script),'--init',str(current),'--tag',tag,'--lambda-cap','330','--lambda-eff-cap','51']
            command+=['--steps','300','--normalization-penalty','0.0001'] if phase=='words' else ['--steps','1000','--ridge','0.001','--step-norm','0.2','--project-constraints']
            log=HERE/'results'/f'{tag}.log';source_sha=hashlib.sha256(script.read_bytes()).hexdigest()
            with log.open('w') as output:subprocess.run(command,cwd=HERE,stdout=output,stderr=subprocess.STDOUT,check=True)
            history=json.loads((HERE/'results'/f'{tag}_history.json').read_text());last=history[-1];minimum=min(history,key=lambda row:row['error2'])
            current=HERE/'parameters'/f'{tag}_last.npz'
            checkpoint=HERE/'parameters'/f'{tag}_best_error.npz' if phase=='words' else current
            record=dict(cycle=cycle,phase=phase,seconds=time.time()-started,source_sha256=source_sha,log=str(log),last=last,best_error=minimum)
            if minimum['error2']<best:
                best=minimum['error2'];destination=HERE/'parameters'/f'{args.tag}_best_error.npz';shutil.copyfile(checkpoint,destination)
                record.update(best_checkpoint=str(destination),best_checkpoint_sha256=hashlib.sha256(destination.read_bytes()).hexdigest())
            records.append(record);print(json.dumps(record),flush=True)
            (HERE/'results'/f'{args.tag}_progress.json').write_text(json.dumps(records,indent=2)+'\n')
            if best<.019995:
                print(json.dumps(dict(status='POTENTIAL_TWO_BODY_PASS_NEEDS_INDEPENDENT_PHYSICS_AND_NATIVE_ORACLES',error2=best)),flush=True);sys.exit(0)
    print(json.dumps(dict(status='FINITE_EXPERIMENT_COMPLETE',best_error2=best)),flush=True)
