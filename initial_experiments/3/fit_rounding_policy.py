"""Screen pairs of partial-row rounding rules, chosen by a lookup bit."""
import itertools,json,os,subprocess
from pathlib import Path
import numpy as np
import binary_centered_invsqrt as m
root=Path(__file__).resolve().parent
configs=[(8,(26,16,8)),(8,(26,17,8)),(7,(27,17,8)),(7,(27,16,9))]
for f,w in configs:
    full=(1<<(21-f))-1
    masks=[(full,full)]+[(full^(1<<k),full) for k in range(21-f)]+[(full,full^(1<<k)) for k in range(21-f)]
    masks+= [(0,full),(full,0),(0,0),(0x5555&full,full),(0xaaaa&full,full)]
    outputs=[];tab=m.table_for(f,w,0)
    for p0,p1 in masks:
        tag=f'policy_f{f}_w'+ '_'.join(map(str,w))+f'_{p0}_{p1}'
        src=root/'candidates'/f'{tag}.in.txt';out=src.with_name(f'{tag}.out.txt');np.savetxt(src,tab,fmt='%d',header=str(len(tab)),comments='')
        r=subprocess.run([str(root/'fit_coefficients'),'binarycentered_reachable_rounded_roundmask_allerrors',str(f),*map(str,w),str(src),str(out),f'{p0},{p1}'],text=True,capture_output=True,env=dict(os.environ,OMP_NUM_THREADS='4'),check=True)
        arr=np.loadtxt(out,skiprows=1);outputs.append(arr)
        print(json.dumps(dict(frac=f,widths=w,masks=[p0,p1],max_error=float(arr[:,3].max()),failed=int((arr[:,3]>=2**-24).sum()))),flush=True)
    errors=np.stack([v[:,3] for v in outputs]);eligible=errors<2**-24
    best=errors.argmin(0);coverage=eligible.any(0)
    pairs=[(i,j) for i,j in itertools.combinations(range(len(masks)),2) if np.all(eligible[i]|eligible[j])]
    print(json.dumps(dict(frac=f,widths=w,uncovered=int((~coverage).sum()),pairs=[(masks[i],masks[j]) for i,j in pairs])),flush=True)
    np.savez(root/'results'/(f'rounding_policy_f{f}_w'+'_'.join(map(str,w))+'.npz'),masks=np.array(masks),fits=np.array(outputs))
