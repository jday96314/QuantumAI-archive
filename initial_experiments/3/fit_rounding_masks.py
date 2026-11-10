"""Fit coefficients while changing which partial rows round upward."""
import json,os,subprocess
from pathlib import Path
import numpy as np
import binary_centered_invsqrt as m
root=Path(__file__).resolve().parent;original=m.table_for;records=[]
configs=[(7,(27,17,w)) for w in [6,7,8,9]]+[(8,(26,16,w)) for w in [7,8,9]]+[(8,(27,15,w)) for w in [7,8,9]]
for f,widths in configs:
    full=(1<<(21-f))-1;patterns=[0,full,0x5555&full,0xaaaa&full]
    for k in range(1,5):patterns.extend([full^((1<<k)-1),full>>k,full^(1<<k)])
    patterns=list(dict.fromkeys(patterns));pairs=[(p,full) for p in patterns]+[(full,p) for p in patterns]+[(p,p) for p in patterns]
    for i,(p0,p1) in enumerate(dict.fromkeys(pairs)):
        tag=f'roundmask_f{f}_w'+ '_'.join(map(str,widths))+f'_{p0}_{p1}'
        tab=original(f,widths,0);src=root/'candidates'/f'{tag}.in.txt';out=src.with_name(f'{tag}.out.txt')
        np.savetxt(src,tab,fmt='%d',header=str(len(tab)),comments='')
        r=subprocess.run([str(root/'fit_coefficients'),'binarycentered_reachable_rounded_roundmask_deep',str(f),*map(str,widths),str(src),str(out),f'{p0},{p1}'],text=True,capture_output=True,env=dict(os.environ,OMP_NUM_THREADS='4'),check=True)
        err=float(r.stdout.strip());rec=dict(tag=tag,frac=f,widths=widths,masks=[p0,p1],error=err)
        if err<2**-24:
            fitted=np.loadtxt(out,dtype=np.int64,skiprows=1);m.table_for=lambda *args:fitted
            xs=np.arange(1,3133453);actual=np.max(abs(m.model(xs,f,widths,rounded=True,rounding_masks=[p0,p1])/2.**(widths[0]+10)*np.sqrt(xs)-1))
            sizes=[int(col.max()).bit_length() for col in fitted.T]
            split=[sizes[j+1]+widths[j]-widths[j+1]-f-1 if j<2 else None for j in range(3)]
            c=m.build(f,widths,rounded=True,rounding_masks=[p0,p1],accum_bits=split,carry_bits=1,reuse_normalized_top=True,output_extra_bits=10)
            c.to_file(root/'candidates'/f'{tag}.kmb',format='kmb');np.save(root/'candidates'/f'{tag}.npy',fitted)
            rec.update(full_pass=bool(actual<2**-24),actual_error=float(actual),Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,estimate=c.num_qubits*(c.operation_statistics().num_toffoli_like+245),split=split)
            m.table_for=original
        records.append(rec);print(json.dumps(rec),flush=True)
        (root/'results/rounding_mask_fit.json').write_text(json.dumps(records,indent=2))
