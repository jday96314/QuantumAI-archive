"""Revisit interpolation widths using partitioned carry arithmetic."""
import json,subprocess
from pathlib import Path
import numpy as np
import binary_centered_invsqrt as p
HERE=Path(__file__).resolve().parent;orig=p.table_for;configs=[]
for f in [6,7,8,9]:
 for wc in [27,28,29]:
  for w1 in range(wc-f-6,wc-f+1):
   for w2 in range(max(2,wc-2*f-9),max(3,wc-2*f+1)):
    widths=(wc,w1,w2);tab=orig(f,widths,0)
    sizes=[int(x.max()).bit_length() for x in tab.T]
    split=[min(sizes[j]-1,max(1,sizes[j+1]+widths[j]-widths[j+1]-f-1)) for j in range(2)]+[None]
    c=p.build(f,widths,accum_bits=split,carry_bits=1);q=c.num_qubits;t=c.operation_statistics().num_toffoli_like
    est=q*(t+245)
    if est<250000:configs.append(dict(frac=f,widths=widths,split=split,Q=q,T=t,estimate=est))
configs.sort(key=lambda r:r['estimate']);path=HERE/'results/split_interpolation_fit.json'
records=json.loads(path.read_text()) if path.exists() else [];done={(r['frac'],tuple(r['widths'])) for r in records}
for rec in configs:
 f=rec['frac'];widths=tuple(rec['widths'])
 if (f,widths) in done:continue
 tag='splitfit_f'+str(f)+'_w'+'_'.join(map(str,widths));base=HERE/'candidates'/tag
 tab=orig(f,widths,0);src=base.with_suffix('.in.txt');dst=base.with_suffix('.out.txt')
 src.write_text(str(len(tab))+'\n'+'\n'.join(' '.join(map(str,r)) for r in tab)+'\n')
 rr=subprocess.run([str(HERE/'fit_coefficients.next'),'binarycentered',str(f),*map(str,widths),str(src),str(dst),'unused'],check=True,capture_output=True,text=True)
 err=float(rr.stdout.strip());rec.update(error=err,full_pass=err<2**-24)
 if rec['full_pass']:
  tab=np.loadtxt(dst,skiprows=1,dtype=np.int64);p.table_for=lambda ff,ww,bias:tab
  try:c=p.build(f,widths,accum_bits=rec['split'],carry_bits=1)
  finally:p.table_for=orig
  c.to_file(str(base.with_suffix('.kmb')),format='kmb');np.save(base.with_suffix('.npy'),tab)
  rec.update(circuit=str(base.with_suffix('.kmb')),table=str(base.with_suffix('.npy')))
 print(json.dumps(rec),flush=True);records.append(rec);path.write_text(json.dumps(records,indent=2)+'\n')
