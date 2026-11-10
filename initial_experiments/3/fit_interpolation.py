"""Search nearby integer polynomial coefficients with a compiled exact model."""
from pathlib import Path
import importlib,json,subprocess
import numpy as np
HERE=Path(__file__).resolve().parent
configs=[]
for mode,name in [('centered','centered_invsqrt'),('binary','binary_invsqrt'),('binarycentered','binary_centered_invsqrt')]:
    p=importlib.import_module(name)
    for f in [6,7,8]:
        for wc in [27,28]:
            for w1 in range(wc-f-4,wc-f):
                for w2 in range(max(2,wc-2*f-5),max(3,wc-2*f)):
                    widths=(wc,w1,w2);c=p.build(f,widths,0)
                    q=c.num_qubits;t=c.operation_statistics().num_toffoli_like
                    est=q*(t+260)
                    if est<310000:
                        configs.append(dict(mode=mode,module=name,frac=f,widths=widths,qubits=q,toffolis=t,estimated_score=est))
configs.sort(key=lambda x:x['estimated_score'])
result_path=HERE/'results/interpolation_fit.json'
records=json.loads(result_path.read_text()) if result_path.exists() else []
done={(r['mode'],r['frac'],tuple(r['widths'])) for r in records}
for rec in configs:
    if (rec['mode'],rec['frac'],tuple(rec['widths'])) in done:continue
    p=importlib.import_module(rec['module']);original=p.table_for
    f=rec['frac'];widths=tuple(rec['widths']);mode=rec['mode']
    tag=mode+'_f'+str(f)+'_w'+'_'.join(map(str,widths));path=HERE/'candidates'/tag
    tab=original(f,widths,0).copy()
    src=path.with_suffix('.in.txt');dst=path.with_suffix('.out.txt')
    src.write_text(str(len(tab))+'\n'+'\n'.join(' '.join(map(str,row)) for row in tab)+'\n')
    result=subprocess.run([str(HERE/'fit_coefficients.next'),mode,str(f),*map(str,widths),str(src),str(dst),'unused'],capture_output=True,text=True,check=True)
    err=float(result.stdout.strip());rec['max_relative_error']=err;rec['full_pass']=err<2.**-24
    if rec['full_pass']:
        tab=np.loadtxt(dst,skiprows=1,dtype=np.int64)
        def fit(ff,ww,bias=0):return tab if ff==f and tuple(ww)==widths else original(ff,ww,bias)
        p.table_for=fit
        try:c=p.build(f,widths,0)
        finally:p.table_for=original
        np.save(path.with_suffix('.npy'),tab);c.to_file(str(path.with_suffix('.kmb')),format='kmb')
        rec.update(circuit=str(path.with_suffix('.kmb')),table=str(path.with_suffix('.npy')))
    print(json.dumps(rec),flush=True);records.append(rec)
    tmp=HERE/'results/interpolation_fit.tmp'
    tmp.write_text(json.dumps(records,indent=2));tmp.replace(result_path)
