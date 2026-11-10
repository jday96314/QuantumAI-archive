"""Compare actual compacted resources for fitted coordinate sizes."""
import hashlib,json
from pathlib import Path
import numpy as np
import binary_centered_invsqrt as p
from optimize_workspace import compact

root=Path(__file__).resolve().parent
rows=[];circuits={}
for frac,lf,widths,split,stem in [
    (8,8,(27,17,7),(17,7),'fine8_w27_17_7'),
    (8,8,(27,17,8),(17,7),'fine8_w27_17_8'),
    (6,7,(27,19,10),(19,11),'coarse6_w27_19_10'),
    (6,7,(27,19,11),(19,11),'coarse6_w27_19_11'),
]:
    logfile=root/f'results/{stem}_fit.log'
    if not logfile.exists() or logfile.read_text().strip().split()[-1]!='1':continue
    table=np.loadtxt(root/f'parameters/{stem}.out.txt',skiprows=1,dtype=np.int64)
    np.save(root/f'candidates/{stem}.npy',table);p.table_for=lambda *unused:table
    for limit in [80,82,84,86,88,90,92,94,96]:
     for keep in [0,2,4,6]:
      if keep>=int(table[:,2].max()).bit_length():continue
      for delay in [0,2,4,6]:
       for product in [16,None]:
        kw=dict(frac=frac,widths=widths,lookup_frac=lf,accum_bits=[*split,None],carry_bits=1,
                rounded=True,booth=True,reuse_normalized_top=True,output_extra_bits=10,
                defer_workspace=True,qubit_limit=limit,curve_keep_bits=keep,defer_curve_phase=True,
                product_limit=product,borrow_integration=True,monomial_delete=True,delay_high_bits=delay)
        try:c=p.build(**kw)
        except ValueError:continue
        T=c.operation_statistics().num_toffoli_like
        # Compaction is needed for every resource comparison: an allocator
        # hole may inflate the largest wire index by several qubits.
        c=compact(c,root/'candidates/quadratic_screen_tmp.kmb')
        cost=c.num_qubits*(T-11+250)
        row=dict(stem=stem,Q=c.num_qubits,T=T,estimate=cost,parameters=kw)
        rows.append(row)
        key=(stem,c.num_qubits,T)
        if key not in circuits and cost<190000:
            path=root/f'candidates/{stem}_q{c.num_qubits}_t{T}.kmb';c.to_file(str(path),format='kmb')
            row.update(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest())
            circuits[key]=row
    print(json.dumps(dict(stem=stem,best=sorted([r for r in rows if r['stem']==stem],key=lambda r:r['estimate'])[:2])),flush=True)
rows.sort(key=lambda r:r['estimate'])
(root/'results/quadratic_fit_resources.json').write_text(json.dumps(rows,indent=2)+'\n')
