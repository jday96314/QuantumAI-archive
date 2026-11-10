"""Resource screen for refitted, fine-coordinate quadratic inverse square root."""
import json
from pathlib import Path
import numpy as np
import binary_centered_invsqrt as p
from optimize_workspace import compact

root=Path(__file__).resolve().parent
rows=[]
for widths in [(27,17,7),(27,17,8)]:
    label='_'.join(map(str,widths))
    table=np.loadtxt(root/f'parameters/fine8_w{label}.out.txt',skiprows=1,dtype=np.int64)
    np.save(root/f'candidates/fine8_w{label}.npy',table)
    p.table_for=lambda *unused:table
    for limit in [82,84,86,88,90,92]:
      for keep in [0,2,4]:
       for delay in [4,6,8]:
        for product in [12,16,None]:
         try:
            kw=dict(frac=8,widths=widths,lookup_frac=8,accum_bits=[17,8,None],carry_bits=1,
                    rounded=True,booth=True,reuse_normalized_top=True,output_extra_bits=10,
                    defer_workspace=True,qubit_limit=limit,curve_keep_bits=keep,defer_curve_phase=True,
                    product_limit=product,borrow_integration=True,monomial_delete=True,delay_high_bits=delay)
            c=p.build(**kw)
            # The allocator can leave holes after a lifetime; the final index
            # compactor is run only on the finalists.
            T=c.operation_statistics().num_toffoli_like
            row=dict(widths=widths,limit=limit,keep=keep,delay=delay,product=product,
                     Q=c.num_qubits,T=T,estimate=c.num_qubits*(T-11+245))
            rows.append(row)
         except ValueError:pass
    print(json.dumps(dict(widths=widths,best=sorted([r for r in rows if tuple(r['widths'])==widths],key=lambda r:r['estimate'])[:3])),flush=True)
ranked=sorted(rows,key=lambda r:r['estimate'])
seen=set()
for row in ranked[:16]:
    widths=tuple(row['widths']);label='_'.join(map(str,widths))
    table=np.load(root/f'candidates/fine8_w{label}.npy');p.table_for=lambda *unused:table
    kw=dict(frac=8,widths=widths,lookup_frac=8,accum_bits=[17,8,None],carry_bits=1,
            rounded=True,booth=True,reuse_normalized_top=True,output_extra_bits=10,
            defer_workspace=True,qubit_limit=row['limit'],curve_keep_bits=row['keep'],defer_curve_phase=True,
            product_limit=row['product'],borrow_integration=True,monomial_delete=True,delay_high_bits=row['delay'])
    tag=f"fine8_w{label}_q{row['limit']}_k{row['keep']}_d{row['delay']}_p{row['product']}"
    path=root/f'candidates/{tag}.kmb';c=compact(p.build(**kw),path)
    c.to_file(str(path),format='kmb');row.update(Q=c.num_qubits,path=str(path),tag=tag,parameters=kw)
    print(json.dumps(row),flush=True)
(root/'results/fine8_resource_screen.json').write_text(json.dumps(ranked,indent=2)+'\n')
