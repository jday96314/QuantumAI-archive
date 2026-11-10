"""Search exact workspace schedules with the verified coefficient table fixed."""
import argparse, json, subprocess, time
from pathlib import Path
import numpy as np
import kickmix as km
import binary_centered_invsqrt as poly

HERE=Path(__file__).resolve().parent
TABLE=np.load(HERE/'candidates/margin_curve8_t5.9e-08.npy')
poly.table_for=lambda *args:TABLE

def build(limit,keep,product,splits):
    return poly.build(7,(27,18,8),lookup_frac=8,accum_bits=[*splits,None],
        carry_bits=1,rounded=True,booth=True,reuse_normalized_top=True,
        output_extra_bits=10,defer_workspace=True,qubit_limit=limit,
        curve_keep_bits=keep,defer_curve_phase=True,product_limit=product,
        borrow_integration=True)

def compact(c,path):
    raw=path.with_suffix('.raw.kmx'); txt=path.with_suffix('.kmx')
    c.to_file(str(raw),format='kmx')
    subprocess.run([str(HERE.parent.parent/'4/codex/compact_wires'),str(raw),str(txt),'22'],check=True,capture_output=True)
    d=km.Circuit.from_file(str(txt));raw.unlink();txt.unlink()
    return d

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',default='results/workspace_search.jsonl');a=p.parse_args()
    best=190711.95;start=time.time()
    with open(a.output,'w') as log:
      for limit in range(78,101,2):
       for keep in [0,2,4,6]:
        for product in [8,12,16,20,None]:
         for splits in [(17,8),(18,9),(19,10),(20,11)]:
          try:
            c=build(limit,keep,product,splits)
            c=compact(c,HERE/'candidates/workspace_temp.kmb')
            T=c.operation_statistics().num_toffoli_like
            # The same numerical phase and rotation set give the same synthesis margin.
            cost=c.num_qubits*(T-11+250.85021)
            row=dict(limit=limit,keep=keep,product=product,splits=splits,Q=c.num_qubits,T=T,predicted_cost=cost)
            log.write(json.dumps(row)+'\n');log.flush()
            if cost<best:
                best=cost;c.to_file(str(HERE/'candidates/workspace_best.kmb'),format='kmb')
                (HERE/'results/workspace_best.json').write_text(json.dumps(row,indent=2)+'\n')
                print(json.dumps(dict(**row,seconds=time.time()-start)),flush=True)
          except ValueError:pass
