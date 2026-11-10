"""Resource screen followed by native smoke checks for borrowed phase bits."""
import hashlib,json,subprocess,sys,time
from pathlib import Path
import numpy as np
import binary_centered_invsqrt as p
from optimize_workspace import compact
from space_increment import install

root=Path(__file__).resolve().parent;install();p.table_for=lambda *args:np.load(root/'candidates/margin_curve8_t5.9e-08.npy')
rows=[];began=time.time();seen=set()
for carry in [1,2]:
 for limit in [84,86,88,90,92,94,96,98]:
  for split in ([(18,9)] if carry==1 else [(14,6),(15,7),(16,7),(16,8),(17,8),(18,9)]):
   for product in [8,12,16,17,None]:
    try:
     c=p.build(7,(27,18,8),lookup_frac=8,accum_bits=[*split,None],carry_bits=carry,rounded=True,booth=True,reuse_normalized_top=True,output_extra_bits=10,defer_workspace=True,qubit_limit=limit,curve_keep_bits=4,defer_curve_phase=True,product_limit=product,borrow_integration=True,monomial_delete=True,borrow_phase_headroom=True,unsigned_carries=True)
     c=compact(c,root/'candidates/headroom_search_temp.kmb');T=c.operation_statistics().num_toffoli_like;estimate=c.num_qubits*(T-11+250.85021)
     row=dict(carry=carry,limit=limit,split=split,product=product,Q=c.num_qubits,T=T,estimate=estimate)
     if estimate<189960:
      name=f'headroom_c{carry}_q{limit}_s{split[0]}_{split[1]}_p{product}'
      path=root/f'candidates/{name}.kmb';c.to_file(str(path),format='kmb');digest=hashlib.sha256(path.read_bytes()).hexdigest()
      if digest in seen:continue
      seen.add(digest);row.update(path=str(path),sha256=digest,tag=name)
     rows.append(row)
    except ValueError:pass
(root/'results/headroom_resource_search.json').write_text(json.dumps(rows,indent=2)+'\n')
ranked=sorted([r for r in rows if 'path' in r],key=lambda r:r['estimate']);print(json.dumps(dict(candidates=len(ranked),best=ranked[:3],seconds=time.time()-began)),flush=True)
valid=[]
for row in ranked[:24]:
 result=subprocess.run([sys.executable,str(root/'../../../from_host/verifiers/invsqrt/verify_submission_invsqrt.py'),row['path'],'--max-x','4096','--max-workers','2'],capture_output=True,text=True)
 try:report=json.loads(result.stdout)
 except ValueError:report=dict(status='ERROR',stderr=result.stderr[-2000:])
 (root/f"results/{row['tag']}_smoke.json").write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps(dict(tag=row['tag'],estimate=row['estimate'],status=report['status'],error=report.get('error'),seconds=time.time()-began)),flush=True)
 if report['status']=='PASS':valid.append(row);break
(root/'results/headroom_smoke_selection.json').write_text(json.dumps(dict(valid=valid,scope='4096 inputs, 64 shots each; not full default verification'),indent=2)+'\n')
