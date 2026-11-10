"""Measure emitted stage counts without changing the approximation model."""
import inspect,json,collections
from pathlib import Path
import kickmix as km
import fused_invsqrt as base
src=inspect.getsource(base.build).replace('b.info_tag(','mark(')
needle=' b=km.CircuitBuilder();';assert needle in src
src=src.replace(needle,' b=km.CircuitBuilder();',1)
line=src.index('\n',src.index(needle));src=src[:line+1]+" def mark(name):\n  b.info_tag(name)\n"+src[line+1:]
ns=dict(vars(base));exec(src,ns);rows=[]
for product in [14,17,20]:
 c=ns['build'](90,product);stage='initial';counts=collections.Counter()
 for line in str(c).splitlines():
  op=line.split(' ',1)[0]
  if op=='INFO_TAG':
   for name in ['normalize_forward','square_forward','coefficient_lookup','curve_forward','release_curve','slope_lookup','linear_forward','affine_forward','phase','affine_reverse','linear_reverse','curve_replay','curve_reverse','square_erase','coefficient_erase','normalize_reverse']:
    if line in [f'INFO_TAG "@enter:{name}"',f'INFO_TAG "@exit:{name}"',f'INFO_TAG "{name}"']:stage=name
  elif op in ['CCX','CCZ','CSWAP']:counts[stage]+=1
 row=dict(product=product,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,stages=dict(counts));rows.append(row);print(json.dumps(row),flush=True)
Path('results/fused_stage_profiles.json').write_text(json.dumps(rows,indent=2)+'\n')
