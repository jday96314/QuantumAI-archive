"""Native checks and profiles for the parameter frame, same small-oracle contract."""
import argparse,inspect,json
from pathlib import Path
import lookahead_jacobi as lookahead
import compact_lookahead_parameters as compact
import check_lookahead_jacobi as component
import check_jacobi_oracle as wrapper
import check_lookahead_oracle as oracle
HERE=Path(__file__).resolve().parent

def install():
 lookahead.parameters=compact.parameters
 source=inspect.getsource(component.model)
 old='p=(m.gmul(a,c)^bg) or 1';new='p=(b if flag else m.gmul(a,c or 1)) or 1'
 assert source.count(old)==1;source=source.replace(old,new);ns=dict(component.__dict__);exec(source,ns);component.model=ns['model']
 (HERE/'parameters/compact_lookahead_model_source.py').write_text(source)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--profile',action='store_true');p.add_argument('--oracle',action='store_true');a=p.parse_args();install()
 if a.profile:
  wrapper.EMIT_EXTRACT=lookahead.extract;wrapper.EMIT_SYNTHESIZE=lookahead.synthesize;wrapper.EMIT_ARITHMETIC=lookahead.Arithmetic
  *_,report=wrapper.emit(70,count_only=True);report['scope']='leading-zero chart with compact in-head parameters; second-pivot interior gaps unsupported; emitted counts'
  (HERE/'results/compact_lookahead_full_profile.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
 elif a.oracle:
  wrapper.EMIT_EXTRACT=lookahead.extract;wrapper.EMIT_SYNTHESIZE=lookahead.synthesize;wrapper.EMIT_ARITHMETIC=lookahead.Arithmetic
  workspace=HERE/'compact_lookahead'
  for sub in ('parameters','candidates','results'):(workspace/sub).mkdir(parents=True,exist_ok=True)
  cf=wrapper.cf;cf.HERE=workspace;cf.TRANSLATOR=HERE/'../../4/codex/translate_prefix';cf.emit=wrapper.emit;cf.model=oracle.model
  cf.test(14,127,4096,'agent20-compact-lookahead-small-a')
 else:
  check=component.check;check.MODE='compact_lookahead';check.FORCE_ZERO_HEAD=True
  check.block.block=lookahead.block;check.block.CubeArithmetic=lookahead.Arithmetic;check.model=component.model;check.pair=component.check_jacobi.pair
  rows=[]
  for n in (6,14,34):
   row=check.check(n);rows.append(row);(HERE/'results/compact_lookahead_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
