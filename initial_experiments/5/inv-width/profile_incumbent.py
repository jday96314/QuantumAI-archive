"""Rebuild selected circuit; annotate allocations and inspect compactor liveness."""
import collections, hashlib, inspect, json, re, sys
from pathlib import Path
import kickmix as km
import zero_head_affine_invsqrt as chosen
from optimize_workspace import compact

HERE=Path(__file__).resolve().parent
real_builder=km.CircuitBuilder
active=None
class TracedBuilder:
 def __init__(self):
  global active
  self.b=real_builder();active=self
 def __getattr__(self,n):return getattr(self.b,n)
 def alloc_qubits(self,n):
  out=self.b.alloc_qubits(n)
  f=sys._getframe(1)
  source=inspect.getframeinfo(f).code_context
  line=(source[0].strip() if source else '')
  name=line.split('=')[0].strip() if '=' in line else f.f_code.co_name
  if f.f_code.co_name=='<listcomp>':name='coefficient_or_counter'
  meta={k:f.f_locals[k] for k in ['j','g','level'] if k in f.f_locals and isinstance(f.f_locals[k],int)}
  label=f'{f.f_code.co_name}:{name}:line{f.f_lineno}:{meta}'
  if n:self.b.info_tag('ALLOC '+json.dumps([label,[str(q) for q in out]]))
  return out
 def create_quantum_register(self,n,name):
  out=self.b.create_quantum_register(n,name);self.b.info_tag('ALLOC '+json.dumps([name,[str(q) for q in out]]));return out

names={'normalize','lookup','horner','shift','restore_affine','release_curve','restore_curve','load_delayed','integrate','accumulate','erase','add'}
def trace(frame,event,arg):
 if active is None or frame.f_code.co_name not in names:return
 if not str(frame.f_code.co_filename).endswith('.py') and frame.f_code.co_filename!='<string>':return
 name=frame.f_code.co_name
 if event=='call':
  fields={k:v for k,v in frame.f_locals.items() if k in ['j','inv','inverse','subtract'] and isinstance(v,(bool,int))}
  active.b.info_tag('@enter:'+name+json.dumps(fields,sort_keys=True))
  if name in ['add','accumulate']:
   state={k:([str(q) for q in v] if isinstance(v,km.array) else str(v)) for k,v in frame.f_locals.items() if k in ['lhs','target','source','coordinate','carry_in','carry_workspace','borrow','auxiliary']}
   active.b.info_tag('STATE '+json.dumps(state))
 elif event=='return':active.b.info_tag('@exit:'+name)

def analyze(path):
 lines=path.read_text().splitlines();ref=lambda s:[int(x) for x in re.findall(r'(?<!\S)q(\d+)',s)] if not s.startswith(('INFO_TAG','REGISTER','APPEND_TO_REGISTER')) else []
 last={q:i for i,s in enumerate(lines) for q in ref(s)}
 live=set(range(22));labels={q:'x' for q in live};peak=[];peak_n=22;depth=0;stack=[];counts=collections.Counter();stage='setup';stage_spans=[];state={}
 for i,s in enumerate(lines):
  if s.startswith('INFO_TAG'):
   tag=s.split(' ',1)[1][1:-1].replace('\\Q','"')
   if tag.startswith('ALLOC '):
    label,qs=json.loads(tag[6:]);labels.update({int(q[1:]):label for q in qs})
   elif tag.startswith('STATE '):state=json.loads(tag[6:])
   elif tag.startswith('@enter:'):
    name=tag[7:];stack.append(name)
    if len(stack)==1:stage=name;stage_spans.append([stage,i,None])
   elif tag.startswith('@exit:'):
    if len(stack)==1:stage_spans[-1][2]=i+1
    if stack:stack.pop()
    if not stack:stage='setup'
  ids=ref(s);live.update(ids)
  if len(live)>peak_n:
   peak_n=len(live);peak=[dict(line=i,stage=stage,stack=stack.copy(),instruction=s,arithmetic_state=state.copy(),live={str(q):labels.get(q,'unknown') for q in sorted(live)})]
  elif len(live)==peak_n and ids and len(peak)<10:peak.append(dict(line=i,stage=stage,stack=stack.copy(),instruction=s,arithmetic_state=state.copy(),live={str(q):labels.get(q,'unknown') for q in sorted(live)}))
  if s.startswith(('CCX ','CCZ ','CSWAP ')):counts[stage]+=1
  reset=depth==0 and s.startswith(('HMR ','R '))
  for q in ids:
   if q>=22 and (reset or last[q]==i):live.discard(q)
  if s.startswith('PUSH_CONDITION'):depth+=1
  if s.startswith('POP_CONDITION'):depth-=1
 return dict(peak=peak_n,peak_examples=peak,stage_static_toffolis=dict(counts),stage_spans=stage_spans)

if __name__=='__main__':
 km.CircuitBuilder=TracedBuilder;sys.setprofile(trace)
 try:c=chosen.build()
 finally:sys.setprofile(None);km.CircuitBuilder=real_builder
 c.to_file('results/incumbent_tagged.kmx',format='kmx')
 report=analyze(HERE/'results/incumbent_tagged.kmx')
 d=compact(chosen.build(),HERE/'candidates/incumbent_rebuilt.kmb');d.to_file('candidates/incumbent_rebuilt.kmb',format='kmb')
 digest=hashlib.sha256(Path('candidates/incumbent_rebuilt.kmb').read_bytes()).hexdigest();expected=json.loads(Path('reference/zero_head_affine_verify.json').read_text())
 report.update(Q=d.num_qubits,static_toffolis=d.operation_statistics().num_toffoli_like,sha256=digest,selected_byte_match=digest==expected['circuit_sha256'])
 Path('results/incumbent_liveness.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:v for k,v in report.items() if k not in ['peak_examples','stage_spans']},indent=2))
 print(json.dumps(collections.Counter(report['peak_examples'][0]['live'].values()),indent=2))
