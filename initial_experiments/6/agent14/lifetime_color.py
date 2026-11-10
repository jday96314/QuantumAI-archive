"""Recolor zero-to-zero scratch lifetimes, splitting at unconditional HMR.

Public input/output lifetimes stay pinned to their register wires. Explicit
zero-input/output ports may be borrowed outside their live segment. A split
is made only at an unconditional HMR, whose poststate is zero. No gate,
classical condition, measurement, or arithmetic dependency is removed.
"""
import heapq,re
import kickmix as km
QRE=re.compile(r'\bq(\d+)\b')


def compact(c,zero_inputs=(),zero_outputs=()):
 headers=[];lines=[]
 for line in str(c).splitlines():
  if line.startswith(('REGISTER','APPEND_TO_REGISTER')):headers.append(line)
  else:lines.append(line)
 end=len(lines)+1;segments={};current={}
 for time,line in enumerate(lines):
  qs=list(dict.fromkeys(int(x) for x in QRE.findall(line)))
  for q in qs:
   if q not in current:
    item=dict(old=q,start=time,end=time+1,color=None);segments.setdefault(q,[]).append(item);current[q]=item
   current[q]['end']=time+1
  if line.startswith('HMR ') and ' if ' not in line:
   for q in qs:current.pop(q,None)
 public={q.id for name,reg in c.register_data.items() for q in reg if str(q).startswith('q')}
 zi={q.id for name in zero_inputs for q in c.register_data[name]};zo={q.id for name in zero_outputs for q in c.register_data[name]};pinned={}
 for q in public:
  ss=segments.get(q,[])
  if not ss:
   if q not in zi or q not in zo:pinned.setdefault(q,[]).append((0,end))
   continue
  if q not in zi:ss[0]['start']=0;ss[0]['color']=q
  if q not in zo:ss[-1]['end']=end;ss[-1]['color']=q
  for s in ss:
   if s['color'] is not None:pinned.setdefault(q,[]).append((s['start'],s['end']))
 allseg=[s for ss in segments.values() for s in ss];allseg.sort(key=lambda s:(s['start'],s['color'] is None,-s['end']))
 busy={};peak_public=max(public,default=-1)+1;max_color=peak_public-1
 def overlaps_pin(color,a,b):return any(x<b and a<y for x,y in pinned.get(color,()))
 for s in allseg:
  a,b=s['start'],s['end'];color=s['color']
  if color is not None:
   assert busy.get(color,-1)<=a,(s,busy.get(color));busy[color]=b;max_color=max(max_color,color);continue
  for color in range(max_color+2):
   if busy.get(color,-1)<=a and not overlaps_pin(color,a,b):break
  s['color']=color;busy[color]=b;max_color=max(max_color,color)
 # Rewrite each gate using the segment that contains that reference.
 pointers={q:0 for q in segments};out=[]
 for time,line in enumerate(lines):
  def replace(match):
   q=int(match[1]);ss=segments[q];i=pointers[q]
   while i+1<len(ss) and time>=ss[i]['end']:i+=1
   pointers[q]=i;s=ss[i];assert s['start']<=time<s['end'],(q,time,s);return f'q{s["color"]}'
  out.append(QRE.sub(replace,line))
 result=km.Circuit('\n'.join(headers)+'\n'+'\n'.join(out)+'\n');report=dict(original_Q=c.num_qubits,new_Q=result.num_qubits,segments=len(allseg),public_qubits=len(public),original_T=c.operation_statistics().num_toffoli_like,new_T=result.operation_statistics().num_toffoli_like)
 assert report['original_T']==report['new_T'];assert result.register_data==c.register_data
 return result,report
